/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"definition.h"
#include"structure.h"
#include"prototype.h"

/* forest processes ********************************************/
void plant_process(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pflx *flux, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double aaa, bbb, emerge, cost_new_leaf;

	/* plant internal carbon */
	/* leaf emergence ****************/
	if(pchar->season == 2){
			
		if(grid->veg_type == 10){
			/* grasslands */
			emerge = 0.04 * mass->nsch_storage; //Tongyu 
			//emerge = 0.04 * mass->nsch_storage; //KBU
		}else if(grid->veg_type == 3){
			emerge = 0.025 * mass->nsch_storage;
		}else if(grid->veg_type == 4){
			emerge = 0.04 * mass->nsch_storage; //TKY (2013/01/24 by A.Ito)
			//emerge = 0.08 * mass->nsch_storage;
			//emerge = 0.03 * mass->nsch_storage;  //GDK
			//emerge = 0.03 * mass->nsch_storage;  //SMF
		}else{
			emerge = 0.05 * mass->nsch_storage;
		}
	
		/* new leaf expansion */
		if(emerge > 0.0){
			/** construction cost of new leaves **/
			cost_new_leaf = emerge * 0.10; /* */
			/* cost_new_leaf = 0.0; */
		
			mass->fol += emerge - cost_new_leaf;
			
			mass->nsch_storage -= emerge;
						
		}else{
			cost_new_leaf = 0.0;
		}
	}else{
		cost_new_leaf = 0.0;
	}
    
    /* planting (crop with stage==1) ****************************/
	if(pchar->crop_stage == 1){
		mass->fol += 0.07;
		mass->stm += 0.05;
		mass->rot += 0.05;
	}
	
	/* leaf shedding ****************/
	if(pchar->season == 3){
		/* litter */
		if(mass->lai>0.05){
			flux->lf = pchar->dcd*mass->fol; /* leaf-shedding */
			pchar->lf =  pchar->dcd;
		}else if(mass->lai<=0.05){
			/* flux->lf = flf(grid, pchar, mass); */
			flux->lf = mass->fol;
			pchar->lf = 1.0;
		}		
	}else{
		flux->lf = f_lf(grid, pchar, mass);
		
		/* drought-induced leaf shedding: 2011/10/25 by A.Ito */
		if(loct->pot_total_h < -148.075 && mass->lai > 0.05 && loct->tmp_2m > -5.0){
			/* KBU */
            //flux->lf += mass->fol*0.009;
	          
            /* Tongyu */
            if(pchar->phototype==3){
                flux->lf += mass->fol*0.012;
            }else if(pchar->phototype==4){
                flux->lf += mass->fol*0.012;
            } /* */
		}
	}
	
	/* stem mortality */
	flux->lc = f_lc(grid, pchar, mass);
	/* root mortality */
	flux->lr = f_lr(grid, pchar, mass);
	
	if(SCI_SCHEME==1){
		/* stable carbon isotope */
		flux->d13c_lf = mass->d13c_fol;
		flux->d13c_lc = mass->d13c_stm;
		flux->d13c_lr = mass->d13c_rot;
	}

	/* litter fall */
	mass->fol -= flux->lf;
	mass->stm -= flux->lc;
	mass->rot -= flux->lr;
	mass->lai = lai_mass(mass, pchar);

	/* photosynthesis, gross primary production *****************/
	/* Monsi-Saeki */
	flux->gpp = f_gpp(grid, loct, pchar, mass);

	/* DePury-Farquhar */
	/* flux->gpp_df97 = f_gpp_df97(grid, loct, pchar, mass, flux); */
	/* flux->gpp = flux->gpp_df97; **************************************************/
	
	if(SCI_SCHEME==1){
		/* stable carbon isotope */
		flux->d13c_gpp = loct->d13C_aCO2-pchar->photo_13c_frac;
	}
	
	/* maintenance respirations */
	flux->rfm = f_rfm(grid, loct, pchar, mass);
	flux->rcm = f_rcm(grid, loct, pchar, mass);
	flux->rrm = f_rrm(grid, loct, pchar, mass);
	flux->rpm = flux->rfm + flux->rcm + flux->rrm;
	
	if(SCI_SCHEME==1){
		/* stable carbon isotope */
		flux->d13c_rfm = mass->d13c_fol;
		flux->d13c_rcm = mass->d13c_stm;
		flux->d13c_rrm = mass->d13c_rot;
	}
	
	/* tentative primary production */	
	flux->epp = flux->gpp - flux->rpm;
	
	/* translocation of photosynthate */
	f_allocation(pchar, mass, flux);
	
	if(SCI_SCHEME==1){
		/* stable carbon isotope */
		flux->d13c_tpf = flux->d13c_gpp;
		flux->d13c_tpc = flux->d13c_gpp;
		flux->d13c_tpr = flux->d13c_gpp;
		flux->d13c_tpp = flux->d13c_gpp;
	}
	
	if(flux->epp > 0.0){
		/* growth construction respiration */
		flux->rfg = f_rfg(grid, pchar, flux) + cost_new_leaf;
		flux->rcg = f_rcg(grid, pchar, flux);
		flux->rrg = f_rrg(grid, pchar, flux);
	}else if(flux->epp <= 0.0){
		flux->rfg = flux->rcg = flux->rrg = 0.0;
	}
	
	if(SCI_SCHEME==1){
		/* stable carbon isotope */
		flux->d13c_rfg = mass->d13c_fol;
		flux->d13c_rcg = mass->d13c_stm;
		flux->d13c_rrg = mass->d13c_rot;
	}
	
	if(SCI_SCHEME==1){
		/* stable carbon isotope */
		if(flux->epp>0.0){
			mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, flux->d13c_tpf, flux->tpf);
			mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, flux->d13c_tpc, flux->tpc);
			mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, flux->d13c_tpr, flux->tpr);
		}else{
			mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, 
				flux->d13c_gpp, pchar->malloc_f*flux->gpp);
			mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, 
				flux->d13c_gpp, pchar->malloc_c*flux->gpp);
			mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, 
				flux->d13c_gpp, pchar->malloc_r*flux->gpp);
		}
	}
	
	/* partitioning of photosynthate *****************************/
	/* C budget of storage: 2011/10/18 by A.Ito */
	/* leaf */
	mass->fol += (flux->tpf - flux->rfg);
	
	/* stem */
	if(flux->tpc > 0.0){
		if(mass->nsch_storage < pchar->max_nsch_storage){
			mass->nsch_storage += flux->tpc;
		}else{
			mass->stm += flux->tpc;
		}
	}else{
		mass->stm += flux->tpc;
	}
	mass->stm -= flux->rcg;
	
	/* root */
	if(flux->tpr > 0.0){
		if(mass->nsch_storage < pchar->max_nsch_storage){
			mass->nsch_storage += flux->tpr;
		}else{
			mass->rot += flux->tpr;
		}
	}else{
		mass->rot += flux->tpr;
	}
	mass->rot -= flux->rrg;
    
    mass->grain += flux->tpg;
	if(pchar->crop_stage == 5){
		flux->hvst = mass->grain;
		mass->grain = 0.0;
	}
    
	/* if(pchar->phenoltype==1){ */
		reallocation_survival(grid, pchar, mass); /* */
	/* } */
	
	/* update LAI */
	mass->lai = lai_mass(mass, pchar);
	
	/* plant growth respiration */
	flux->rpg = flux->rfg + flux->rcg + flux->rrg;
	/* plant respiration */
	flux->rp =  flux->rpm + flux->rpg;
	/* surplus production */
	flux->spp = flux->gpp - flux->rfm - flux->rfg;
	/* net primary production */
	flux->npp = flux->gpp -flux->rp;

	/* plant litterfall */
	flux->lL = flux->lf + flux->lc + flux->lr;
	
	if(SCI_SCHEME==1){
		flux->d13c_lL = d13c_addition3(
			flux->d13c_lf, flux->lf, flux->d13c_lc, flux->lc, flux->d13c_lr, flux->lr);
	}
	
	/* daily plant mass */	
	mass->plant = mass->fol + mass->stm + mass->rot;
	
	/* nitrogen cycle */
	if(N_CYCLE ==1){
		/* n allocation */
		f_n_alloc(grid, loct, pchar, mass, flux);

		/* n reallocation */
		f_n_realloc(grid, loct, pchar, mass, flux);
		
		/* abandon by litterfall */
		f_n_abandon_salvage(grid, loct, pchar, mass, flux);
	}

}

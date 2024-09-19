/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 21 by Akihiko Ito: Version 0.2.0.
	modify: 2007 / 05 / 21 by Akihiko Ito: Version 0.2.4.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* plant processes ***************************************************/
void plant_process(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pflx *flux, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float emerge, cost_new_leaf;

	/* plant internal carbon ******************************************/
	/* leaf emergence */
	if(mass->season == 2){
		
		if(grid->veg_type == 10){
			/* grasslands */
			emerge = 0.08 * mass->nsch_storage;
		}else if(grid->veg_type == 12 || grid->veg_type == 14){
			/* crop lands */
			emerge = 0.08 * mass->nsch_storage;
		}else if(grid->veg_type == 3){
			emerge = 0.05 * mass->nsch_storage;
		}else if(grid->veg_type == 4){
			emerge = 0.05 * mass->nsch_storage;
		}else{
			emerge = 0.05 * mass->nsch_storage;
		}
        emerge *= 1.0 / (double)DSTEP;  /* hourly rate */
		
		/* new leaf expansion */
		if(emerge > 0.0){
			/** construction cost of new leaves **/
			cost_new_leaf = emerge*0.10; /* */
			/* cost_new_leaf = 0.0; */
			
			mass->fol += emerge - cost_new_leaf;
			
			mass->nsch_storage -= emerge;
			
		}else{
			cost_new_leaf = 0.0;
		}
	}else{
		cost_new_leaf = 0.0;
	}
	
	/* leaf shedding ************************************************/
	if(mass->season == 3){
		/* litter */
		if(mass->lai > 0.05){
			flux->lf = pchar->dcd*mass->fol; /* leaf-shedding */
			pchar->lf =  pchar->dcd;
		}else if(mass->lai <= 0.05){
			/* flux->lf = flf(grid, pchar, mass); */
			flux->lf = mass->fol;
			pchar->lf = 1.0;
		}		
	}else{
		flux->lf = f_lf(grid, pchar, mass);
		
		/* drought-induced leaf shedding: 2011/10/25 by A.Ito */
		if(loct->pot_total_h < -148.075 && mass->lai > 0.1 && loct->tmp_2m > -5.0){
			/* flux->lf += mass->fol*0.001; */
			//flux->lf += mass->fol*0.0005;  /* 2012/05/23 by A.Ito */
			//flux->lf += mass->fol*0.0006;  /* 2013/09/06 by A.Ito */
			flux->lf += mass->fol*0.0001;  /* 2013/10/15 by A.Ito */
		}
	}
	
    /* stem litter */
	flux->lc = f_lc(grid, pchar, mass);
    /* root litter */
	flux->lr = f_lr(grid, pchar, mass);
	
	if(SCI_CALC==1){
		/* stable carbon isotope */
		flux->d13c_lf = mass->d13c_fol;
		flux->d13c_lc = mass->d13c_stm;
		flux->d13c_lr = mass->d13c_rot;
	}

	/* litter fall */
	mass->fol -= flux->lf;
	mass->lai = lai_mass(mass, pchar); /* LAI update */
	mass->stm -= flux->lc;
	mass->rot -= flux->lr;  

	/* photosynthesis, gross primary production *******************************/
	if(mass->season == -1){
		// flux->gpp = f_gpp(grid, loct, pchar, mass) * 0.2;	/* 0.2: assumption of urban greenness */
	}else{
		// flux->gpp = f_gpp(grid, loct, pchar, mass);
	}
    
    /**/
    if(FIX_LAI == 1){
        //pchar->lai = 1.0;
        //pchar->lai = 2.0;
        //pchar->lai = 3.0;
        pchar->lai = 4.0;
    }else{
        pchar->lai = mass->lai;
    }
    pchar->tmp = loct->tmp_2m;
    pchar->vpd_a = loct->vpd;
    
    /***************************************** initialization */
    pchar->ci_sn = 0.7 * loct->aCO2 / 10.0;
    pchar->ci_sd = 0.7 * loct->aCO2 / 10.0;
       
    /* leaf parameters */
    f_leaf_prop(grid, loct, pchar);
    
    /* canopy radiation parameters */
    f_apar_sunshade(grid, loct, pchar);
    
    /* canopy photosynthesis */
    f_gpp_sunshade(mass->season, grid, loct, pchar);
    
    /* unit cobversion */
	flux->gpp = lTs * (pchar->p_sn + pchar->p_sd);
	//flux->gpp = (pchar->photocap_n + pchar->photocap_n);
    
    if(SCI_CALC == 1){
		/* stable carbon isotope */
		flux->d13c_gpp = loct->d13C_aCO2 - pchar->photo_13c_frac;
	}
	
	/* respiuration ******************************************/
	/* maintenance respirations */
    /* leaf */
	flux->rfm = f_rfm(grid, loct, pchar, mass);
    /* stem */
	flux->rcm = f_rcm(grid, loct, pchar, mass);
    /* root */
	flux->rrm = f_rrm(grid, loct, pchar, mass);
    
    /* for debugging */
    //flux->rfm = flux->rfm;
    //flux->rcm = flux->gpp;
    //flux->rrm = flux->gpp;
    
    /* total */
	flux->rpm = flux->rfm + flux->rcm + flux->rrm;

	if(SCI_CALC == 1){
		/* stable carbon isotope */
		flux->d13c_rfm = mass->d13c_fol;
		flux->d13c_rcm = mass->d13c_stm;
		flux->d13c_rrm = mass->d13c_rot;
	}
	
	/* EPP: effective primary production */	
	flux->epp = flux->gpp - flux->rpm;
	
	/* translocation of photosynthate */
    f_allocation(pchar, mass, flux);
    //flux->tpf = flux->tpc = flux->tpr = 0.0;
	
	if(SCI_CALC == 1){
		/* stable carbon isotopic ratio of translocation */
		flux->d13c_tpf = flux->d13c_gpp;
		flux->d13c_tpc = flux->d13c_gpp;
		flux->d13c_tpr = flux->d13c_gpp;
		flux->d13c_tpp = flux->d13c_gpp;
	}
	
	if(flux->epp > 0.0){
		/* growth construction respiration */
        /* leaf */
		flux->rfg = f_rfg(grid, pchar, flux);
        /* stem */
		flux->rcg = f_rcg(grid, pchar, flux);
        /* root */
		flux->rrg = f_rrg(grid, pchar, flux);
	}else if(flux->epp <= 0.0){
        /* null/negative growth: no growth respiration */
		flux->rfg = flux->rcg = flux->rrg = 0.0;
	}
    //flux->rfg = flux->rcg = flux->rrg = 0.0;
	
	if(SCI_CALC == 1){
		/* stable carbon isotope */
		flux->d13c_rfg = mass->d13c_fol;
		flux->d13c_rcg = mass->d13c_stm;
		flux->d13c_rrg = mass->d13c_rot;
	}
	
	if(SCI_CALC == 1){
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
	
	/* partitioning of photosynthate */
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
	
	/* if(pchar->phenoltype==1){ */
	//	reallocation_survival(grid, pchar, mass);
	/* } */
	
	/* update LAI */
	mass->lai = lai_mass(mass, pchar);
    
    /* summation */
	flux->rpg = flux->rfg + flux->rcg + flux->rrg;
	flux->rp = flux->rpm + flux->rpg;
	flux->spp = flux->gpp - flux->rfm - flux->rfg;
	flux->npp = flux->gpp - flux->rp;
	//flux->npp = flux->rpm;
    /* total litterfall */
	flux->lL = flux->lf + flux->lc + flux->lr;
	
	if(SCI_CALC==1){
		flux->d13c_lL = d13c_addition3(
			flux->d13c_lf, flux->lf, flux->d13c_lc, flux->lc, flux->d13c_lr, flux->lr);
	}
	
	/* monthly total biomass */	
	mass->plant = mass->fol + mass->stm + mass->rot;
	//mass->plant = mass->fol;
    
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

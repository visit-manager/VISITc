/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* allocation of photosynthate **********************************/
void f_allocation(
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	double aaa, bbb, cc1, ccc, ddd;
	double alloc_f, alloc_c, alloc_r, alloc_g; /* allocation ratios */
	
	/* during growing-period ****************************************/	
	if(flux->epp > 0.0){ 
		/* EPP = (GPP - maintenance respiration) is positive */
		/* allocation ratios of EPP */
		if(mass->lai > pchar->opt_lai){ 
			/* if holding LAI is greater than the optimam one */
			//ccc = pchar->alloc_ass*1.0;
			ccc = 0.0;
			alloc_f = ccc;
			alloc_c = (1.0-ccc)*pchar->alloc_abg;
			alloc_r = (1.0-ccc)*(1.0-pchar->alloc_abg);
		}else if(mass->lai <= pchar->opt_lai){ 
			/* if holding LAI is smaller than the optimam one */
			aaa = (pchar->opt_lai - mass->lai)*100.0*2.0/2.2/pchar->sla;
			bbb = flux->epp*pchar->alloc_ass;
			
			/* allocate photosyntahte to foliage to attain the optimum one */
			if(aaa<=bbb){
				ccc = pchar->alloc_ass;
				
			}else if(aaa>bbb){
				cc1 = aaa/bbb;
				
				/* maximum allocation ratio to foliage is 5 % */
				ccc = ((pchar->alloc_ass*cc1)<0.05)?(pchar->alloc_ass*cc1):0.05;
				
                //ccc = ((pchar->alloc_ass*cc1)<0.12)?(pchar->alloc_ass*cc1):0.12;
			}

			/* ccc = pchar->alloc_ass*1.0; */
			
			if(pchar->season==0){
				ddd = 0.7;
			}else{
				ddd = 1.0;
			}
			
			alloc_f = ccc*ddd;
			alloc_c = (1.0-ccc*ddd)*pchar->alloc_abg;
			alloc_r = (1.0-ccc*ddd)*(1.0 - pchar->alloc_abg);
		}
        
        /* allocation to grain **/
		if(pchar->crop_stage == 4){
			alloc_g = 0.7;
		}else{
			alloc_g = 0.0;
		}
				
		/* monthly translocation fluxes */
		flux->tpp = (alloc_g + (1.0 - alloc_g) * (alloc_f + alloc_c + alloc_r))*flux->epp;
		flux->tpf = (1.0 - alloc_g)*alloc_f*flux->epp;
		flux->tpc = (1.0 - alloc_g)*alloc_c*flux->epp;
		flux->tpr = (1.0 - alloc_g)*alloc_r*flux->epp;
        
		/**/
		flux->tpg = alloc_g * flux->epp;
		
	}else if(flux->epp<=0.0){ /* during NON growing-period */
		/* allocation ratios of GPP, not EPP */
		alloc_f = pchar->alloc_ass;
		alloc_c = (1.0 - pchar->alloc_ass)*pchar->alloc_abg;
		alloc_r = (1.0 - pchar->alloc_ass)*(1.0 - pchar->alloc_abg);
				
		/* daily translocation fluxes */
		flux->tpp = (alloc_f + alloc_c + alloc_r)*flux->epp;
		flux->tpf = alloc_f*flux->gpp - flux->rfm;
		flux->tpc = alloc_c*flux->gpp - flux->rcm;
		flux->tpr = alloc_r*flux->gpp - flux->rrm;

		//flux->tpf = alloc_f*flux->epp;
		//flux->tpc = alloc_c*flux->epp;
		//flux->tpr = alloc_r*flux->epp;
        
        /* if(flux->tpc < 0.0){
            printf("### %lf %lf %lf %lf\n", flux->tpc, alloc_c, flux->gpp, flux->rcm);
        } */
        
        alloc_g = 0.0;
		flux->tpg = 0.0;
	}
	
	/* monthly partitioning ratios */
	pchar->malloc_f = alloc_f;
	pchar->malloc_c = alloc_c;
	pchar->malloc_r = alloc_r;
    pchar->malloc_g = alloc_g;
}

/* plant-internal mass re-adjustment **********************************/
void reallocation_survival(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass,
    struct Pflx *flux
){
	double crit_lai;
	double ral_cap_stf, ral_cap_rtf;
	double aaa, bbb, ccc;
	double ral_stf, ral_rtf;
	
	switch(grid->veg_type){
		case 1: case 2: case 3: case 4: case 5: case 6: case 7: 
			crit_lai = 0.1;
			break;
		case 8: case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16: 
			crit_lai = 0.05;
			break;
        default:
            crit_lai = 0.05;
	}
	
	ral_cap_stf = 0.05;
	ral_cap_rtf = 0.1;	
		
	/* to fliage: critical reallocation for survival */
	aaa = crit_lai * 100.0 * 2.0/2.2/pchar->sla;
	if(mass->fol < aaa){
		bbb = mass->stm * ral_cap_stf;
		ccc = mass->rot * ral_cap_rtf;
		
		ral_stf = aaa * pchar->alloc_abg * (bbb/aaa)/(0.5 + (bbb/aaa));
		ral_rtf = aaa * (1.0 - pchar->alloc_abg) * (ccc/aaa)/(0.5 + (ccc/aaa));

		mass->fol += ral_stf + ral_rtf;
		mass->stm -= ral_stf;
		mass->rot -= ral_rtf;
        
        flux->rtpc += ral_stf;
        flux->rtpr += ral_rtf;
	}
	
	/* to fliage, improve production */
	/* aaa = (pchar->opt_lai-mass->lai)*100.0*2.0/2.2/pchar->sla;
	if(aaa > 0.0){
		bbb = mass->stm*ral_cap_stf;
		ccc = mass->rot*ral_cap_rtf;
		
		ral_stf = aaa * pchar->alloc_abg*(bbb/aaa)/(1.5+(bbb/aaa));
		ral_rtf = aaa * (1.0-pchar->alloc_abg)*(ccc/aaa)/(1.5+(ccc/aaa));

		mass->fol += ral_stf+ral_rtf;
		mass->stm -= ral_stf;
		mass->rot -= ral_rtf;
	} */
	
	mass->lai = lai_mass(mass, pchar);
	
	/* to stem and branch */
	if(mass->stm <= 0.0){
		if(mass->rot > 0.0){
			aaa = mass->rot * pchar->alloc_abg;
			mass->rot -= aaa;
			mass->stm += aaa;
		}
	}
	
	/* to root */
	if(mass->rot <= 0.0){
		if(mass->stm > 0.0){
			aaa = mass->stm*(1.0 - pchar->alloc_abg);
			mass->rot += aaa;
			mass->stm -= aaa;
		}
	}
}


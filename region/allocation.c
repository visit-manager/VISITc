/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in December 20, 2007						*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2009 / 04 / 25 by Akihiko Ito: Version 1.1.0.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* allocation of photosynthate **********************************/
void allocation(
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	float aaa, bbb, cc1, ccc, ddd;
	float alloc_f, alloc_c, alloc_r; /* allocation ratios */
	
	/* during growing-period ****************************************/	
	if(flux->epp > 0.0){ 
		/* EPP = (GPP - maintenance respiration) is positive */
		/* allocation ratios of EPP */
		if(mass->lai > pchar->opt_lai){ 
			/* if holding LAI is greater than the optimam one */
			/**** avoid too high LAI under good condition: 2009/04/24 A.Ito ****/
			aaa = mass->lai - pchar->opt_lai;
			bbb = 1.0 - 0.5*aaa;
			bbb = (bbb>0.0)?bbb:0.0;
			
			/* formar: bbb = 1.0; */
			
			ccc = pchar->alloc_ass*bbb;
			alloc_f = ccc;
			alloc_c = (1.0 - ccc)*pchar->alloc_abg;
			alloc_r = (1.0 - ccc)*(1.0 - pchar->alloc_abg);
		}else{
			/* if holding LAI is smaller than the optimam one */
			aaa = (pchar->opt_lai - mass->lai) * 100.0*2.0/2.2/pchar->sla;
			bbb = flux->epp * pchar->alloc_ass;
			
			/* allocate photosynthate to foliage to attain the optimum one */
			if(aaa <= bbb){
				ccc = pchar->alloc_ass;
				
			}else{
				cc1 = aaa / bbb;
				
				/* maximum allocation ratio to foliage is 10 % */
				ccc = ((pchar->alloc_ass*cc1)<0.1)?(pchar->alloc_ass*cc1):0.10;
			}

			/* ccc = pchar->alloc_ass*1.0; */
			
			if(mass->season == 0){
				ddd = 0.7;
			}else{
				ddd = 1.0;
			}
			
			alloc_f = ccc*ddd;
			alloc_c = (1.0 - ccc * ddd) * pchar->alloc_abg;
			alloc_r = (1.0 - ccc * ddd) * (1.0 - pchar->alloc_abg);
		}
				
		/* monthly translocation fluxes */
		flux->tpp = (alloc_f + alloc_c + alloc_r)*flux->epp;
		flux->tpf = alloc_f * flux->epp;
		flux->tpc = alloc_c * flux->epp;
		flux->tpr = alloc_r * flux->epp;
		
	}else{ /* during NON growing-period */
		/* allocation ratios of GPP, not EPP */
		alloc_f = pchar->alloc_ass;
		alloc_c = (1.0 - pchar->alloc_ass) * pchar->alloc_abg;
		alloc_r = (1.0 - pchar->alloc_ass) * (1.0 - pchar->alloc_abg);
				
		/* translocation fluxes */
		flux->tpp = (alloc_f + alloc_c + alloc_r) * flux->epp;
		flux->tpf = alloc_f*flux->gpp - flux->rfm;
		flux->tpc = alloc_c*flux->gpp - flux->rcm;
		flux->tpr = alloc_r*flux->gpp - flux->rrm;
	}
	
	/* monthly partitioning ratios */
	pchar->malloc_f = alloc_f;
	pchar->malloc_c = alloc_c;
	pchar->malloc_r = alloc_r;
}

/* plant-internal mass re-adjustment **********************************/
void reallocation_survival(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float crit_lai;
	float ral_cap_stf, ral_cap_rtf;
	float aaa, bbb, ccc;
	float ral_stf, ral_rtf;
	
	crit_lai = 0.1;
	ral_cap_stf = 0.05;
	ral_cap_rtf = 0.1;	
		
	/* to fliage: critical reallocation for survival */
	aaa = crit_lai * 100.0*2.0/2.2/pchar->sla;
	if(mass->fol < aaa){
		bbb = mass->stm * ral_cap_stf;
		ccc = mass->rot * ral_cap_rtf;
		
		ral_stf = aaa * pchar->alloc_abg * (bbb/aaa)/(0.5+(bbb/aaa));
		ral_rtf = aaa * (1.0 - pchar->alloc_abg) * (ccc/aaa)/(0.5+(ccc/aaa));

		mass->fol += ral_stf + ral_rtf;
		mass->stm -= ral_stf;
		mass->rot -= ral_rtf;
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
			aaa = mass->stm * (1.0 - pchar->alloc_abg);
			mass->rot += aaa;
			mass->stm -= aaa;
		}
	}
}


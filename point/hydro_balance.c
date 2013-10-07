/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/*   initialize   hydrology  */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* additional message for hydrological scheme */
#define VERB_HYD 0

/* water condition *****************************************************/
void f_hydrology(
	struct Grid *grid,
	struct Loct *loct, 
	struct Echar *echar,
	struct Mass *mass
){
	double snow, rain, thaw;
	double gain, capt_water, baseflow, retran;
	double dry_inx;
	double aa, bb, cc;
	double cap_incep, lqd_in_soil;

	/* clear ******************/
	loct->pm_incep_tree = 0.0;
	loct->pm_incep_c3 = 0.0;
	loct->pm_incep_c4 = 0.0;
	loct->incep = 0.0;
	loct->incep_tree = 0.0;
	loct->incep_c3 = 0.0;
	loct->incep_c4 = 0.0;
	loct->pm_evpr = 0.0; 		
	loct->evpr = 0.0;			
	loct->pm_trnsp_tree = 0.0;	
	loct->pm_trnsp_c3 = 0.0; 	
	loct->pm_trnsp_c4 = 0.0; 	
	loct->trnsp = 0.0;
	loct->trnsp_tree = 0.0; 	
	loct->trnsp_c3 = 0.0; 	
	loct->trnsp_c4 = 0.0; 	
	loct->thaw = 0.0;
	loct->ro1 = 0.0;
	loct->ro2 = 0.0;

	/* PRECIPITATION: snow/rain portions  **************************/
	/* based on data in Mizu-Kankyo-no-Kisyogaku by J. Kondo */
	loct->snp = 1.0/(1.0+exp(0.75*(loct->tmp_2m-2.0)));
	snow = loct->snp*loct->prate_sfc;
	rain = (1.0 - loct->snp)*loct->prate_sfc;
	
	/* SNOW BUDGET ***********************************************/
	/** snow thaw water  **/
	if(mass->snwa>0.1){
		thaw = (1.0/11.0)/(1.0+exp(-0.5*(loct->tmp_2m-4.0)));	
		/* thaw=1.0/(1.0+exp(-0.375*grid->tmp_sfc));  */
		aa = (1.0 + 10.0/(0.05*mass->snwa+1.0));
		loct->thaw = thaw*aa*mass->snwa;
	}else{
		loct->thaw = mass->snwa;
	}
	
	mass->snwa += snow - loct->thaw;
	if(NOTICE==1 && VERB_HYD ==1 && mass->snwa<0.0){
		printf("!! ERROR: negative snwa (hydrology.c): %lf = %lf - %lf\n", 
			mass->snwa, snow, loct->thaw);
	}
	
	/* interception *************************************************/
	/**  tree canopy interception  **/
	loct->pm_incep_tree = pm_incep(grid, loct, 1);  /*  potential  */
	cap_incep = (mass->tree).lai*0.25;
	capt_water = (rain<cap_incep)?rain:cap_incep;
	aa = 0.85;
	bb = capt_water+loct->pm_incep_tree;
	cc = capt_water*loct->pm_incep_tree;
	loct->incep_tree = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);  /*  actual */
	
	/**  C3 herbaceous canopy interception  **/
	loct->pm_incep_c3 = pm_incep(grid, loct, 2);  /*  potential  */
	cap_incep = loct->funder_c3*(mass->c3).lai*0.25;
	aa = loct->funder_c3*(rain - loct->incep_tree);
	capt_water = (aa<cap_incep)?aa:cap_incep;
	aa = 0.85;
	bb = capt_water+loct->pm_incep_c3;
	cc = capt_water*loct->pm_incep_c3;
	loct->incep_c3 = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);  /*  actual */

	/**  c3 herbaceous canopy interception  **/
	loct->pm_incep_c4 = pm_incep(grid, loct, 3);  /*  potential  */
	cap_incep = loct->funder_c4*(mass->c4).lai*0.25;
	aa = loct->funder_c4*(rain - loct->incep_tree);
	capt_water = (aa<cap_incep)?aa:cap_incep;
	aa = 0.85;
	bb = capt_water+loct->pm_incep_c4;
	cc = capt_water*loct->pm_incep_c4;
	loct->incep_c4 = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);  /*  actual */
	
	loct->incep = loct->incep_tree + loct->incep_c3 + loct->incep_c4;
	
	/* liquid water input to soil */
	lqd_in_soil = (rain - loct->incep) + loct->thaw;
	
	/**  runoff 1, estimated by the bucket model  **/
	dry_inx = grid->fieldcap30 - mass->sw30;
	gain = lqd_in_soil;
	loct->ro1 = pow(gain, 3.0) + pow(dry_inx, 3.0);
	loct->ro1 = (loct->ro1>=0.0)?loct->ro1:0.0;
	loct->ro1 = pow(cc, 0.33333) - dry_inx;
	loct->ro1 = (loct->ro1>=0.0)?loct->ro1:0.0;

	/************  water balance for fast flows   ****************/
	mass->sw30 += lqd_in_soil - loct->ro1;
	if(NOTICE==1 && VERB_HYD ==1 && mass->sw30<0.0){
		printf("!! ERROR: negative sw30 (hydrology.c): %lf\n", mass->sw30);
	}

	/**  soil evaporation  **/
	loct->pm_evpr = pm_evap(grid, loct);  /*  potential  */
	aa = 0.85;
	bb = (mass->sw30)+loct->pm_evpr;
	cc = (mass->sw30)*loct->pm_evpr;
	loct->evpr = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);  /*  actual */
	loct->evpr = (loct->evpr>0.0)?loct->evpr:0.0;
	
	/**  extraction of evaporation **/
	mass->sw30 -= loct->evpr;
	if(NOTICE==1 && VERB_HYD ==1 && mass->sw30<0.0){
		printf("!! ERROR: negative sw30 (hydrology.c): %lf\n", mass->sw30);
	}
	
	/**  transpiration from C3 herbaceous using upper soil water  **/
	loct->gc_c3 = (echar->c3).gc * loct->funder_c3;
	loct->pm_trnsp_c3 = pm_transp(grid, loct, 2);  /*  potential  */
	aa = 0.85;
	bb = (mass->sw30)+loct->pm_trnsp_c3;
	cc = (mass->sw30)*loct->pm_trnsp_c3;
	loct->trnsp_c3 = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);  /*  actual */
	if(NOTICE==1 && VERB_HYD ==1 && loct->trnsp_c3<0.0){
		printf("!! ERROR: negative trnsp_f (hydrology.c): %lf\n", loct->trnsp_c3);
	}

	/**  transpiration from c3 herbaceous using upper soil water  **/
	loct->gc_c4 = (echar->c4).gc * loct->funder_c4;
	loct->pm_trnsp_c4 = pm_transp(grid, loct, 3);  /*  potential  */
	aa = 0.85;
	bb = (mass->sw30)+loct->pm_trnsp_c4;
	cc = (mass->sw30)*loct->pm_trnsp_c4;
	loct->trnsp_c4 = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);  /*  actual */
	if(NOTICE==1 && VERB_HYD ==1 && loct->trnsp_c4<0.0){
		printf("!! ERROR: negative trnsp_f (hydrology.c): %lf\n", loct->trnsp_c3);
	}

	/************  water balance for slow flows   ****************/
	mass->sw30 -= (loct->trnsp_c3 + loct->trnsp_c4);	
	if(NOTICE==1 && VERB_HYD ==1 && mass->sw30<0.0){
		printf("!! ERROR: negative sw30 (hydrology.c): %lf\n", mass->sw30);
	}

	/**  transpiration from tree using lower soil water  **/
	loct->gc_tree = (echar->tree).gc;
	loct->pm_trnsp_tree = pm_transp(grid, loct, 1);  /*  potential  */
	aa = 0.85;
	bb = (mass->sww)+loct->pm_trnsp_tree;
	cc = (mass->sww)*loct->pm_trnsp_tree;
	loct->trnsp_tree = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);  /*  actual */
	if(NOTICE==1 && VERB_HYD ==1 && loct->trnsp_tree<0.0){
		printf("!! ERROR: negative trnsp_t (hydrology.c): %lf\n", loct->trnsp_tree);
	}
	
	loct->trnsp = loct->trnsp_c3 + loct->trnsp_c4 + loct->trnsp_tree;
	
	/* RUNOFF ************************************************/
	/**  runoff, baseflow  **/
	if(loct->tmp200_soil > 0.0){
		baseflow = 0.001 * mass->sww;   /* after Sellers et al. (1996) */
		if((strcmp(grid->site_id, "QHB")==0)){
			baseflow = 0.003 * mass->sww;
		}
	}else{
		baseflow = 0.0;
	}
	mass->sww -= baseflow;

	/**  runoff 2, estimated by the bucket model  **/
	dry_inx = grid->fieldcap - mass->sww;
	gain = loct->ro1;
	cc = pow(gain, 3.0) + pow(dry_inx, 3.0);
	cc = (cc>=0.0)?cc:0.0;
	loct->ro2 = pow(cc,0.33333) - dry_inx;
	loct->ro2 = (loct->ro2>=0.0)?loct->ro2:0.0;
	
	loct->ro2 += baseflow;

	/** re-translocation ***********************************/
	/**  penetration capacity from SW30 to SWW  **/
	aa = grid->hyd_cond*1000.0*3600.0*24.0;

	retran = (mass->sww*grid->fieldcap30/grid->fieldcap - mass->sw30)/(1.0+grid->fieldcap30/grid->fieldcap);
	if(retran>0.0){
		retran *= 0.5;
		if(retran>aa){
			retran = aa;
		}
	}else{
		if(retran<(-aa)){
			retran = -aa;
		}
	}
	mass->sw30 += retran;
	mass->sww -= retran;  
	
	if(mass->sw30<0.0){
		mass->sw30 = 0.0;
	}
	if(mass->sww<0.0){
		mass->sww = 0.0;
	}

	/************  water balance for slow flows   ****************/
	mass->sww += loct->ro1 - loct->ro2 - loct->trnsp_tree;  
	if(NOTICE==1 && VERB_HYD ==1 && mass->sww<0.0){
		printf("!! ERROR: negative sww (hydrology.c): %lf\n", mass->sww);
	}
	
	loct->aet = loct->incep + loct->trnsp + loct->evpr;
	loct->pet = loct->pm_evpr + loct->pm_incep_tree + loct->pm_incep_c3 + loct->pm_incep_c4 + 
				loct->pm_trnsp_tree + loct->pm_trnsp_c3 + loct->pm_trnsp_c4;
}

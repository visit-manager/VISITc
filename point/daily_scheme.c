/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/*************************************************************************************/
void daily_scheme(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	/* nitrogen input ************************/
	if(N_CYCLE ==1){
		/* biological N2 fixation */
		f_biolfix(grid, loct, mass, flux);
	
		/* plant N uptake */
		f_n_uptake(grid, loct, echar, mass, flux);
	}
	
	/***** vegetation processes *****/
	f_growth_period(grid, loct, &(echar->tree));
	f_growth_period(grid, loct, &(echar->c3));
	f_growth_period(grid, loct, &(echar->c4));
    
    if(FIX_PHENOLOGY == 1 && loct->adyear>=1990){
        (echar->tree).season = fixlp_tree_season[loct->doy];
        (echar->c3).season = fixlp_c3_season[loct->doy];
        (echar->c4).season = fixlp_c4_season[loct->doy];
        
        (echar->tree).day_flush = fixlp_tree_dayflush[loct->doy];
        (echar->c3).day_flush = fixlp_c3_dayflush[loct->doy];
        (echar->c4).day_flush = fixlp_c4_dayflush[loct->doy];
        
        (echar->tree).day_shed = fixlp_tree_dayshed[loct->doy];
        (echar->c3).day_shed = fixlp_c3_dayshed[loct->doy];
        (echar->c4).day_shed = fixlp_c4_dayshed[loct->doy];
    }

	/* plant process *********/
	plant_process(grid, loct, &(flux->tree), &(echar->tree), &(mass->tree));
	plant_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
	plant_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
	
	/* VOC emission */
	f_voc_emit_guenther97(grid, loct, echar, mass, flux);
	/* plant CH4 emission */
	f_ch4emit_plant(grid, loct, echar, mass, flux);

	/**/
	(flux->soil).li_tf = (flux->tree).lf;			
	(flux->soil).li_tc = (flux->tree).lc;			
	(flux->soil).li_tr = (flux->tree).lr;			
	(flux->soil).li_gf = loct->funder_c3 * (flux->c3).lf + loct->funder_c4 * (flux->c4).lf;	
	(flux->soil).li_gc = loct->funder_c3 * (flux->c3).lc + loct->funder_c4 * (flux->c4).lc;	
	(flux->soil).li_gr = loct->funder_c3 * (flux->c3).lr + loct->funder_c4 * (flux->c4).lr;	
	/**/
	if(SCI_SCHEME == 1){
		(flux->soil).d13c_li_f = d13c_addition((flux->tree).d13c_lf,(flux->tree).lf, 
				(flux->c3).d13c_lf,(flux->c3).lf);
		(flux->soil).d13c_li_c = d13c_addition((flux->tree).d13c_lc,(flux->tree).lc, 
				(flux->c3).d13c_lc,(flux->c3).lc);
		(flux->soil).d13c_li_r = d13c_addition((flux->tree).d13c_lr,(flux->tree).lr, 
				(flux->c3).d13c_lr,(flux->c3).lr);
	}
	
	/* soil processes *********/
	f_cycle_soil(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));
	
	/* dissolved organic carbon */
	f_doc_boyer(grid, loct, mass, flux);
	if(loct->time > 100){
		(mass->soil).msl_a -= (flux->soil).doc_boyer / 1000000.0; 
	}

	/* ecosystem mass balance *********************/	
	flux->gpp = (flux->tree).gpp + loct->funder_c3 * (flux->c3).gpp + loct->funder_c4 * (flux->c4).gpp;
	flux->npp = (flux->tree).npp + loct->funder_c3 * (flux->c3).npp + loct->funder_c4 * (flux->c4).npp;
	/* net ecosystem production */
	flux->nep = ((flux->tree).npp + loct->funder_c3 * (flux->c3).npp + 
				loct->funder_c4 * (flux->c4).npp) - (flux->soil).hr;
	/* total ecosystem carbon storage */
	mass->total_c = (mass->tree).plant + loct->funder_c3 * (mass->c3).plant + 
				loct->funder_c4 * (mass->c4).plant + (mass->soil).soil;
	/** net carbon balance taking crop harvest into account **/
	flux->necb = flux->nep + ((flux->tree).hvst + loct->funder_c3 * (flux->c3).hvst + 
				loct->funder_c4 * (flux->c4).hvst);
	
	flux->rr = (flux->tree).rrm + (flux->tree).rrg + 
				loct->funder_c3 * ((flux->c3).rrm + (flux->c3).rrg) + 
				loct->funder_c4 * ((flux->c4).rrm + (flux->c4).rrg);
	flux->sr = flux->rr + (flux->soil).hr;
	flux->er = (flux->tree).rp + loct->funder_c3 * (flux->c3).rp + 
				loct->funder_c4 * (flux->c4).rp + (flux->soil).hr;
	
	/* stable carbon isotope */
	if(SCI_SCHEME==1){
		d13c_efflux(grid, loct, flux);
	}
	
	loct->vmc = mass->sw30/300.0;
	loct->d_c_min = (flux->soil).mr_tf + (flux->soil).mr_tc + (flux->soil).mr_tr + 
				(flux->soil).mr_gf + (flux->soil).mr_gc + (flux->soil).mr_gr + 
				(flux->soil).mr_a + (flux->soil).mr_i + (flux->soil).mr_p;

	/* CASA moisture ******************************/
	f_casa_moisture(grid, loct);
	
	/* GHG budget ****************************************/
	if(GHG_CALC==1 || N_CYCLE ==1){
		/* N2O by NGAS */
		f_n2o_emit_ngas(grid, loct, mass, flux);
		f_n2o_emit_ngas_2(grid, loct, mass, flux);
		
		/* CH4 oxy by Del Grosso scheme */
		f_ch4oxy_delgrosso(grid, loct, flux);
		
		/* CH4 oxy by Ridgewell scheme */
		f_ch4oxy_ridgewell(grid, loct, mass, flux);
		/* CH4 oxy by Curry scheme */
		f_ch4oxy_curry(grid, loct, mass, flux);
		/* CH4 oxy by CASA scheme */
		f_ch4oxy_casa(grid, loct, flux);
		
		/* CH4 emission by Walter & Heimann (2000) */
		if(WH_CH4==1){	/* wetland */
			if(strcmp(grid->site_id, "MSE")==0){
                
                /* revised: 2013/03/17 by A.Ito */
                if(loct->adyear>=2001 && loct->adyear<=2001){
                    loct->water_table_depth = paddy_wtd[loct->adyear-2001][loct->doy];
                }else{
                    loct->water_table_depth = paddy_wtd[0][loct->doy];
                }
			}else{
				loct->water_table_depth = 0.2;
			}
            /* scheme => ch4_emit.c */
			f_ch4emit_walter(grid, loct, flux);
            
		}else{
			loct->water_table_depth = -2.0;
			(flux->soil).ch4_wh_plant = 0.0;
			(flux->soil).ch4_wh_ebull = 0.0;
			(flux->soil).ch4_wh_diff = 0.0;
		}
		
		/* N2O by CASA */
		f_n2o_emit_casa(grid, loct, mass, flux);
	}
	
	/* nitrogen budget */
	if(N_CYCLE ==1){
		n_budget(grid, loct, mass, flux);
	}
	
	/* store data as previous condition in the next step */
	loct->m_casa_pre = loct->m_casa;
	loct->vmc_pre = loct->vmc;
    loct->water_table_depth_pre = loct->water_table_depth;
}

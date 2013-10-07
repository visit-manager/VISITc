/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* N budget ******************************************/
void n_budget(
	struct Grid *grid, 
	struct Loct *loct,
	struct Mass *mass, 
	struct Flux *flux
){
	double abdn, uptake_no3, uptake_nh4;

	abdn = (flux->tree).n_abdn_cnpy + (flux->tree).n_abdn_strg + 
			(flux->c3).n_abdn_cnpy + (flux->c3).n_abdn_strg + 
			(flux->c4).n_abdn_cnpy + (flux->c4).n_abdn_strg;
	uptake_no3 = (flux->tree).uptake_no3 + (flux->c3).uptake_no3 + (flux->c4).uptake_no3;
	uptake_nh4 = (flux->tree).uptake_nh4 + (flux->c3).uptake_nh4 + (flux->c4).uptake_nh4;

	/* plant N ******************************************/
	/* tree canopy */
	(mass->tree).n_canopy += (flux->tree).n_alloc_cnpy 
						+ (flux->tree).n_realloc 
						- (flux->tree).n_salvage 
						- (flux->tree).n_abdn_cnpy;
	if((mass->tree).n_canopy < 0.0){
		(mass->tree).n_canopy = 0.0;
	}

	/* C3 harb canopy */
	(mass->c3).n_canopy += (flux->c3).n_alloc_cnpy 
						+ (flux->c3).n_realloc 
						- (flux->c3).n_salvage 
						- (flux->c3).n_abdn_cnpy;
	if((mass->c3).n_canopy < 0.0){
		(mass->c3).n_canopy = 0.0;
	}

	/* C4 harb canopy */
	(mass->c4).n_canopy += (flux->c4).n_alloc_cnpy 
						+ (flux->c4).n_realloc 
						- (flux->c4).n_salvage 
						- (flux->c4).n_abdn_cnpy;
	if((mass->c4).n_canopy < 0.0){
		(mass->c4).n_canopy = 0.0;
	}
	
	/* storage */
	(mass->tree).n_storage += (flux->tree).n_alloc_strg 
						- (flux->tree).n_realloc 
						+ (flux->tree).n_salvage 
						- (flux->tree).n_abdn_strg;
	if((mass->tree).n_storage < 0.0){
		(mass->tree).n_storage = 0.0;
	}
	
	(mass->c3).n_storage += (flux->c3).n_alloc_strg 
						- (flux->c3).n_realloc 
						+ (flux->c3).n_salvage 
						- (flux->c3).n_abdn_strg;
	if((mass->c3).n_storage < 0.0){
		(mass->c3).n_storage = 0.0;
	}
	
	(mass->c4).n_storage += (flux->c4).n_alloc_strg 
						- (flux->c4).n_realloc 
						+ (flux->c4).n_salvage 
						- (flux->c4).n_abdn_strg;
	if((mass->c4).n_storage < 0.0){
		(mass->c4).n_storage = 0.0;
	}
	
	/* soil N *******************************************/
	/* microbe */
	(mass->soil).n_mcrb += (flux->soil).n_immbl
						- (flux->soil).n_mcrb_abdn;
	if((mass->soil).n_mcrb < 0.0){
		(mass->soil).n_mcrb = 0.0;
	}
		
	/* litter */
	(mass->soil).n_lttr += abdn - (flux->soil).n_minerlz_lttr;
	if((mass->soil).n_lttr < 0.0){
		(mass->soil).n_lttr = 0.0;
	}

	/* humus */
	(mass->soil).n_hums += (flux->soil).n_mcrb_abdn
						- (flux->soil).n_minerlz_hums;
	if((mass->soil).n_hums < 0.0){
		(mass->soil).n_hums = 0.0;
	}

	/* nitrate */
	(mass->soil).n_no3 += loct->depo_no3 + loct->nfert_no3 
						+ (flux->soil).n_nitrif
						- uptake_no3 
						- (flux->soil).d_n2o_dnt_ngas*28.0/44.0
						- (flux->soil).d_n2_ngas
						- (flux->soil).n_leach;
	if((mass->soil).n_no3 < 0.0){
		(mass->soil).n_no3 = 0.0;
	}

	/* ammonium */
	(mass->soil).n_nh4 += loct->depo_nh4 + loct->nfert_nh4 
						+ ((flux->soil).n_minerlz_lttr +(flux->soil).n_minerlz_hums
							- (flux->soil).n_immbl)
						- uptake_nh4 
						- (flux->soil).n_nitrif
						- (flux->soil).d_n2o_ntr_ngas*28.0/44.0
						- (flux->soil).n_nh3vlt;
	if((mass->soil).n_nh4 < 0.0){
		(mass->soil).n_nh4 = 0.0;
	}
}

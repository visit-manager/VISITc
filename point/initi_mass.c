/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/* header files */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* initialize the tree biomass to a juvenile stage *************/
void initTree(
	struct Mass *mass
){
	/* C, Mg C ha-1 */
	(mass->tree).fol = INIT_C; 
	(mass->tree).stm = INIT_C; 
	(mass->tree).rot = INIT_C;
    (mass->tree).nsch_storage = 0.1*INIT_C;
	(mass->tree).plant = (mass->tree).fol + (mass->tree).stm + (mass->tree).rot;
    (mass->tree).grain = 0.0;
    (mass->tree).rooting_depth = 0.1;

	/* N, kg N ha-1 */
	(mass->tree).n_canopy = INIT_N; 
	(mass->tree).n_storage = INIT_N; 
	(mass->tree).n_plant = (mass->tree).n_canopy + (mass->tree).n_storage;
}

/* initialize the floor C3/C4 plant biomass to a juvenile stage *******/
void initFloor(
	struct Mass *mass
){
	/* C, Mg C ha-1 */
	(mass->c3).fol = INIT_C; 
	(mass->c3).stm = INIT_C; 
	(mass->c3).rot = INIT_C;
    (mass->c3).nsch_storage = 0.1*INIT_C;
	(mass->c3).plant = (mass->c3).fol + (mass->c3).stm + (mass->c3).rot;
    (mass->c3).grain = 0.0;
    (mass->c3).rooting_depth = 0.1;

	(mass->c4).fol = INIT_C; 
	(mass->c4).stm = INIT_C; 
	(mass->c4).rot = INIT_C;
    (mass->c4).nsch_storage = 0.1*INIT_C;
	(mass->c4).plant = (mass->c4).fol + (mass->c4).stm + (mass->c4).rot;
    (mass->c4).grain = 0.0;
    (mass->c4).rooting_depth = 0.1;

	/* N, kg N ha-1 */
	(mass->c3).n_canopy = INIT_N; 
	(mass->c3).n_storage = INIT_N; 
	(mass->c3).n_plant = (mass->c3).n_canopy + (mass->c3).n_storage;

	(mass->c4).n_canopy = INIT_N; 
	(mass->c4).n_storage = INIT_N; 
	(mass->c4).n_plant = (mass->c4).n_canopy + (mass->c4).n_storage;
}

/* initialize the soil organic carbon **********/
void initSoil(
	struct Mass *mass
){
	/* C, Mg C ha-1 */
	(mass->soil).ltr_tf = INIT_C;
	(mass->soil).ltr_tc = INIT_C;
	(mass->soil).ltr_tr = INIT_C; 
	(mass->soil).ltr_gf = INIT_C;
	(mass->soil).ltr_gc = INIT_C;
	(mass->soil).ltr_gr = INIT_C; 
	(mass->soil).msl_a = INIT_C;
	(mass->soil).msl_i = INIT_C;
	(mass->soil).msl_p = INIT_C;
	
	(mass->soil).mcrb = INIT_C;
	(mass->soil).doc = 0.0;
	
	(mass->soil).soil = (mass->soil).ltr_tf + (mass->soil).ltr_tc + (mass->soil).ltr_tr + 
		(mass->soil).ltr_gf + (mass->soil).ltr_gc + (mass->soil).ltr_gr + 
		(mass->soil).msl_a + (mass->soil).msl_i + (mass->soil).msl_p;
	
	/* N, kg N ha-1 */
	(mass->soil).n_lttr = INIT_N;
	(mass->soil).n_hums = INIT_N;

	(mass->soil).n_nh4 = INIT_N;
	(mass->soil).n_no3 = INIT_N;
	
	(mass->soil).n_mcrb = INIT_N;

	(mass->soil).n_soil = (mass->soil).n_lttr + (mass->soil).n_hums + (mass->soil).n_mcrb + 
		(mass->soil).n_nh4 + (mass->soil).n_no3;
}

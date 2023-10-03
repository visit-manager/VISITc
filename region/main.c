/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & JAMSTEC							*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version VISITc for 1-hourly run                                     */

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito
	modify: 2013 / 01 / 04 by Akihiko Ito
*/

/* header files */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

#ifdef _OPENMP
#include<omp.h>
#endif

/* global variables */
#include"global_var.h"

/*** main simulation roop ***********************************************/
int main(void){
	FILE *fp_r[N_OFILE];
	
	/* *******************************************************************/
	/** setting configuration **/
	f_setting(&grid[0]);
    srand(1393);
	
	/* *******************************************************************/
	/** initialization & open files **/
	f_initialize(grid, &loct, &echar, echar_type, mass, &flux, fp_r);
	
	/* *******************************************************************/
	/** spin-up **/
    if(USE_RESTART == 0){
        f_spinup(grid, &loct, &echar, echar_type, mass, &flux, fp_r[0]);
    }
	
	/** experiment **/
	f_experiment(grid, &loct, &echar, echar_type, mass, &flux, fp_r); /* */
	
	/* *******************************************************************/
	/** close files **/
	f_terminate(fp_r);
	
	return 0;
}

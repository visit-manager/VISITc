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
#include"definition.h"			/* definition of constants */
#include"structure.h"			/* definition of structures */
#include"prototype.h"			/* prototypes of functions */
#include"global_variables.h"	/* definition of global variables */

/* main simulation roop *************************************************/
int main(int argc, char *argv[]){
	
	/* output files */
	FILE *fp_r[NFILE];
    /* experiment number using LARS-generated weathewr data */
    short rnum;
	
	/* rnum = (short)atol(argv[1]); */ /* used in ensemble simulations */
	rnum = 0; /* */
	
	/********************************************************************/	
	/* optional: parameter perturbation for sensitivity analysis */
	f_parameter_perturbation(rnum);
	
	/********************************************************************/	
	/** setting configuration: =>[setting.c] **/
	f_setting(&grid[0][0]);
	
	/** open output files: =>[io_file.c] **/
	f_open_files(&grid[0][0], rnum, fp_r);

	/********************************************************************/	
	/** initialization (open output files): =>[initialize.c] **/
	f_initialize(grid, &loct, &echar, echar_type, mass, &flux, fp_r);
	
	/* printf("DEBUG:  %lf\n", (echar.soil).sr_hp); */
	
	/********************************************************************/	
	/** spin-up: =>[spinup.c] **/
	f_spinup(grid, &loct, &echar, echar_type, mass, &flux, fp_r[0]);
	
	/** experiment: =>[experiment.c] **/
	f_experiment(grid, &loct, &echar, echar_type, mass, &flux, fp_r);
	
	/********************************************************************/	
	/** close output files: =>[io_file.c] **/
	f_close_files(fp_r);
	
	return 0;
}

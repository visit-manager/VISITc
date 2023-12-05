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
#include<time.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

#ifdef _OPENMP
#include<omp.h>
#endif

/* global variables */
#include"global_var.h"

/* main simulation roop *************************************************/
int main(void){
    struct Grid *grid2;
    struct Mass *mass2;
    struct Loct *loct2;
    struct Echar *echar2;
    struct Flux *flux2;
    
    grid2 = (struct Grid *)malloc(sizeof(struct Grid) * NROW * NCOL);
    mass2 = (struct Mass *)malloc(sizeof(struct Mass) * NROW * NCOL);
    loct2 = (struct Loct *)malloc(sizeof(struct Loct) * NROW * NCOL);
    echar2 = (struct Echar *)malloc(sizeof(struct Echar) * NROW * NCOL);
    flux2 = (struct Flux *)malloc(sizeof(struct Flux) * NROW * NCOL);

	/* *******************************************************************/
	/** setting configuration **/
	f_setting(&grid2[0]);
    srand( time(NULL)%(100000) );
	
	/* *******************************************************************/
	/** initialization & open files **/
	f_initialize(grid2, loct2, echar2, mass2, flux2);
	
	/* *******************************************************************/
	/** spin-up **/
    if(USE_RESTART == 0 || USE_RESTART == 2){
        f_spinup(grid2, loct2, echar2, mass2, flux2);
    }
	
	/** experiment **/
	f_experiment(grid2, loct2, echar2, mass2, flux2); /* */
	
	/* *******************************************************************/
    
    free(grid2);
    free(mass2);
    free(loct2);
    free(echar2);
    free(flux2);

	return 0;
}

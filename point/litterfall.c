/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/* functions of litterfall */
#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* from foliage **********************************/
double f_lf(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double lf;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lf = pchar->lf*mass->fol;
    
    if(pchar->crop_stage==5){
        lf = 0.7*mass->fol;
    }
		
	return (lf);
}

/* from stem and branch ***************************/		
double f_lc(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double lc;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lc = pchar->lc*mass->stm;
    
    if(pchar->crop_stage==5){
        lc = 0.7*mass->stm;
    }
	
	return (lc);
}

/* from root system ******************************/
double f_lr(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double lr;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lr = pchar->lr*mass->rot;
    
    if(pchar->crop_stage==5){
        lr = 0.7*mass->rot;
    }
	
	return (lr);
}

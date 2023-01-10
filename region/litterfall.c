/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
*/

/* functions of litterfall */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* from foliage **********************************/
float f_lf(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float lf;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lf = pchar->lf * mass->fol / (double)DSTEP;
	
	return (lf);
}

/* from stem and branch ***********************/		
float f_lc(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float lc;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lc = pchar->lc * mass->stm / (double)DSTEP;
	
	return (lc);
}

/* from root system ******************************/
float f_lr(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float lr;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lr = pchar->lr * mass->rot / (double)DSTEP;
	
	return (lr);
}

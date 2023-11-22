/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2007 / 01 / 26 by Akihiko Ito: Version 0.2.0.
*/

/* functions of decomposition */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* soil respiration coefficient of litter layer ****************************/
float frl(
	struct Grid *grid,
	struct Loct *loct,
	struct Schar *schar
){
	float ftl, fwl, fal, fsl;
	
	/* temperature effect, exponential */
	/* ftl=exp(log(soil->qTl)/10.0*(grid->tmp10_soil-to)); */
	/* after Lloyd & Taylor (1994), Func. Ecol. */
	if(loct->tmp10_soil > -20.0){
		/* ftl=0.1+0.9*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil+46.02))); */
		ftl = 0.01 + exp(308.56*(1.0/56.02 - 1.0/(loct->tmp10_soil + 46.02)));
	}else{
		ftl = 0.01;
	} 
	
	if(NOTICE == 1 && ftl<0.0){
		printf("FARTAL ERROR: negative ftl in frl (decomposition.c): %f\n", ftl);
	}

	/* soil moisture effect, saturating */
	fwl = 0.8 * loct->soilwtr_l/(schar->kml*grid->fieldcap30 + loct->soilwtr_l)+0.2;
	/* soil apparence effect */
	fal = 0.4 * loct->soilappr_l*(1.0*schar->kmsl)/(schar->kmsl + loct->soilappr_l)+0.6;
	
    fsl = (fwl>fal)?fal:fwl;
	
	if(NOTICE == 1 && fwl<0.0){
		printf("FARTAL ERROR: negative fsl in frl (decomposition.c): %f %f %f\n", fsl, fwl, fal);
	}

	return (ftl*fsl);
}

/* soil respiration coefficient of humus layer *********************************/
float frh(
	struct Grid *grid,
	struct Loct *loct, 
	struct Schar *schar
){
	float fth, fwh, fah, fsh;
	
	/* temperature effect, exponential */
	/* fth=exp(log(soil->qTh)/10.0*(grid->tmp200_soil-to)); */
	/* after Lloyd & Taylor (1994), Func. Ecol. */
	if(loct->tmp200_soil > -20.0){
		/* fth=0.1+0.9*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil+46.02))); */
		fth = 0.01 + exp(308.56*(1.0/56.02 - 1.0/(loct->tmp200_soil + 46.02)));
	}else{
		fth = 0.01;
	} 

	if(NOTICE==1 && fth<0.0){
		printf("FARTAL ERROR: negative fth in frh (decomposition.c): %f\n", fth);
	}

	/* soil moisture effect, saturating */
	fwh = 0.8*loct->soilwtr_h/(schar->kmh*grid->fieldcap + loct->soilwtr_h)+0.2;
	/* soil apparence effect */
	fah = 0.4*loct->soilappr_w*(1.0*schar->kmsh)/(schar->kmsh + loct->soilappr_w)+0.6;
	
    fsh = (fwh>fah)?fah:fwh;
		
	if(NOTICE==1 && fsh<0.0){
		printf("FARTAL ERROR: negative fsh in frh (decomposition.c): %f %f %f\n", fsh, fwh, fah);
	}

	return (fth*fsh);
}

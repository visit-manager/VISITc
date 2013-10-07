/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/* functions of decomposition */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* soil respiration coefficient of litter layer ****************************/
double frl(
	struct Grid *grid,
	struct Loct *loct,
	struct Schar *schar
){
	double ftl, fwl, fal, fsm;
	double satmp_var, sasw_var;
	
	if(SA_PARA==6){
		satmp_var = SA_PARA_VAR;
	}else{
		satmp_var = 1.0;
	}

	if(SA_PARA==7){
		sasw_var = SA_PARA_VAR;
	}else{
		sasw_var = 1.0;
	}
	
	/* temperature effect, exponential */
	/* ftl=exp(log(soil->qTl)/10.0*(grid->tmp10_soil-to)); */
	/* after Lloyd & Taylor (1994), Func. Ecol. */
	if(loct->tmp10_soil > -20.0){
		/* ftl=0.1+0.9*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil+46.02))); */
		ftl = 0.01 + exp(308.56*(1.0/56.02-1.0/(loct->tmp10_soil+46.02*satmp_var)));
	}else{
		ftl = 0.01;
	} 
	
	if(NOTICE==1 && ftl<0.0){
		printf("FARTAL ERROR: negative ftl (soil_proc.c): %lf\n", ftl);
	}

	/* soil moisture effect, saturating */
	fwl = 0.8*loct->soilwtr_l/(schar->kml*sasw_var*grid->fieldcap30 + loct->soilwtr_l)+0.2;
	/* soil apparence effect */
	fal = 0.4*loct->soilappr_l*(1.0*schar->kmsl)/(schar->kmsl + loct->soilappr_l)+0.6;
	
	fsm = (fwl>fal)?fal:fwl;
	
	if(NOTICE==1 && fwl<0.0){
		printf("FARTAL ERROR: negative fsm in l (soil_proc.c): %lf\n", grid->fieldcap30);
	}

	return (ftl*fsm);
}

/* soil respiration coefficient of humus layer *********************************/
double frh(
	struct Grid *grid,
	struct Loct *loct, 
	struct Schar *schar
){
	double fth, fwh, fah, fsm;	
	double satmp_var, sasw_var;
	
	if(SA_PARA==6){
		satmp_var = SA_PARA_VAR;
	}else{
		satmp_var = 1.0;
	}

	if(SA_PARA==7){
		sasw_var = SA_PARA_VAR;
	}else{
		sasw_var = 1.0;
	}
	
	/* temperature effect, exponential */
	/* fth=exp(log(soil->qTh)/10.0*(grid->tmp200_soil-to)); */
	/* after Lloyd & Taylor (1994), Func. Ecol. */
	if(loct->tmp200_soil > -20.0){
		/* fth=0.1+0.9*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil+46.02))); */
		fth = 0.01 + exp(308.56*(1.0/56.02-1.0/(loct->tmp200_soil+46.02*satmp_var)));
	}else{
		fth = 0.01;
	} 

	if(NOTICE==1 && fth<0.0){
		printf("FARTAL ERROR: negative fth (soil_proc.c): %lf\n", fth);
	}

	/* soil moisture effect, saturating */
	fwh = 0.8*loct->soilwtr_h/(schar->kmh*sasw_var*grid->fieldcap + loct->soilwtr_h)+0.2;
	/* soil apparence effect */
	fah = 0.4*loct->soilappr_w*(1.0*schar->kmsh)/(schar->kmsh + loct->soilappr_w)+0.6;
	
	fsm = (fwh>fah)?fah:fwh;
		
	if(NOTICE==1 && fsm<0.0){
		printf("FARTAL ERROR: negative fsm in h (soil_proc.c): %lf\n", grid->fieldcap);
	}

	return (fth*fsm);
}

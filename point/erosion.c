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

/* soil erosion estmated with RUSLE scheme **************************************/
/*
 Ito, A. 2007. Simulated impacts of climate and land-cover change on soil erosion and 
 implication for the carbon cycle, 1901 to 2100. Geophysical Research Letters 
 34:10.1029/2007GL029342.
*/
void f_erosion_rusle(
	struct Grid *grid,  
	struct Loct *loct,  
	struct Mass *mass, 
	struct Flux *flux
){
	double ff, mm, ss, dg, silt, soc;

	grid->slope_angle = 10.0;		/* angle, degree */
	grid->slope_length = 100.0;		/* length, m */
	soc = 1.0;	/* % weight */
	
	/* R: rain *****************************************/
	if(loct->prate_ann<=850.0){
		loct->f_rain = 0.0483 * pow(loct->prate_ann, 1.610);
	}else{
		loct->f_rain = 587.8 - 1.219*loct->prate_ann 
				+ 0.004105 * loct->prate_ann * loct->prate_ann;
	}
	if(loct->f_rain<0.0){
		loct->f_rain = 0.0;
	}
	
	/* LS: slope *****************************************/
	if(tan(grid->slope_angle*dTr)<0.09){
		ss = 10.8 * sin(grid->slope_angle*dTr) + 0.03;
	}else{
		ss = 16.8 * sin(grid->slope_angle*dTr) - 0.50;
	}
	ff = sin(grid->slope_angle*dTr)/0.0896/(3.0*pow(sin(grid->slope_angle*dTr), 0.8) + 0.56);
	mm = ff/(1.0 + ff);
	loct->f_slope = ss*pow( (grid->slope_length/22.13), mm);
	
	/* K: erodibility *****************************************/
	if(grid->clay_frac > 0.0){
		silt = 1.0 - grid->sand_frac - grid->clay_frac;
		dg = -3.5*grid->sand_frac -2.0*silt -0.5*grid->clay_frac;
		loct->f_erodibility = 0.0293 * (0.65 - dg + 0.24*dg*dg)*
				exp(-0.0021*soc/grid->clay_frac - 0.00037*(soc/grid->clay_frac)*(soc/grid->clay_frac) 
				- 4.02*grid->clay_frac + 1.72*grid->clay_frac*grid->clay_frac);
	}else{
		loct->f_erodibility = 0.0;
	}
	
	/* C: vegetation cover *****************************************/
    if(loct->veg_state == 0){
        loct->f_cover = 0.001;	/* forest */
        //loct->f_cover = 0.01;	/* forest */
    }else if(loct->veg_state == 1){
        loct->f_cover = 0.5;	/* cropland */
    }
	
	/* P: soil conservation *****************************************/
	if(loct->veg_state == 0){
        loct->f_conservation = 1.0;	/* forest */
    }else if(loct->veg_state == 1){
        loct->f_conservation = 0.5;	/* cropland */
    }
	
	/* erosion ********************************************/
	/* t/ha/yr */
	(flux->soil).erosion_soil = loct->f_rain * loct->f_slope * loct->f_erodibility 
			* loct->f_cover * loct->f_conservation;
			
	(flux->soil).erosion_orgmat = (flux->soil).erosion_soil * soc/100.0;
	(flux->soil).erosion_carbon = (flux->soil).erosion_orgmat * 0.45;
}

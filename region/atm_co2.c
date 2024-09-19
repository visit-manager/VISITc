/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in December 20, 2007						*/

/*  [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2007 / 12 / 25 by Akihiko Ito: Version 0.2.0.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* GHG scenario ***************************************/
/* source: http://crga.atmos.uiuc.edu/research/post-sres.html
 M.E.Schlesinger and S.Malyshev			*/
/* atmospheric CO2, ppmv */
extern float		atm_co2[4][GHG_TS];	/* CO2 scenarios */
/* atmospheric CH4, pptv*/
extern float		atm_ch4[4][GHG_TS];	/* CH4 scenarios */
/* atmospheric N2O, pptv*/
extern float		atm_n2o[4][GHG_TS];	/* N2O scenarios */

/* atmospheric CO2  ****************************************************************/
void atmco2_trend(
	struct  Grid *grid, 
	struct  Loct *loct
){
	float base, inc, lgrd, season;
	float amplitude;
	float time, aa0, aa1, aa2, aa3, aa4, aa5;
	
	/**  time (AD)  **/
	time = (float)(loct->CO2y);
	if(time<1900.0 || time>2100.0){
		printf("!! Wrong CO2 year: %.1lf\n", time);
	}
	
	/** BASE: annual mean concentration, ppmv **/
	/* regression from Mauna Loa data */
	aa0 = 1904299.0;
	aa1 = -3322.4242 * pow(time, 1.0);
	aa2 = 1.6541596 * pow(time, 2.0);
	aa3 = 1.3362655 * pow(time, 3.0)/10000.0;
	aa4 = -3.0828809 * pow(time, 4.0)/10000000.0;
	aa5 = 6.2121261 * pow(time, 5.0)/100000000000.0;
	base = aa0+aa1+aa2+aa3+aa4+aa5;

	/** INC: increment rate, ppmv yr-1  **/
	inc = 0.0; /*  (-52.11+0.026984*time)*((float)(grid->doy)-5.5)/12.0;  */

	/** LGRD: north-south gradient, ppmv  **/
	lgrd = 1.6*(grid->lat / 85.0);

	/** SEASON: seasonal oscillation, ppmv  **/
	amplitude = exp(0.04 * grid->lat);
    season = 0.0;
	if(grid->lat>=0.0){
		season = amplitude/2.0*sin(((float)(loct->doy) - 0.0)/365.0*2.0*PI);
	}else if(grid->lat<0.0){
		season = amplitude/2.0*sin(((float)(loct->doy) + 180.0)/365.0*2.0*PI);
	}
	/* background (free atmosphere) CO2 concentration, ppmv  */
	//loct->bCO2 = base + lgrd + season + inc;
    
    /* loct->bCO2 = atm_co2_a1[loct->CO2y - 1750]; */
    loct->bCO2 = atm_co2[2][loct->CO2y - FDY_GHG];
	
	/* future */
	if(loct->CO2y > EYR){
		loct->bCO2 = atm_co2[1][loct->CO2y - FDY_GHG];
	}
	
	/* stable carbon isotope composition (d13C), permille  *****/
	/* regression from Mauna Loa data */
	/*  aa0 = -93986.5526;
	aa1 = pow(time, 1.0)*(255.1757597);
	aa2 = pow(time, 2.0)*(-0.2766839721);
	aa3 = pow(time, 3.0)*(1.497414535/10000.0);
	aa4 = pow(time, 4.0)*(-4.04456997/100000000.0);
	aa5 = pow(time, 5.0)*(4.361334303/1000000000000.0);
	grid->d13C_bCO2 = aa0+aa1+aa2+aa3+aa4+aa5;  */
	/* inverse relationship between concentration and d13C */
	loct->d13C_bCO2 = -6.0 + (-0.02 * (base-280.0)) + (0.05*season);
}

/* intra-canopy CO2-d13C  ********************************************************/
void f_co2_in_canopy(
	struct  Grid *grid, 
	struct  Loct *loct, 
	struct  Mass *mass, 
	struct  Flux *flux
){	
	/* ambient CO2 concentration, ppmv */
	loct->aCO2 = loct->bCO2;
	/* this will be modified with some microclimatic models */
	
	/* ambient and background d13C-CO2 */
	loct->d13C_aCO2 = loct->d13C_bCO2;
}

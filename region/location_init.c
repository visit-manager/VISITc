/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in November 02, 2007						*/

/*  [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
*/

/*   initialize environmental locations  */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

#define TERM_HYD  0.1

/* climatology: *_d[] means the time-series during 1948 to 2004- */	
extern float tmp_sfc_d[PERIOD][366];		/* ground surface temperature, degree Celcius */
extern float tmp_2m_d[PERIOD][366];		/* 2m air temperature, degree Celcius */
extern float tmp10_soil_d[PERIOD][366];	/* soil temperature at 10 cm depth, degree Celcius */
extern float tmp200_soil_d[PERIOD][366];	/* soil temperature at 200 cm depth, degree Celcius */
extern float dswrf_sfc_d[PERIOD][366];		/* downward shortwave radiation at the surface, W m-2 */
extern float tcdc_clm_d[PERIOD][366];		/* total cloudiness, fraction */
extern float prate_sfc_d[PERIOD][366];		/* precipitation, mm mon-1 */
extern float spfh_2m_d[PERIOD][366];		/* specific humidity, kg kg-1 */
extern float wind_10m_d[PERIOD][366];		/* u-wind velocity, m s-1 */
extern float vpd_d[PERIOD][366];   		/* VPD, hPa */
extern float tmp_2m_ncep_dav[366];			/* 2m air temperature, degree Celcius */
extern float prate_sfc_ncep_mav[12];		/* precipitation, mm mon-1 */
extern float tmp_sfc_dav[366];	
extern float tmp_2m_dav[366];	
extern float tmp10_soil_dav[366];	
extern float tmp200_soil_dav[366];	
extern float prate_sfc_dav[366];	
extern float dswrf_sfc_dav[366];	
extern float spfh_2m_dav[366];	
extern float tcdc_clm_dav[366];	

extern long	WMODE;
/* GHG scenario ***************************************/
/* source: http://crga.atmos.uiuc.edu/research/post-sres.html
 M.E.Schlesinger and S.Malyshev			*/
/* atmospheric CO2, ppmv */
extern float		atm_co2_a1[553];	/* SRES A1 */
extern float		atm_co2_a2[553];	/* SRES A2 */
extern float		atm_co2_b1[553];	/* SRES B1 */
extern float		atm_co2_b2[553];	/* SRES B2 */
/* atmospheric CH4, pptv*/
extern float		atm_ch4_a1[553];	/* SRES A1 */
extern float		atm_ch4_a2[553];	/* SRES A2 */
extern float		atm_ch4_b1[553];	/* SRES B1 */
extern float		atm_ch4_b2[553];	/* SRES B2 */
/* atmospheric N2O, pptv*/
extern float		atm_n2o_a1[553];	/* SRES A1 */
extern float		atm_n2o_a2[553];	/* SRES A2 */
extern float		atm_n2o_b1[553];	/* SRES B1 */
extern float		atm_n2o_b2[553];	/* SRES B2 */

/* location conditions derived from the primary data (Secondary data1) **********/
void f_init_loct(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	short h;
	float atmp2m_a[366];
	float stmp10cm_a[366];
	float prate_sfc_a[366];
    
    loct->phase = 0;
	
	loct->npp_max = 1.0;
	
	/* initial atmospheric CO2 level at Year *******************/
	loct->CO2y = BYR - 1; 
	loct->doy = 364;
	atmco2_trend(grid, loct);
	
	/* long-term average climate *******************************/
	for(h=0;h<366;h++){
		atmp2m_a[h] = stmp10cm_a[h] = prate_sfc_a[h] = 0.0;
	}
	
    /* initialize annual metrics */
    grid->tmp2m_ann = grid->stmp10cm_av = 0.0;
    grid->tmp2m_max = 0.0;
    grid->tmp2m_min = 0.0;
    grid->gp_atem = grid->prec_ann = 0.0;
    loct->tsoil_annav = 0.0;
    loct->tair_dayav = loct->tair_dayav_c = 0.0;
    grid->n_clim = grid->tmp_clim = grid->prec_clim = grid->aet_clim = grid->pet_clim = 0.0;

	/* initialization of vegetation condition ********************************/
	(mass->tree).lai = lai_mass(&(mass->tree), &(echar->tree));
	(mass->c3).lai = lai_mass(&(mass->c3), &(echar->c3));
	(mass->c4).lai = lai_mass(&(mass->c4), &(echar->c4));
    
    loct->laimax_meta = 3.0;

	/* initialization of water condition *************************************/
	mass->sww = 0.5 * grid->fieldcap;     
	mass->sw30 = 0.5 * grid->fieldcap30;
	mass->snwa = 0.0;
	
	mass->sww = (mass->sww>0.0)?mass->sww:0.0;
	mass->sw30 = (mass->sw30>0.0)?mass->sw30:0.0;
	mass->snwa = (mass->snwa>0.0)?mass->snwa:0.0;
	
	/* initialize soil CH4 concentration */
	for(h=0;h<=SOILWET_LAYER;h++){
		loct->prof_ch4[h] = atm_ch4_a1[loct->CO2y - 1750]/1000.0 
			* loct->air_prsr / (UGC*(loct->tmp10_soil + ZAT));
	}
}

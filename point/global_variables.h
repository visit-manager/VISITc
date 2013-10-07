/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include"definition.h"

/* global variables ****************************/
/* configuration */
long WMODE;				/* working mode, site or region */
long WROW;				/* working row of 1st layer */
long WCOL;				/* working column of 1st layer */
long NROW;				/* number of row */
long NCOL;				/* number of column */
long YLDIST;			/* year of the last disturbance */
long SPUPT;				/* spin-up time (years) */

long BYR;				/* begin year */
long EYR;				/* end year */
long PERIOD;			/* period */

/* structures */
struct Grid grid[DROW][DCOL];
struct Mass mass[DROW][DCOL];
struct Loct loct;
struct Echar echar;
struct Flux flux;
struct Echar echar_type[NFILE];		/* default parameters */

/* climatology: *_d[] means the time-series during 1948 to 2005 */	
long	month_day[12];
double	tmp_sfc_d[HCLIMD][YSTEP];		/* ground surface temperature, degree Celcius */
double	tmp_2m_d[HCLIMD][YSTEP];		/* 2m air temperature, degree Celcius */
double	tmp10_soil_d[HCLIMD][YSTEP];	/* soil temperature at 10 cm depth, degree Celcius */
double	tmp200_soil_d[HCLIMD][YSTEP];	/* soil temperature at 200 cm depth, degree Celcius */
double	dswrf_sfc_d[HCLIMD][YSTEP];	/* downward shortwave radiation at the surface, W m-2 */
double	tcdc_clm_d[HCLIMD][YSTEP];	/* total cloudiness, fraction */
double	prate_sfc_d[HCLIMD][YSTEP];	/* precipitation, mm mon-1 */
double	spfh_2m_d[HCLIMD][YSTEP];		/* specific humidity, kg kg-1 */
double	wind_10m_d[HCLIMD][YSTEP];	/* wind velocity, m s-1 */
double	vpd_d[HCLIMD][YSTEP]; 		/* VPD, hPa */

double	tmp_2m_ncep_dav[YSTEP];		/* 2m air temperature, degree Celcius */
double	prate_sfc_ncep_mav[12];		/* precipitation, mm mon-1 */

double	tmp_sfc_dav[YSTEP];	
double	tmp_2m_dav[YSTEP];	
double	tmp10_soil_dav[YSTEP];	
double	tmp200_soil_dav[YSTEP];	
double	prate_sfc_dav[YSTEP];	
double	dswrf_sfc_dav[YSTEP];	
double	spfh_2m_dav[YSTEP];	
double	tcdc_clm_dav[YSTEP];	

double	paddy_wtd[5][YSTEP];

double prec_altdata[HCLIMD][YSTEP];

/* GCM climate projection data */
double gcm_temp[121][YSTEP];
double gcm_prec[121][YSTEP];
double gcm_ahmd[121][YSTEP];
double gcm_dswrf[121][YSTEP];
double gcm_temp_av[YSTEP];
double gcm_prec_av[YSTEP];
double gcm_ahmd_av[YSTEP];
double gcm_dswrf_av[YSTEP];

/* LARS-generated dummy weather data */
double lars_tmin[63][YSTEP];
double lars_tmax[63][YSTEP];
double lars_prec[63][YSTEP];
double lars_srad[63][YSTEP];

/* sensitivity analysis */
short SENSANS_DIST = 0;
/* 1: no disturbance */
short SENSANS_ACO2 = 0;
/* 1: no atmospheric CO2 rise */
short SENSANS_TEMP = 0;
/* 1: no temperature change */
/* 2: use LARS-generated dammy data */
short SENSANS_PREC = 0;
/* 1: no precipitation change */
/* 2: use LARS-generated dammy data */
short SENSANS_SRAD = 0;
/* 1: no solar radiation change */
/* 2: use LARS-generated dammy data */
short N_LARS = 0;

/* parameter sensitivity analysis */
short SA_PARA = 0;
/* 1: photosyntehtic optimal temperature */
/* 2: photosyntehtic CO2 - Michaelis constant */
/* 3: photosynthetic soil water dependence */
/* 4: photosyntehtic light use efficiency */
/* 5: maintenance respiration temperature */
/* 6: decomposition temperature dependence */
/* 7: decomposition soil water dependence */
/* 8: disturbance intensity */

/* ensemble number */
short SA_PARA_EN = 0;
double SA_PARA_VAR = 0.0;


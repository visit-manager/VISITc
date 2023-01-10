/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
		create: 2009 / 06 / 15 by Akihiko Ito
*/

/*********************************/
/* declaration of structures	*/
/*********************************/

/* configuration */
long WMODE;				/* working mode, site or region */
long WGRIDS;			/* working row of 1st layer */
long YLDIST;			/* year of the last disturbance */
long SPUPT;				/* spin-up time (years) */
long NBIOME;            /* number of biomes */

/* climatology: *_d[] means the time-series during 1948 to 2005 */	
long	month_day[12];
float	tmp_sfc_d[PERIOD][366];		/* ground surface temperature, degree Celcius */
float	tmp_2m_d[PERIOD][366];		/* 2m air temperature, degree Celcius */
float	tmp10_soil_d[PERIOD][366];	/* soil temperature at 10 cm depth, degree Celcius */
float	tmp200_soil_d[PERIOD][366];	/* soil temperature at 200 cm depth, degree Celcius */
float	dswrf_sfc_d[PERIOD][366];	/* downward shortwave radiation at the surface, W m-2 */
float	tcdc_clm_d[PERIOD][366];	/* total cloudiness, fraction */
float	prate_sfc_d[PERIOD][366];	/* precipitation, mm mon-1 */
float	spfh_2m_d[PERIOD][366];		/* specific humidity, kg kg-1 */
float	wind_10m_d[PERIOD][366];	/* wind velocity, m s-1 */
float	vpd_d[PERIOD][366]; 		/* VPD, hPa */

float	tmp_2m_ncep_dav[366];		/* 2m air temperature, degree Celcius */
float	prate_sfc_ncep_mav[12];		/* precipitation, mm mon-1 */

float	tmp_sfc_dav[366];	
float	tmp_2m_dav[366];	
float	tmp10_soil_dav[366];	
float	tmp200_soil_dav[366];	
float	prate_sfc_dav[366];	
float	dswrf_sfc_dav[366];	
float	spfh_2m_dav[366];	
float	tcdc_clm_dav[366];	

float	paddy_wtd[366];

/* GHG scenario ***************************************/
/* source: http://crga.atmos.uiuc.edu/research/post-sres.html
 M.E.Schlesinger and S.Malyshev			*/
/* atmospheric CO2, ppmv */
float		atm_co2_a1[553];	/* SRES A1 */
float		atm_co2_a2[553];	/* SRES A2 */
float		atm_co2_b1[553];	/* SRES B1 */
float		atm_co2_b2[553];	/* SRES B2 */
/* atmospheric CH4, pptv*/
float		atm_ch4_a1[553];	/* SRES A1 */
float		atm_ch4_a2[553];	/* SRES A2 */
float		atm_ch4_b1[553];	/* SRES B1 */
float		atm_ch4_b2[553];	/* SRES B2 */
/* atmospheric N2O, pptv*/
float		atm_n2o_a1[553];	/* SRES A1 */
float		atm_n2o_a2[553];	/* SRES A2 */
float		atm_n2o_b1[553];	/* SRES B1 */
float		atm_n2o_b2[553];	/* SRES B2 */

/* global variables */
/* structures */
struct Grid grid[NROW * NCOL];
struct Mass mass[NROW * NCOL];
struct Loct loct;
struct Echar echar;
struct Flux flux;
struct Echar echar_type[MAX_BIOME];		/* default parameters */

struct Loct loct2[NROW*NCOL];
struct Echar echar2[NROW*NCOL];
struct Flux flux2[NROW*NCOL];

float fdat[NROW * NCOL], fdatc[NROW * NCOL];
unsigned char cdat[NROW * NCOL];
int idat[NROW * NCOL];

#if OUT_DAY==1
float outdat01[NROW*NCOL];
float outdat02[NROW*NCOL];
float outdat03[NROW*NCOL];
float outdat04[NROW*NCOL];
#endif

float out_a[N_ANNRES][NROW*NCOL];



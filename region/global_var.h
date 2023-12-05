/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
		create: 2009 / 06 / 15 by Akihiko Ito
*/

/* ********************************/
/* declaration of structures	*/
/* ********************************/

/* configuration */
long WMODE;				/* working mode, site or region */
long WGRIDS;			/* working row of 1st layer */
long YLDIST;			/* year of the last disturbance */
long SPUPT;				/* spin-up time (years) */
long NBIOME;            /* number of biomes */

/* climatology: *_d[] means the time-series during 1948 to 2005 */	
long	month_day[12];

float	tmp_2m_ncep_dav[366];		/* 2m air temperature, degree Celcius */
float	prate_sfc_ncep_mav[12];		/* precipitation, mm mon-1 */

float	tmp_sfc_dav[366];	
float	tmp_2m_dav[366];	
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
/* struct Grid grid2[NROW * NCOL];
struct Mass mass2[NROW * NCOL];
struct Loct loct2[NROW * NCOL];
struct Echar echar2[NROW * NCOL];
struct Flux flux2[NROW * NCOL]; */
struct Echar echar_type[MAX_BIOME];		/* default parameters */


float fdat[NROW * NCOL], fdatc[NROW * NCOL];
unsigned char cdat[NROW * NCOL];
int idat[NROW * NCOL];

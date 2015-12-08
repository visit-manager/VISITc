/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/**********************************************************/
/* UNIVERSAL CONSTANTS :DO NOT EDIT!!! */

#define YES 1
#define NO 0

/** pai **/
#define PI 3.141592653	

/* angle conversion, from degree to radian */
#define dTr 0.0174533	

/* angle conversion, from radian to degree */
#define rTd 57.29577951 

/* from CO2-base to Crabon-base */
#define cdTc 0.272727	

/* from dry-matter-base to Carbon base */
#define dmTc 2.2		

/* from dry-matter-base to Carbon base */
#define cTdm 0.4545	

/* from micro-mol m-2 s-1 to Mg C ha-1 day-1*/
#define lTs (3600.0*12.0/100000000.0) 

/* zero degree centigrade in absolute temperature */
#define ZAT 273.15		

/* universal gas constant */
#define UGC  8.314		

/*  latent heat of water, in W h kg-1  */
#define LHT 695.0		

/* standard stable carbon isotope ratio */
#define STCIR (0.0111/0.9889) 

/* Stefan-Boltzman constant (longwave radiation) */
#define SB_CONST  (5.6703/100000000.0)


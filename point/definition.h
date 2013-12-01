/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/* calculation step */
#define DSTEP 48		/* dayly 48 steps: i.e.,  30-minute interval */
#define YSTEP 366		/* dayly 48 steps: i.e.,  30-minute interval */
/* number of soil layers */
#define SOIL_LAYER 20

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

/**********************************************************/
/* CONDITIONAL CONSTANTS :EDIT ON DEMAND!!! */

/* number of rows */
#define DROW 1
/* number of columns */
#define DCOL 1
/* number of outputfiles */
#define NFILE 9

/* number of climate data length (years) */
#define HCLIMD 65 /* NCEP/NCAR 1948-2012 */
//#define HCLIMD 110 /* CEAMIP 1901-2010 */

/* control debugging messages */
#define NOTICE  1
/* 1: with debug message */

/** initial mass, Mg C ha-1 **/
#define INIT_C  0.01
#define INIT_N  0.001

/* nitrogen cycle scheme */
#define N_CYCLE 1
/* 0: skip, i.e. not calculate */
/* 1: calculate */

/* 30-min canopy scheme */
#define FLUX_SCHEME  0 /***************************************/
/* 0: skip, i.e. not calculate */
/* 1: calculate */

/* stable carbon isotope */
#define SCI_SCHEME  0

/* restart */
#define USE_RESTART  0
/* 0: none */

/* choice of algorithm to solve leaf gas exchange */
#define SOLV_LEAF  0
/* 0: Newton */
/* 1: analytical */

/* CH4 & N2O schemes */
#define GHG_CALC 1
#define GHG_OUT 0
/* 0: none */

/* CH4 emission with a multi-layer model */
#define WH_CH4 0
/* 0: off */
/* 1: on (activate Walter & Heimann CH4 scheme) */

/* control future projection */
#define PREDICT 0
/* 0: no prediction */
/* 1: MIROC-high A1B */
/* 10: MPI-ECHAM A1B */
/* 20: CCCma-T47 A1B */
/* 21: CCCma-T47 A2 */
/* 22: CCCma-T47 B1 */
/* 100: MIROC-high A1B -2100 */
/* 101: MPI-ECHAM A1B -2050 */
/* 102: CCCma-T47 A1B -2050 */
/* 103: CCCma-T47 A2 -2050 */
/* 104: CCCma-T47 B1 -2050 */

/*  soil degradation after land-use conversion */
#define DGSOIL 0
/* 0: off */
/* 1: on */

/*  land-cover conversion */
#define LCCONV 0
/* 0: off */
/* 1: on */

/* flags for seasonal change in leaf properties */
#define AG_VMX 1		/* Vcmax */
#define AG_JVR 1		/* Jmax/Vcmax */
#define AG_LMA 1		/* LMA */
#define AG_RD 1			/* Rd */

#define PARA_CHOICE 0
/* 0: site-specific parameters */
/* 1: default (globally-used) parameters */

#define F_DOWNREG 0
/* 0: no downregulation */
/* 1: down regulation by N dilution */

#define FIX_CLIM 0
/* 0: using normal climate data */
/* 1: using fixed climate data */

#define FIX_GSCO2 0
/* 0: off */
/* 1: on (fix CO2 level in stomatal schemes) */

/* mesophyll conductance */
#define G_MESOPHYL 0
/* 0: off */
/* 1: on */

/* perturbation for multiple maparemets */
#define PERTURB_MPARA 0
/* 0: off */
/* 1: on */

#define EX_N_ECOPHYS 0
/* 0: control */
/* 1: kn=>0.5 */

/* fix LAI phenology */
#define FIX_PHENOLOGY 0
/* 0: off */
/* 1: on */

/* SASU implimentation: 2013/10/29 */
#define USE_SASU 0
/* 0: not */
/* 1: use for soil */
/* 2: use for soil + vegetation */

/* randomization of climate data in spin-up */
#define RANDCLIM_SU 1
/* 0: off */
/* 1: on */

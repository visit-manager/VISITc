/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & JAMSTEC							*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version diurnal-CFSR run                                            */

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito
	modify: 2013 / 01 / 04 by Akihiko Ito
*/

#define N_ANSDAY 32
#define N_ANSBM 32
#define N_CLIMD 7

#define N_OUTANN 12
#define N_OUTHR 12

/* **********************************************/
/* initial mass, Mg C ha-1 ****/
#define INIT_C 1.0
#define INIT_N 0.1

/* *********************************************************/
/* experimental area ********/

/* ERA5: 2021/03/30 by A.Ito */
/* ERA5: 2023/10/20 by A.Ito */
#define NROW 720
#define NCOL 1440

/* CFSR: 2012/01/02 by A.Ito */
//#define NROW 576
//#define NCOL 1152

/* Bamiyan: 2017/07/31 by A.Ito */
//#define NROW 483
//#define NCOL 844

/* East Asia: 2017/10/30 by A.Ito */
//#define NROW 160
//#define NCOL 320

/* JAPAN: 2021/12/27 by A.Ito */
//#define NROW 3120
//#define NCOL 3840

/* BHUTAN-BANGRADESH: 2022/04/22 by A.Ito */
//#define NROW 1200
//#define NCOL 1200

/* JAPAN center: 2022/09/22 by A.Ito */
//#define NROW 300
//#define NCOL 300

/* PAWCs center: 2023/06/13 by A.Ito */
//#define NROW 200
//#define NCOL 1490

/* experimental period ***********************/

/* ERA5: 2021/03/30 by A.Ito */
/* ERA5: 2023/10/20 by A.Ito */
#define BYR 2019
#define EYR 2021

/* CFSR: 2012/01/02 by A.Ito */
//#define BYR 1979
/* #define EYR 2010 */
//#define EYR 2015

/* Bamiyan: 2017/07/31 by A.Ito */
//#define BYR 1979
/* #define EYR 2010 */
/* #define EYR 2016 */
//#define EYR 2055

/* East Asia: 2017/10/30 by A.Ito */
//#define BYR 1979
//#define EYR 2015

/* ERA5: 2021/03/30 by A.Ito */
//#define BYR 2020
//#define EYR 2020

/* PAWCs: 2023/06/13 by A.Ito */
//#define BYR 2003
//#define EYR 2020

#define PERIOD (EYR-BYR+1)

/* ******************************************/
/* cauculate only every Nth grids */
#define RUNFAST 1
/* 1: every grid */

/* show message */
#define NOTICE 1
/* 0: no */
/* 1: with debug message */

/* Semi-analytical spi-up */
#define USE_SASU 1
/* 0: no */
/* 1: on */

#define P_MONI 0
// 0 all area
// 345988 Manaus
// 199160 TKY

/* leap-year: 2017/08/15 by A.Ito */
#define LEAP_YR 0

/* restart file */
#define USE_RESTART 0
/* 0: off */
/* 1: use restart file (skip spin-up) */
/* 2: use restart file (continue spin-up) */

/* **********************************************/
/* number of MAX biomes */
#define MAX_BIOME 128
/* number of analysis variables */
#define MAX_ANSIS 256

/* number of land-use data */
#define N_LUDATA 34 /* UNH harmonized data */

/* nitrogen cycle scheme */
#define N_CYCLE 0
/* 0: skip, i.e. not calculate */
/* 1: calculate */

/* first year of GHG data */
#define FDY_GHG 1750

/* *********************************************/
/* solve equilibrium leaf gas exchange */
#define SOLV_LEAF 0
/* 0: Discrete Newton */
/* 1: analytical */

/* CH4 & N2O schemes */
#define GHG_CALC 0
#define GHG_OUT 0

/* CH4 emission with a multi-layer model */
#define WH_CH4 0
/* 0: off */
/* 1: on */

/* stable carbon isotope */
#define SCI_CALC 0
/* 0: neglect */
/* 1: include */

/* flags for seasonal change in leaf properties */
#define AG_VMX 0		/* Vcmax */
#define AG_JVR 0		/* Jmax/Vcmax */
#define AG_LMA 0		/* LMA */
#define AG_RD 0			/* Rd */

/* calculation using potential vegetation */
#define VEG_POT 0
/* 0: off */
/* 1: on */

#define EX_URBAN_FOREST 0
/* 0: off */
/* 1: replaced by evergreen broadleaved forest */

/* experiment with fixed LAI */
#define FIX_LAI 0
/* 0: off */
/* 1: on */
/* 2: use fix-co2 LAI: fix LAI */
/* 3: use fix-co2 LAI: fix canopy N */

#define FIX_CO2 0
/* 0: off */
/* 1: on */
/* 2: on for stomata */

#define FIX_NTOP 0
/* 0: default */
/* 1: +10% */

#define FIX_PPFD 0
/* 0: default */
/* 1: all diffuse */
/* 2: all direct */
/* 3: diffuse:direct = 50:50 */

#define CONST_KN 1
/* 0: default */
/* 1: de Pury & Farquhar default (KN = 0.713 = Kdiff) */
/* 2: high (KN = 1.22: meadian meta-analysis: 2015/01/26 by A.Ito) */
/* 3: meta-analysis-based KN: 2015/04/17 by A.Ito */

/* constrain LAI max by Iio et al. (2014): 2016/08/15 by A.Ito */
#define CONST_LAIMAX 0
/* 0: off */
/* 1: on */

/* experiment of sensitivity to diffuse fraction: 2017/03/21 by A.Ito */
/* revised 2017/09/14 by A.Ito */
#define EX_PAR 0
/* 0: off */
/* 1: increase +10% diffuse fraction */
/* 2: increase -10% diffuse fraction */
/* 3: increase +10% cloudiness */
/* 4: increase -10% cloudiness */

/* experiment of sensitivity to cloudiness: 2017/05/30 by A.Ito */
#define EX_CLOUD 0
/* 0: off */
/* 1: increase cloudiness */

#define EX_BAMIYAN 0
/* 0: historocal only */
/* 1: projection for rcp2.6 */
/* 2: projection for rcp8.5 */

#define EX_BAMIYAN_S 0
/* 0: no grazing */
/* 1: constant grazing */
/* 2: population-proportional grazing */
/* 3: -2% grazing */
/* 4: -3% grazing */
/* 5: -4% grazing */

/* land use setting */
#define LANDUSE 0
/* 0: natural vegetation */
/* 1: no land-use change since 1901 */
/* 2: no land-use change since 1990 */
/* 3: linear land-use change since 1990 */
/* 4: high future land-use change, coupling */
/* 5: IMAGE2 land-use change since 1990 */
/* 6: EOS-WEBSTER Hurtt land-use change, 1700-2000 */
/* 7: Ramankutty land-use change, 1700-2007 */
/* 8: Hurtt harmonized land-use change, 1700-2010 (revised 2012/02/28) */

/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2007 / 05 / 17 by Akihiko Ito: Version 0.2.4.
*/

/* initialize hydrology */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

extern long	month_day[12];

/* additional message for hydrological scheme */
#define VERB_HYD 0

/* air density, kg m-3 ******************************************/
float f_airdens(
	struct Grid *grid, 
	struct Loct *loct
){
	float aaa, bbb, ccc, air_density;
	
	aaa = ZAT/(loct->tmp_2m + ZAT);
	bbb = loct->air_prsr/1013.25;
	ccc = 1.0 - 0.378 * loct->vp/loct->air_prsr;
	air_density = 1.293 * aaa * bbb * ccc;
	
	return (air_density);
}

/* saturation vapour pressure (in hPa) ******************************/
float f_vap_pre_sat(
	struct Grid *grid, 
	struct Loct *loct
){
	float vps;
	
	/* Tetens's equation */
	/* Mizu-Kankyo-no-Kisyogaku pp.26 Eq 2.14 */
	if(loct->tmp_2m > 0.0){ /* at water surface */
		vps = 6.1078 * pow(10.0,(7.5*loct->tmp_2m)/(237.3 + loct->tmp_2m));
	}else{ /* at ice surface */
		vps = 6.1078 * pow(10.0,(9.5*loct->tmp_2m)/(265.3 + loct->tmp_2m));
	}
	vps = (vps>=0.0)?vps:0.0;
	
	return (vps);
}

/* slope of the saturation vapour pressure curve (hPa K-1) ***************/
float f_slope_vps(
	struct Grid *grid,
	struct Loct *loct
){
	float slope, aaa, bbb, ccc;
	
	if(loct->tmp_2m > 0.0){ /* at water surface */
		/* Mizu-Kankyo-no-Kisyogaku pp26 Eq 2.15 */
		aaa = 6.1078 * (2500.0 - 2.4*loct->tmp_2m);
		bbb = 0.4615 * (ZAT + loct->tmp_2m)*(ZAT + loct->tmp_2m);
		ccc = pow(10.0, (7.5*loct->tmp_2m)/(237.3 + loct->tmp_2m));
	}else{ /* at ice surface */
		/* Mizu-Kankyo-no-Kisyogaku pp26 Eq 2.16 */
		aaa = 6.1078 * 2834.0;
		bbb = 0.4615*(ZAT + loct->tmp_2m)*(ZAT + loct->tmp_2m);
		ccc = pow(10.0, (9.5*loct->tmp_2m)/(265.3 + loct->tmp_2m));
	}
	slope = (aaa/bbb) * ccc;
	
	return (slope);
}

/* Aerodynamic Resistance ********************************************/
float f_r_aero(
	struct Grid *grid, 
	struct Loct *loct
){
	float k_con;
	float r_aero;
	
	k_con = 0.41; /* Karman's constant */
	
	loct->wnd_10m = (loct->wnd_10m>=0.1)?loct->wnd_10m:0.1;
	
	r_aero = (log(10.0)*log(10.0))/(k_con * k_con * loct->wnd_10m);
	r_aero = (r_aero>=0.1)?r_aero:0.1;
	r_aero = (r_aero<=59.5)?r_aero:59.5;
	
	return (r_aero);
}

/* Penman-Monteith Transpiration ****************************************/
float pm_incep(
	struct Grid *grid, 
	struct Loct *loct, 
	long canopy_layer
){
	float cp, psycon, spwt, eta, net_rad;
	float rc_p, interception;
	float aaa, bbb;
	
	spwt = loct->air_dns; /** density of air, in kg m-3 **/
	cp = 0.2813; /** specific heat of air, in W h kg-1 K-1 **/
	psycon = 0.667; /** psychlometer constant, in hPa K-1 **/
	eta = 0.0224*1.0/1000.0; /** unit conversion of conductance from mmol H2O m-2 s-1 to m s-1 **/
	
	/** no vegetation resistance **/
	rc_p = 0.0;
	
    net_rad = 0.0;
	if(canopy_layer == 1){
		net_rad = loct->rn_tree;
	}else if(canopy_layer == 2){
		net_rad = loct->rn_c3;
	}else if(canopy_layer == 3){
		net_rad = loct->rn_c4;
	}
	
	aaa = (loct->slope_vps * net_rad) + (cp * spwt * loct->vpd/loct->r_aero);
	bbb = loct->slope_vps + psycon*(1.0 + rc_p/loct->r_aero);	
	/*  interception = loct->daylen*aaa/bbb/LHT; */
	interception = aaa / bbb / LHT;
	interception = (interception>=0.0)?interception:0.0;
	
	return (interception);
}

/* Penman-Monteith Evaporation ***************************************/
float pm_evap(
	struct Grid *grid, 
	struct Loct *loct
){
	float cp, psycon, spwt, eta;
	float ggc, evaporation;
	float aaa, bbb;
	
	spwt = loct->air_dns; /** density of air, in kg m-3 **/
	cp = 0.2813; /** specific heat of air, in W h kg-1 K-1 **/
	/* cp = 1014.0; */ /** specific heat of air, in J kg-1 K-1 **/
	psycon = 0.667; /** psychlometer constant, in hPa K-1 **/
	eta = 0.0224*1.0/1000.0; /** unit conversion of conductance from mmol H2O m-2 s-1 to m s-1 **/
	
	/* ground resistance */
	ggc = 1000.0*(1.0 - (grid->fieldcap30 - loct->soilwtr_l)/grid->fieldcap30) + 100.0;
	loct->rsst_soil = 1.0/(ggc*eta);
	
	aaa = (loct->slope_vps * loct->rn_ground) + (cp * spwt * loct->vpd/loct->r_aero);
	bbb = loct->slope_vps + psycon * (1.0 + loct->rsst_soil/loct->r_aero);
	
	/* evaporation = loct->daylen*aaa/bbb/LHT; */
	evaporation = aaa / bbb / LHT;
	evaporation = (evaporation>=0.0)?evaporation:0.0;
	
	return (evaporation);
}

/* Penman-Monteith Transpiration ****************************************/
float pm_transp(
	struct Grid *grid, 
	struct Loct *loct, 
	long canopy_layer
){
	float cp, psycon, spwt, eta;
	float canopy_con, rad_net;
	float rc_p, transpiration;
	float aaa, bbb;
	
	spwt = loct->air_dns;	/** density of air, in kg m-3 **/
	cp = 0.2813; 		/** specific heat of air, in W h kg-1 K-1 **/
	psycon = 0.667; /** psychlometer constant, in hPa K-1 **/
	eta = 0.0224*1.0/1000.0; /** unit conversion of conductance from mmol H2O m-2 s-1 to m s-1 **/
	
    canopy_con = rad_net = 0.0;
	if(canopy_layer == 1){
		canopy_con = loct->gc_tree;
		rad_net = loct->rn_tree;
	}else if(canopy_layer == 2){
		canopy_con = loct->gc_c3;
		rad_net = loct->rn_c3;
	}else if(canopy_layer == 3){
		canopy_con = loct->gc_c4;
		rad_net = loct->rn_c4;
	}
	
	if(canopy_con>0.0 && rad_net>0.0){
		/** canopy resistance **/
		rc_p = 1.0/(canopy_con*eta);
		
		aaa = (loct->slope_vps * rad_net) + (cp * spwt * loct->vpd/loct->r_aero);
		bbb = loct->slope_vps + psycon*(1.0 + rc_p/loct->r_aero);	
		/* transpiration = loct->daylen*aaa/bbb/LHT; */
		transpiration = aaa / bbb / LHT;
	}else{
		transpiration = 0.0;
	}
	transpiration = (transpiration>=0.0)?transpiration:0.0;
	
	return (transpiration);
}

/* Daily CASA soil moisture submodel ****************************/
void f_casa_moisture(
	struct Grid *grid, 
	struct Loct *loct
){
	float fc;
	float pc;
	
	fc = grid->field_cap/300.0;
	pc = grid->pore_cap/300.0;
	
	loct->rdr_casa = (1.0 + grid->a_sw)/(1.0 + grid->a_sw * pow(loct->vmc, grid->b_sw));
	
	if(loct->prate_sfc >= loct->pet){
		loct->m_casa = loct->m_casa_pre + ((loct->prate_sfc - loct->pet)/1000.0);
	}else if(loct->prate_sfc < loct->pet){
		loct->m_casa = loct->m_casa_pre - ((loct->pet - loct->prate_sfc)/1000.0) * loct->rdr_casa;
	}
	
	if(loct->m_casa >= fc){
		loct->m_casa = fc;
	}
	if(loct->m_casa < 0.0){
		loct->m_casa = 0.0;
	}
	
	loct->e_casa = (((loct->prate_sfc - loct->pet)/1000.0) - 
		(fc - loct->m_casa_pre)); /* /(float)month_day[loct->month] */
	
	/* moinsture scalar */
	loct->iw_casa = 0.0;
	if(loct->e_casa > 0.0){
		/* Eq.3a in Potter et al. (1996) */
		loct->iw_casa = ((loct->e_casa + fc) / pc) * 100.0;
	}else if(loct->e_casa == 0.0){
		/* Eq.3b in Potter et al. (1996) */
		loct->iw_casa = (loct->m_casa / pc) * 100.0;
	}else{
		/* Eq.3b in Potter et al. (1996) */
		/*i_w = 0.0;*/
		loct->iw_casa = (loct->m_casa / pc) * 100.0;
	}
	
	loct->iw_casa = (loct->iw_casa<0.0)?0.0:loct->iw_casa;
	loct->iw_casa = (loct->iw_casa>=100.0)?100.0:loct->iw_casa;
}

/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* from foliage, maintenance respiration **********************************/
double f_rfm(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double rfmt0, t0, rfm, ft;
	
	/* specific rate, at 15 deg C */
	rfmt0 = pchar->rmf/1000.0; 
	t0 = 15.0; /* base temperature */
	/* temperature dependence, exponential */
	ft = exp(log(pchar->qTf)/10.0*(loct->tmp_sfc-t0));
	
	if(mass->fol >= 0.0){	
		rfm = mass->fol*rfmt0*ft;
	}else{
		rfm = 0.0;
	}
	
	return (rfm);
}

/***** from foliage, growth / construction respiration ****/
double f_rfg(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pflx *flux
){
	double rfg;
	
	/* construction cost is proportional to biomass growth */
	rfg = pchar->rgf * flux->tpf; 
	/* flux->tpf: carbon translocation */
	
	return (rfg);
}

/**** from stem and branch, maintenance respiration *****************/
double f_rcm(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double rfmt0, t0, rfm, ft;
	
	/* specific rate, at 15 deg C */
	rfmt0 = pchar->rmc/1000.0; 
	t0 = 15.0; 
	/* temperature dependence, exponential */
	ft = exp(log(pchar->qTc)/10.0*(loct->tmp_sfc - t0));
	
	if(mass->stm >= 0.0){	
		rfm = mass->stm * rfmt0 * ft;
	}else{
		rfm = 0.0;
	}
	
	return (rfm);
}

/***** from stem and branch, growth and construction respiration ****/
double f_rcg(
	struct Grid *grid, 
	struct Pchar *pchar,
	struct Pflx *flux
){
	double rcg;
	
	/* construction cost is proportional to biomass growth */
	rcg = pchar->rgc * flux->tpc; 
	
	return (rcg);
}

/**** from root, maintenance respiration ***************************/
double f_rrm(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double rfmt0, t0, rfm, ft;
	
	/* specific rate, at 15 deg C */
	rfmt0 = pchar->rmr/1000.0; 
	t0 = 15.0; 

	/* temperature dependence, exponential */
	/* ft=exp(log(pchar->qTr)/10.0*(grid->tmp_sfc-t0)); */
	ft = exp(log(pchar->qTr)/10.0*(loct->tmp10_soil-t0));
	
	if(mass->rot>=0.0){	
		rfm = mass->rot * rfmt0 * ft;
	}else{
		rfm = 0.0;
	}
	
	return (rfm);
}

/***** from root, growth and construction respiration ******/
double f_rrg(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pflx *flux
){
	double rrg;
	
	/* construction cost is proportional to biomass growth */
	rrg = pchar->rgr * flux->tpr; 
	
	return (rrg);
}

/***** 30-min step CO2 emission ************************************************************/
void f_resp_h(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	double ftl, fwl, fal, fsm;
	double rf_t1, rf_t2, rc_t, rc_f;
	double rf_f, rr_t, rr_f;
	double hr_tf, hr_tc, hr_tr;
	double hr_gf, hr_gc, hr_gr;
	double hr_a, hr_i, hr_p;
	double conv_dh;
	
	/* units conversion ***********************************/
	conv_dh = 100.0/12.0*1000000.0/48.0/1800.0;
	/*	100.0:		Mg C ha-1		 to	g C m-2			*/
	/*	1.0/12.0:		g C m-2		 to	mol C m-2			*/
	/*	1000000.0:	mol C		 to	micro mol C		*/
	/*	1.0/48.0:		per 48 time-steps in a day				*/
	/*	1.0/1800.0:	per second (30 minutes)				*/
		
	/* leaf dark respiration rate for canopy ***********************************/
	rf_t1 = (echar->overs1).lai*(echar->overs1).rd0*(echar->overs1).ft_rd;
	rf_t2 = (echar->overs2).lai*(echar->overs2).rd0*(echar->overs2).ft_rd;
	rf_f = (echar->unders).lai*(echar->unders).rd0*(echar->unders).ft_rd;
	
	(echar->overs1).rd = rf_t1; /* unit: micro mol CO2 m-2 s-1 */
	(echar->overs2).rd = rf_t2;
	(echar->unders).rd = rf_f;
	
	/*************************************/
	ftl = exp(log((echar->tree).qTc)/10.0*((echar->overs1).tmp-15.0)); /* */
	/* ftl = exp(log(2.0)/10.0*((echar->overs1).tmp-15.0)); */	/* Q10 = 2.0 */
	rc_t = (mass->tree).stm*(echar->tree).rmc/1000.0*ftl + (flux->tree).rcg;
	rc_f = (mass->c3).stm*(echar->c3).rmc/1000.0*ftl + (flux->c3).rcg;
	rr_t = (mass->tree).rot*(echar->tree).rmr/1000.0*ftl + (flux->tree).rrg;
	rr_f = (mass->c3).rot*(echar->c3).rmr/1000.0*ftl + (flux->c3).rrg;
	rc_t *= conv_dh;
	rc_f *= conv_dh;
	rr_t *= conv_dh;
	rr_f *= conv_dh;
	
	/* unit: micro mol CO2 m-2 s-1 */
	(flux->tree).ar_h = rf_t1*loct->comp_over1 + rf_t2*loct->comp_over2 + rc_t + rr_t;
	(flux->c3).ar_h = rf_f + rc_f + rr_f;
	
	/* temperature effect, exponential */
	if((echar->soil).tmp>-20.0){
		ftl = 0.01+exp(308.56*(1.0/56.02-1.0/((echar->soil).tmp+46.02)));
	}else{
		ftl = 0.01;
	} 
	
	/* soil moisture effect, saturating */
	fwl = 0.8*loct->soilwtr_l/((echar->soil).kml*grid->fieldcap30 + loct->soilwtr_l) + 0.2;
	/* soil apparence effect */
	fal = 0.4*loct->soilappr_l*(1.0*(echar->soil).kmsl)/((echar->soil).kmsl + loct->soilappr_l) + 0.6;
	/* co-limited */
	fsm = (fwl>fal)?fal:fwl;
	
	hr_tf = (mass->soil).ltr_tf*(echar->soil).sr_lf/1000.0*ftl*fsm * conv_dh;
	hr_tc = (mass->soil).ltr_tc*(echar->soil).sr_lc/1000.0*ftl*fsm * conv_dh;
	hr_tr = (mass->soil).ltr_tr*(echar->soil).sr_lr/1000.0*ftl*fsm * conv_dh;
	
	hr_gf = (mass->soil).ltr_gf*(echar->soil).sr_lf/1000.0*ftl*fsm * conv_dh;
	hr_gc = (mass->soil).ltr_gc*(echar->soil).sr_lc/1000.0*ftl*fsm * conv_dh;
	hr_gr = (mass->soil).ltr_gr*(echar->soil).sr_lr/1000.0*ftl*fsm * conv_dh;
	
	hr_a = (flux->soil).mr_a * conv_dh;
	hr_i = (flux->soil).mr_i * conv_dh;
	hr_p = (flux->soil).mr_p * conv_dh;
	
	/* unit: micro mol CO2 m-2 s-1 */
	(flux->soil).hr_h = hr_tf + hr_tc + hr_tr + hr_gf + hr_gc + hr_gr + hr_a + hr_i + hr_p;
}

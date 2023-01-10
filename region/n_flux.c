/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/*
	created 2007 / 09 / 20	by A.Ito
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

extern long	WMODE;

/** NH3 volatilization **********************************************/
void f_nh3_volatilization(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Smas *mass, 
	struct Sflx *flux
){
	float f_ph, f_tmp, f_sw;
	float swp, nh4_soil;
	
	/* g N ha-1 */
	/*  nh4_soil = grid->total_n_1m * 0.3 * 0.01*10000.0 *0.5;  */
	/*  nh4_soil = (mass->soil).n_nh4 /1000.0 /10000.0;  */
	nh4_soil = mass->n_nh4;
	
	/* pH: Lin et al. (2000) Eq.(16) */
	f_ph = 5.8/30.0 * pow(10.0, grid->soil_ph - 10.0);
	if(f_ph<0.0){
		f_ph = 0.0;
	}
	
	/* temperature: Thornley (1998) Eq.(3.11a) */
	f_tmp = pow((loct->tmp10_soil - 0.0), 2.0)*(45.0 - loct->tmp10_soil) 
			/ (pow((20.0 - 0.0), 2.0)*(45.0 - 20.0));
	if(f_tmp<0.0){
		f_tmp = 0.0;
	}
			
	/* soil water: Thornley (1998) Eq.(6.7a) */
	/* soil water potential Eq.(6.2g) */
	if(loct->soilwtr_l > 0.0){
		swp = -10.0 * pow((grid->fieldcap30/loct->soilwtr_l), 5.0);
		f_sw = exp((18.0 * swp)/(8314.0*(loct->tmp10_soil + ZAT)));
	}else{
		f_sw = 0.0;
	}
	if(f_sw<0.0){
		f_sw = 0.0;
	}
	
	/* Thornley (1998) Eq.(5.4i) */
	/* g NH3 ha-1 day-1 */
	flux->n_nh3vlt = nh4_soil * 0.02 * f_ph * f_tmp * f_sw 
			* 17.0/14.0;
}

/* biological N2 fixation model by Cleveland **********************************/
/*
Cleveland, C.C., Townsend, A.R., Schimel, D.S., Fisher, H., Howarth, 
R.W., Hedin, L.O., Perakis, S.S., Latty, E.F., Von Fischer, J.C., 
Elseroad, A. and Wasson, M.F., 1999. Global patterns of terrestrial 
biological nitrogen (N2) fixation in natural ecosystems. 
Global Biogeochemical Cycles 13, 623-645.
*/
/**  biological N2 fixation, both symbiotic and non-symbiotic  **/
void f_biolfix(
	struct Grid *grid, 
	struct Loct *loct,
	struct Flux *flux
){
	float aet, total_nbiofix;
	
	/* annual actual evapotranspiration */
	aet = loct->aet;
	aet *= 100.0; /* 100 mm day-1 */
	
	/* Cleveland et al. (1999) Figure 1, Central */
	total_nbiofix = 0.234 * aet - 0.172; /* kg N ha-1 yr-1 */
	/* total_nbiofix *= 1000.0; */ /* g N ha-1 yr-1 */
	if(total_nbiofix < 0.0){
		total_nbiofix = 0.0;
	}
	
	(flux->tree).n_biofix = total_nbiofix * 0.5;
	(flux->c3).n_biofix = total_nbiofix * 0.5 * loct->funder_c3;
	(flux->c4).n_biofix = total_nbiofix * 0.5 * loct->funder_c4;
}

/* leaching of NO3 ************************************************/
void f_n_leaching(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Smas *mass, 
	struct Sflx *flux
){
	float ntr_conc, aa;

	/* g N / ha */
	/* kg H2O / m2 */
	if(loct->soilwtr_l > 0.0){
		ntr_conc = 0.005 * (mass->n_no3/10000.0) / loct->soilwtr_l;
	}else{
		ntr_conc = 0.0;
	}
	/* g N / kg H2O */
	
	aa = loct->ro2 * ntr_conc;
	if(aa > mass->n_no3){
		aa = mass->n_no3;
	}
	
	/* g NO3-N m-2 month-1 */
	flux->n_leach = aa*10000.0;
	/* g NO3-N ha-1 month-1 */
}

/* N deposition ********************************************************/
/* 
Dentener, F. J. 2006. Global Maps of Atmospheric Nitrogen Deposition, 
1860, 1993, and 2050. Data set. Available on-line from Oak Ridge National 
Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, U.S.A. 
*/
void f_n_deposit(
	struct Grid *grid, 
	struct Loct *loct
){
	float pre_ann, ndepo_dry, ndepo_wet, aa;
	float f_no3, f_nh4;
		
	pre_ann =grid->prec_ann + 1.0;

    f_no3 = f_nh4 = 0.0;
	if(WMODE == 1){
		/* kg N ha-1 yr-1 */
		if(strcmp(grid->site_id, "TKY")==0){
			aa = 28.7;		
			f_no3 = 0.7297;
			f_nh4 = 1.0 - f_no3;
		}else{
			f_no3 = 0.5;
			f_nh4 = 1.0 - f_no3;
			aa = 10.0;
		}
	
	/*	switch(grid->site_id){
			case TKY:		
			case 2:		
				f_no3 = 0.5;
				f_nh4 = 1.0 - f_no3;
				aa = 10.0;		
				break;
			case 3:		
				f_no3 = 0.5;
				f_nh4 = 1.0 - f_no3;
				aa = 10.0;		
				break;
			case 4:		
				aa = 16.1;		
				f_no3 = 0.8070;
				f_nh4 = 1.0 - f_no3;
				break;
			default:	
				f_no3 = 0.5;
				f_nh4 = 1.0 - f_no3;
				aa = 10.0;
		} */

		ndepo_dry = 0.5 * aa /365.0;
		ndepo_wet = 0.5 * aa * loct->prate_sfc/pre_ann;
	}else{
		/* kg N ha-1 yr-1 */
		loct->ndepo[0] = 10.0;
		loct->ndepo[1] = 10.0;
		loct->ndepo[2] = 10.0;
		
        ndepo_dry = ndepo_wet = 0.0;
		if(loct->climy<=1850){
			ndepo_dry = 0.5 * loct->ndepo[0] /365.0;
			ndepo_wet = 0.5 * loct->ndepo[0] * (loct->prate_sfc + 1.0)/pre_ann;
		}else if(loct->climy>1850 && loct->climy<=1993){
			aa = loct->ndepo[0] + (loct->ndepo[1] - loct->ndepo[0])*(float)(loct->climy-1850)/143.0;
		
			ndepo_dry = 0.5 * aa /365.0;
			ndepo_wet = 0.5 * aa * (loct->prate_sfc + 1.0)/pre_ann;
		}else if(loct->climy>1993 && loct->climy<=2050){
			aa = loct->ndepo[1] + (loct->ndepo[2] - loct->ndepo[1])*(float)(loct->climy-1993)/57.0;
		
			ndepo_dry = 0.5 * aa /365.0;
			ndepo_wet = 0.5 * aa * (loct->prate_sfc + 1.0)/pre_ann;
		}else if(loct->climy>2050){
			ndepo_dry = loct->ndepo[2] * 0.5 /365.0;
			ndepo_wet = loct->ndepo[2] * 0.5 * (loct->prate_sfc + 1.0)/pre_ann;
		}
		
		if(ndepo_dry<0.0){
			ndepo_dry = 0.0;
		}
		if(ndepo_wet<0.0){
			ndepo_wet = 0.0;
		}
	}
		
	/* original unit: mg N m-2 yr-1 */
	/* g N ha-1 day-1 */
	loct->depo_no3 = f_no3 * (ndepo_dry + ndepo_wet) * 1000.0;
	loct->depo_nh4 = f_nh4 * (ndepo_dry + ndepo_wet) * 1000.0;
}

/* plant N uptake *******************************************/
void f_n_uptake(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	float f_temp;
	float n_max_no3 = 1.9;
	float n_max_nh4 = 1.4;
	float ks, navil;
	float uptake_no3, uptake_nh4;

	/* Raich et al. (1991): Appendix  Eq.1.16 */
	
	/*****
	Effect of N allocatiom to root ?
	*****/
	
	/* temperature factor */
	f_temp = exp(0.0693 * loct->tmp10_soil);
	/* soil diffusion parameter */
	ks = 0.90 * pow(loct->soilwtr_l/grid->fieldcap30, 3.0) + 0.1;
	
	/* NO3 uptake */
	navil = (mass->soil).n_no3/10000.0;
	uptake_no3 = navil * n_max_no3 * ks / (90.0 + ks*navil) * f_temp;
	/* g N ha-1 month-1 */
	(flux->tree).uptake_no3 = uptake_no3 * 0.5 * 10000.0;
	(flux->c3).uptake_no3 = uptake_no3 * 0.5 * loct->funder_c3 * 10000.0;
	(flux->c4).uptake_no3 = uptake_no3 * 0.5 * loct->funder_c4 * 10000.0;
	
	/* NH4 uptake */
	navil = (mass->soil).n_nh4/10000.0;
	uptake_nh4 = navil * n_max_nh4 * ks / (90.0 + ks*navil) * f_temp;
	/* g N ha-1 month-1 */
	(flux->tree).uptake_nh4 = uptake_nh4 * 0.5 * 10000.0;
	(flux->c3).uptake_nh4 = uptake_nh4 * loct->funder_c3 * 10000.0;
	(flux->c4).uptake_nh4 = uptake_nh4 * loct->funder_c4 * 10000.0;
}

/* N abandoned as litter *********************************************/
void f_n_abandon_salvage(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	pchar->n_salvage = 0.5;
	pchar->cn0_strg = 0.5;

	/* canopy *******/
	/* abandon */
	flux->n_abdn_cnpy = (1.0 - pchar->n_salvage) * pchar->lf * 
					mass->n_cnpy;
	/* salvage (draw back to storage) */
	flux->n_salvage = pchar->n_salvage * pchar->lf * 
					mass->n_cnpy;

	/* storage *******/
	flux->n_abdn_strg = (flux->lc + flux->lr)
				* 1000000.0 / pchar->cn0_strg;
	if(flux->n_abdn_strg > mass->n_strg*0.9){
		flux->n_abdn_strg = mass->n_strg*0.9;
	}
}

/* N immobilization by microbes **********************************/
void f_n_immoblz(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Schar *schar, 
	struct Smas *mass, 
	struct Sflx *flux
){
	/* flux->n_immbl = 0.07 * flux->n_minerlz_lttr + 1.2 * flux->n_minerlz_hums + 
		(0.004 * mass->n_no3 + 0.03 * mass->n_nh4); */

	flux->n_immbl = 0.07 * flux->n_minerlz_lttr + 0.8 * flux->n_minerlz_hums + 
		(0.02 * mass->n_no3 + 0.03 * mass->n_nh4); /* TMK */

	/* flux->n_immbl = 0.070 * flux->n_minerlz_lttr + 1.2 * flux->n_minerlz_hums + 
		(0.004 * mass->n_no3 + 0.03 * mass->n_nh4);   FJY */
}

/* N abandoned from microbes ******************/
void f_n_mcrb_abdn(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Schar *schar, 
	struct Smas *mass, 
	struct Sflx *flux
){
	float f_temp;
	
	f_temp = exp(log(2.0)/10.0 * (loct->tmp10_soil - 10.0));

	flux->n_mcrb_abdn = 1.1 * f_temp * mass->n_mcrb;  /*  TKY */
	/*  flux->n_mcrb_abdn = 1.1 * f_temp * mass->n_mcrb;   FJY */
}

/* allocation of uptake N ******************************/
void f_n_alloc(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	float n_obtain;
	float n_leaf_conc;
	float amax;	/* maximum photosynthetic rate at high N */
	float kn;		/* half-saturation constant */
	float rd_nsp;	/* specific dark respiration rate for N */
	float a_nlmt;
	float n_demand1;
	
	/* nitrogen obtained by root uptake and symbiotic fixation */
	n_obtain = flux->n_biofix + flux->uptake_no3 + flux->uptake_nh4;

	amax = pchar->amax_nphoto = 30.0;
	kn = pchar->kn_nphoto = 100.0;
	rd_nsp = pchar->rd_n = 0.05;
	
	/* optimal N for nitrogen utilization */
	pchar->n_opt_leaf = sqrt(amax*kn/rd_nsp) - kn;
	
	/* leaf N concentration, mmol N m-2 */
	n_leaf_conc = pchar->n_conc_larea;
	a_nlmt = amax * n_leaf_conc / (kn + n_leaf_conc);
	
	if(pchar->n_opt_leaf > n_leaf_conc){
		n_demand1 = (pchar->n_opt_leaf - n_leaf_conc) * 14.0 / 1000.0 * 10000.0;
		
		if(n_obtain > n_demand1){
			flux->n_alloc_cnpy = n_demand1;
		}else{
			flux->n_alloc_cnpy = n_obtain;
		}
	}else if(pchar->n_opt_leaf <= n_leaf_conc){
		
		/* n_demand2 = 1.0 * (amax - a_nlmt);
		
		if(n_obtain > n_demand2){
			flux->n_alloc_cnpy = n_demand2;
		}else{
			flux->n_alloc_cnpy = n_obtain;
		} */
		
		flux->n_alloc_cnpy = 0.0;
	}
	
	flux->n_alloc_strg = n_obtain - flux->n_alloc_cnpy;
}

/* N reallocation from storage to canopy **********/
void f_n_realloc(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	float amax;	/* maximum photosynthetic rate at high N */
	float kn;		/* half-saturation constant */
	float rd_nsp;	/* specific dark respiration rate for N */
	float n_leaf_conc;
	float n_demand1;
	float n_stock;
	float a_nlmt;
	
	/* nitrogen stock for N reallocation */
	n_stock = mass->n_strg;

	amax = pchar->amax_nphoto = 30.0;
	kn = pchar->kn_nphoto = 100.0;
	rd_nsp = pchar->rd_n = 0.05;

	/* optimal N for nitrogen utilization */
	pchar->n_opt_leaf = sqrt(amax*kn/rd_nsp) - kn;
	
	/* leaf N concentration, mmol N m-2 */
	n_leaf_conc = pchar->n_conc_larea;
	a_nlmt = amax * n_leaf_conc / (kn + n_leaf_conc);
	
	if(pchar->n_opt_leaf > n_leaf_conc){
		/* g N ha-1 month-1 */
		n_demand1 = (pchar->n_opt_leaf - n_leaf_conc) * 14.0 / 1000.0 * 10000.0;
		
		if(n_stock > n_demand1){
			flux->n_realloc = n_demand1;
		}else{
			flux->n_realloc = n_stock;
		}
	}else if(pchar->n_opt_leaf <= n_leaf_conc){
		
		/* n_demand2 = 1.0 * (amax - a_nlmt);
		
		if(n_stock > n_demand2){
			flux->n_realloc = n_demand2;
		}else{
			flux->n_realloc = n_stock;
		} */
		
		flux->n_alloc_cnpy = 0.0;
	}
}

/*************************************************/
void f_n_mineralz(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Smas *mass, 
	struct Sflx *flux
){
	float f_c_min;
	float snmass, decom;
	
	/** litter **/
	snmass = mass->ltr_tf + mass->ltr_tc + mass->ltr_tr + 
			mass->ltr_gf + mass->ltr_gc + mass->ltr_gr;
	decom = flux->mr_tf + flux->mr_tc + flux->mr_tr + 
			flux->mr_gf + flux->mr_gc + flux->mr_gr;
	if(snmass > 0.0){
		f_c_min = decom / snmass;
	}else{
		f_c_min = 0.0;
	}
	
	/* g N ha-1 month-1 */
	flux->n_minerlz_lttr = mass->n_lttr * f_c_min * 1.0;  /*  */
	/* flux->n_minerlz_lttr = mass->n_lttr * f_c_min * 0.9; FJY */

	/** humus **/
	snmass = mass->msl_a + mass->msl_i + mass->msl_p;
	decom = flux->mr_a + flux->mr_i + flux->mr_p;
	if(snmass > 0.0){
		f_c_min = decom / snmass * (1.0 - 0.75); 
		/*  f_c_min = decom / snmass * (1.0 - 0.90);   FJY */
	}else{
		f_c_min = 0.0;
	}
	
	/* g N ha-1 month-1 */
	flux->n_minerlz_hums = mass->n_hums * f_c_min;
}


/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/*  [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* ecophysiological vegetation processes ******************************/
void f_ecophysiology(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar, 
	struct Pmas *mass
){
	long g;	
	float aaa, bbb;
	
	/** initial intercellular CO2 concentration, ppmv ****************/
	pchar->ci = loct->aCO2*0.7;  /* typically, 70% of ambient */
	
	f_n_leaf_conc(grid, pchar, mass);
	
	/* stabilization of single-leaf processes ***************************/
	/* should be modified to Newton-like iteration method ? */
	for(g=0; g<6; g++){			
		 /* quantum yield , mol CO2 mol photon-1*/
		f_photo_qy(loct, pchar);	
		
		/* light-saturated photosynthesis rate, micromol m-2 s-1 */
		f_pc_sat(grid, loct, pchar);

		/* canopy-top photosynthetic rate*/
		aaa = pchar->psat*pchar->lue*loct->ppfd_h;
		bbb = pchar->psat+pchar->lue*loct->ppfd_h;
		if(bbb>0.0){
			pchar->ptop = aaa/bbb;
		}else{
			pchar->ptop = 0.0;
		}

		/* stomatal conductance, mmol H2O m-2 s-1 */
		f_stom_cond(grid, loct, pchar); 

		/* intercellular CO2 concentration, ppmv */
		f_incelco2(loct, pchar); 
	}
	/** photosynthetic 13C discrimination, permille **/
	if(SCI_CALC == 1){
		photo_13c_frac(grid, loct, pchar);
	}
	
	/** canopy conductance **/	
	pchar->gc = f_canopy_cond(grid, loct, pchar, mass);

	/** plant respiration **/
	f_q10_ar(loct, pchar);	/* Q10 */
	f_spcfc_resp(pchar, mass);  /* woody specific respiration rate */

	/** litterfall of plant respiration **/
	f_mortality(grid, pchar);

	/** optimum leaf area index **/	
	f_opt_lai(grid, loct, pchar);
	
	/* maximum non-structural carbohydrate storage in plant body */
	switch(grid->veg_type){
		case 1: case 2: case 3: case 4: case 5: case 6: case 7: 
			pchar->max_nsch_storage = 0.1*mass->stm_sp + 0.3*mass->rot_fn;
			break;
		case 8: case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16: 
			if(pchar->phototype == 3){
				pchar->max_nsch_storage = 0.4*mass->stm_sp + 0.8*mass->rot_fn;
			}else if(pchar->phototype == 4){
				pchar->max_nsch_storage = 0.2*mass->stm_sp + 0.4*mass->rot_fn;
			}
			break;
	}
}

/* leaf area index *****************************/
float lai_mass(
	struct Pmas *mass, 
	struct Pchar *pchar
){
	float sla, lai_est;
	
	/** specific leaf area as a function of... what? **/
	sla = pchar->sla;
	
	lai_est = sla*mass->fol * dmTc /100.0 /2.0;
	lai_est = (lai_est>=0.0)?lai_est:0.0;
	
	/* dmTc: dry-matter to carbon */
	/* 100.0: cm2 g dm-1 to Mg ha-1 base */
	/* 2.0: single-sided leaf area */	
	
	/*  lai_est = 1.0;  */
    
	return (lai_est);
}

/* give irradiance attenuation coefficient ****************************/
float irr_attn(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	float aaa, bbb;
	
	/* a function of solar hight angle */
	aaa = sin(loct->solhgt_h*dTr); 
	aaa = (aaa <= 1.0)?aaa:1.0; 
	aaa = (aaa >= 0.3)?aaa:0.3; /* to avoid extreme values*/
	bbb = pchar->eK0/aaa;
	
	return (bbb);
}

/*  intercellular CO2 concentration  **********************************/
void f_incelco2(
	struct Loct *loct,
	struct Pchar *pchar
){
	float ci, gs_co2;
	/* give intercellular CO2 concentration, as a function of ambient CO2 level and stomatal conductance */
	
	gs_co2 = pchar->gs/1.56; 
	/* 1.56: conversion from H2O to CO2 conductance */
	
	ci = loct->aCO2 - (pchar->ptop/(gs_co2/1000.0));
	/* 1000.0: conbert from mmol to micro mol */
	
	ci = (ci>=0.0)?ci:0.0;
	ci = (ci<=loct->aCO2)?ci:loct->aCO2;
		
	pchar->ci = ci; 
}

/*  quantum yield C3 and C4 plants ***********************************/
void f_photo_qy(
	struct Loct *loct,
	struct Pchar *pchar
){
	float eftem, efci;
	
	/* basic data: Osmond et al. (1980) 
	Physiological processes in plant ecology, 
	toward a synthesis with Atriplex, Springer  */
	/** plant_type: 3=C3,  4=C4, (5=CAM)  **/
	if(pchar->phototype == 3){
		/* temperature dependence */
		eftem = (52.0 - loct->tmp_sfc)/(3.5 + 0.75*(52.0 - loct->tmp_sfc)); 
		/* CO2 dependence */
		efci = pchar->ci/(90.0 + 0.6*pchar->ci); 
		/* 3.5, 52.0, etc.: empirical parameters */
	}else{  /* pchar->phototype == 4 */
		/* insensitive QE of c3 species */
		eftem = 1.0;
		efci = 1.0;
	}
	
	/* eftem = 1.0; */
	
	/* give quantum yield */
	pchar->lue = pchar->lue0 * eftem * efci;
}

/* stomatal conductance ************************************/
void f_stom_cond(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	float b1d, cc;
	
	/* stomatal conductance model by Ball, Woodraw, and Berry (1987) */
	/* modified by Leuning (1995) */
	b1d = pchar->gs_b1 / ((loct->aCO2 - pchar->cmpcd)*(1.0 + loct->vpd/pchar->gs_b2));
	
	/** add soil water factor **/
	cc = 1.0; /* not defined yet */

	if(pchar->ptop>0.0){
		pchar->gs = pchar->gs_b0 + b1d*pchar->ptop*cc;
	}else{
		pchar->gs = pchar->gs_b0; 
	}
}

/* canopy conductance ********************************/
float f_canopy_cond(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float gpp, gc;
	float cc1, cc2, cc3, bb, f_co2, f_vpd;
	
	/* give daily GPP based on the Monsi-Saeki theory, or Kuroiwa's equation */
	if(pchar->psat > 0.0){
		cc1 = 2.0 * pchar->psat / pchar->eK;
		bb = pchar->eK * pchar->lue * pchar->ppfd_top/pchar->psat;	/* check pchar->ppfd_top : 2009/04/19  *****/
		cc2 = 1.0 + sqrt(1.0 + bb);
		cc3 = 1.0 + sqrt(1.0 + bb*exp(-1.0*pchar->eK*mass->lai));
		gpp = cc1 * log(cc2/cc3);

	/*	cc1 = 2.0*pchar->psat*loct->daylen[grid->doy]/pchar->eK; 
		bb = pchar->eK*pchar->lue*loct->par_h[grid->hour]/pchar->psat;
		cc2 = 1.0+sqrt(1.0+bb);
		cc3 = 1.0+sqrt(1.0+bb*exp(-1.0*pchar->eK*mass->lai));
		gpp = cc1*log(cc2/cc3)/24.0;  */
	}else{
		gpp = 0.0;
	}

	f_co2 = 1.0/(loct->aCO2 - pchar->cmpcd);
	f_vpd = 1.0/(1.0 + loct->vpd/pchar->gs_b2);
	bb = pchar->gs_b1 * (f_co2*f_vpd);

	gc = pchar->gs_b0 * mass->lai + bb*gpp;
	
	return (gc);
}

/* optimum LAI by Kuroiwa (1966) **************************/
void f_opt_lai(
	struct Grid *grid,
	struct Loct *loct, 
	struct Pchar *pchar
){
	float aaa, bbb, ccc, cc3, ddd, eee;
	float arm, arg, ar;	
	
	aaa = 1.0 / pchar->eK;
	bbb = pchar->eK * pchar->lue * pchar->ppfd_top;
	
	/* daily respiratory cost */
	eee = log(pchar->qTc)/10.0 * (loct->tmp_sfc - 15.0);
	arm = pchar->rmf * exp(eee) / 1000.0*dmTc*10000.0/(pchar->sla);
	arg = pchar->lf * dmTc * 10000.0 / (pchar->sla)*(1.0 + pchar->rgf);
	ar = arm + arg;
	
	if(pchar->psat > 0.0){
		cc3 = (pchar->psat * loct->daylen)/
				(pchar->psat * loct->daylen - ar*24.0);
		ccc = pchar->psat*(cc3 - 1.0);

		ddd = bbb / ccc;
		ddd = (ddd > 1.0)?ddd:1.0;
	
		pchar->opt_lai = aaa * log(ddd);
	}else{
		pchar->opt_lai = 0.0;
	}
}

/* Q10 of autotrophic respiration **************************/
void f_q10_ar(
	struct Loct *loct, 
	struct Pchar *pchar
){
	float aaa;
	
	/* seasonal acclimation in Q10 ****************/
	/* larger at cool and smaller at warm */
	/* ex. Yokota and Hagihara (1996) Tree Physiology */
	aaa = exp(-0.009*(loct->tmp_sfc - 15.0)); 
	
	pchar->qTf = pchar->qTf0 * aaa;
	pchar->qTc = pchar->qTc0 * aaa;
	pchar->qTr = pchar->qTr0 * aaa;
}

/* lowering specific respiration rate of massive biomass with heart wood tissue *******/
void f_spcfc_resp(
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float powstm, powrot;
	
	pchar->rmf = pchar->rmf0;
	
	/* specific respiration increasing in a power of 2/3 manner */
	/* size-dependence: Yokota and Hagihara (1998) Tree Physiology */
	powstm = 1.0 - 0.33334 * mass->stm/(pchar->f_sz_s + mass->stm);
	powrot = 1.0 - 0.33334 * mass->rot/(pchar->f_sz_r + mass->rot);
	
	mass->stm_sp = pow(mass->stm, powstm); /* sapwood mass in stem */
	if(mass->stm_sp > mass->stm){
		mass->stm_sp = mass->stm;
	}
	mass->stm_ht = mass->stm - mass->stm_sp; /* heartwood mass in stem */
	
	mass->rot_fn = pow(mass->rot, powrot); /* fine root */
	if(mass->rot_fn > mass->rot){
		mass->rot_fn = mass->rot;
	}
	mass->rot_tp = mass->rot - mass->rot_fn; /* tap[coarse] root */
		
	pchar->rmc = (pchar->rmc_s*mass->stm_sp + pchar->rmc_h*mass->stm_ht)/(mass->stm + 0.00001);
	pchar->rmr = (pchar->rmr_s*mass->rot_fn + pchar->rmr_h*mass->rot_tp)/(mass->rot + 0.00001);
}

/* mortarity of plant organs *************************/
void f_mortality(
	struct Grid *grid, 
	struct Pchar *pchar
){
	float bbb;
	
	/* larger at warm, and smaller at cool */
	bbb = 1.0; /* +(grid->tmp_sfc-grid->tmp_sfc_am)/100.0; */
	/* !! pending 2005/Feb/21  */
	
	pchar->lf = pchar->lf0 * bbb; /* leaf */
	pchar->lc = pchar->lc0 * bbb; /* stem */
	pchar->lr = pchar->lr0 * bbb; /* root */
}

/* leaf N concentration ****************************/
void f_n_leaf_conc(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	/* mmol N m-2 leaf area */
	if(mass->lai > 0.01){
		pchar->n_conc_larea = mass->n_cnpy / 10000.0 / 14.0 * 1000.0 / mass->lai;
	}else{
		pchar->n_conc_larea = 1.0;
	}
		
	/* mmol N g-1 leaf weight */
	if(mass->fol > 0.01){
		pchar->n_conc_lmass = mass->n_cnpy / 14.0 *1000.0 / (mass->fol * 1000000.0);
	}else{
		pchar->n_conc_lmass = 1.0;
	}
    
    /* De Pury & Farquhar canopy N *******************/
    pchar->photocap_n = pchar->photocap_n0;
}



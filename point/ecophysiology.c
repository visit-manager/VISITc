/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* ecophysiological vegetation processes ******************************/
void f_ecophysiology(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar, 
	struct Pmas *mass
){
	long g;	
	double aaa, bbb, ccc, ddd;
	
	/** initialize intercellular CO2 concentration, ppmv ****************/
	pchar->ci = loct->aCO2*0.7;  /* typically, 70% of ambient */
	
	/* stabilization of single-leaf processes ***************************/
	/* should be modified to Newton-like iteration method ? */
	for(g=0; g<6; g++){			
		 /* quantum yield , mol CO2 mol photon-1*/
		f_photo_qy(loct, pchar);	
		
		/* light-saturated photosynthesis rate, micromol m-2 s-1 */
		f_pc_sat(grid, loct, pchar);

		/* canopy-top photosynthetic rate*/
		aaa = pchar->psat*pchar->lue*loct->ppfd_h[loct->hour];
		bbb = pchar->psat+pchar->lue*loct->ppfd_h[loct->hour];
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
	if(SCI_SCHEME==1){
		photo_13c_frac(grid, loct, pchar);
	}
	
	/* canopy properties *************************************************/	
	/** canopy conductance **/	
	pchar->gc = f_canopy_cond(grid, loct, pchar, mass);

	/* respiration properties **********************************************/	
	/** plant respiration **/
	f_q10_ar(loct, pchar);	/* Q10 */
	f_spcfc_resp(pchar, mass);  /* woody specific respiration rate */

	/* litterfall properties **********************************************/	
	/** litterfall of plant respiration **/
	f_mortality(grid, loct, pchar);

	/* allocation properties **********************************************/	
	/** optimum leaf area index **/	
	f_opt_lai(grid, loct, pchar);
	
	/* nitrogen allocation properties */
	f_n_conc(grid, pchar, mass);
	
	/* maximum non-structural carbohydrate storage in plant body */
	switch(grid->veg_type){
		case 1: case 2: case 3: case 4: case 5: case 6: case 7: 
			pchar->max_nsch_storage = 0.1*mass->stm_sp + 0.3*mass->rot_fn;
			break;
		case 8: case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16: 
			if(pchar->phototype == 3){
				pchar->max_nsch_storage = 0.1*mass->stm_sp + 0.4*mass->rot_fn;
			}else if(pchar->phototype == 4){
				pchar->max_nsch_storage = 0.1*mass->stm_sp + 0.4*mass->rot_fn;
			}
			break;
	}
    
    /* rooting depth */
	g = 1;
	ccc = 0.0;
	do{
		aaa = 0.01*(double)g;
		
		ddd = (1.0 - 0.5*(exp(-pchar->root_dpt_a*aaa) + exp(-pchar->root_dpt_b*aaa)));
		bbb = (mass->rot*1000000.0/10000.0) * (ddd - ccc);
		
		ccc = ddd;
		g++;
	}while(bbb>0.5 && aaa<0.6);
	
	mass->rooting_depth = aaa;
}

/* leaf area index *****************************************/
double lai_mass(
	struct Pmas *mass, 
	struct Pchar *pchar
){
	double sla, lai_est;
	
	/** specific leaf area as a function of... what? **/
	sla = pchar->sla;
	
	lai_est = sla*mass->fol*dmTc/100.0/2.0;
	lai_est = (lai_est>=0.0)?lai_est:0.0;
	
	/* dmTc: dry-matter to carbon */
	/* 100.0: cm2 g dm-1 to Mg ha-1 base */
	/* 2.0: single-sided leaf area */	
	
	/*  lai_est = 1.0;  */
	
	return (lai_est);
}

/* give irradiance attenuation coefficient *******************************/
double irr_attn(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	double aaa, bbb;
	
	/* a function of solar hight angle */
	aaa = sin(loct->solhgt_h[loct->doy][loct->hour]*dTr); 
	aaa = (aaa <= 1.0)?aaa:1.0; 
	aaa = (aaa >= 0.3)?aaa:0.3; /* to avoid extreme values*/
	bbb = pchar->eK0/aaa;
	
	return (bbb);
}

/*  intercellular CO2 concentration  *******************************/
void f_incelco2(
	struct Loct *loct,
	struct Pchar *pchar
){
	double ci, gs_co2;
	/* give intercellular CO2 concentration, as a function of ambient CO2 
	level and stomatal conductance */
	
	gs_co2 = pchar->gs/1.56; 
	/* 1.56: conversion from H2O to CO2 conductance */
	
	ci = loct->aCO2-(pchar->ptop/(gs_co2/1000.0));
	/* 1000.0: conbert from mmol to micro mol */
	
	ci = (ci>=0.0)?ci:0.0;
	ci = (ci<=loct->aCO2)?ci:loct->aCO2;
		
	pchar->ci = ci; 
}

/*  quantum yield C3 and c3 plants *********************************/
void f_photo_qy(
	struct Loct *loct,
	struct Pchar *pchar
){
	double eftem, efci;
	
	/* basic data: Osmond et al. (1980) 
	Osmond, C. B., O. Bjorkaman, and D. J. Anderson. 1980. 
	Physiological processes in plant ecology, toward a synthesis 
	with Atriplex. Springer-Verlag. */
	/** plant_type: 3=C3,  4=c3, (5=CAM)  **/
	if(pchar->phototype==3){
		/* temperature dependence */
		eftem = (52.0-loct->tmp_sfc)/(3.5+0.75*(52.0-loct->tmp_sfc)); 
		/* CO2 dependence */
		efci = pchar->ci/(90.0+0.6*pchar->ci); 
		/* 3.5, 52.0, etc.: empirical parameters */
	}else if(pchar->phototype==4){ 
		/* insensitive QE of c3 species */
		eftem = 1.0;
		efci = 1.0;
	}
    
	/* give quantum yield */
	pchar->lue = pchar->lue0*eftem*efci;
}

/* stomatal conductance ************************************/
void f_stom_cond(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	double b1d, cc, aco2;
	
	aco2 = loct->aCO2;
	/* flx stomatal CO2 level: 2011/12/30 by A.Ito */
	if(FIX_GSCO2 == 1){
		aco2 = 350.0;
	}
	
	/* stomatal conductance model by Ball, Woodraw, and Berry (1987) */
	/* modified by Leuning (1995) */
	/*
	Leuning, R. 1995. A critical appraisal of a combined stomatal-
	photosynthesis model for C3 plants. Plant, Cell and Environment 
	18:339-355.
	*/
	b1d = pchar->gs_b1/((aco2 - pchar->cmpcd)*(1.0 + loct->vpd/pchar->gs_b2));
	
	/** add soil water factor **/
	cc = 1.0; /* not defined yet */

	if(pchar->ptop>0.0){
		pchar->gs = pchar->gs_b0+b1d*pchar->ptop*cc;
	}else{
		pchar->gs = pchar->gs_b0; 
	}
}

/* canopy conductance ********************************/
double f_canopy_cond(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double gpp, gc;
	double cc1, cc2, cc3, bb, f_co2, f_vpd;
	
	/* give daily GPP based on the Monsi-Saeki theory, or Kuroiwa's equation */
	if(pchar->psat > 0.0){
		cc1 = 2.0*pchar->psat/pchar->eK;  
		bb = pchar->eK*pchar->lue*pchar->ppfd_t/pchar->psat;
		cc2 = 1.0 + sqrt(1.0 + bb);
		cc3 = 1.0 + sqrt(1.0 + bb*exp(-1.0*pchar->eK*mass->lai));
		gpp = cc1*log(cc2/cc3);

	/*	cc1 = 2.0*pchar->psat*loct->daylen[grid->doy]/pchar->eK; 
		bb = pchar->eK*pchar->lue*loct->par_h[grid->hour]/pchar->psat;
		cc2 = 1.0+sqrt(1.0+bb);
		cc3 = 1.0+sqrt(1.0+bb*exp(-1.0*pchar->eK*mass->lai));
		gpp = cc1*log(cc2/cc3)/24.0;  */
	}else{
		gpp = 0.0;
	}

	f_co2 = 1.0/(loct->aCO2 - pchar->cmpcd);
	/* flx stomatal CO2 level: 2011/12/30 by A.Ito */
	if(FIX_GSCO2 == 1){
		f_co2 = 1.0/(350.0 - pchar->cmpcd);
	}
	
	f_vpd = 1.0/(1.0 + loct->vpd/pchar->gs_b2);
	
	bb = pchar->gs_b1*(f_co2*f_vpd);

	gc = pchar->gs_b0 * mass->lai + bb*gpp;
	
	return (gc);
}

/* optimum LAI by Kuroiwa (1966) **************************/
void f_opt_lai(
	struct Grid *grid,
	struct Loct *loct, 
	struct Pchar *pchar
){
	double aaa, bbb, ccc, cc3, ddd, eee;
	double arm, arg, ar;	
	
	aaa = 1.0/pchar->eK;
	bbb = pchar->eK*pchar->lue*pchar->ppfd_t; 
	
	/* daily respiratory cost */
	eee = log(pchar->qTc)/10.0*(loct->tmp_sfc-15.0);
	arm = pchar->rmf*exp(eee)/1000.0*dmTc*10000.0/(pchar->sla);
	arg = pchar->lf*dmTc*10000.0/(pchar->sla)*(1.0+pchar->rgf);
	ar = arm + arg;
	
	if(pchar->psat>0.0){
		cc3 = (pchar->psat * loct->daylen[loct->doy])/
				(pchar->psat * loct->daylen[loct->doy]-ar*24.0);
		ccc = pchar->psat*(cc3-1.0);

		ddd = bbb/ccc;
		ddd = (ddd>1.0)?ddd:1.0;
	
		pchar->opt_lai = aaa*log(ddd);
	}else{
		pchar->opt_lai = 0.0;
	}
}

/* Q10 of autotrophic respiration **************************/
void f_q10_ar(
	struct Loct *loct, 
	struct Pchar *pchar
){
	double aaa;
	
	/* seasonal acclimation in Q10 ****************/
	/* larger at cool and smaller at warm */
	/* ex. Yokota and Hagihara (1996) Tree Physiology */
	/*
	Yokota, T., and A. Hagihara. 1996. Seasonal change in the temperature 
	coefficient Q10 for respiration of field-grown hinoki cypress 
	(Chamaecyparis obtusa) trees. Journal of Forest Research 1:165-168.
	*/
	aaa = exp(-0.009*(loct->tmp_sfc-15.0)); 
	
	pchar->qTf = pchar->qTf0*aaa; 
	pchar->qTc = pchar->qTc0*aaa;
	pchar->qTr = pchar->qTr0*aaa;
}

/* lowering specific respiration rate of massive biomass with heart 
											wood tissue *******/
void f_spcfc_resp(
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double powstm, powrot;
	double ncon, f_ncon;
	
	pchar->rmf = pchar->rmf0;
	
	/* specific respiration increasing in a power of 2/3 manner */
	/* size-dependence: Yokota and Hagihara (1998) Tree Physiology */
	/*
	Yokota, T., and A. Hagihara. 1998. Changes in the relationship between 
	tree size and aboveground respiration in field-grown hinoki cypress 
	(Chamaecyparis obtusa) trees over three years. Tree Physiology 18:37-43.
	*/
	powstm = 1.0 - 0.33334*mass->stm/(pchar->f_sz_s + mass->stm);
	powrot = 1.0 - 0.33334*mass->rot/(pchar->f_sz_r + mass->rot);
	
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
	
	/*
	ncon = mass->n_storage / ((mass->stm + mass->rot)*1000000.0);
	f_ncon = 0.1 + ncon * 1000.0;
	*/
	f_ncon = 1.0;
	
	pchar->rmc = (pchar->rmc_s*mass->stm_sp + pchar->rmc_h*mass->stm_ht)/
					(mass->stm+0.00001) * f_ncon;
	pchar->rmr = (pchar->rmr_s*mass->rot_fn + pchar->rmr_h*mass->rot_tp)/
					(mass->rot+0.00001) * f_ncon;
}

/* mortarity of plant organs *************************/
void f_mortality(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar
){
	double bbb;
	
	/* larger at warm, and smaller at cool */
	bbb = 1.0; /* +(grid->tmp_sfc-grid->tmp_sfc_am)/100.0; */
	/* !! pending 2005/Feb/21  */
	
	pchar->lf = pchar->lf0*bbb; /* leaf */
	pchar->lc = pchar->lc0*bbb; /* stem */
	pchar->lr = pchar->lr0*bbb; /* root */
}

/* leaf N concentration ****************************/
void f_n_conc(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	/* leaf ********************************/
	/* mmol N m-2 leaf area */
	if(mass->lai > 0.01){
		pchar->n_conc_larea = mass->n_canopy / 10000.0 / 14.0 * 1000.0 / mass->lai;
	}else{
		pchar->n_conc_larea = 1.0;
	}
		
	/* mmol N g-1 leaf weight */
	if(mass->fol > 0.01){
		pchar->n_conc_lmass = mass->n_canopy / 14.0 *1000.0 / (mass->fol * 1000000.0);
	}else{
		pchar->n_conc_lmass = 1.0;
	}
	
	/* storage (g N/ha) ***********************/
	mass->n_struct = mass->stm*1000000.0/250.0 + mass->rot*1000000.0/150.0;
	if(mass->n_struct > (mass->n_storage*0.95)){
		mass->n_struct = mass->n_storage*0.95;
	}
	
	/* capacity of mobile nitrogen stock in stems and roots */
	pchar->n_max_mobile = 0.025 * (mass->stm_sp + mass->rot_fn)*1000000.0;
	if(pchar->n_max_mobile < 1.0){
		pchar->n_max_mobile = 1.0;
	}
}

/* optimal leaf area index estimated by Hikosaka & Anten (2012) */
/*
Hikosaka K & Anten NPR (2012) An evolutionary game of leaf dynamics and 
its consequences for canopy structure. Func.Ecol. 26:1024–1032.
*/
void f_opt_lai_hikosaka_anten(
    struct Grid *grid,
    struct Loct *loct,
    struct Cchar *cchar,
    struct Pmas *mass
){
    long h, i;
    double spp_t, spp_n, lai_t, spp_opt, lai_opt;
    struct Loct *cloct;
    struct Cchar *c_target, *c_neighbor;
    
    /* same properties */
    
    cloct = loct;
    c_target = cchar;
    c_neighbor = cchar;
 
    c_target->photocap_n =  c_neighbor->photocap_n = cchar->photocap_n0 * cchar->vcmx_b;
    c_target->lai =  c_neighbor->lai = cchar->lai;

    /* initialize air conditions */
    c_target->tmp =  c_neighbor->tmp = loct->tmp_2m;
    c_target->vpd_a =  c_neighbor->vpd_a = loct->vpd;
    
    spp_opt = lai_opt = 0.0;
    for(i=0;i<20;i++){
        
        lai_t = 10.0 - 0.5 * (double)i;
        
        c_target->lai = lai_t;
        c_neighbor->lai = lai_t;
    
        spp_t = spp_n = 0.0;
        for(h=0;h<DSTEP;h++){
            cloct->hour = h;
            cloct->hangle = -180.0+(double)h*7.5;
        
            /* canopy-top PPFD */
            c_target->ppfdb_top =  c_neighbor->ppfdb_top = loct->ppfdb_h[h];
            c_target->ppfdd_top =  c_neighbor->ppfdd_top = loct->ppfdd_h[h];
            
            c_target->ci_sn =  c_neighbor->ci_sn = 0.7 * loct->aCO2 / 10.0;
            c_target->ci_sd =  c_neighbor->ci_sd = 0.7 * loct->aCO2 / 10.0;
            
            /* leaf properties */
            f_leaf_prop(grid, cloct, c_target);
            f_leaf_prop(grid, cloct, c_neighbor);
            
            /* sun/shade canopy radiation absorption */
            f_apar_sunshade(grid, cloct, c_target);
            f_apar_sunshade(grid, cloct, c_neighbor);
            
            /* sun/shade carbon assimilation */
            f_gpp_sunshade(1, grid, cloct, c_target);
            f_gpp_sunshade(1, grid, cloct, c_neighbor);
            
            //gpp_t += c_target->photocap_n;
            //gpp_n += c_neighbor->photocap_n;

            spp_t += (c_target->a_sn + c_target->a_sd) * 12.0*1800.0/1000000.0;
            spp_n += (c_neighbor->a_sn + c_neighbor->a_sd) * 12.0*1800.0/1000000.0;
        }
        
        if(spp_t > spp_opt){
            spp_opt = spp_t;
            lai_opt = lai_t;
        }
        //printf("T: %6.1lf %8.2lf  N: %6.1lf %8.2lf\n", lai_t, spp_t, lai_t,spp_n);
    }
    
    //printf("%ld %ld Opt: %6.1lf %8.2lf\n", loct->climy, loct->doy, lai_opt, spp_opt);
}

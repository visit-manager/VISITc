/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/** This file gives description of functions releted to photosynthesis, or GPP **/
/*
De Pury, D.G.G. and Farquhar, G.D., 1997. Simple scaling of photosynthesis from 
leaves to canopies without the errors of big-leaf models. 
Plant, Cell and Environment, 20:537-557.
*/
#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* leaf boundary conductance */
#define GB 250.0

/* leaf biochemical properties *************************************************************/
void f_leaf_prop(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Cchar *cchar
){
	/* canopy CO2 concentration */
	cchar->co2_a = loct->aCO2/10.0;	/* ppmv -> Pa */
	
	/* O2 concentration, 21% of atmospheric pressure, Pa */
	cchar->o2_i = 0.21 * loct->air_prsr*100.0; 
	
	/* CO2 compensation point, Pa */
	/* Brooks & Farquhar (1985) */
	/* Eq.9 in De Pury and Farquhar (1997) */
	cchar->comp_co2 = 3.69+0.188*(cchar->tmp-25.0)+0.0036*(cchar->tmp-25.0)*(cchar->tmp-25.0);
		
	/* Michaelis constant of carboxylation and oxygenation rates */
	/* Eq.8 in De Pury and Farquhar (1997) */
	/* Rubisco carboxylation */
	cchar->kc = cchar->kc0*exp(cchar->acen_kc*(cchar->tmp-25.0)/(298.15*UGC*(cchar->tmp+ZAT)));
	/* Rubisco oxygenation */
	cchar->ko = cchar->ko0*exp(cchar->acen_ko*(cchar->tmp-25.0)/(298.15*UGC*(cchar->tmp+ZAT)));
	/* effective Michaelis constant of Rubisco */
	cchar->k_effc = cchar->kc*(1.0+cchar->o2_i/cchar->ko);
	
	/* temperature coefficient of electron transport */
	/* Eq.10 in De Pury and Farquhar (1997) */
	cchar->ft_et = exp(((cchar->tmp+ZAT)-298.15)*cchar->acen_et/(UGC*(cchar->tmp+ZAT)*298.15))*
		(1.0+exp((cchar->str_et*298.15-cchar->hc_et)/(UGC*298.15)))/
		(1.0+exp((cchar->str_et*(cchar->tmp+ZAT)-cchar->hc_et)/(UGC*(cchar->tmp+ZAT))));		
	
	/* light attenuation coefficient */
	if(loct->sinb_h[loct->doy][loct->hour]>0.0){
		/* beam radiation */
		cchar->ke_b1 = 0.5/loct->sinb_h[loct->doy][loct->hour];
		cchar->ke_b2 = 0.46/loct->sinb_h[loct->doy][loct->hour];
		
		/* cchar->ke_b2*=0.9; */	/* for sensitivity experiment */
	}else{
		/* dark condition */
		cchar->ke_b1 = 0.9;
		cchar->ke_b2 = 0.9;
	}
	
	/* temperature coefficient of dark respiration */
	cchar->ft_rd = exp(66400.0*(cchar->tmp-25.0)/(UGC*298.0*(273.15+cchar->tmp)));
	/* 66400.0:	activation energy */
}

/* absorbed PAR *********************************************************************************/
void f_apar_sunshade(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Cchar *cchar
){
	double irr_b;
	double appfd_sn1, appfd_sn2, appfd_sn3;
	double appfd_sd1, appfd_sd2;
	
	/** sunny and shaded LAI **/
	/* cchar->lai_sn = (1.0-exp(-cchar->ke_b1 * cchar->lai))/loct->sinb_h[grid->doy][grid->hour]; */
	 /* Eq.18 of DF97 */
	cchar->lai_sn = (1.0-exp(-cchar->ke_b1 * cchar->lai))/cchar->ke_b1;
	if(cchar->lai_sn > cchar->lai){
		cchar->lai_sn = cchar->lai;			/* sunny */
	}
	 /* Eq.21 of DF97 */
	cchar->lai_sd = cchar->lai - cchar->lai_sn;	/* shade */
	
	/* daytime **************************************************/
	if(loct->sinb_h[loct->doy][loct->hour]>0.0){
		
		/* reflection */
		irr_b = (1.0-sqrt(1.0-cchar->scttr))/(1.0+sqrt(1.0-cchar->scttr));
		cchar->rfl_b = 1.0-exp(-2.0*irr_b*cchar->ke_b1/(1.0+cchar->ke_b1));
		
		/* absorbed PPFD by canopy, total */
		 /* Eq.13 of DF97 */
		cchar->appfd = (1.0-cchar->rfl_b)*cchar->ppfdb_top*(1.0-exp(-cchar->ke_b2*cchar->lai))+
					(1.0-cchar->rfl_d)*cchar->ppfdd_top*(1.0-exp(-cchar->ke_d*cchar->lai));
		
		/* absorbed PPFD by sunny leaves */
		/* beam: Eq.20b of DF97 */
		appfd_sn1 = cchar->ppfdb_top*(1.0-cchar->scttr)*(1.0-exp(-cchar->ke_b1*cchar->lai));
		/* diffused: Eq.20c of DF97 */
		appfd_sn2 = cchar->ppfdd_top*(1.0-cchar->rfl_d)*(1.0-exp(-(cchar->ke_d+cchar->ke_b1)*cchar->lai))*
			cchar->ke_d/(cchar->ke_d+cchar->ke_b1);
		/* beam scattered within canopy: Eq.20d of DF97 */
		appfd_sn3 = cchar->ppfdb_top*((1.0-cchar->rfl_b)*(1.0-exp(-(cchar->ke_b1+cchar->ke_b2)*cchar->lai))*
			cchar->ke_b2/(cchar->ke_b1+cchar->ke_b2)-(1.0-cchar->scttr)*(1.0-exp(-2.0*cchar->ke_b1*cchar->lai))/2.0);
		/* absorbed PPFD by sunny leaves: Eq.20a of DF97 */
		cchar->appfd_sn = appfd_sn1+appfd_sn2+appfd_sn3;
		/* absorbed PPFD by sunny leaves PS II: Eq.6 of DF97 */
		/* see p.35 of Von Caemmerer (2000) */
		cchar->apar_sn = cchar->appfd_sn*(1.0-cchar->spect)/2.0;
		/* cchar->spect : spectral coefficinet for chlorophyll absorption */
		/* 2.0 : half separation between PS I and PS II */
		
		/* absorbed PPFD by shaded leaves */
		/* diffused: Eq.A26b of DF97 */
		appfd_sd1 = cchar->ppfdd_top*(1.0-cchar->rfl_d)*(1.0-exp(-cchar->ke_d*cchar->lai)-
			(1.0-exp(-(cchar->ke_d+cchar->ke_b1)*cchar->lai))*cchar->ke_d/(cchar->ke_d+cchar->ke_b1));
		/* beam scattered within canopy: Eq.A26c of DF97 */
		appfd_sd2 = cchar->ppfdb_top*((1.0-cchar->rfl_b)*(1.0-exp(-cchar->ke_b2*cchar->lai)-
			(1.0-exp(-(cchar->ke_b1+cchar->ke_b2)*cchar->lai))*cchar->ke_b2/(cchar->ke_b1+cchar->ke_b2))-
			(1.0-cchar->scttr)*(1.0-exp(-cchar->ke_b1*cchar->lai)-(1.0-exp(-2.0*cchar->ke_b1*cchar->lai))/2.0));
		/* absorbed PPFD by shaded leaves: Eq.A26a of DF97 */
		cchar->appfd_sd = appfd_sd1+appfd_sd2;
		/* absorbed PPFD by shaded leaves PS II */
		cchar->apar_sd = cchar->appfd_sd*(1.0-cchar->spect)/2.0; /* Eq.6 of DF97 */
		
		/* canopy-absorbed PAR */
		cchar->apar = cchar->apar_sn+cchar->apar_sd;
		
		cchar->ppfdb_btm = cchar->ppfdb_top-appfd_sn1-appfd_sn3-appfd_sd2;
		cchar->ppfdd_btm = cchar->ppfdd_top-appfd_sn2-appfd_sd1;
	}else{
		cchar->appfd = 0.0;
		cchar->appfd_sn = 0.0;
		cchar->appfd_sd = 0.0;
		cchar->apar = 0.0;
		cchar->apar_sn = 0.0;
		cchar->apar_sd = 0.0;
		cchar->ppfdb_btm = 0.0;
		cchar->ppfdd_btm = 0.0;
	}
}

/* canopy GPP by sun/shade model **********************************************/
void f_gpp_sunshade(
	long phenol_stage, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Cchar *cchar
){
	long nn;
	double crit, step, kn, nadj, nc_con, nc_ex, n_top;
	
	step = 0.1;
    
    if(EX_N_ECOPHYS==1){
        kn = 0.4;
        nadj = 1.0;
    }else if(EX_N_ECOPHYS==2){
        
        kn = 0.68;
   
        nc_con = cchar->lai * ((cchar->nc_ct-cchar->nc_np) * exp(-cchar->ke_n) / cchar->ke_n + cchar->nc_np);
        nc_ex = cchar->lai * ((cchar->nc_ct-cchar->nc_np) * exp(-kn) / kn + cchar->nc_np);
        
        if(nc_ex > 0.0){
            nadj = nc_con / nc_ex;
        }else{
            nadj = 0.0;
        }
        
    }else{
        kn = cchar->ke_n;
        nadj = 1.0;
    }
    
    n_top = nadj * cchar->nc_ct;
	
	/****** sunny ******/
	/* maximum carboxylation rate at 25 deg-C: Eq.22 in DF97 */
	cchar->vcmx_sn0 = cchar->lai*cchar->photocap_n*(n_top - cchar->nc_np)*
		(1.0-exp(-(kn + cchar->ke_b1*cchar->lai)))/(kn + cchar->ke_b1*cchar->lai);
	/* maximum carboxylation rate: Eq.8 in DF97 */
	cchar->vcmx_sn = cchar->vcmx_sn0*exp(cchar->acen_vc*(cchar->tmp-25.0)/(298.15*UGC*(cchar->tmp+ZAT)));
    
	/* maximum electron transport rate: proportional Jmax/Vcmax relationship */
	cchar->jmx_sn = cchar->jvr*cchar->vcmx_sn0*cchar->ft_et;

    if(G_MESOPHYL==1){
        cchar->vcmx_sn = 0.95*cchar->vcmx_sn + 13.87;
        cchar->jmx_sn = 1.03*cchar->jmx_sn + 15.46;
    }
    
    /* stabilization */
	if(cchar->apar_sn>0.0 && loct->sinb_h[loct->doy][loct->hour]>0.0 && phenol_stage !=0){
		/* electron transport rate: Eq.5 in DF97 */
		cchar->j_sn = ((cchar->apar_sn+cchar->jmx_sn)-sqrt((cchar->apar_sn+cchar->jmx_sn)*
			(cchar->apar_sn+cchar->jmx_sn)-4.0*cchar->cv_j*cchar->apar_sn*cchar->jmx_sn))/(2.0*cchar->cv_j);
		
		if(SOLV_LEAF==0){
			/** stabilization with Discrete Newton method **/
			nn = 0;
			crit = 10.0;
			while(crit>0.1&&nn<20){
				cchar->ci_sn = cchar->ci_sn - f_gas_sun(cchar, cchar->ci_sn)*step/
						(f_gas_sun(cchar, cchar->ci_sn+step) - f_gas_sun(cchar, cchar->ci_sn));
				crit = fabs(f_gas_sun(cchar, cchar->ci_sn));
				nn++;
				
				if(NOTICE==1 && cchar->ci_sn<0.0){
					printf("FATAL ERROR: negative Ci of sun-leaf: %lf\n", cchar->ci_sn);
				}
			} 
		}else if(SOLV_LEAF==1){
			/** analytical solution by Baldocchi (1994) **/
			leaf_ansolv(0, cchar);
		}
	}else{
		/* dark condition */
		cchar->j_sn = 0.0;
		cchar->pc_sn = 0.0;
		cchar->pj_sn = 0.0;
		cchar->p_sn = 0.0;
		cchar->ci_sn = cchar->co2_a;
		
		f_gas_sun(cchar, cchar->ci_sn);
	}

	/****** shade ******/
	/* maximum carboxylation rate at 25 deg-C: Eq.A27 in DF97 */
	/* Vcmax(shade) = Vcmax - Vcmax(sun) */
	cchar->vcmx_sd0 = cchar->lai*cchar->photocap_n*(n_top - cchar->nc_np)*((1.0-exp(-kn))/
		kn - (1.0-exp(-(kn + cchar->ke_b1*cchar->lai)))/(kn + cchar->ke_b1*cchar->lai));
	/* maximum carboxylation rate: Eq.8 in DF97 */
	cchar->vcmx_sd = cchar->vcmx_sd0*exp(cchar->acen_vc*(cchar->tmp-25.0)/(298.15*UGC*(cchar->tmp+ZAT)));
	/* maximum electron transport rate */
	cchar->jmx_sd = (cchar->jvr*cchar->vcmx_sd0)*cchar->ft_et; 

    if(G_MESOPHYL==1){
        cchar->vcmx_sd = 0.95*cchar->vcmx_sd + 13.87;
        cchar->jmx_sd = 1.03*cchar->jmx_sd + 15.46;
    }

    /* stabilization */
	if(cchar->apar_sd>0.0 && loct->sinb_h[loct->doy][loct->hour]>0.0){
		/* electron transport rate: Eq.5 in DF97 */
		cchar->j_sd = ((cchar->apar_sd+cchar->jmx_sd)-sqrt((cchar->apar_sd+cchar->jmx_sd)*(cchar->apar_sd+
			cchar->jmx_sd)-4.0*cchar->cv_j*cchar->apar_sd*cchar->jmx_sd))/(2.0*cchar->cv_j);
		
		if(SOLV_LEAF==0){
			/** stabilization with Discrete Newton methid **/
			nn=0;
			crit=10.0;
			while(crit>0.1 && nn<20){
				cchar->ci_sd = cchar->ci_sd - f_gas_shade(cchar, cchar->ci_sd)*step/
						(f_gas_shade(cchar, cchar->ci_sd+step) - f_gas_shade(cchar, cchar->ci_sd));
				crit = fabs(f_gas_shade(cchar, cchar->ci_sd));
				nn++;

				if(NOTICE==1 && cchar->ci_sd<0.0){
					printf("FATAL ERROR: negative Ci of shade-leaf: %lf\n", cchar->ci_sd);
				}
			} 
		}else if(SOLV_LEAF==1){
			/** analytical solution by Baldocchi (1994) **/
			leaf_ansolv(1, cchar); 
		}
	}else{
		/* dark condition */
		cchar->j_sd = 0.0;
		cchar->pc_sd = 0.0;
		cchar->pj_sd = 0.0;
		cchar->p_sd = 0.0;
		cchar->ci_sd = cchar->co2_a;
		
		f_gas_shade(cchar, cchar->ci_sd);
	}
}

/* objective function of Discrete Newton method, sun-leaf ***************************/
double f_gas_sun(
	struct Cchar *cchar,
	double ci
){	
	double fgc, gb, cs, p_co2;
    
	/* leaf-boundary conductance of sun leaf */
	if(cchar->lai_sn>1.0){
		fgc = cchar->lai_sn;
	}else{
		fgc = 1.0;
	}
	gb = GB*fgc/100.0;
	cchar->gb = gb;
    
    /* dark respiration */
	cchar->rd_sn = cchar->rd0*cchar->ft_rd*fgc;
    
	/* CO2 assimilation rate */
	if(ci > cchar->comp_co2){
		/* Rubisco-limited */
		cchar->pc_sn = cchar->vcmx_sn*(ci-cchar->comp_co2)/(ci+cchar->k_effc);
		/* RuBP-regeneration-limited */
		cchar->pj_sn = cchar->j_sn*(ci-cchar->comp_co2)/(4.0*(ci+2.0*cchar->comp_co2));
	}else{
		cchar->pc_sn = 0.0;
		cchar->pj_sn = 0.0;
	}
	
	/* combination of limitations ***********/
	cchar->p_sn = (cchar->pc_sn<cchar->pj_sn)?cchar->pc_sn:cchar->pj_sn; /*	*/
	/* cchar->p_sn = ((cchar->pc_sn+cchar->pj_sn)-sqrt((cchar->pc_sn+cchar->pj_sn)*(cchar->pc_sn+cchar->pj_sn)-
		4.0*cchar->cv_p*cchar->pc_sn*cchar->pj_sn))/(2.0*cchar->cv_p); */
	/* cv_p, convexity: 0.877 () */
	
	if(NOTICE==1 && cchar->p_sn<0.0){
		printf("FATAL ERROR: negative P of sun-leaf: %lf\n", cchar->p_sn);
	}
	
	cchar->a_sn = cchar->p_sn - cchar->rd_sn;
	
	cs = cchar->co2_a - cchar->a_sn/gb;
	/* flx stomatal CO2 level: 2011/12/30 by A.Ito */
	if(FIX_GSCO2 == 1){
		cs = 35.0 - cchar->a_sn/gb;
	}
	
    /* stomatal conductance */
	if(cchar->a_sn>0.0){
		cchar->gs_sn = cchar->gss_b0/100.0 + cchar->gss_b1/100.0/10.0*cchar->a_sn/
					(cs*(1.0+cchar->vpd_a/cchar->gss_b2));
		cchar->gs_sn *= fgc;
	}else{
		cchar->gs_sn = cchar->gss_b0/100.0*fgc;
	}
	
    /* mesophyll conductance: 2012/01/21 (A.Ito) */
    cchar->g_msp = cchar->gs_sn; /* !!!! assumption */
	
	if(NOTICE==1 && cchar->gs_sn<0.0){
		/* printf("FARTAL ERROR: negative gs of sun-leaf: %lf %lf\n", cchar->a_sn, cchar->gs_sn); */
	}
    
    if(G_MESOPHYL == 1){
        /* considering mesophyll diffusion */
        return (ci + cchar->a_sn/gb + cchar->a_sn/cchar->gs_sn + cchar->a_sn/cchar->g_msp - cchar->co2_a);
    }else{
        /* default */
        return (ci + cchar->a_sn/gb + cchar->a_sn/cchar->gs_sn - cchar->co2_a);
    }
}

/**** objective function of Discrete Newton method, shade-leaf ********************************/
double f_gas_shade(
	struct Cchar *cchar,
	double ci
){	
	double fgc, gb, cs;

	/* leaf-boundary conductance of shade leaf */
	if(cchar->lai_sd>1.0){
		fgc = cchar->lai_sd;
	}else{
		fgc = 1.0;
	}
	gb = GB*fgc/100.0;
	cchar->gb = gb;
    
    /* dark respiration */
	cchar->rd_sd = cchar->rd0*cchar->ft_rd*fgc;
	
	/* CO2 assimilation rate */
	if(ci > cchar->comp_co2){
		/* Rubisco-limited */
		cchar->pc_sd = cchar->vcmx_sd*(ci-cchar->comp_co2)/(ci+cchar->k_effc);
		/* RuBP-regeneration-limited */
		cchar->pj_sd = cchar->j_sd*(ci-cchar->comp_co2)/(4.0*(ci+2.0*cchar->comp_co2));
	}else{
		cchar->pc_sd = 0.0;
		cchar->pj_sd = 0.0;
	}
	/* combination */
	cchar->p_sd = (cchar->pc_sd<cchar->pj_sd)?cchar->pc_sd:cchar->pj_sd; /*	*/
	/* cchar->p_sd = ((cchar->pc_sd+cchar->pj_sd)-sqrt((cchar->pc_sd+cchar->pj_sd)*(cchar->pc_sd+cchar->pj_sd)-
		4.0*cchar->cv_p*cchar->pc_sd*cchar->pj_sd))/(2.0*cchar->cv_p); */
		
	if(NOTICE==1 && cchar->p_sd<0.0){
		printf("FATAL ERROR: negative P of shade-leaf: %lf %lf %lf %lf\n", 
			ci, cchar->p_sd, cchar->pc_sd, cchar->pj_sd);
	}
	
	/* leaf-level net assimilation rate */
	cchar->a_sd = cchar->p_sd - cchar->rd_sd;
	
	/* sub stomatal CO2 partial pressure */
	cs = cchar->co2_a - cchar->a_sd/gb;
	/* flx stomatal CO2 level: 2011/12/30 by A.Ito */
	if(FIX_GSCO2 == 1){
		cs = 35.0 - cchar->a_sd/gb;
	}
	
    /* stomatal conductance */
	if(cchar->a_sd>0.0){
		cchar->gs_sd = cchar->gss_b0/100.0 + cchar->gss_b1/100.0/10.0*cchar->a_sd/
				(cs*(1.0+cchar->vpd_a/cchar->gss_b2));
		/* cchar->gs_sd = cchar->gss_b0 + cchar->gss_b1*cchar->a_sd/
				((ci-cchar->comp_co2)*10.0*(1.0+cchar->vpd_a/cchar->gss_b2)); */
		cchar->gs_sd *= fgc;
	}else{
		cchar->gs_sd = cchar->gss_b0*fgc;
	}

    /* mesophyll conductance: 2012/01/21 (A.Ito) */
    cchar->g_msp = cchar->gs_sd; /* !!!! assumption */
	
	if(NOTICE==1 && cchar->gs_sd<0.0){
		/* printf("FARTAL ERROR: negative gs of shade-leaf: %lf %lf\n", cchar->a_sd, cchar->gs_sd); */
	}
	
    if(G_MESOPHYL == 1){
        /* considering mesophyll diffusion */
        return (ci + cchar->a_sd/gb + cchar->a_sd/cchar->gs_sd + cchar->a_sd/cchar->g_msp - cchar->co2_a);
    }else{
        /* default */
        return (ci + cchar->a_sd/gb + cchar->a_sd/cchar->gs_sd - cchar->co2_a);
    }
}

/********************************************************/
 void leaf_ansolv(
 	long frac, 
 	struct Cchar *cchar
 ){
 	double aa, bb, dd, ee, fgc;
 	double ca, cs, gb, bd, rh, m, cs_gs;
 	double alpha, beta, gamma, theta;
 	double pp, qq, rr, qqq, rrr;
 	double xx1, xx2, xx3;
 	double f_pp, f_qq, f_cube;
 	 	
 	/**************************************************************************************/
 	gb = GB;						/* leaf boundary layer conductance*/
 	bd = cchar->gss_b0;				/* minimum conductance */
 	rh = 1.0/(1.0+cchar->vpd_a/cchar->gss_b2);	/* relative humidity */
 	m = cchar->gss_b1/10.0;			/* coefficient of slope */
 	ca = cchar->co2_a;				/* atmospheric CO2, Pa */
 	
 	if(frac==0){ 
		 /* for sun leaf */
		if(cchar->lai_sn>1.0){
			fgc = cchar->lai_sn;
		}else{
			fgc = 1.0;
		}
		cchar->rd_sn = cchar->rd0*cchar->ft_rd*fgc;
	}else if(frac==1){
	 	/* for shade leaf */
		if(cchar->lai_sd>1.0){
			fgc = cchar->lai_sd;
		}else{
			fgc = 1.0;
		}
		cchar->rd_sd = cchar->rd0*cchar->ft_rd*fgc;
	}
	gb *= fgc/100.0;
	bd *= fgc/100.0;
	m *= fgc/100.0;
	
	cchar->gb = gb;
 	
 	/* Rubisco-limited ************************************************************************/
 	if(frac==0){
		aa = cchar->vcmx_sn;
	}else if(frac==1){
		aa = cchar->vcmx_sd;
	}
	bb = cchar->k_effc;
	dd = cchar->comp_co2;
	ee = 1.0;
	
	/******************************************************************/
 	alpha = 1.0+bd/gb - m*rh;
 	beta = ca * (gb * m * rh - 2.0*bd - gb);
 	gamma = ca*ca * bd * gb;
 	theta = gb * m * rh - bd;
 	
 	/******************************************************************/
 	if(frac==0){
		pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*cchar->rd_sn)/(ee * alpha);
	 	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*cchar->rd_sn*beta 
	 		+ cchar->rd_sn*bb*theta)/(ee * alpha);
	 	rr = (-aa*gamma + aa*dd*gamma/ca + ee*cchar->rd_sn*gamma 
	 		+ cchar->rd_sn*bb*gamma/ca)/(ee * alpha);
 	}else if(frac==1){
		pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*cchar->rd_sd)/(ee * alpha);
	 	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*cchar->rd_sd*beta 
	 		+ cchar->rd_sd*bb*theta)/(ee * alpha);
	 	rr = (-aa*gamma + aa*dd*gamma/ca + ee*cchar->rd_sd*gamma 
	 		+ cchar->rd_sd*bb*gamma/ca)/(ee * alpha);
 	}
 	f_pp = (pp/3.0)*(pp/3.0) - qq/3.0;
 	f_qq = ((pp/3.0)*qq - 2.0*(pp/3.0)*(pp/3.0)*(pp/3.0) - rr)/2.0;
 	f_cube = f_qq*f_qq - f_pp*f_pp*f_pp;
 	
 	if(f_cube<0.0){
	 	/******************************************************************/
		qqq = (pp*pp - 3.0*qq)/9.0;
		rrr = (2.0*pp*pp*pp - 9.0*pp*qq + 27.0*rr)/54.0;
		theta = acos(rrr/sqrt(qqq*qqq*qqq));
		 	
	 	xx1 = -2.0*sqrt(qqq) * cos(theta/3.0) - pp/3.0;
	 	xx2 = -2.0*sqrt(qqq) * cos((theta+2.0*PI)/3.0) - pp/3.0;
	 	xx3 = -2.0*sqrt(qqq) * cos((theta+4.0*PI)/3.0) - pp/3.0; 	/**/
	 		 	
	 	if(frac==0){
	 		cchar->pc_sn = xx3 + cchar->rd_sn;
	 	}else if(frac==1){
	 		cchar->pc_sd = xx3 + cchar->rd_sd;
	 	}
	 }else{
	 	if(frac==0){
	 		cchar->pc_sn = 0.0;
	 	}else if(frac==1){
	 		cchar->pc_sd = 0.0;
	 	}
	 }

 	/* RUBP-limited ************************************************************************/
 	if(frac==0){
 		aa = cchar->j_sn;
 	}else if(frac==1){
 		aa = cchar->j_sd;
 	}
 	bb = 8.0*cchar->comp_co2;
 	dd = cchar->comp_co2;
 	ee = 4.0;
 	
	/******************************************************************/
 	alpha = 1.0+bd/gb - m*rh;
 	beta = ca * (gb * m * rh - 2.0*bd - gb);
 	gamma = ca*ca * bd * gb;
 	theta = gb * m * rh - bd;
 	
 	/******************************************************************/
 	if(frac==0){
		pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*cchar->rd_sn)/(ee * alpha);
	 	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*cchar->rd_sn*beta 
	 		+ cchar->rd_sn*bb*theta)/(ee * alpha);
	 	rr = (-aa*gamma + aa*dd*gamma/ca + ee*cchar->rd_sn*gamma 
	 		+ cchar->rd_sn*bb*gamma/ca)/(ee * alpha);
 	}else if(frac==1){
		pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*cchar->rd_sd)/(ee * alpha);
	 	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*cchar->rd_sd*beta 
	 		+ cchar->rd_sd*bb*theta)/(ee * alpha);
	 	rr = (-aa*gamma + aa*dd*gamma/ca + ee*cchar->rd_sd*gamma 
	 		+ cchar->rd_sd*bb*gamma/ca)/(ee * alpha);
 	}
 	
 	f_pp = (pp/3.0)*(pp/3.0) - qq/3.0;
 	f_qq = ((pp/3.0)*qq - 2.0*(pp/3.0)*(pp/3.0)*(pp/3.0) - rr)/2.0;
 	f_cube = f_qq*f_qq - f_pp*f_pp*f_pp;
 	
 	if(f_cube<0.0){
		/******************************************************************/
		qqq = (pp*pp - 3.0*qq)/9.0;
		rrr = (2.0*pp*pp*pp - 9.0*pp*qq + 27.0*rr)/54.0;
		theta = acos(rrr/sqrt(qqq*qqq*qqq));
		 	
	 	/******************************************************************/
	 	xx1 = -2.0*sqrt(qqq) * cos(theta/3.0) - pp/3.0;
	 	xx2 = -2.0*sqrt(qqq) * cos((theta+2.0*PI)/3.0) - pp/3.0;
	 	xx3 = -2.0*sqrt(qqq) * cos((theta+4.0*PI)/3.0) - pp/3.0;
	 	
	 	if(frac==0){
	 		cchar->pj_sn = xx3 + cchar->rd_sn;
	 	}else if(frac==1){
	 		cchar->pj_sd = xx3 + cchar->rd_sd;
	 	}
	 }else{
	 	if(frac==0){
	 		cchar->pj_sn = 0.0;
	 	}else if(frac==1){
	 		cchar->pj_sd = 0.0;
	 	}
	 }

 	/* Solution ************************************************************************/
 	if(frac==0){ /* sunny */
 		/* cchar->p_sn = ((cchar->pc_sn+cchar->pj_sn)-sqrt((cchar->pc_sn+cchar->pj_sn)*(cchar->pc_sn+cchar->pj_sn)-
					4.0*cchar->cv_p*cchar->pc_sn*cchar->pj_sn))/(2.0*cchar->cv_p); */
					
		cchar->p_sn = (cchar->pc_sn<cchar->pj_sn)?cchar->pc_sn:cchar->pj_sn;
		
		if(cchar->p_sn>=0.0 && cchar->p_sn<=100.0){
			;
		}else{
			/* printf("!! Warning: [p_sn, %lf] in [leaf_ansolv()] !!\n", cchar->p_sn); */
			/* printf("!! Warning: [p_sn, %lf %lf %lf] in [leaf_ansolv()] !!\n", cchar->p_sn, cchar->pc_sn, cchar->pj_sn); */
		}
		
		cchar->a_sn = cchar->p_sn - cchar->rd_sn;
		
		/* substomatal cavity CO2 */
		cs = cs_gs = ca - cchar->a_sn/gb;
		/* flx stomatal CO2 level: 2011/12/30 by A.Ito */
		if(FIX_GSCO2 == 1){
			cs_gs = 35.0 - cchar->a_sn/gb;
		}
		
		cchar->gs_sn = bd + m*cchar->a_sn*rh/cs_gs;		
		cchar->ci_sn = cs - cchar->a_sn/cchar->gs_sn;
	}else if(frac==1){ /* shade */
		/* cchar->p_sd = ((cchar->pc_sd+cchar->pj_sd)-sqrt((cchar->pc_sd+cchar->pj_sd)*(cchar->pc_sd+cchar->pj_sd)-
					4.0*cchar->cv_p*cchar->pc_sd*cchar->pj_sd))/(2.0*cchar->cv_p); */
					
		cchar->p_sd = (cchar->pc_sd<cchar->pj_sd)?cchar->pc_sd:cchar->pj_sd;
		
		if(cchar->p_sd>=0.0 && cchar->p_sd<=100.0){
			;
		}else{
			/* printf("!! Warning: [p_sd, %lf] in [leaf_ansolv()] !!\n", cchar->p_sd); */
			/* printf("!! Warning: [p_sd, %lf %lf %lf] in [leaf_ansolv()] !!\n", cchar->p_sd, cchar->pc_sd, cchar->pj_sd); */
		}
		
		cchar->a_sd = cchar->p_sd - cchar->rd_sd;
		
		/* substomatal cavity CO2 */
 		cs = cs_gs = ca - cchar->a_sd/gb;
		/* flx stomatal CO2 level: 2011/12/30 by A.Ito */
		if(FIX_GSCO2 == 1){
			cs_gs = 35.0 - cchar->a_sd/gb;
		}
		
		cchar->gs_sd = bd + m*cchar->a_sd*rh/cs_gs;		
		cchar->ci_sd = cs - cchar->a_sd/cchar->gs_sd;
	}
}

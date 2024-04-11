/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
*/

/** This file gives description of functions releted to photosynthesis, or GPP **/
/*
De Pury, D.G.G. and Farquhar, G.D., 1997. Simple scaling of photosynthesis from 
leaves to canopies without the errors of big-leaf models. 
Plant, Cell and Environment, 20:537-557.
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* leaf boundary conductance */
#define GB 250.0

/* leaf biochemical properties *****************************************************************/
void f_leaf_prop(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	/* canopy CO2 concentration */
	pchar->co2_a = loct->aCO2/10.0;	/* ppmv -> Pa */
	
	/* O2 concentration, 21% of atmospheric pressure, Pa */
	pchar->o2_i = 0.21 * loct->air_prsr * 100.0;
	
	/* CO2 compensation point, Pa */
	/* Brooks & Farquhar (1985) */
	/* Eq.9 in De Pury and Farquhar (1997) */
	pchar->comp_co2 = 3.69 + 0.188*(pchar->tmp-25.0)+0.0036*(pchar->tmp - 25.0)*(pchar->tmp - 25.0);
		
	/* Michaelis constant of carboxylation and oxygenation rates */
	/* Eq.8 in De Pury and Farquhar (1997) */
	/* Rubisco carboxylation */
	pchar->kc = pchar->kc0 * exp(pchar->acen_kc*(pchar->tmp-25.0)/(298.15*UGC*(pchar->tmp + ZAT)));
	/* Rubisco oxygenation */
	pchar->ko = pchar->ko0 * exp(pchar->acen_ko*(pchar->tmp-25.0)/(298.15*UGC*(pchar->tmp + ZAT)));
	/* effective Michaelis constant of Rubisco */
	pchar->k_effc = pchar->kc * (1.0 + pchar->o2_i/pchar->ko);
	
	/* temperature coefficient of electron transport */
	/* Eq.10 in De Pury and Farquhar (1997) */
	pchar->ft_et = exp(((pchar->tmp + ZAT)-298.15)*pchar->acen_et/(UGC*(pchar->tmp + ZAT)*298.15))*
		(1.0 + exp((pchar->str_et*298.15-pchar->hc_et)/(UGC*298.15)))/
		(1.0 + exp((pchar->str_et*(pchar->tmp + ZAT)-pchar->hc_et)/(UGC*(pchar->tmp + ZAT))));		
	
	/* light attenuation coefficient */
	if(loct->sinb_h > 0.0){
		/* beam radiation */
		pchar->ke_b1 = 0.5/loct->sinb_h;
		pchar->ke_b2 = 0.46/loct->sinb_h;
		
		/* pchar->ke_b2*=0.9; */	/* for sensitivity experiment */
	}else{
		/* dark condition */
		pchar->ke_b1 = 0.9;
		pchar->ke_b2 = 0.9;
	}
	
	/* temperature coefficient of dark respiration */
	pchar->ft_rd = exp(66400.0*(pchar->tmp-25.0)/(UGC*298.0*(273.15 + pchar->tmp)));
	/* 66400.0:	activation energy */
}

/* absorbed PAR *********************************************************************************/
void f_apar_sunshade(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	float irr_b;
	float appfd_sn1, appfd_sn2, appfd_sn3;
	float appfd_sd1, appfd_sd2;
	
	/** sunny and shaded LAI **/
	/* pchar->lai_sn = (1.0-exp(-pchar->ke_b1 * pchar->lai))/loct->sinb_h[grid->doy][grid->hour]; */
	 /* Eq.18 of DF97 */
	pchar->lai_sn = (1.0-exp(-pchar->ke_b1 * pchar->lai))/pchar->ke_b1;
	if(pchar->lai_sn > pchar->lai){
		pchar->lai_sn = pchar->lai;			/* sunny */
	}
	 /* Eq.21 of DF97 */
	pchar->lai_sd = pchar->lai - pchar->lai_sn;	/* shade */
	
	/* daytime **************************************************/
	if(loct->sinb_h > 0.0){
		
		/* reflection */
		irr_b = (1.0-sqrt(1.0-pchar->scttr))/(1.0+sqrt(1.0 - pchar->scttr));
		pchar->rfl_b = 1.0-exp(-2.0*irr_b*pchar->ke_b1/(1.0 + pchar->ke_b1));
		
		/* absorbed PPFD by canopy, total */
		 /* Eq.13 of DF97 */
		pchar->appfd = (1.0-pchar->rfl_b)*pchar->ppfdb_top*(1.0-exp(-pchar->ke_b2*pchar->lai))+
					(1.0-pchar->rfl_d)*pchar->ppfdd_top*(1.0-exp(-pchar->ke_d*pchar->lai));
		
		/* absorbed PPFD by sunny leaves */
		/* beam: Eq.20b of DF97 */
		appfd_sn1 = pchar->ppfdb_top*(1.0-pchar->scttr)*(1.0-exp(-pchar->ke_b1*pchar->lai));
		/* diffused: Eq.20c of DF97 */
		appfd_sn2 = pchar->ppfdd_top*(1.0-pchar->rfl_d)*(1.0-exp(-(pchar->ke_d+pchar->ke_b1)*pchar->lai))*
			pchar->ke_d/(pchar->ke_d+pchar->ke_b1);
		/* beam scattered within canopy: Eq.20d of DF97 */
		appfd_sn3 = pchar->ppfdb_top*((1.0-pchar->rfl_b)*(1.0-exp(-(pchar->ke_b1+pchar->ke_b2)*pchar->lai))*
			pchar->ke_b2/(pchar->ke_b1+pchar->ke_b2)-(1.0-pchar->scttr)*(1.0-exp(-2.0*pchar->ke_b1*pchar->lai))/2.0);
		/* absorbed PPFD by sunny leaves: Eq.20a of DF97 */
		pchar->appfd_sn = appfd_sn1 + appfd_sn2 + appfd_sn3;
		/* absorbed PPFD by sunny leaves PS II: Eq.6 of DF97 */
		/* see p.35 of Von Caemmerer (2000) */
		pchar->apar_sn = pchar->appfd_sn*(1.0 - pchar->spect)/2.0;
		/* pchar->spect : spectral coefficinet for chlorophyll absorption */
		/* 2.0 : half separation between PS I and PS II */
		
		/* absorbed PPFD by shaded leaves */
		/* diffused: Eq.A26b of DF97 */
		appfd_sd1 = pchar->ppfdd_top*(1.0-pchar->rfl_d)*(1.0 - exp(-pchar->ke_d*pchar->lai)-
			(1.0-exp(-(pchar->ke_d + pchar->ke_b1)*pchar->lai))*pchar->ke_d/(pchar->ke_d+pchar->ke_b1));
		/* beam scattered within canopy: Eq.A26c of DF97 */
		appfd_sd2 = pchar->ppfdb_top*((1.0-pchar->rfl_b)*(1.0 - exp(-pchar->ke_b2*pchar->lai)-
			(1.0-exp(-(pchar->ke_b1 + pchar->ke_b2)*pchar->lai))*pchar->ke_b2/(pchar->ke_b1+pchar->ke_b2))-
			(1.0-pchar->scttr)*(1.0-exp(-pchar->ke_b1*pchar->lai)-(1.0-exp(-2.0*pchar->ke_b1*pchar->lai))/2.0));
		/* absorbed PPFD by shaded leaves: Eq.A26a of DF97 */
		pchar->appfd_sd = appfd_sd1 + appfd_sd2;
		/* absorbed PPFD by shaded leaves PS II */
		pchar->apar_sd = pchar->appfd_sd*(1.0 - pchar->spect)/2.0; /* Eq.6 of DF97 */
		
		/* canopy-absorbed PAR */
		pchar->apar = pchar->apar_sn + pchar->apar_sd;
		
		pchar->ppfdb_btm = pchar->ppfdb_top - appfd_sn1 - appfd_sn3 - appfd_sd2;
		pchar->ppfdd_btm = pchar->ppfdd_top - appfd_sn2 - appfd_sd1;
	}else{
		pchar->appfd = 0.0;
		pchar->appfd_sn = 0.0;
		pchar->appfd_sd = 0.0;
		pchar->apar = 0.0;
		pchar->apar_sn = 0.0;
		pchar->apar_sd = 0.0;
		pchar->ppfdb_btm = 0.0;
		pchar->ppfdd_btm = 0.0;
	}
}

/* canopy GPP by sun/shade model **********************************************/
void f_gpp_sunshade(
	long phenol_stage, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	long nn;
	float crit, step, ksw, nc0, nc_ct1, kn, kb, kl, lat;
	
	step = 0.1;
	
    /* soil water effect on photosynthetic capacity: 2013/05/28 by A.Ito */
    ksw = grid->fieldcap * pchar->km_nstl;
    if(ksw > 120.0){
        ksw = 120.0;
    }
    
    pchar->f_sw = 0.99 * loct->soilwtr_h/(loct->soilwtr_h + ksw) + 0.01;
    if(pchar->f_sw < 0.0){
        pchar->f_sw = 0.0;
    }
    if(pchar->f_sw > 1.0){
        pchar->f_sw = 1.0;
    }
    
    /* sensitivity to KN: 2014/09/22 by A.Ito ************************/
    if(CONST_KN == 0){
        kn = pchar->ke_n;
    }else if(CONST_KN == 1){
        /* de Pury & Farquhar (1997) default value */
        kn = 0.713;
    }else if(CONST_KN == 2){
        /* meta-analysis median value */
        kn = 1.22;
    }else if(CONST_KN == 3){
        /* meta-analysis median value */
        
        kl = 0.6404; /* mean non-wheat */
        lat = fabs(grid->lat);
        kb = 0.156 - 0.026404 * pchar->lai + 0.396881 * kl + 0.001723 * lat - 0.000134 * grid->topo;
        
        kn = pchar->lai * kb;
        
        if(kn > 1.0){
            kn  = 1.5;
        }
        if(kn < 0.1){
            kn  = 0.1;
        }
    }
    
    /* fix canopy-top N concentration to fix-CO2 level: 2014/08/17 by A.Ito */
    nc_ct1 = pchar->nc_ct;
    if(FIX_LAI == 3){
        if(pchar->lai > 0.0 && loct->lai > 0.0 && grid->p_lai > 0.0){
            
            nc0 = grid->p_lai * ((pchar->nc_ct - pchar->nc_np)*(1.0 - exp(-1.0*kn))/kn +pchar->nc_np);
            
            nc_ct1 = pchar->nc_np + ((nc0 / loct->lai - pchar->nc_np) * kn) / (1.0 - exp(-1.0*kn));
            
            if(nc_ct1 > pchar->nc_ct * 5.0){
                nc_ct1 = pchar->nc_ct * 5.0;
            }
            if(nc_ct1 < pchar->nc_ct * 0.2){
                nc_ct1 = pchar->nc_ct * 0.2;
            }
            
        }else{
            nc_ct1 = pchar->nc_ct;
        }
    }
    
    /* sensitivity to canopy-top N conc.: 2014/12/10 by A.Ito */
    if(FIX_NTOP == 1){
        nc_ct1 *= 1.1;
    }
    
	/****** sunny ******/
	/* maximum carboxylation rate at 25 deg-C: Eq.22 in DF97 */
	pchar->vcmx_sn0 = pchar->lai * pchar->photocap_n*(nc_ct1 - pchar->nc_np)*
		(1.0-exp(-(kn + pchar->ke_b1*pchar->lai)))/(kn + pchar->ke_b1*pchar->lai);
    //pchar->vcmx_sn0 = pchar->lai*pchar->photocap_n*(nc_ct1 - pchar->nc_np);
    
	/* maximum carboxylation rate: Eq.8 in DF97 */
	pchar->vcmx_sn = pchar->vcmx_sn0 * exp(pchar->acen_vc*(pchar->tmp - 25.0)/(298.15*UGC*(pchar->tmp + ZAT)));
    
    pchar->vcmx_sn *= pchar->f_sw;
    
	/* maximum electron transport rate: proportional Jmax/Vcmax relationship */
	pchar->jmx_sn = pchar->jvr * pchar->vcmx_sn0 * pchar->ft_et;
    
    pchar->jmx_sn *= pchar->f_sw;
    
	/* stabilization */
	if(pchar->apar_sn>0.0 && loct->sinb_h>0.0 && phenol_stage !=0){
		/* electron transport rate: Eq.5 in DF97 */
		pchar->j_sn = ((pchar->apar_sn + pchar->jmx_sn)-sqrt((pchar->apar_sn + pchar->jmx_sn)*
			(pchar->apar_sn+pchar->jmx_sn) - 4.0*pchar->cv_j*pchar->apar_sn*pchar->jmx_sn))/(2.0*pchar->cv_j);
		
        if(pchar->j_sn >= 0.0 && pchar->j_sn <= 10000.0){
            ;
        }else{
            pchar->j_sn = 0.0;
        }
        
		if(SOLV_LEAF == 0){
			/** stabilization with Discrete Newton method **/
			nn = 0;
			crit = 10.0;
			while(crit>0.1 && nn<20){
				pchar->ci_sn = pchar->ci_sn - f_gas_sun(pchar, pchar->ci_sn)*step/
						(f_gas_sun(pchar, pchar->ci_sn+step) - f_gas_sun(pchar, pchar->ci_sn));
				crit = fabs(f_gas_sun(pchar, pchar->ci_sn));
				nn++;
				
				if(NOTICE==1 && pchar->ci_sn<0.0){
					//printf("FARTAL ERROR: negative Ci of sun-leaf: %f %f %f\n", pchar->ci_sn, pchar->p_sn, pchar->gs_sn);
				}
			} 
		}else if(SOLV_LEAF == 1){
			/** analytical solution by Baldocchi (1994) **/
			leaf_ansolv(0, pchar);
		}
	}else{
		/* dark condition */
		pchar->j_sn = 0.0;
		pchar->pc_sn = 0.0;
		pchar->pj_sn = 0.0;
		pchar->p_sn = 0.0;
		pchar->ci_sn = pchar->co2_a;
		
		f_gas_sun(pchar, pchar->ci_sn);
	}

	/****** shade ******/
	/* maximum carboxylation rate at 25 deg-C: Eq.A27 in DF97 */
	/* Vcmax(shade) = Vcmax - Vcmax(sun) */
	pchar->vcmx_sd0 = pchar->lai * pchar->photocap_n*(nc_ct1 - pchar->nc_np)*((1.0 - exp(-kn))/
		kn - (1.0 - exp(-(kn + pchar->ke_b1*pchar->lai)))/(kn + pchar->ke_b1*pchar->lai));
    //pchar->vcmx_sd0 = pchar->lai*pchar->photocap_n*(nc_ct1 - pchar->nc_np);
    
	/* maximum carboxylation rate: Eq.8 in DF97 */
	pchar->vcmx_sd = pchar->vcmx_sd0 * exp(pchar->acen_vc*(pchar->tmp - 25.0)/(298.15*UGC*(pchar->tmp + ZAT)));
    
    pchar->vcmx_sd *= pchar->f_sw;
    
	/* maximum electron transport rate */
	pchar->jmx_sd = (pchar->jvr * pchar->vcmx_sd0) * pchar->ft_et;
    
    pchar->jmx_sd *= pchar->f_sw;
    
	/* stabilization */
	if(pchar->apar_sd>0.0 && loct->sinb_h>0.0){
		/* electron transport rate: Eq.5 in DF97 */
		pchar->j_sd = ((pchar->apar_sd + pchar->jmx_sd)-sqrt((pchar->apar_sd + pchar->jmx_sd)*(pchar->apar_sd+
			pchar->jmx_sd) - 4.0*pchar->cv_j * pchar->apar_sd*pchar->jmx_sd))/(2.0*pchar->cv_j);
		
        if(pchar->j_sd >= 0.0 && pchar->j_sd <= 10000.0){
            ;
        }else{
            pchar->j_sd = 0.0;
        }
        
		if(SOLV_LEAF == 0){
			/** stabilization with Discrete Newton method **/
			nn=0;
			crit=10.0;
			while(crit>0.1 && nn<20){
				pchar->ci_sd = pchar->ci_sd - f_gas_shade(pchar, pchar->ci_sd)*step/
						(f_gas_shade(pchar, pchar->ci_sd+step) - f_gas_shade(pchar, pchar->ci_sd));
				crit = fabs(f_gas_shade(pchar, pchar->ci_sd));
				nn++;

				if(NOTICE==1 && pchar->ci_sd<0.0){
					//printf("FARTAL ERROR: negative Ci of shade-leaf: %f\n", pchar->ci_sd);
				}
			} 
		}else if(SOLV_LEAF == 1){
			/** analytical solution by Baldocchi (1994) **/
			leaf_ansolv(1, pchar); 
		}
	}else{
		/* dark condition */
		pchar->j_sd = 0.0;
		pchar->pc_sd = 0.0;
		pchar->pj_sd = 0.0;
		pchar->p_sd = 0.0;
		pchar->ci_sd = pchar->co2_a;
		
		f_gas_shade(pchar, pchar->ci_sd);
	}
    
    pchar->vcmx = pchar->vcmx_sn + pchar->vcmx_sd;
}

/* objective function of Discrete Newton method, sun-leaf ***************************/
float f_gas_sun(
	struct Pchar *pchar,
	float ci
){	
	float fgc, gb, cs;
	
	/* stomatal conductance of sun leaf */
	if(pchar->lai_sn>1.0){
		fgc = pchar->lai_sn;
	}else{
		fgc = 1.0;
	}
	gb = GB*fgc/100.0;
	pchar->gb = gb;
	pchar->rd_sn = pchar->rd0*pchar->ft_rd*fgc;
	
	/* CO2 assimilation rate */
	if(ci > pchar->comp_co2){
		/* Rubisco-limited */
		pchar->pc_sn = pchar->vcmx_sn*(ci - pchar->comp_co2)/(ci + pchar->k_effc);
		/* RuBP-regeneration-limited */
		pchar->pj_sn = pchar->j_sn*(ci - pchar->comp_co2)/(4.0*(ci + 2.0*pchar->comp_co2));
	}else{
		pchar->pc_sn = 0.0;
		pchar->pj_sn = 0.0;
	}
	
	/* combination of limitations ***********/
	pchar->p_sn = (pchar->pc_sn<pchar->pj_sn)?pchar->pc_sn:pchar->pj_sn; /*	*/
	/* pchar->p_sn = ((pchar->pc_sn+pchar->pj_sn)-sqrt((pchar->pc_sn+pchar->pj_sn)*(pchar->pc_sn+pchar->pj_sn)-
		4.0*pchar->cv_p*pchar->pc_sn*pchar->pj_sn))/(2.0*pchar->cv_p); */
	/* cv_p, convexity: 0.877 () */
	
	if(NOTICE==1 && pchar->p_sn<0.0){
		printf("FARTAL ERROR: negative P of sun-leaf: %f\n", pchar->p_sn);
	}
	
	pchar->a_sn = pchar->p_sn - pchar->rd_sn;
	
	cs = pchar->co2_a - pchar->a_sn/gb;
    
    /* fix CO2 for stomata: 2014/09/19 */
    if(FIX_CO2 == 2){
        cs = 33.0;
    }
	
	if(pchar->a_sn>0.0){
		pchar->gs_sn = pchar->gss_b0/100.0 + pchar->gss_b1/100.0/10.0*pchar->a_sn/
					(cs*(1.0 + pchar->vpd_a/pchar->gss_b2));
		pchar->gs_sn *= fgc;
	}else{
		pchar->gs_sn = pchar->gss_b0/100.0*fgc;
	}
	
	if(NOTICE==1 && pchar->gs_sn<0.0){
		/* printf("FARTAL ERROR: negative gs of sun-leaf: %f %f\n", pchar->a_sn, pchar->gs_sn); */
	}

	return (ci + pchar->a_sn/gb + pchar->a_sn/pchar->gs_sn - pchar->co2_a);
}

/**** objective function of Discrete Newton method, shade-leaf ********************************/
float f_gas_shade(
	struct Pchar *pchar,
	float ci
){	
	float fgc, gb, cs;

	/* stomatal conductance of shade leaf */
	if(pchar->lai_sd>1.0){
		fgc = pchar->lai_sd;
	}else{
		fgc = 1.0;
	}
	gb = GB*fgc/100.0;
	pchar->gb = gb;
	pchar->rd_sd = pchar->rd0*pchar->ft_rd*fgc;
	
	/* CO2 assimilation rate */
	if(ci > pchar->comp_co2){
		/* Rubisco-limited */
		pchar->pc_sd = pchar->vcmx_sd*(ci - pchar->comp_co2)/(ci + pchar->k_effc);
		/* RuBP-regeneration-limited */
		pchar->pj_sd = pchar->j_sd*(ci - pchar->comp_co2)/(4.0*(ci + 2.0*pchar->comp_co2));
	}else{
		pchar->pc_sd = 0.0;
		pchar->pj_sd = 0.0;
	}
	/* combination */
	pchar->p_sd = (pchar->pc_sd<pchar->pj_sd)?pchar->pc_sd:pchar->pj_sd; /*	*/
	/* pchar->p_sd = ((pchar->pc_sd+pchar->pj_sd)-sqrt((pchar->pc_sd+pchar->pj_sd)*(pchar->pc_sd+pchar->pj_sd)-
		4.0*pchar->cv_p*pchar->pc_sd*pchar->pj_sd))/(2.0*pchar->cv_p); */
		
	if(NOTICE==1 && pchar->p_sd<0.0){
		printf("FARTAL ERROR: negative P of shade-leaf: %f %f %f %f\n", 
			ci, pchar->p_sd, pchar->pc_sd, pchar->pj_sd);
	}
	
	/* leaf-level net assimilation rate */
	pchar->a_sd = pchar->p_sd - pchar->rd_sd;
	
	/* sub stomatal CO2 partial pressure */
	cs = pchar->co2_a - pchar->a_sd/gb;
	
    /* fix CO2 for stomata: 2014/09/19 */
    if(FIX_CO2 == 2){
        cs = 33.0;
    }
	
	if(pchar->a_sd > 0.0){
		pchar->gs_sd = pchar->gss_b0/100.0 + pchar->gss_b1/100.0/10.0*pchar->a_sd/
				(cs*(1.0 + pchar->vpd_a/pchar->gss_b2));
		/* pchar->gs_sd = pchar->gss_b0 + pchar->gss_b1*pchar->a_sd/
				((ci-pchar->comp_co2)*10.0*(1.0+pchar->vpd_a/pchar->gss_b2)); */
		pchar->gs_sd *= fgc;
	}else{
		pchar->gs_sd = pchar->gss_b0*fgc;
	}

	if(NOTICE==1 && pchar->gs_sd<0.0){
		/* printf("FARTAL ERROR: negative gs of shade-leaf: %f %f\n", pchar->a_sd, pchar->gs_sd); */
	}
	
	return (ci + pchar->a_sd/gb + pchar->a_sd/pchar->gs_sd - pchar->co2_a);
}

/* *******************************************************/
/*
Baldocchi D. (1994). An analytical solution for coupled leaf photosynthesis and 
stomatal conductance models. Tree Physiology, 14, 1069-1079.
*/
 void leaf_ansolv(
 	long frac, 
 	struct Pchar *pchar
 ){
 	float aa, bb, dd, ee, fgc;
 	float ca, cs, gb, bd, rh, m;
 	float alpha, beta, gamma, theta;
 	float pp, qq, rr, qqq, rrr;
 	float xx1, xx2, xx3;
 	float f_pp, f_qq, f_cube;
 	 	
 	/**************************************************************************************/
 	gb = GB;						/* leaf boundary layer conductance*/
 	bd = pchar->gss_b0;				/* minimum conductance */
 	rh = 1.0/(1.0+pchar->vpd_a/pchar->gss_b2);	/* relative humidity */
 	m = pchar->gss_b1/10.0;			/* coefficient of slope */
 	ca = pchar->co2_a;				/* atmospheric CO2, Pa */
 	
    fgc = 1.0;
 	if(frac == 0){
		 /* for sun leaf */
		if(pchar->lai_sn > 1.0){
			fgc = pchar->lai_sn;
		}else{
			fgc = 1.0;
		}
		pchar->rd_sn = pchar->rd0 * pchar->ft_rd * fgc;
	}else if(frac == 1){
	 	/* for shade leaf */
		if(pchar->lai_sd>1.0){
			fgc = pchar->lai_sd;
		}else{
			fgc = 1.0;
		}
		pchar->rd_sd = pchar->rd0 * pchar->ft_rd * fgc;
	}
	gb *= fgc/100.0;
	bd *= fgc/100.0;
	m *= fgc/100.0;
	
	pchar->gb = gb;
 	
 	/* Rubisco-limited ************************************************************************/
    aa = 0.0;
    if(frac==0){
		aa = pchar->vcmx_sn;
	}else if(frac==1){
		aa = pchar->vcmx_sd;
	}
	bb = pchar->k_effc;
	dd = pchar->comp_co2;
	ee = 1.0;
	
	/******************************************************************/
 	alpha = 1.0+bd/gb - m*rh;
 	beta = ca * (gb * m * rh - 2.0*bd - gb);
 	gamma = ca*ca * bd * gb;
 	theta = gb * m * rh - bd;
 	
 	/******************************************************************/
 	if(frac == 0){
		pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*pchar->rd_sn)/(ee * alpha);
	 	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*pchar->rd_sn*beta 
	 		+ pchar->rd_sn*bb*theta)/(ee * alpha);
	 	rr = (-aa*gamma + aa*dd*gamma/ca + ee*pchar->rd_sn*gamma 
	 		+ pchar->rd_sn*bb*gamma/ca)/(ee * alpha);
 	}else if(frac == 1){
		pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*pchar->rd_sd)/(ee * alpha);
	 	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*pchar->rd_sd*beta 
	 		+ pchar->rd_sd*bb*theta)/(ee * alpha);
	 	rr = (-aa*gamma + aa*dd*gamma/ca + ee*pchar->rd_sd*gamma 
	 		+ pchar->rd_sd*bb*gamma/ca)/(ee * alpha);
 	}else{
        printf("Wrong frac: %ld\n", frac);
        pp = qq = rr = 1.0;
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
	 		 	
	 	if(frac == 0){
	 		pchar->pc_sn = xx3 + pchar->rd_sn;
	 	}else if(frac == 1){
	 		pchar->pc_sd = xx3 + pchar->rd_sd;
	 	}
	 }else{
	 	if(frac == 0){
	 		pchar->pc_sn = 0.0;
	 	}else if(frac == 1){
	 		pchar->pc_sd = 0.0;
	 	}
	 }

 	/* RUBP-limited ************************************************************************/
 	if(frac == 0){
 		aa = pchar->j_sn;
 	}else if(frac == 1){
 		aa = pchar->j_sd;
 	}
 	bb = 8.0*pchar->comp_co2;
 	dd = pchar->comp_co2;
 	ee = 4.0;
 	
	/******************************************************************/
 	alpha = 1.0+bd/gb - m*rh;
 	beta = ca * (gb * m * rh - 2.0*bd - gb);
 	gamma = ca*ca * bd * gb;
 	theta = gb * m * rh - bd;
 	
 	/******************************************************************/
 	if(frac==0){
		pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*pchar->rd_sn)/(ee * alpha);
	 	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*pchar->rd_sn*beta 
	 		+ pchar->rd_sn*bb*theta)/(ee * alpha);
	 	rr = (-aa*gamma + aa*dd*gamma/ca + ee*pchar->rd_sn*gamma 
	 		+ pchar->rd_sn*bb*gamma/ca)/(ee * alpha);
 	}else if(frac==1){
		pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*pchar->rd_sd)/(ee * alpha);
	 	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*pchar->rd_sd*beta 
	 		+ pchar->rd_sd*bb*theta)/(ee * alpha);
	 	rr = (-aa*gamma + aa*dd*gamma/ca + ee*pchar->rd_sd*gamma 
	 		+ pchar->rd_sd*bb*gamma/ca)/(ee * alpha);
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
	 		pchar->pj_sn = xx3 + pchar->rd_sn;
	 	}else if(frac==1){
	 		pchar->pj_sd = xx3 + pchar->rd_sd;
	 	}
	 }else{
	 	if(frac==0){
	 		pchar->pj_sn = 0.0;
	 	}else if(frac==1){
	 		pchar->pj_sd = 0.0;
	 	}
	 }

 	/* Solution ************************************************************************/
 	if(frac==0){ /* sunny */
 		/* pchar->p_sn = ((pchar->pc_sn+pchar->pj_sn)-sqrt((pchar->pc_sn+pchar->pj_sn)*(pchar->pc_sn+pchar->pj_sn)-
					4.0*pchar->cv_p*pchar->pc_sn*pchar->pj_sn))/(2.0*pchar->cv_p); */
					
		pchar->p_sn = (pchar->pc_sn<pchar->pj_sn)?pchar->pc_sn:pchar->pj_sn;
		
		if(pchar->p_sn>=0.0 && pchar->p_sn<=100.0){
			;
		}else{
			/* printf("!! Warning: [p_sn, %f] in [leaf_ansolv()] !!\n", pchar->p_sn); */
			/* printf("!! Warning: [p_sn, %f %f %f] in [leaf_ansolv()] !!\n", pchar->p_sn, pchar->pc_sn, pchar->pj_sn); */
		}
		
		pchar->a_sn = pchar->p_sn - pchar->rd_sn;
		
		/* substomatal cavity CO2 */
		cs = ca - pchar->a_sn/gb;
		
		pchar->gs_sn = bd + m*pchar->a_sn*rh/cs;		
		pchar->ci_sn = cs - pchar->a_sn/pchar->gs_sn;
	}else if(frac==1){ /* shade */
		/* pchar->p_sd = ((pchar->pc_sd+pchar->pj_sd)-sqrt((pchar->pc_sd+pchar->pj_sd)*(pchar->pc_sd+pchar->pj_sd)-
					4.0*pchar->cv_p*pchar->pc_sd*pchar->pj_sd))/(2.0*pchar->cv_p); */
					
		pchar->p_sd = (pchar->pc_sd<pchar->pj_sd)?pchar->pc_sd:pchar->pj_sd;
		
		if(pchar->p_sd>=0.0 && pchar->p_sd<=100.0){
			;
		}else{
			/* printf("!! Warning: [p_sd, %f] in [leaf_ansolv()] !!\n", pchar->p_sd); */
			/* printf("!! Warning: [p_sd, %f %f %f] in [leaf_ansolv()] !!\n", pchar->p_sd, pchar->pc_sd, pchar->pj_sd); */
		}
		
		pchar->a_sd = pchar->p_sd - pchar->rd_sd;
		
		/* substomatal cavity CO2 */
 		cs = ca - pchar->a_sd/gb;
				
		pchar->gs_sd = bd + m*pchar->a_sd*rh/cs;		
		pchar->ci_sd = cs - pchar->a_sd/pchar->gs_sd;
	}
}

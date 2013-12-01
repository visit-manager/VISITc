/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/*  functions  of  soil  process  */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* daily prcesses in soil organic matter  **********************************/
void f_cycle_soil(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Schar *schar, 
	struct Smas *mass,  
	struct Sflx *flux
){
	double degrade_tf, degrade_tc, degrade_tr;
	double degrade_gf, degrade_gc, degrade_gr;
	double degrade_ha, degrade_hi, degrade_hp;
	double hf_a, hf_i, hf_p;
	double d13c_hf_a,  d13c_hf_i,  d13c_hf_p;
	double d13c_mrl, d13c_mrh;
	double d13c_ltr, d13c_msl;
		
	/** environmental scalars **/
	schar->f_tm_l = frl(grid, loct, schar);		/* litter */
	schar->f_tm_h = frh(grid, loct, schar);		/* humus */
	
	/* degraded soil organic carbon */
	degrade_tf = mass->ltr_tf * schar->sr_lf/1000.0 * schar->f_tm_l;
	degrade_tc = mass->ltr_tc * schar->sr_lc/1000.0 * schar->f_tm_l;
	degrade_tr = mass->ltr_tr * schar->sr_lr/1000.0 * schar->f_tm_l;	
	degrade_gf = mass->ltr_gf * schar->sr_lf/1000.0 * schar->f_tm_l;
	degrade_gc = mass->ltr_gc * schar->sr_lc/1000.0 * schar->f_tm_l;
	degrade_gr = mass->ltr_gr * schar->sr_lr/1000.0 * schar->f_tm_l;	
	degrade_ha = mass->msl_a * schar->sr_ha/1000.0 * schar->f_tm_h;
	degrade_hi = mass->msl_i * schar->sr_hi/1000.0 * schar->f_tm_h;
	degrade_hp = mass->msl_p * schar->sr_hp/1000.0 * schar->f_tm_h;

	/* microbial respiration */
	flux->mr_tf = degrade_tf * schar->f_co2_lf;
	flux->mr_tc = degrade_tc * schar->f_co2_lc;
	flux->mr_tr = degrade_tr * schar->f_co2_lr;
	flux->mr_gf = degrade_gf * schar->f_co2_lf;
	flux->mr_gc = degrade_gc * schar->f_co2_lc;
	flux->mr_gr = degrade_gr * schar->f_co2_lr;
	flux->mr_a = degrade_ha;
	flux->mr_i = degrade_hi;
	flux->mr_p = degrade_hp;
	
	/* humus formation */
	flux->hf_tfa = degrade_tf*(1.0-schar->f_co2_lf)*schar->f_hm_a;	/* leaf litter -> active */
	flux->hf_tfi = degrade_tf*(1.0-schar->f_co2_lf)*schar->f_hm_i;	/* leaf litter -> intermediate */
	flux->hf_tfp = degrade_tf*(1.0-schar->f_co2_lf)*schar->f_hm_p;	/* leaf litter -> passive */
	flux->hf_tca = degrade_tc*(1.0-schar->f_co2_lc)*schar->f_hm_a;	/* stem litter -> active */
	flux->hf_tci = degrade_tc*(1.0-schar->f_co2_lc)*schar->f_hm_i;	/* stem litter -> intermediate */
	flux->hf_tcp = degrade_tc*(1.0-schar->f_co2_lc)*schar->f_hm_p;	/* stem litter -> passive */
	flux->hf_tra = degrade_tr*(1.0-schar->f_co2_lr)*schar->f_hm_a;	/* root litter -> active */
	flux->hf_tri = degrade_tr*(1.0-schar->f_co2_lr)*schar->f_hm_i;	/* root litter -> intermediate */
	flux->hf_trp = degrade_tr*(1.0-schar->f_co2_lr)*schar->f_hm_p;	/* root litter -> passive */

	flux->hf_gfa = degrade_gf*(1.0-schar->f_co2_lf)*schar->f_hm_a;	/* leaf litter -> active */
	flux->hf_gfi = degrade_gf*(1.0-schar->f_co2_lf)*schar->f_hm_i;	/* leaf litter -> intermediate */
	flux->hf_gfp = degrade_gf*(1.0-schar->f_co2_lf)*schar->f_hm_p;	/* leaf litter -> passive */
	flux->hf_gca = degrade_gc*(1.0-schar->f_co2_lc)*schar->f_hm_a;	/* stem litter -> active */
	flux->hf_gci = degrade_gc*(1.0-schar->f_co2_lc)*schar->f_hm_i;	/* stem litter -> intermediate */
	flux->hf_gcp = degrade_gc*(1.0-schar->f_co2_lc)*schar->f_hm_p;	/* stem litter -> passive */
	flux->hf_gra = degrade_gr*(1.0-schar->f_co2_lr)*schar->f_hm_a;	/* root litter -> active */
	flux->hf_gri = degrade_gr*(1.0-schar->f_co2_lr)*schar->f_hm_i;	/* root litter -> intermediate */
	flux->hf_grp = degrade_gr*(1.0-schar->f_co2_lr)*schar->f_hm_p;	/* root litter -> passive */

	hf_a = flux->hf_tfa + flux->hf_tca + flux->hf_tra 
		+ flux->hf_gfa + flux->hf_gca + flux->hf_gra;	/* input to active humus */
	hf_i = flux->hf_tfi + flux->hf_tci + flux->hf_tri 
		+ flux->hf_gfi + flux->hf_gci + flux->hf_gri;		/* input to intermediate humus */
	hf_p = flux->hf_tfp + flux->hf_tcp + flux->hf_trp 
		+ flux->hf_gfp + flux->hf_gcp + flux->hf_grp;	/* input to passive humus */

	/*  mass  balance  */
	/* litter */
	mass->ltr_tf += flux->li_tf - degrade_tf; 
	mass->ltr_tc += flux->li_tc - degrade_tc; 
	mass->ltr_tr += flux->li_tr - degrade_tr; 
	mass->ltr_gf += flux->li_gf - degrade_gf; 
	mass->ltr_gc += flux->li_gc - degrade_gc; 
	mass->ltr_gr += flux->li_gr - degrade_gr; 
	/* humus: mineral soil */
	mass->msl_a += hf_a - degrade_ha;
	mass->msl_i += hf_i - degrade_hi;
	mass->msl_p += hf_p - degrade_hp;
	mass->soil = mass->ltr_tf+mass->ltr_tc+mass->ltr_tr + 
			mass->ltr_gf+mass->ltr_gc+mass->ltr_gr + 
			mass->msl_a+mass->msl_i+mass->msl_p;
	
	/* error massages for negative mass balance */
	if(mass->ltr_tf<0.0){		
		printf("!!!---NEGATIVE    TREE  LEAF  LITTER---!!!\n");	
		mass->ltr_tf = 0.0;		
	}
	if(mass->ltr_tc<0.0){		printf("!!!---NEGATIVE    TREE  STEM  LITTER---!!!\n");	mass->ltr_tc = 0.0;		}
	if(mass->ltr_tr<0.0){		printf("!!!---NEGATIVE    TREE  ROOT  LITTER---!!!\n");	mass->ltr_tr = 0.0;		}
	if(mass->ltr_gf<0.0){		printf("!!!---NEGATIVE    GRASS  LEAF  LITTER---!!!\n");	mass->ltr_gf = 0.0;		}
	if(mass->ltr_gc<0.0){		printf("!!!---NEGATIVE     GRASS  STEM  LITTER---!!!\n");	mass->ltr_gc = 0.0;		}
	if(mass->ltr_gr<0.0){		printf("!!!---NEGATIVE     GRASS  ROOT  LITTER---!!!\n");	mass->ltr_gr = 0.0;		}
	if(mass->msl_a<0.0){		printf("!!!---NEGATIVE    ACTIVE  HUMUS---!!!\n");	mass->msl_a = 0.0;	}
	if(mass->msl_i<0.0){		printf("!!!---NEGATIVE    INT.  HUMUS---!!!\n");	mass->msl_i = 0.0;		}
	if(mass->msl_p<0.0){		printf("!!!---NEGATIVE    PASSIVE  HUMUS---!!!\n");	mass->msl_p = 0.0;	}  
	
	/*  total microbial respiration  */
	flux->hr = flux->mr_tf+flux->mr_tc+flux->mr_tr  +  
			flux->mr_gf+flux->mr_gc+flux->mr_gr  +  
			flux->mr_a+flux->mr_i+flux->mr_p;
	
	/* stable carbon isotope ******************************************************************/	
	if(SCI_SCHEME==1){
		flux->d13c_mr_f = mass->d13c_ltr_f;
		flux->d13c_mr_c = mass->d13c_ltr_c;
		flux->d13c_mr_r = mass->d13c_ltr_r;
		flux->d13c_mr_a = mass->d13c_msl_a;
		flux->d13c_mr_i = mass->d13c_msl_i;
		flux->d13c_mr_p = mass->d13c_msl_p;

		flux->d13c_hf_f = mass->d13c_ltr_f;
		flux->d13c_hf_c = mass->d13c_ltr_c;
		flux->d13c_hf_r = mass->d13c_ltr_r;

		d13c_hf_a = d13c_addition3(mass->d13c_ltr_f, flux->hf_tfa,
								mass->d13c_ltr_c, flux->hf_tca,
								mass->d13c_ltr_r, flux->hf_tra);
		d13c_hf_i = d13c_addition3(mass->d13c_ltr_f, flux->hf_tfi,
								mass->d13c_ltr_c, flux->hf_tci,
								mass->d13c_ltr_r, flux->hf_tri);
		d13c_hf_p = d13c_addition3(mass->d13c_ltr_f, flux->hf_tfp,
								mass->d13c_ltr_c, flux->hf_tcp,
								mass->d13c_ltr_r, flux->hf_trp);

		mass->d13c_ltr_f = d13c_addition(flux->d13c_li_f, flux->li_tf, mass->d13c_ltr_f,  mass->ltr_tf-degrade_tf);
		mass->d13c_ltr_c = d13c_addition(flux->d13c_li_c, flux->li_tc, mass->d13c_ltr_c,  mass->ltr_tc-degrade_tc);
		mass->d13c_ltr_r = d13c_addition(flux->d13c_li_r, flux->li_tr, mass->d13c_ltr_r,  mass->ltr_tr-degrade_tr);
		mass->d13c_msl_a = d13c_addition(d13c_hf_a, hf_a,  mass->d13c_msl_a, mass->msl_a-degrade_ha);
		mass->d13c_msl_i = d13c_addition(d13c_hf_i, hf_i,  mass->d13c_msl_i, mass->msl_i-degrade_hi);
		mass->d13c_msl_p = d13c_addition(d13c_hf_p, hf_p,  mass->d13c_msl_p, mass->msl_p-degrade_hp);

		/* stable carbon isotope */
		d13c_mrl = d13c_addition3(flux->d13c_mr_f, flux->mr_tf, 
								flux->d13c_mr_c, flux->mr_tc, 
								flux->d13c_mr_r, flux->mr_tr);
		d13c_mrh = d13c_addition3(flux->d13c_mr_a, flux->mr_a, 
								flux->d13c_mr_i, flux->mr_i, 
								flux->d13c_mr_p, flux->mr_p);
		flux->d13c_hr = d13c_addition(d13c_mrl, flux->mr_tf+flux->mr_tc+flux->mr_tr,
								d13c_mrh, flux->mr_a+flux->mr_i+flux->mr_p);
		
		d13c_ltr = d13c_addition3(mass->d13c_ltr_f, mass->ltr_tf, 
								mass->d13c_ltr_c, mass->ltr_tc, 
								mass->d13c_ltr_r, mass->ltr_tr);
		d13c_msl = d13c_addition3(mass->d13c_msl_a, mass->msl_a, 
								mass->d13c_msl_i, mass->msl_i, 
								mass->d13c_msl_p, mass->msl_p);
		mass->d13c_soil = d13c_addition(d13c_ltr, mass->ltr_tf+mass->ltr_tc+mass->ltr_tr, 
									d13c_msl, mass->msl_a+mass->msl_i+mass->msl_p);
	}
	
	if(N_CYCLE ==1){
		/* nitrogen dynamics *******************************************/
		/* mineralization of organic N */
		f_n_mineralz(grid, loct, mass, flux);
		
		/* NH3 volatilization */
		f_nh3_volatilization(grid, loct, mass, flux);
		
		/* NO3- leaching */
		f_n_leaching(grid, loct, mass, flux);
		
		/* microbial immobilization */
		f_n_immoblz(grid, loct, schar, mass, flux);
		
		/* microbial abandonment */
		f_n_mcrb_abdn(grid, loct, schar, mass, flux);
	}
}


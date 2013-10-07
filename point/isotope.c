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

/* initialization ************************************************/
void f_init_d13c(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	/* mass */
	(mass->tree).d13c_fol = (mass->c3).d13c_fol = loct->d13C_bCO2; 	
	(mass->tree).d13c_stm = (mass->c3).d13c_stm = loct->d13C_bCO2; 	
	(mass->tree).d13c_rot = (mass->c3).d13c_rot = loct->d13C_bCO2; 	
	(mass->soil).d13c_ltr_f = loct->d13C_bCO2; 	
	(mass->soil).d13c_ltr_c = loct->d13C_bCO2; 	
	(mass->soil).d13c_ltr_r = loct->d13C_bCO2; 	
	(mass->soil).d13c_msl_a = loct->d13C_bCO2; 	
	(mass->soil).d13c_msl_i = loct->d13C_bCO2; 	
	(mass->soil).d13c_msl_p = loct->d13C_bCO2; 	

	(mass->tree).d13c_plant = (mass->c3).d13c_plant = loct->d13C_bCO2; 	
	(mass->soil).d13c_soil = (mass->soil).d13c_soil = (mass->soil).d13c_soil = loct->d13C_bCO2; 	
	mass->d13c_total = loct->d13C_bCO2; 	
	
	flux->d13c_er = loct->d13C_bCO2; 

	/* ecophysiology */
	(echar->tree).photo_13c_frac = 0.0;
	(echar->c3).photo_13c_frac = 0.0;

	/* flux */
	(flux->tree).d13c_gpp=(flux->c3).d13c_gpp=loct->d13C_bCO2;
	(flux->tree).d13c_epp=(flux->c3).d13c_epp=loct->d13C_bCO2;
	(flux->tree).d13c_spp=(flux->c3).d13c_spp=loct->d13C_bCO2;
	(flux->tree).d13c_npp=(flux->c3).d13c_npp=loct->d13C_bCO2;
	
	(flux->tree).d13c_rfm=(flux->c3).d13c_rfm=loct->d13C_bCO2;
	(flux->tree).d13c_rcm=(flux->c3).d13c_rcm=loct->d13C_bCO2;
	(flux->tree).d13c_rrm=(flux->c3).d13c_rrm=loct->d13C_bCO2;
	(flux->tree).d13c_rpm=(flux->c3).d13c_rpm=loct->d13C_bCO2;
	(flux->tree).d13c_rfg=(flux->c3).d13c_rfg=loct->d13C_bCO2;
	(flux->tree).d13c_rcg=(flux->c3).d13c_rcg=loct->d13C_bCO2;
	(flux->tree).d13c_rrg=(flux->c3).d13c_rrg=loct->d13C_bCO2;
	(flux->tree).d13c_rpg=(flux->c3).d13c_rpg=loct->d13C_bCO2; 
	(flux->tree).d13c_rp=(flux->c3).d13c_rp=loct->d13C_bCO2; 
			
	(flux->tree).d13c_lf=(flux->c3).d13c_lf=loct->d13C_bCO2; 
	(flux->tree).d13c_lc=(flux->c3).d13c_lc=loct->d13C_bCO2; 
	(flux->tree).d13c_lr=(flux->c3).d13c_lr=loct->d13C_bCO2; 
	(flux->tree).d13c_lL=(flux->c3).d13c_lL=loct->d13C_bCO2; 
	
	(flux->tree).d13c_lf_c=(flux->c3).d13c_lf_c=loct->d13C_bCO2; 

	(flux->tree).d13c_tpp=(flux->c3).d13c_tpp=loct->d13C_bCO2; 
	(flux->tree).d13c_tpf=(flux->c3).d13c_tpf=loct->d13C_bCO2; 
	(flux->tree).d13c_tpc=(flux->c3).d13c_tpc=loct->d13C_bCO2; 
	(flux->tree).d13c_tpr=(flux->c3).d13c_tpr=loct->d13C_bCO2; 
	
	(flux->tree).d13c_hvst=(flux->c3).d13c_hvst=loct->d13C_bCO2; 	
	
	(flux->soil).d13c_li_f=(flux->soil).d13c_li_c=(flux->soil).d13c_li_r=loct->d13C_bCO2; 
	(flux->soil).d13c_hf_f=(flux->soil).d13c_hf_c=(flux->soil).d13c_hf_r=loct->d13C_bCO2; 
	(flux->soil).d13c_mr_f=(flux->soil).d13c_mr_c=(flux->soil).d13c_mr_r=loct->d13C_bCO2; 
	(flux->soil).d13c_mr_a=(flux->soil).d13c_mr_i=(flux->soil).d13c_mr_p=loct->d13C_bCO2; 
	(flux->soil).d13c_hr=loct->d13C_bCO2; 
	
	flux->d13c_nep=loct->d13C_bCO2; 
	flux->d13c_ncb=loct->d13C_bCO2; 
}

/***** photosynthetic carbon isotope descrimination *****/
void photo_13c_frac(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *plant
){
	double aaa, bbb, bb3, bb4, fai;
	
	aaa = 4.4;
	bbb = 27.0;
	bb3 = 30.0;
	bb4 = -5.7;
	fai = 0.21;

	if(plant->phototype == 3){
		/* c3 */
		plant->photo_13c_frac = aaa+(bbb-aaa)*plant->ci/loct->aCO2;
	}else if(plant->phototype == 4){
		/* c3 */
		plant->photo_13c_frac = aaa+(bb4+bb3*fai-aaa)*plant->ci/loct->aCO2;
	}else{
		printf("!! Bad photosynthetic type\n");
		exit(1);
	}
}

/***** conversion from d13C to 12C/13C ratio *****/
double deltaTratio(
	double delta
){
	double ratio;
	
	ratio = (1.0 + delta/1000.0)*STCIR;

	return ratio;
}

/***** conversion from 12C/13C ratio to d13C *****/
double ratioTdelta(
	double ratio
){
	double delta;
	
	delta = (ratio/STCIR - 1.0)*1000.0;
	
	return delta;
}

/*******  addition of two isotopically different substances, a and b  *******/
double d13c_addition(
	double d13c_a, 
	double mass_a, 
	double d13c_b, 
	double mass_b
){
	double  ratio_a,  ratio_b;
	double  c12_a, c13_a,  c12_b,  c13_b,  c12_pro,  c13_pro;
	double  d13c_product;
	
	ratio_a = deltaTratio(d13c_a); 
	ratio_b = deltaTratio(d13c_b); 
	
	c12_a = mass_a/(1.0+ratio_a);  c13_a=mass_a-c12_a;
	c12_b = mass_b/(1.0+ratio_b);  c13_b=mass_b-c12_b;
	
	c12_pro = c12_a+c12_b;
	c13_pro = c13_a+c13_b;
	
	if(c12_pro>0.0){
		d13c_product = (c13_pro/c12_pro/STCIR-1.0)*1000.0;
	}else{
		d13c_product = d13c_a;
	}
	
	return(d13c_product);
}

/*******  addition of two isotopically different substances, a and b  *******/
double d13c_addition3(
	double d13c_a, 
	double mass_a, 
	double d13c_b, 
	double mass_b, 
	double d13c_c, 
	double mass_c
){
	double  d13c_ab,  mass_ab, d13c_abc;

	mass_ab = mass_a+mass_b;
	d13c_ab = d13c_addition(d13c_a, mass_a, d13c_b, mass_b);
	
	d13c_abc = d13c_addition(d13c_ab, mass_ab, d13c_c, mass_c);
		
	return(d13c_abc);
}

/******  d13c of efflux CO2 from terrestrial ecosystems  *******/
void d13c_efflux(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Flux *flux
){
	double mr_ltr, mr_msl;
	double d13c_mr_ltr, d13c_mr_msl;
	
	/*  total  ecosystem  respiratory  efflux  of CO2 */
	flux->er = (flux->tree).rp+(flux->c3).rp+(flux->soil).hr;

	(flux->tree).d13c_rpm = d13c_addition3((flux->tree).d13c_rfm, (flux->tree).rfm,   
			(flux->tree).d13c_rcm, (flux->tree).rcm,   (flux->tree).d13c_rrm, (flux->tree).rrm);
	(flux->tree).d13c_rpg = d13c_addition3((flux->tree).d13c_rfg, (flux->tree).rfg,   
			(flux->tree).d13c_rcg, (flux->tree).rcg,   (flux->tree).d13c_rrg, (flux->tree).rrg);
	(flux->tree).d13c_rp = d13c_addition((flux->tree).d13c_rpm,    (flux->tree).rpm, 
			(flux->tree).d13c_rpg, (flux->tree).rpg);

	(flux->c3).d13c_rpm = d13c_addition3((flux->c3).d13c_rfm, (flux->c3).rfm,   
			(flux->c3).d13c_rcm, (flux->c3).rcm,   (flux->c3).d13c_rrm, (flux->c3).rrm);
	(flux->c3).d13c_rpg = d13c_addition3((flux->c3).d13c_rfg, (flux->c3).rfg,   
			(flux->c3).d13c_rcg, (flux->c3).rcg,   (flux->c3).d13c_rrg, (flux->c3).rrg);
	(flux->c3).d13c_rp = d13c_addition((flux->c3).d13c_rpm,    
			(flux->c3).rpm, (flux->c3).d13c_rpg, (flux->c3).rpg);
	
	mr_ltr = (flux->soil).mr_tf+(flux->soil).mr_tc+(flux->soil).mr_tr;
	mr_msl = (flux->soil).mr_a+(flux->soil).mr_i+(flux->soil).mr_p;
	d13c_mr_ltr = d13c_addition3(	(flux->soil).d13c_mr_f, (flux->soil).mr_tf, 
							(flux->soil).d13c_mr_c, (flux->soil).mr_tc, 
							(flux->soil).d13c_mr_r, (flux->soil).mr_tr);
	d13c_mr_msl = d13c_addition3(	(flux->soil).d13c_mr_a, (flux->soil).mr_a, 
							(flux->soil).d13c_mr_i, (flux->soil).mr_i, 
							(flux->soil).d13c_mr_p, (flux->soil).mr_p);
	
	(flux->soil).d13c_hr = d13c_addition(d13c_mr_ltr, mr_ltr,  d13c_mr_msl, mr_msl);
	
	/*  total  ecosystem  respiratory  efflux  of  d13C - CO2 */
	flux->d13c_er = d13c_addition(
		d13c_addition((flux->tree).d13c_rp, (flux->tree).rp, (flux->c3).d13c_rp, 
		(flux->c3).rp), (flux->tree).rp+(flux->c3).rp,   
		(flux->soil).d13c_hr, (flux->soil).hr);
}

/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/*  Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 14, 2007							*/

/*  [history]
	create: 2006 / 07 / 09 by Akihiko Ito
	modify: GHG schemes addition by Motoko Inatomi
	modify: 2007 / 05 / 30 by Akihiko Ito: Version 0.2.4.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* formula of  DAILY N2O emission  **********************************/
/* NGAS by Parton et al. (1996)  */
/*
Parton WJ, Mosier AR, Ojima DS, Valentine DW, Schimel DS, 
Weier K, Kulmala AE (1996) Generalized model for N2 and N2O 
production from nitrification and denitrification. 
Global Biogeochemical Cycles 10:401-412
*/
void f_n2o_emit_ngas(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Mass *mass,  
	struct Flux *flux
){
	float aa, bb, cc;
	float n_h2o, n_t, n_ph, n_nh4, day_n_n2o;
	float fd_wfps, fd_no3, fd_co2, dt;
	float fr_wfps, fr_no3, fr_co2;
	float day_d_n2, day_d_n2o;
	float nh4_soil;		/* micro g g-1*/ /* assumed */
	float no3_soil;						/* assumed */
	float kmax = 3.8;		/* g N ha-1 day-1 */ /* assumed */
	float nmax = 30.0;							/* assumed */
	float r_n2_n2o;
	
/*	if(WSITE == TKY){
		nh4_soil = 20.0;
		no3_soil = 8.8167;
	}else if(WSITE == QHB){
		nh4_soil = 20.0;
		no3_soil = 10.0;
	}else if(WSITE == FJY){
		nh4_soil = 5.0;
		no3_soil = 5.0;
	}else{
		nh4_soil = 20.0;
		no3_soil = 8.8167;
	}
	
	nh4_soil = 5.0;
	no3_soil = 2.0;  */
	
	/* prognostic soil inorganic N */
	nh4_soil = (mass->soil).n_nh4/10000.0*1000000.0 / (grid->bulkdens*1000.0*1000.0);
	no3_soil = (mass->soil).n_no3/10000.0*1000000.0 / (grid->bulkdens*1000.0*1000.0);
	
	/* NITRIFICATION *******************************************/
	/* water factor */
	/* Fig.(2a) in Parton et al. (1996) */
    n_h2o = 0.0;
	if(grid->stexture==0){
		/* sandy */
		aa = (loct->wfps - 1.70) / (0.55 - 1.70);
		bb = 3.22 * ((1.70 - 0.55) / (0.55 + 0.007));
		cc = (loct->wfps + 0.007) / (0.55 + 0.007);
		n_h2o = pow(aa, bb) * pow(cc, 3.22);
	}else if(grid->stexture==1||grid->stexture==2){
		/* medium and fine */
		aa = (loct->wfps - 1.27) / (0.60 - 1.27);
		bb = 2.84 * ((1.27 - 0.60) / (0.60 - 0.0012));
		cc = (loct->wfps - 0.0012) / (0.60 - 0.0012);
		n_h2o = pow(aa, bb) * pow(cc, 2.84);
	}
	
	/* temperature factor */
	/* Fig.(2b) in Parton et al. (1996) */
	n_t = -0.06 + 0.13 * exp(0.07*loct->tmp10_soil);
	
	/* pH factor */
	/* Fig.(2c) in Parton et al. (1996) */
	n_ph = 0.56 + atan(PI * 0.45 * (-5.0 + grid->soil_ph))/PI;
	
	/* ammonium factor*/
	/* Fig.(2d) in Parton et al. (1996) */
	n_nh4 = 1.0 - exp(-0.0105 * nh4_soil);
	
	/* N2O emission */
	/* Eq.(1) in Parton et al. (1996) */
	day_n_n2o = n_h2o * n_ph * n_t * (kmax + nmax*n_nh4);
	
	/* DENITRIFICATION **************************************/
	/* Fig.(3a) in Parton et al. (1996) */
    fd_wfps = 0.0;
	if(grid->stexture == 0){
		/* sandy */
		fd_wfps = 1.56 / pow(12.0, (16.0 / pow(12.0, 2.01*loct->wfps)));
	}else if(grid->stexture == 1){
		/* medium */
		fd_wfps = 4.82 / pow(14.0, (16.0 / pow(14.0, 1.39*loct->wfps)));
	}else if(grid->stexture == 2){
		/* fine */
		fd_wfps = 60.0 / pow(18.0, (22.0 / pow(18.0, 1.06*loct->wfps)));
	}
	
	/* Fig.(3b) in Parton et al. (1996) */
	fd_no3 = 11000.0 + (40000.0 * atan(PI * 0.002 * (no3_soil - 180.0)))/PI;
	/* Fig.(3c) in Parton et al. (1996) */
	fd_co2 = (24000.0 / (1.0 + 200.0 / exp(0.35 * (flux->soil).hr*1000.0))) - 100.0;
	
	/* Gas flux, Eq.(2) in Parton et al. (1996) */
	/* g N ha-1 day-1 */
	dt = fd_wfps * ((fd_no3>fd_co2)?fd_co2:fd_no3);
	
	/* Fig.(5a) in Parton et al. (1996) */
	fr_wfps = 1.4 / pow(13.0, 17.0/pow(13.0, 2.2*loct->wfps));
	/* Fig.(5b) in Parton et al. (1996) */
	fr_no3 = (1.0 - (0.5 + (1.0 * atan(PI * 0.01 * (no3_soil - 190.0))) / PI)) * 25.0;
	/* Fig.(5c) in Parton et al. (1996) */
	fr_co2 = 13.0 + (30.78 * atan(PI * 0.07 * ((flux->soil).hr*1000.0 - 13.0))) / PI;
	/*  change from flux->sr to (flux->soil).hr :070530 */
	
	r_n2_n2o = fr_wfps * ((fr_no3>fr_co2)?fr_co2:fr_no3);
	
	/* N2O emission */
	/* Eqs.(3+4) in Parton et al. (1996) */
	day_d_n2o = dt / (1.0 + r_n2_n2o);
	/* Eqs.(3+5) in Parton et al. (1996) */
	day_d_n2 = dt / (1.0 + 1.0/r_n2_n2o);
	
	/* total ***************************************************/
	/* g N ha-1 day-1 */
	/* mg N20 m-2 day-1 */
	(flux->soil).d_n2o_ntr_ngas = day_n_n2o*44.0/28.0*0.1;
	(flux->soil).d_n2o_dnt_ngas = day_d_n2o*44.0/28.0*0.1;
	(flux->soil).d_n2o_ngas = (flux->soil).d_n2o_ntr_ngas + (flux->soil).d_n2o_dnt_ngas;
	
	/* mg N2 m-2 day-1 */
	(flux->soil).d_n2_ngas = day_d_n2*0.1;

	/* nitrification */
	(flux->soil).n_nitrif = day_n_n2o * 75.0;
}

/*
Parton, W.J., Holland, E.A., Del Grosso, S.J., Hartman, M.D., Martin, R.E., Mosier, A.R., 
Ojima, D.S. and Schimel, D.S., 2001. Generalized model for NOx and N2O emissions from soils. 
Journal of Geophysical Research 106, 17403-17419.
*/
/*
Del Grosso, S.J., Parton, W.J., Mosier, A.R., Ojima, D.S., Kulmala, A.E. and Phongpan, S., 2000. 
General model for N2O and N2 gas emissions from soils due to denitrification. 
Global Biogeochemical Cycles 14, 1045-1060.
*/
void f_n2o_emit_ngas_2(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Mass *mass,  
	struct Flux *flux
){
	float aa, dd, ee, ff, gg, hh, ii;
	float n_ph, day_n_n2o;
	float fd_wfps, fd_no3, fd_co2, dt, l_m, x_i, co2_soil;
	float fr_wfps, fr_no3_co2, k1;
	float day_d_n2, day_d_n2o;
	float co2_emit;
	float nh4_soil;		/* micro g g-1*/ /* assumed */
	float no3_soil;						/* assumed */
	float r_n2_n2o;
	float d_n_min;		/* nitrogen gross mineralization, g N ha-1 day-1 */
	float fno3;			/* soil nitrification rate, g N m-2 d-1 */
	float n_k1 = 0.2;		/* fraction of d_n_min from Parton et al. (2001) */
	float f_nh4;			/* model-derived soil ammonium conc., g N m-2 */
	float kmax = 0.1;		/* maxmum fraction of nh4; nitrified, d-1 from Parton et al. (2001) */
	float f_t;			/* effect of soil temperature on nitrification form Fig.(3a) in Parton et al. (2001) */
	float f_wfps;			/* the effect of SWC and soil texture on nitrificcation form Fig.(3b) in Parton et al. (2001) */
	float fn2o;			/* N2O flux from nitrification, g N m-2 d-1*/
	float n_k2 = 0.02;		/* fraction of nitrified N lost as N2O flux */
	float r_nox;			/* ratio of NOx to N2O fluxes */
	float nox;			/* potential soil gas glux of NOx, g N ha-1 d-1 */
	float st_min, st_max, st_opt;
	float d_d0;		/* normalized diffusivity in aggregated media */
	float ss;			/* solid phase of the aggregated medium, m3 m-3 bed volume */
	float a2;			/* intra-aggregate pore space, m3 m-3 bed volume */
	float pp;		/* inter-aggregate pore space, m3 m-3 bed volume */
	float th_a;		/* volume of water per unit bed volume contained in intra-arregate pore space, m3 m-3 */
	float th_p;		/* volume of water per unit bed volume contained in inter-arregate pore space, m3 m-3 */
	float s_wa, s_wp;		/* fractional liquid saturation of aa and pp components of the total void volume */
	float xx, yy, zz;		/* exponent terms */
	float fc;		/* field capacity */
	float aaa, bbb, ccc, ddd, eee, fff, ggg, hhh, iii, jjj;
	/* float kmax = 3.8; */		/* g N ha-1 day-1 */ /* assumed */
	/* float k_reac;
	float norg_total; */
	
	/* assumed */
	/* if(WSITE == TKY){
		nh4_soil = 20.0;		
		no3_soil = 8.8167;	
	}else if(WSITE == QHB){
		nh4_soil = 20.0;	
		no3_soil = 10.0;	
	}else if(WSITE == FJY){
		nh4_soil = 5.0;		
		no3_soil = 5.0;		
	}else{
		nh4_soil = 20.0;	
		no3_soil = 8.8167;	
	}
	nh4_soil = 5.0;		
	no3_soil = 2.0;		*/
	
	/* prognostic soil inorganic N */
	/* micro g N/ g soil */
	nh4_soil = (mass->soil).n_nh4/10000.0*1000000.0 / (grid->bulkdens*1000.0*1000.0);
	no3_soil = (mass->soil).n_no3/10000.0*1000000.0 / (grid->bulkdens*1000.0*1000.0);
	
	/* CO2 conc. in soil, temporary value 360.0 (atmospheric CO2 conc.) */
	co2_soil = 360.0;		/* ppm */ 
						
	/* from CHEM96_Potter *************************************************/		
	fc = a2 = grid->field_cap/300.0;
	pp = grid->pore_cap/300.0 - a2;			/* inter-aggregate pore space */
	ss = 1.0 - grid->pore_cap/300.0;
	
	th_a = loct->m_casa;
    th_p = 0.0;
	if(loct->m_casa <= a2){
		th_p = 0.0;
	}else if(loct->m_casa>a2 && pp>loct->e_casa){
		th_p = loct->e_casa;
	}else if(loct->m_casa>a2 && pp<=loct->e_casa){
		th_p = pp;
	}
	
	s_wa = th_a / a2;
	s_wp = th_p / pp;
	
	aaa = pp - th_p;
	bbb = (a2 - th_a) / (a2 + ss);
	
	/* Eq. 4a in Potter (1996) */
	xx = 0.477 * pow(pp, 3.0) - 0.596 * pow(pp, 2.0) + 0.437 * pp + 0.564;
	/* Eq. 4b in Potter (1996) */
	yy =  0.477 * pow(aaa, 3.0) - 0.596 * pow(aaa, 2.0) + 0.437 * aaa + 0.564;
	/* Eq. 4c in Potter (1996) */
	zz = 0.477 * pow(bbb, 3.0) - 0.596 * pow(bbb, 2.0) + 0.437 * aaa + 0.564;
	
    aa = 0.0;
 
	/* Eq. 3 in Potter (1996) */
	fff = pow((1.0 - s_wa), 2.0);
	jjj = (a2 - th_a / aa + ss);
	ggg = pow(jjj, 2.0 * zz);
	hhh = (1.0 - pow(pp, 2.0 * xx));
	iii = (aaa - pow(aaa, 2.0 * yy));
	ccc = fff * ggg * hhh * iii;
	ddd = pow((1.0- s_wa), 2.0) * pow(jjj, 2.0) * (1.0 - pow(pp, 2.0 * xx)) + aaa - pow(aaa, 2.0 * yy);
	eee = pow((1.0- s_wp), 2.0) * pow(aaa, 2.0 * yy);
    d_d0 = 0.0;
	if(jjj>0){
		d_d0 = ccc / ddd + eee;
	}else if(jjj<=0){
		d_d0 = eee;
	}
	
	/* N2O from NITRIFICATION *******************************************/
	
	/* N mineralization rate ******/
	/* g N m-2 day-1 */
	
	/*  norg_total = (mass->soil).n_lttr + (mass->soil).n_hums;
	k_reac = 200000000.0 * exp(-16000.0/(1.987*(loct->tmp10_soil + ZAT)));
	d_n_min = norg_total * (1.0 - exp(-k_reac*1.0)) *10000.0 * pow(100.0, -2.0);  */
	
	d_n_min = ((flux->soil).n_minerlz_hums + (flux->soil).n_minerlz_lttr 
				- (flux->soil).n_immbl) / 10000.0;
	if(d_n_min < 0.0){
		d_n_min =0.0;
	}

	/* WFPS factor ******/
	/* Fig.(3b) in Parton et al. (2001) and Doran et al. (1988)*/
    f_wfps = 0.0;
	if(grid->stexture == 0){
		/* sandy */
		f_wfps = 5.44 * loct->wfps -5.03 * pow(loct->wfps, 2.0) - 0.492;
	}else if(grid->stexture==1 || grid->stexture==2){
		/* medium and fine */
		f_wfps = 5.63 * loct->wfps -4.64 * pow(loct->wfps, 2.0) - 0.745;
	}else if(grid->stexture == 3){
		/* volcanic */
		f_wfps = 5.76 * loct->wfps -3.70 * pow(loct->wfps, 2.0) - 1.280;
	}
	if(f_wfps < 0.0){
		f_wfps = 0.0;
	}
	
	/* temperature factor ******/
	/* Fig.(3a) in Parton et al. (2001)  and Fig.(2) in Malhi et al. (1982)*/
	if(grid->stmp10cm_av > 0.0){
		st_min = 9.98 * log(grid->stmp10cm_av) - 13.251;
		st_max = 6.4642 * log(grid->stmp10cm_av) + 34.215;
		st_opt = 6.6804 * log(grid->stmp10cm_av) + 13.926;
		hh = (loct->tmp10_soil - st_min) * (loct->tmp10_soil - st_max);
		ii = loct->tmp10_soil - st_opt;
		/* Fig.(A3) in Raich et al. (1991) */
		f_t = hh/(hh - pow(ii, 2.0));
	}else{
		f_t = 0.0;
	}
	if(f_t < 0.0){
		f_t =0.0;
	}
	
	/* pH factor ******/
	/* Fig.(2c) in Parton et al. (1996) */
	n_ph = 0.56 + atan(PI * 0.45 * (-5.0 + grid->soil_ph))/PI;
	if(n_ph < 0.0){
		n_ph =0.0;
	}	
	
	/* ammonium factor ******/
	f_nh4 = nh4_soil * pow(10.0, -6.0) * grid->bulkdens * 3.0 * (1.0/(pow(0.025, 2.0) * PI));
	if(f_nh4 < 0.0){
		f_nh4 =0.0;
	}	
	
	/* Eq.(2) in Parton et al. (2001) */
	/* g N /m2/day */
	fno3 = d_n_min * n_k1 * kmax * f_nh4 * f_t * f_wfps * n_ph;
	
	/* Eq.(3) in Parton et al. (2001) */
	fn2o = n_k2 * fno3;
	
	/* Eq.(5) in Parton et al. (2001) */
	r_nox = 15.2 + (35.5 * atan(0.68 * PI * (10.0 * d_d0 - 1.86))) / PI;
	
	/* N2O emission */
	/* Eq.(1) in Parton et al. (1996) */
	/* day_n_n2o = n_h2o * n_ph * n_t * (kmax + nmax*n_nh4); */
	day_n_n2o = fn2o * pow(100.0, 2.0);
	
	/* DENITRIFICATION **************************************/	
	/* WFPS factor */
	/* Fig. (6a) in DelGrosso et al. (2000) */
	l_m = ((0.113 < d_d0)?0.113:d_d0) * -3.05 + 0.36;
	/* Eq. (3) in DelGrosso et al. (2000) */
	x_i = 0.90 - l_m;
	/* Eq. (1) in Parton et al. (2001) */
	fd_wfps = 0.5* atan(0.6 * PI * (10.0 * loct->wfps - x_i))/PI;
	if(fd_wfps < 0.0){
		fd_wfps = 0.0;
	}
	
	/* nitrate factor */
	/* Fig. (1c) in Delgrosso et al. (2000) */
	fd_no3 = 1.15 * pow(no3_soil, 0.57);
	
	/* CO2 factor */
	/* Fig. (1d) in DelGrosso et al. (2000) */	
	/* micro g C / g soil / day */
	co2_emit = (flux->soil).hr *1000000.0/10000.0 *1000000.0 / (grid->bulkdens*1000.0*1000.0);
	fd_co2 = 0.1 * pow(co2_emit, 1.3);
	
	/* Gas flux, Eq.(2) in Parton et al. (1996) & Eq.(2) in Del Grosso et al. (2000)*/
	/* g N ha-1 day-1 */
	dt = fd_wfps * ((fd_no3>fd_co2)?fd_co2:fd_no3);

	/* Fig. (5h) in DelGrosso et al. (2000) */
	dd = 0.015 * loct->wfps*100.0 - 0.32;
	fr_wfps = (dd > 0.1)?dd:0.1;
	
	/**/
	/* Fig. (6b) in DelGrosso et al. (2000) */
	ee = 38.4 - 350.0 * d_d0;
	k1 = (ee > 1.7)?ee:1.7;


	
	
		
	
	
	/* Fig. (7c) in DelGrosso et al. (2000) */
	ff = 0.16 * k1;
	gg = k1 * exp(-0.8 * (no3_soil/((flux->soil).hr * 1000.0)));
	fr_no3_co2 = (ff>gg)?gg:ff;
	
	/*  change from flux->sr to (flux->soil).hr :070530   */
	
	/*r_n2_n2o = fr_wfps * ((fr_no3>fr_co2)?fr_co2:fr_no3); */
	
	/* Eq. (4) in DelGrosso et al. (2000) */
	r_n2_n2o = fr_wfps * fr_no3_co2;
	
	/* N2O emission */
	/* Eqs.(3+4) in Parton et al. (1996) */
	day_d_n2o = dt / (1.0 + r_n2_n2o);
	/* Eqs.(3+5) in Parton et al. (1996) */
	day_d_n2 = dt / (1.0 + 1.0/r_n2_n2o);
	
	/* g N ha-1 d-1 */
	nox = r_nox * day_d_n2o + r_nox * day_n_n2o /* * pls */;
	
	/* chemodenitrification from Table.5 in Li et al. (2000) ***********************/
	/*dndc_ft = 0.03 * loct->tmp10_soil + 0.2;
	dndc_fph = 2236 * exp(-2.5 * grid-> ph_soil);
	dndc_rn = fno3 * pow(10, -3) * pow(100, 2); */
	/* kg N ha-1 d-1 */
	/*dndc_rchem = dndc_a * dndc_rn * dndc_ft * dndc_fph; */
	
	/* total ***************************************************/
	/* g N ha-1 day-1 */
	/* mg N20 m-2 day-1 */
	(flux->soil).d_n2o_ntr_ngas_2 = day_n_n2o*44.0/28.0*0.1;
	(flux->soil).d_n2o_dnt_ngas_2 = day_d_n2o*44.0/28.0*0.1;
	(flux->soil).d_n2o_ngas_2 = (flux->soil).d_n2o_ntr_ngas_2 + (flux->soil).d_n2o_dnt_ngas_2;
	
	/* mg N2 m-2 day-1 */
	(flux->soil).d_n2_ngas_2 = day_d_n2*0.1;
}

/* Daily step CASA nitrogen trace gas emission from soil ****************************/
/* Potter et al. (1996) */
/*
Potter CS, Matson PA, Vitousek PM, Davidson EA (1996) 
Process modeling of controls on nitrogen trace gas emissions from soils 
worldwide. Journal of Geophysical Research 101:1361-1377
*/
void f_n2o_emit_casa(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Mass *mass,  
	struct Flux *flux
){
	float d_n_min;		/* nitrogen gross mineralization, g N ha-1 day-1 */
	float d_no, d_n2o, d_n2;		/* soil trace gas emission, g N ha-1 day-1 */
	float f_emit;		/* fraction of gas emission per mineralization */
	float n_total;
	
	n_total = 300.0;
	f_emit = 0.01;
	
	d_no = d_n2o = 0.0;
		
	/* N mineralization rate */
	/* g N ha-1 day-1 */
	/*  d_n_min = ((loct->d_c_min / (mass->soil).soil) * (mass->soil).n_soil) * pow(10, 6);  */
	
	/*  k_reac = 200000000.0 * exp(-16000.0/(1.987*(loct->tmp10_soil + ZAT)));
	d_n_min = n_total * (1.0 - exp(-k_reac*1.0)) *10000.0;  */
	
	d_n_min = (flux->soil).n_minerlz_hums + (flux->soil).n_minerlz_lttr 
				- (flux->soil).n_immbl;
				
	if(d_n_min < 0.0){
		d_n_min =0.0;
	}
	
	/* Fig.3 of Potter et al. (1996) */
    d_no = d_n2o = d_n2 = 0.0;
	if(loct->iw_casa < 10.0){
		d_no = 0.1 * loct->iw_casa * d_n_min * f_emit;
		d_n2o = 0.0;
		d_n2 = 0.0;
	}else if(10.0<=loct->iw_casa && loct->iw_casa<30.0){
		d_no = (-0.005 * loct->iw_casa + 1.05) * d_n_min * f_emit;
		d_n2o = (0.005 * loct->iw_casa - 0.05) * d_n_min * f_emit;
		d_n2 = 0.0;
	}else if(30.0<=loct->iw_casa && loct->iw_casa<60.0){
		d_no = (-0.01 * loct->iw_casa + 1.2) * d_n_min * f_emit;
		d_n2o = (0.01 * loct->iw_casa - 0.2) * d_n_min * f_emit;
		d_n2 = 0.0;
	}else if(60.0<=loct->iw_casa && loct->iw_casa<70.0){
		d_no = (-0.04 * loct->iw_casa + 3.0) * d_n_min * f_emit;
		d_n2o = (0.04 * loct->iw_casa - 2.0) * d_n_min * f_emit;
		d_n2 = 0.0;
	}else if(70.0<=loct->iw_casa && loct->iw_casa<80.0){
		d_no = (-0.02 * loct->iw_casa + 1.6) * d_n_min * f_emit;
		d_n2o = (0.02 * loct->iw_casa - 0.6) * d_n_min * f_emit;
		d_n2 = 0.0;
	}else if(80.0<=loct->iw_casa && loct->iw_casa<90.0){
		d_no = 0.0;
		d_n2o = (-0.02 * loct->iw_casa + 2.6) * d_n_min * f_emit;
		d_n2 = (0.02 * loct->iw_casa - 1.6) * d_n_min * f_emit;
	}else if(90.0<=loct->iw_casa){
		d_no = 0.0;
		d_n2o = (-0.08 * loct->iw_casa + 8.0) * d_n_min * f_emit;
		d_n2 = (0.08 * loct->iw_casa - 7.0) * d_n_min * f_emit;
	}
	
	/* mg NO m-2 day-1 */
	(flux->soil).d_no_casa = d_no*30.0/14.0*0.1;
	/* mg N2 m-2 day-1 */
	(flux->soil).d_n2_casa = d_n2*0.1;

	/* mg N20 m-2 day-1 */
	(flux->soil).d_n2o_casa = d_n2o*44.0/28.0*0.1;
}

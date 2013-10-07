/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/** This file gives description of functions releted to photosynthesis, or GPP **/
#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

#define DF_STEP 24
/* leaf boundary conductance */
#define GB 250.0

/* formula of hourly gross primary production *******************/
double f_gpp_df97(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *pmas, 
	struct Pflx *pflx
){
	short h;
	double aa, bb, cc, jj, hangle, ge, dtc, co2_a, o2_i;
	double f_cloud, lai_t, kt, hd, e2p_d, e2p_b;
	double h_sinh[DF_STEP], toprad[DF_STEP], sfcrad[DF_STEP];
	double par_b[DF_STEP], par_d[DF_STEP], ppfd_b[DF_STEP], ppfd_d[DF_STEP];
	double ke_b1, ke_b2, ke_d, lai_sn, lai_sd, scttr, irr_b, rfl_b, rfl_d;
	double appfd, appfd_sn1, appfd_sn2, appfd_sn3, appfd_sd1, appfd_sd2;
	double appfd_sn, appfd_sd, ps2_appfd_sn[DF_STEP], ps2_appfd_sd[DF_STEP], spect;
	double comp_co2, temp, kc, ko, k_effc, ft_et, ft_rd;
	double n_photocap, n_top, n_nonphoto, vcmax0, vcmax, jmax;
	double j_sn, j_sd, assim_sn, rdark_sn, stocon_sn, inco2_sn;
	double assim_sd, rdark_sd, stocon_sd, inco2_sd;
	double debug1, debug2, debug3, f_ds, f_ds1, f_ds2, root1, root2, gpp_df97;
	double monitor1[DF_STEP], monitor2[DF_STEP], monitor3[DF_STEP], monitor4[DF_STEP], monitor5[DF_STEP];
	
	ge = 2.0*PI/365.0*loct->doy;
	aa = 0.034221*cos(ge) + 0.00128*sin(ge);
	bb = 0.000719*cos(2.0*ge) + 0.000077*sin(2.0*ge);
	dtc = 1.00011 + aa + bb;
	
	f_cloud = loct->tcdc_clm;
	jj = 0.803 - 0.34*f_cloud - 0.458*f_cloud*f_cloud;
	jj = (jj<=1.0)?jj:1.0; 
	jj = (jj>=0.0)?jj:0.0;
	
	lai_t = pmas->lai;
	scttr = 0.15;
	rfl_d = 0.036;
	ke_d = 0.719;
	spect = 0.15;
	temp = loct->tmp_sfc;
	n_photocap = 1.16;
	
	/*  if(pchar->c_layer == 0){
		n_top = 70.0;
	}else if(pchar->c_layer == 1){
		n_top = 50.0;
	} */
/*	n_top = pchar->n_conc_larea; ******************************************************************************/
	n_top = 50.0;
	n_nonphoto = 25.0;
		
	/* canopy CO2 concentration */
	co2_a = loct->aCO2/10.0;	/* ppmv -> Pa */
	/* O2 concentration, 21% of atmospheric pressure, Pa */
	o2_i = 0.21 * loct->air_prsr*100.0;
	
	/* CO2 compensation point, Pa */
	/* Brooks & Farquhar (1985) */
	/* Eq.9 in De Pury and Farquhar (1997) */
	comp_co2 = 3.69 + 0.188*(temp-25.0) + 0.0036*(temp-25.0)*(temp-25.0);
	
	/* Michaelis constant of carboxylation and oxygenation rates */
	/* Eq.8 in De Pury and Farquhar (1997) */
	/* Rubisco carboxylation */
	kc = 40.4*exp(59400.0*(temp-25.0)/(298.15*UGC*(temp+ZAT)));
	/* Rubisco oxygenation */
	ko = 24800.0*exp(36000.0*(temp-25.0)/(298.15*UGC*(temp+ZAT)));
	/* effective Michaelis constant of Rubisco */
	k_effc = kc*(1.0 + o2_i / ko);
	
	/* temperature coefficient of electron transport */
	/* Eq.10 in De Pury and Farquhar (1997) */
	ft_et = exp(((temp+ZAT)-298.15)*37000.0/(UGC*(temp+ZAT)*298.15))*
		(1.0 + exp((710.0*298.15 - 220000.0)/(UGC*298.15)))/
		(1.0 + exp((710.0*(temp+ZAT) - 220000.0)/(UGC*(temp+ZAT))));
	
	/* temperature coefficient of dark respiration */
	ft_rd = exp(66400.0*(temp-25.0)/(UGC*298.0*(ZAT+temp)));
	
	/* drought stress scholar *************/
	root1 = 1.0 - 0.5 * (exp(-0.3 * pchar->root_dpt_a) + exp(-0.3 * pchar->root_dpt_b));
	root2 = 1.0 - root1;
	
	if(pchar->phototype==3){			/* C3 plants */
		f_ds1 = 0.95*loct->soilwtr_l / (loct->soilwtr_l + grid->fieldcap30*0.5) + 0.05;
		f_ds2 = 0.95*loct->soilwtr_h / (loct->soilwtr_h + grid->field_cap*0.5) + 0.05;
	}else if(pchar->phototype==4){		/* C4 plants */
		f_ds1 = 0.86*loct->soilwtr_l / (loct->soilwtr_l + grid->fieldcap30*0.5) + 0.14;
		f_ds2 = 0.86*loct->soilwtr_h / (loct->soilwtr_h + grid->field_cap*0.5) + 0.14;
	}
	f_ds = root1*f_ds1 + root2*f_ds2;
	f_ds = (f_ds<=1.0)?f_ds:1.0; 
	f_ds = (f_ds>=0.0)?f_ds:0.0;
	
	gpp_df97 = 0.0;
	for(h=0;h<DF_STEP;h++){
		hangle = -180.0+(double)h* (360.0/(double)DF_STEP);
		
		/* solar angle */
		h_sinh[h] = sin(grid->lat*dTr)*sin(loct->soldec[loct->doy]*dTr) 
			+cos(grid->lat*dTr)*cos(loct->soldec[loct->doy]*dTr)*cos(hangle*dTr);
		h_sinh[h] = (h_sinh[h]<=1.0)?h_sinh[h]:1.0; 
		h_sinh[h] = (h_sinh[h]>=-1.0)?h_sinh[h]:-1.0;
		
		/* hourly atmosphere-top radiation, W m-2 */
		toprad[h] = 1367.0 * dtc * h_sinh[h];
		toprad[h] = (toprad[h]>0.0)?toprad[h]:0.0;
		
		/* surface radiation, W m-2 */
		sfcrad[h] = jj * toprad[h];
		
		if(sfcrad[h]>0.0){
			kt = sfcrad[h] / toprad[h];
			hd = sfcrad[h] * (0.958 - 0.982 * kt);
			hd = (hd>0.0)?hd:0.0;
			
			/* PAR, W m-2 */
			par_b[h] = 0.57 * hd;
			par_d[h] = 0.43 * (sfcrad[h] - hd);
			
			/* PAR => PPFD by Dye parameterization */
			cc = 0.958 - 0.982 * kt;
			cc = (cc>0.01)?cc:0.01;
			e2p_d = 4.5886 * cc / (0.010773 + cc);
			e2p_d = (e2p_d>4.2)?e2p_d:4.2;
			e2p_b = 4.576 - 0.033144*cc;
			ppfd_d[h] = e2p_d * par_d[h];
			ppfd_b[h] = e2p_b * par_b[h];
			
			if(pchar->c_layer == 1){
				ppfd_d[h] *= loct->grdrad_h[h*2];
				ppfd_b[h] *= loct->grdrad_h[h*2];
			}
			
			/* sun/shade canopy */
			ke_b1 = 0.5 / h_sinh[h];
			ke_b2 = 0.46 / h_sinh[h];
			
			/* Eq.18 of DF97 */
			lai_sn = (1.0 - exp(-ke_b1 * lai_t))/ke_b1;
			if(lai_sn > lai_t){
				lai_sn = lai_t;
			}
			/* Eq.21 of DF97 */
			lai_sd =lai_t - lai_sn;
			
			irr_b = (1.0 - sqrt(1.0 - scttr))/(1.0 + sqrt(1.0 - scttr));
			rfl_b = 1.0-exp(-2.0*irr_b*ke_b1 / (1.0 + ke_b1));
			
			appfd = (1.0 - rfl_b)*ppfd_b[h] * (1.0-exp(-ke_b2 * lai_t))+
					(1.0 - rfl_d)*ppfd_d[h] * (1.0-exp(-ke_d * lai_t));
			
			/* absorbed PPFD by sunny leaves */
			/* beam: Eq.20b of DF97 */
			appfd_sn1 = ppfd_b[h]*(1.0 - scttr)*(1.0-exp(-ke_b1 * lai_t));
			/* diffused: Eq.20c of DF97 */
			appfd_sn2 = ppfd_d[h]*(1.0 - rfl_d)*(1.0-exp(-(ke_d + ke_b1)*lai_t))*
				ke_d/(ke_d + ke_b1);
			/* beam scattered within canopy: Eq.20d of DF97 */
			appfd_sn3 = ppfd_b[h]*((1.0 - rfl_b)*(1.0 - exp(-(ke_b1 + ke_b2)*lai_t))*
				ke_b2/(ke_b1 + ke_b2)-(1.0 - scttr)*(1.0 - exp(-2.0*ke_b1*lai_t))/2.0);
			/* absorbed PPFD by sunny leaves: Eq.20a of DF97 */
			appfd_sn = appfd_sn1 + appfd_sn2 + appfd_sn3;
			/* absorbed PPFD by sunny leaves PS II: Eq.6 of DF97 */
			/* see p.35 of Von Caemmerer (2000) */
			ps2_appfd_sn[h] = appfd_sn*(1.0 - spect)/2.0;
			/* cchar->spect : spectral coefficinet for chlorophyll absorption */
			/* 2.0 : half separation between PS I and PS II */
			
			/* absorbed PPFD by shaded leaves */
			/* diffused: Eq.A26b of DF97 */
			appfd_sd1 = ppfd_d[h]*(1.0 - rfl_d)*(1.0-exp(-ke_d * lai_t)-
				(1.0-exp(-(ke_d + ke_b1)*lai_t))*ke_d/(ke_d + ke_b1));
			/* beam scattered within canopy: Eq.A26c of DF97 */
			appfd_sd2 = ppfd_b[h]*((1.0 - rfl_b)*(1.0-exp(-ke_b2 * lai_t)-
				(1.0 - exp(-(ke_b1 + ke_b2)*lai_t))*ke_b2/(ke_b1 + ke_b2))-
				(1.0 - scttr)*(1.0-exp(-ke_b1*lai_t)-(1.0 - exp(-2.0 * ke_b1 * lai_t))/2.0));
			/* absorbed PPFD by shaded leaves: Eq.A26a of DF97 */
			appfd_sd = appfd_sd1 + appfd_sd2;
			/* absorbed PPFD by shaded leaves PS II */
			ps2_appfd_sd[h] = appfd_sd*(1.0 - spect)/2.0; /* Eq.6 of DF97 */
			
			/* sunny leaves ***************/
			/* maximum carboxylation rate at 25 deg-C: Eq.22 in DF97 */
			vcmax0 = lai_t*n_photocap*(n_top - n_nonphoto)*
				(1.0-exp(-(0.713 + ke_b1*lai_t)))/(0.713 + ke_b1*lai_t);
			/* maximum carboxylation rate: Eq.8 in DF97 */
			vcmax = vcmax0 * exp(64800.0*(temp-25.0)/(298.15*UGC*(temp+ZAT)));
			/* maximum electron transport rate: proportional Jmax/Vcmax relationship */
			jmax =2.1 * vcmax0 * ft_et;
			
			/* drought stress */
			vcmax *= f_ds;
			jmax *= f_ds;
			
			/* electron transport rate: Eq.5 in DF97 */
			bb = ps2_appfd_sn[h] + jmax;
			j_sn = (bb - sqrt(bb*bb - 4.0*0.7*ps2_appfd_sn[h]*jmax))/(2.0*0.7);
			
			leaf_ansolv_df97(loct->vpd, co2_a, lai_sn, ft_rd, 
				vcmax, k_effc, comp_co2, j_sn, 
				&assim_sn, &rdark_sn, &stocon_sn, &inco2_sn, 
				&debug1, &debug2, &debug3);
			
			monitor1[h] = assim_sn;
			monitor2[h] = stocon_sn;
			monitor3[h] = inco2_sn;
			monitor4[h] = debug1;
			monitor5[h] = debug2;

			/* shade leaves ****************/
			/* maximum carboxylation rate at 25 deg-C: Eq.A27 in DF97 */
			/* Vcmax(shade) = Vcmax - Vcmax(sun) */
			vcmax0 = lai_t*n_photocap*(n_top - n_nonphoto)*((1.0-exp(-0.713))/
				0.713 - (1.0-exp(-(0.713 + ke_b1*lai_t)))/(0.713 + ke_b1*lai_t));
			/* maximum carboxylation rate: Eq.8 in DF97 */
			vcmax = vcmax0*exp(64800.0*(temp-25.0)/(298.15*UGC*(temp+ZAT)));
			/* maximum electron transport rate */
			jmax = 2.1 * vcmax0 * ft_et; 

			/* drought stress */
			vcmax *= f_ds;
			jmax *= f_ds;

			/* electron transport rate: Eq.5 in DF97 */
			bb = ps2_appfd_sd[h] + jmax;
			j_sd = (bb - sqrt(bb*bb - 4.0*0.7*ps2_appfd_sd[h]*jmax))/(2.0*0.7);
			
			leaf_ansolv_df97(loct->vpd, co2_a, lai_sd, ft_rd, 
				vcmax, k_effc, comp_co2, j_sd, 
				&assim_sd, &rdark_sd, &stocon_sd, &inco2_sd, 
				&debug1, &debug2, &debug3);
			
			/* total CO2 assimilation */
			gpp_df97 += (assim_sn + assim_sd) * 3600.0 * 12.0 / 100000000.0;
			
		/*	monitor1[h] += vcmax;
			monitor2[h] += jmax;
			monitor3[h] += lai_sd; */
		}else{
			gpp_df97 += 0.0;
			monitor1[h] = monitor2[h] = monitor3[h] = monitor4[h] = monitor5[h] = 0.0;
		}
	}
	
	loct->x1 = monitor1[12];
	loct->x2 = monitor2[12];
	loct->x3 = monitor3[12];
	
	return (gpp_df97);
}

/**************************************************************************/
/*
Baldocchi, D. (1994), An analytical solution for coupled leaf photosynthesis 
and stomatal conductance models, Tree Physiology, 14, 1069-1079.
*/
void leaf_ansolv_df97(
	double vpd, 
	double aco2, 
	double lai, 
	double ft_rd, 
	double vcmax, 
	double k_effc, 
	double comp_co2, 
	double j_et, 
	double *assim,
	double *rdark,
	double *stocon,
	double *inco2, 
	double *debug1, 
	double *debug2,
	double *debug3
){
 	float aa, bb, dd, ee, fgc;
 	float ca, cs, gb, bd, rh, m;
 	float alpha, beta, gamma, theta;
 	float pp, qq, rr, qqq, rrr;
 	float xx1, xx2, xx3;
 	float f_pp, f_qq, f_cube, rd, ci, gs, pc, pj, p_gross, a_net;
 	 	
 	/**********************************************************************/
 	gb = GB;					/* leaf boundary layer conductance*/
 	bd = 5.0;					/* minimum conductance */
 	rh = 1.0/(1.0 + vpd/5.0);	/* relative humidity */
 	m = 5000.0/10.0;			/* coefficient of slope */
 	ca = aco2;					/* atmospheric CO2, Pa */
 	
	/* for sun leaf */
	if(lai > 1.0){
		fgc = lai;
	}else{
		fgc = 1.0;
	}
	rd = (0.0089*vcmax)*ft_rd*fgc;

	gb *= fgc/100.0;
	bd *= fgc/100.0;
	m *= fgc/100.0;
 	
 	/* Rubisco-limited ***************************************************/
	aa = vcmax;
	bb = k_effc;
	dd = comp_co2;
	ee = 1.0;
	
	/******************************************************************/
 	alpha = 1.0+bd/gb - m*rh;
 	beta = ca * (gb * m * rh - 2.0*bd - gb);
 	gamma = ca*ca * bd * gb;
 	theta = gb * m * rh - bd;
 	
 	/******************************************************************/
	pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*rd)/(ee * alpha);
	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*rd*beta 
		+ rd*bb*theta)/(ee * alpha);
	rr = (-aa*gamma + aa*dd*gamma/ca + ee*rd*gamma 
		+ rd*bb*gamma/ca)/(ee * alpha);
 	f_pp = (pp/3.0)*(pp/3.0) - qq/3.0;
 	f_qq = ((pp/3.0)*qq - 2.0*(pp/3.0)*(pp/3.0)*(pp/3.0) - rr)/2.0;
 	f_cube = f_qq*f_qq - f_pp*f_pp*f_pp;
 	
 	if(f_cube<0.0){
	 	/***************************************************************/
		qqq = (pp*pp - 3.0*qq)/9.0;
		rrr = (2.0*pp*pp*pp - 9.0*pp*qq + 27.0*rr)/54.0;
		theta = acos(rrr/sqrt(qqq*qqq*qqq));
		 	
	 	xx1 = -2.0*sqrt(qqq) * cos(theta/3.0) - pp/3.0;
	 	xx2 = -2.0*sqrt(qqq) * cos((theta+2.0*PI)/3.0) - pp/3.0;
	 	xx3 = -2.0*sqrt(qqq) * cos((theta+4.0*PI)/3.0) - pp/3.0; 	/**/
	 		 	
		pc = xx3 + rd;
	 }else{
		pc = 0.0;
	 }

 	/* RUBP-limited ***************************************************/
	aa = j_et;
 	bb = 8.0*comp_co2;
 	dd = comp_co2;
 	ee = 4.0;
 	
	/******************************************************************/
 	alpha = 1.0+bd/gb - m*rh;
 	beta = ca * (gb * m * rh - 2.0*bd - gb);
 	gamma = ca*ca * bd * gb;
 	theta = gb * m * rh - bd;
 	
 	/******************************************************************/
	pp = (ee*beta + bb*theta - aa*alpha + ee*alpha*rd)/(ee * alpha);
	qq = (ee*gamma + bb*gamma/ca - aa*beta + aa*dd*theta + ee*rd*beta 
		+ rd*bb*theta)/(ee * alpha);
	rr = (-aa*gamma + aa*dd*gamma/ca + ee*rd*gamma 
		+ rd*bb*gamma/ca)/(ee * alpha);

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
	 	
		pj = xx3 + rd;
	 }else{
		pj = 0.0;
	 }

 	/* Solution ***********************************************************/
	/* cchar->p_sn = ((cchar->pc_sn+cchar->pj_sn)-sqrt((cchar->pc_sn+cchar->pj_sn)*(cchar->pc_sn+cchar->pj_sn)-
				4.0*cchar->cv_p*cchar->pc_sn*cchar->pj_sn))/(2.0*cchar->cv_p); */
				
	p_gross = (pc < pj)?pc:pj;
	
	if(p_gross>=0.0 && p_gross<=100.0){
		;
	}else{
		p_gross = 0.0;
		/* printf("!! Warning: [p_sn, %f] in [leaf_ansolv()] !!\n", cchar->p_sn); */
		/* printf("!! Warning: [p_sn, %f %f %f] in [leaf_ansolv()] !!\n", cchar->p_sn, cchar->pc_sn, cchar->pj_sn); */
	}
	
	a_net = p_gross - rd;
	
	/* substomatal cavity CO2 */
	cs = ca - a_net/gb;
	
	gs = bd + m*a_net*rh/cs;		
	ci = cs - a_net/gs;
	
	*assim = p_gross;
	*rdark = rd;
	*stocon = gs;
	*inco2 = ci;
	
	*debug1 = pc;
	*debug2 = pj;
}


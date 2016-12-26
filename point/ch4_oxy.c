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

/* CH4 oxydation by Ridgewell et al. (1999) **********************************/
/*
Ridgwell AJ, Marshall SJ, Gregson K (1999) 
Consumption of atmospheric methane by soils: A process-based model. 
Global Biogeochemical Cycles 13:59-70
*/
void f_ch4oxy_ridgewell(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Mass *mass,  
	struct Flux *flux
){
	double r_sm;
	double r_t;
	double r_n;
	double k_d;
	double b;		
	double co_ch4 ;	/* ppmv */
	double zd = 30.0;		/* depth, cm */
	double f = 616.9;		/* mg CH4 ppmv-1 cm-1 */
	double pd = 2.65;		/* particle density, mg cm-3 */
	double phi;			/* total pore volume, cm3 cm-3 */
	double e;				/* air-filled porosity, cm3 cm-3 */
	double g_soil;
	double g_t;
	double d0_ch4 = 0.196;
	double d_ch4;			/* diffusivity of CH4 in free air, STP of 0.196 cm2 s-1 */
	double k0 = 0.00087;	/* base oxydation rate constant for uncultivated moist soil at 0 deg C */
	double lcult = 0.0;		/* temporary, fractional intensity of cultivation by JCAM83_Matthews */
	double d_j_ch4;			/*  */
	
	/* empirical constant by Saxton et al. (1986) */
	/* b = -5.311792308 */
	b = grid->b_sw;
	
	phi = 1.0 - (grid->bulkdens/pd);		/* total pore volume, cm3 cm-3 */
	if(NOTICE==1 && phi<=0.0){
		printf("!! ERROR: negative phi (ch4_oxy.c): %lf\n", phi);
	}
	
	/* co_ch4 = 1.72; */ /* default */
	co_ch4 = loct->atm_ch4_a1[loct->CO2y-1750] / 1000.0;
	
	/* moisture factor */
	if(((loct->prate_sfc + mass->sw30)/loct->pet) > 1.0){
		/* Eq. (11a) in Ridgewell et al. (1999) */
		r_sm = 1.0;
	}else{
		/* Eq. (11b) in Ridgewell et al. (1999) */
		r_sm = (loct->prate_sfc + mass->sw30)/loct->pet;
	}
	/* temperature factor */
	if(loct->tmp10_soil < 0.0){
		/* Eq. (9a) in Ridgewell et al. (1999) */
		r_t = 0.0;
	}else{
		/* Eq. (9b) in Ridgewell et al. (1999) */
		r_t = exp(0.0693 * loct->tmp10_soil - 8.56 * pow(10.0, -7.0) * pow(loct->tmp10_soil, 4.0));
	}
	/* nitrogen factor */
	/* Eq. (10) in Ridgewell et al. (1999) */
	r_n = 1.0 - (0.75 * lcult);
	
	/* Eq. (8) in Ridgewell et al. (1999) */
	k_d = k0 * r_sm * r_t * r_n;

	e = (phi - ((mass->sw30/10.0)/30.0));
	if(e<0.0){
		e = 0.0;
	}
	/* Eq. (6) in Ridgewell et al. (1999) */
	g_soil = pow(phi, 4.0/3.0) * pow(e/phi, (1.5+3.0/b));
	/* Eq. (7) in Ridgewell et al. (1999) */
	g_t = 1.0 + 0.0055 * loct->tmp10_soil;
	
	/* Eq. (5) in Ridgewell et al. (1999) */
	d_ch4 = g_soil * g_t * d0_ch4;
	/* Eq. (4) in Ridgewell et al. (1999) */
	if((d_ch4 + k_d*zd)>0.0){
		d_j_ch4 = ((co_ch4 * d_ch4)/ zd) * (1.0 - (d_ch4/(d_ch4 + k_d*zd))) * f;
	}else{
		d_j_ch4 = ((co_ch4 * d_ch4)/ zd) * f;
	}
	
	/* mg CH4 m-2 day-1 */
	(flux->soil).ch4oxy_ridg = d_j_ch4;
}

/* CH4 oxydation by Curry (2007) **********************************/
/*
	Curry, C.L., 2007. Modeling the soil consumption of atmospheric methane at the global scale. 
	Global Biogeochemical Cycles 21, 10.1029/2006GB002818.
*/
void f_ch4oxy_curry(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass,  
	struct Flux *flux
){
	double f_i;					/* fractional inundation */
	double r_w;					/* fraction of wetland */
	double g_t;					/* function of soil temperature */
	double phi;					/* total porosity, cm3 cm-3 in R99*/
	double pd = 2.65;			/* particle density, mg cm-3 */
	double phi_air;				/* air-filled porosity, cm3 cm-3 in R99*/
	double b;
	double g_soil;				/*function of soil moisuture content */
	double d_soil;				/* diffusion coefficient in soil, cm2 s-1 */
	double d_air = 0.196;		/* cm2 s-1 (STP) */
	double r_t;					/* temperature factor */
	double ps_sat;				/* saturation soil water potential, kPa **********/
	double theta;
	double frac_water;			/* fractional water content */
	double frac_ice;			/* fractional ice content **************/
	double ps;					/* soil water potential, absolute value ************/
	double beta = 0.8;			/* arbitary constant, see Sec 2.3 */
	double r_sm;				/* soil moisture factor */
	double k;					/* first-order oxidation rate constant, s-1 */
	double k_0 = 5.0 * pow(10, -5);		/* s-1, see Sec 2.3 */
	double j_0;					/* flux at the surface, mg CH4 m-2 d-1 */
	double g_0 = 586.7;			/* a factor, mg CH4 ppmv-1 sd-1 m2 cm */
	double c_0 = 1.753;			/* CH4 concentraton at surface, ppmv */
	double r_c = 1.0;			/* fraction of agricultural regions, non-agricultural = 1.0 */
	
	double aaa, bbb;
	
	c_0 = loct->atm_ch4_a1[loct->CO2y-1750]/1000.0;

	/* upland fraction */
	f_i = 0.0;			/* inundation */
	r_w = 1.0 - f_i;	/* non-wetland */
	frac_water = 1.0;
	frac_ice = 0.0;
	
	/***************************************************************/
	/* temperature factor */
	/* eq. 5 in Curry (2007) */
	g_t = 1.0 + 0.0055 * loct->tmp10_soil;
	
	/* total porosity, cm3 cm-3 */
	phi = 1.0 - (grid->bulkdens/pd);		
	/* theta = frac_water + frac_ice;	*/
	theta = (mass->sw30/10.0)/30.0;
	phi_air = phi - theta;
	if(phi_air<0.0){
		phi_air = 0.0;
	}
	/* soil water content - potential curve parameter */
	b = 15.9 * grid->clay_frac + 2.91;
	/* eq. 5 in Curry (2007) */
	g_soil = pow(phi, 4.0/3.0) * pow((phi_air / phi), (1.5 + 3.0/b));
	
	/* diffusion coefficient in soil, cm2 s-1 ****/
	/* eq. 4 in Curry (2007) */
	d_soil = d_air * g_t * g_soil;
	
	/***************************************************************/
	/* temperature factor */
	/* eq. 7 in Curry (2007) */
	if(loct->tmp10_soil < -10.0 || loct->tmp10_soil >= 43.3){
		r_t = 0.0;
	}else if(-10.0 <= loct->tmp10_soil && loct->tmp10_soil < 0.0){
		r_t = pow((0.1 * loct->tmp10_soil +1.0), 2);
	}else if(0.0 <= loct->tmp10_soil && loct->tmp10_soil < 43.3){
		r_t = exp(0.0693 * loct->tmp10_soil - 8.65 * pow(10, -7) * pow(loct->tmp10_soil, 4));
	}
	
	/* saturation water potential, kPa */
	ps_sat = pow(10.0, -1.31 * grid->sand_frac + 1.88) / 10.0;
	/* eq. 8 in Curry (2007) */
	/* b is replaced by Saxton-derived one ? */
	ps = ps_sat * pow((frac_water / phi), -b);
	
	/* eq. 9 in Curry (2007) */
	if(ps < 0.2){
		r_sm = 1.0;
	}else if (0.2 >= ps && ps <= 100.0){
		aaa = log10(ps) - log10(0.2);
		bbb = log10(100.0) - log10(0.2);
		r_sm = pow((1.0 - aaa/bbb), beta);
	}else if(100.0 < ps){
		r_sm = 0.0;
	}
	
	/* first-order oxidation constant */
	/* eq. 6 in Curry (2007) */
	k = k_0 * r_t * r_sm;
	
	/* surface Ch4 flux */
	/* eq. 10 in Curry (2007) */
	j_0 = g_0 * c_0 * r_c * r_w * pow((d_soil * k), 0.5);
	/* mg CH4 m-2 day-1 */
	
	(flux->soil).ch4oxy_curry = j_0;
}

/* CH4 oxydation by Potter et al. (1996) **********************************/
/*  
Potter CS, Davidson EA, Verchot LV (1996) Estimation of global 
biogeochemical controls and seasonality in soil methane consumption. 
Chemosphere 32:2219-2246
*/
void f_ch4oxy_casa(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Flux *flux
){
	double j_ch4;		/* methane consumption flux, mg CH4 m-2 day-1 */
	double d_ch4;		/* diffusion coefficient of CH4 in the topsoil(0-0.3m) */
	double c_ch4;		/* concentration gradient of CH4 in the topsoil profile, ppmv cm-1 */
	double d0_ch4;		/* diffusion coefficient of CH4 in air */
	double d_d0;		/* normalized diffusivity in aggregated media */
	double ss;			/* solid phase of the aggregated medium, m3 m-3 bed volume */
	double aa;			/* intra-aggregate pore space, m3 m-3 bed volume */
	double pp;			/* inter-aggregate pore space, m3 m-3 bed volume */
	double th_a;		/* volume of water per unit bed volume contained in intra-arregate pore space, m3 m-3 */
	double th_p;		/* volume of water per unit bed volume contained in inter-arregate pore space, m3 m-3 */
	double s_wa, s_wp;		/* fractional liquid saturation of aa and pp components of the total void volume */
	double xx, yy, zz;		/* exponent terms */
	double fc;		/* field capacity */
	double aaa, bbb, ccc, ddd, eee, fff, ggg, hhh, iii, jjj;
	
	/* from CHEM96_Potter */
	c_ch4 = 0.04;		
	fc = aa = grid->field_cap/300.0;
	pp = grid->pore_cap/300.0 - aa;			/* inter-aggregate pore space */
	ss = 1.0 - grid->pore_cap/300.0;
	
	th_a = loct->m_casa;
	if(loct->m_casa <= aa){
		th_p = 0.0;
	}else if(loct->m_casa>aa && pp>loct->e_casa){
		th_p = loct->e_casa;
	}else if(loct->m_casa>aa && pp<=loct->e_casa){
		th_p = pp;
	}
	
	s_wa = th_a / aa;
	s_wp = th_p / pp;
	
	aaa = pp - th_p;
	bbb = (aa - th_a) / (aa + ss);
	
	/* Eq. 4a in Potter (1996) */
	xx = 0.477 * pow(pp, 3.0) - 0.596 * pow(pp, 2.0) + 0.437 * pp + 0.564;
	/* Eq. 4b in Potter (1996) */
	yy =  0.477 * pow(aaa, 3.0) - 0.596 * pow(aaa, 2.0) + 0.437 * aaa + 0.564;
	/* Eq. 4c in Potter (1996) */
	zz = 0.477 * pow(bbb, 3.0) - 0.596 * pow(bbb, 2.0) + 0.437 * aaa + 0.564;
	
	/* Eq. 3 in Potter (1996) */
	fff = pow((1.0 - s_wa), 2.0);
	jjj = (aa - th_a / aa + ss);
	ggg = pow(jjj, 2.0 * zz);
	hhh = (1.0 - pow(pp, 2.0 * xx));
	iii = (aaa - pow(aaa, 2.0 * yy));
	ccc = fff * ggg * hhh * iii;
	ddd = pow((1.0- s_wa), 2.0) * pow(jjj, 2.0) *  (1.0 - pow(pp, 2.0 * xx)) + aaa - pow(aaa, 2.0 * yy);
	eee = pow((1.0- s_wp), 2.0) * pow(aaa, 2.0 * yy);
	if(jjj>0){
		d_d0 = ccc / ddd + eee;
	}else if(jjj<=0){
		d_d0 = eee;
	}
	
	/* Eq. 2 in Potter (1996) */
	d0_ch4 = (0.194 * (0.9734 + 0.0055 * loct->tmp10_soil));
	
	d_ch4 = d0_ch4 * d_d0;
	
	if(loct->tmp10_soil>0){
		/* Eq. 1 in Potter (1996) */
		j_ch4 = d_ch4 * c_ch4 * (2.7 * pow(10, -4)) * 30 * 24 * 3600;  /* * n[ff]*/
	}else{
		j_ch4 = 0.0;
	}
	
	/* mg CH4 m-2 day-1 */
	(flux->soil).ch4oxy_casa = j_ch4;
}

/* CH4 oxydation by Del Grosso et al. (2000) **********************************/
/*
Del Grosso, S. J., et al. (2000), General CH4 oxidation model and 
comparisons of CH4 oxidation in natural and managed systems, 
Global Biogeochemical Cycles, 14, 999-1019.
*/
void f_ch4oxy_delgrosso(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Flux *flux
){
	double max_ch4oxy;
	double f_tsoil, f_dopt, f_swater;
	double vac, dopt, rfc;
	double wmin, wopt, wmax, vmc;
	double aa, bb, cc, dd, ee, ff;
	double tp1, tp2, tp3, tp4, tp5, tp6, tp7, tp8;
	double pfc, vfrac, theta_V, theta_P, theta_A, sw_p, s_wat;
	
    /* 2016/10/28 revised by A.Ito, thanks to Ueyama-san */
    rfc = grid->field_cap / 300.0;

	if(grid->veg_type==3||grid->veg_type==4){
		/* deciduous forests */
		
		/*  bulk_dens = 1.3;  */
		max_ch4oxy = 40.0 - 18.3 * grid->bulkdens;
		
		aa = (0.1 * loct->wfps*100.0 - 0.5)/(1.84 - 0.5);
		bb = (0.1 * loct->wfps*100.0 - 55.0)/(1.84 - 55.0);
		cc = 0.13*(55.0 - 1.84)/(1.84 - 0.5);
		f_swater = pow(aa, 0.13) * pow(bb, cc);
		if(f_swater < 0.1){
			f_swater = 0.1;
		}
		
		f_tsoil = 0.0209 * loct->tmp10_soil + 0.845;
		if(loct->tmp10_soil < 0.0){
			f_tsoil = 0.0;
		}
		
		/* g C ha-1 day-1 */
		(flux->soil).ch4oxy_delgrosso = max_ch4oxy * f_swater * f_tsoil;
	}else{
		/* grassland, tropical and coniferous forests */
		
		/* inter-aggregate pore (air) space */
		vac = grid->pore_cap - grid->field_cap;	
		/* percent water fraction, % */
		pfc = loct->wfps*100.0 / (grid->field_cap / grid->pore_cap); 
		/* volumetric fraction */
		vfrac = (pfc > 100.0)?((loct->wfps*grid->pore_cap - grid->field_cap)/vac):0.0;
		
		aa = (vfrac>1.0)?1.0:vfrac;
		theta_V = (pfc < 100.0)? (grid->field_cap + aa * vac) : (pfc/100.0 * grid->field_cap);
		
		theta_P = (theta_V < grid->field_cap)? 0.0 : (theta_V - grid->field_cap);
		
		theta_A = (theta_V > grid->field_cap)? grid->field_cap : theta_V;
		
		s_wat = ((theta_V / grid->field_cap)<1.0)?(theta_V / grid->field_cap):1.0;
		
		sw_p =  ((theta_P / vac)<1.0)?(theta_P / vac):1.0;
		
		if((1.0 - s_wat)>0.0){
			tp1 = pow((1.0 - s_wat), 2.0);
		}else{
			tp1 = 0.0;
		}
		
		tp2 = (grid->field_cap - theta_A) / (grid->field_cap + (1.0 - grid->pore_cap));
		
		if(tp2>0.0){
			tp3 = pow(tp2, (0.5 * tp2 + 1.16));
		}else{
			tp3 = 0.0;
		}
		
		tp4 = 1.0 - pow(vac, (0.5 * vac + 1.16));
		
		tp5 = vac - theta_P;
		
		if(tp5>0.0){
			tp6 = pow(tp5, (0.5 * tp5 + 1.16));
		}else{
			tp6 = 0.0;
		}
		
		tp7 = pow((1.0 - sw_p), 2.0);
		
		tp8 = tp1 * tp3 * tp4 * (tp5 - tp6) / (0.000001 + (tp1 * tp3 * tp4) + tp5 - tp6) * 10000000.0;
		if(tp8<0.0){
			tp8 = 0.0;
		}
		
		dopt = tp8/10000000.0 + tp7 * tp6;
		dopt = (dopt>0.0)?dopt:0.0;
		dopt = (dopt<1.0)?dopt:1.0;
		
		/* in Fig. 5a */
		max_ch4oxy = 53.8 * dopt + 0.58;
		
		/* Fig.5a */
		wmin = 3.0 * rfc - 0.28;   /* Wmin */
		wopt = 6.3 * rfc - 0.58;   /* Wopt */
		wmax = 10.6 * rfc + 1.9;   /* Wmax */
		
		vmc = loct->soilwtr_l/300.0;
		
		if(vmc*10.0<wmin ||vmc*10.0>wmax){
			f_swater = 0.1;
		}else{
			dd = (10.0 * vmc - wmin)/(wopt - wmin);
			ee = (10.0 * vmc - wmax)/(wopt - wmax);
			ff = 0.4 * (wmax - wopt)/(wopt - wmin);
			f_swater = pow(dd, 0.4) * pow(ee, ff);
		}
		
		/* Fig.5b */
		aa = (dopt>0.11)?dopt:0.11;
		f_tsoil = loct->tmp10_soil * aa * 0.095 + 0.9;
		
		/* agriculture */
		if(grid->veg_type == 12 || grid->veg_type == 13 || grid->veg_type == 14){
			/* Fig. 5c */
			if(dopt<0.1){
				f_dopt = 0.9;
			}else if(dopt>0.28){
				f_dopt = 0.28;
			}else{
				f_dopt = -4.6 * dopt + 1.6;
			}
		}else{
			f_dopt = 1.0;
		}
		
		/* g C ha-1 day-1 */
		(flux->soil).ch4oxy_delgrosso = max_ch4oxy * f_dopt  * f_tsoil * f_swater;  /*    */
	}
	
	/* mg CH4 m-2 day-1 */
	(flux->soil).ch4oxy_delgrosso *= 1000.0/10000.0 * 16.0/12.0;
}

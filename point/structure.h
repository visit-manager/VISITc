/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/*********************************/
/* declaration of structures	*/
/*********************************/

#include"definition.h"
#include"setting.h"

#define N_ANSIS 256

/* --suffix-- */
/* f: foliage */
/* c: stem and branch */
/* r: root */
/* l: litter */
/* h: humus */

/* grid conditions, derived from data files ****************************************/
struct Grid{ 
		short		flag_datavl;		/* data availability, =1 if avalable */

		/* site information *********************************************/
		char		site_id[16];		/* site ID, e.g. TKY */
		char		file_name[128];		/* file name identifier */
		char		file_clim[128];		/* NCEP file name for site simulation */
		char		file_site[128];		/* site file name */
		char		file_para[128];		/* parameter file name */
		
		double		lat;				/* latitude, degree */
		double		lon;				/* longitude, degree */
		short		dg_row;				/* row number in the area */
		short		dg_col;				/* column number in the area */
		
		long		veg_type;			/* vegetation type */
		long		veg_glc;			/* vegetation type by GLC2000 */
		long		veg_nlc;			/* vegetation type by new land cover */
		long		veg_mod12;			/* vegetation type by MOD12: primary IGBP type */
		long		veg_mod12a;			/* vegetation type by MOD12: UMD class */
		long		veg_mod12b;			/* vegetation type by MOD12: fAPAR/LAI class */
		long		veg_mod12c;			/* vegetation type by MOD12: BGC class */
		long		veg_mod12d;			/* vegetation type by MOD12: PFT class */
		long		veg_mod12e;			/* vegetation type by MOD12: secondary IGBP type */
		double		area;				/* grid area, ha */
		double		agrPt;				/* fraction of agricultural area */
		
		/* climatology: *_a[] means the average during 1948 to 2005 ***********/
		double		temp_wclim[12];		/* WorldClim monthly temperature, deg C */
		double		prec_wclim[12];		/* WorldClim monthly precipitation, mm */
		double		atmp2m_av;			/* mean temperature, degree C */
		double		atmp2m_max;			/* maximum temperature, degree C */
		double		atmp2m_min;			/* minimum temperature, degree C */
		double		stmp10cm_av;		/* average soil temperature, deg C */
		double		gp_atem;			/* average air temperature during the growing period, degree C */
		double		gp_tem;				/* average temperature during the growing period, degree C */
		double		gp_pre;				/* precipitation during the prowing period, mm */
		double		prate_sfc_ann;		/* precipitation, mm year-1 */
		
		/* soil ****************************************************************/
		long		stexture;			/* soil texture */
		long		soil_usda;			/* soil type by USDA map */
		long		soil_type;			/* soil type (simplified) */
		double		topo;				/* topography, orology, and altitude, m above MSL */
		double		fieldcap;			/* soil water holding capacity, mm */
		double		wiltpoint;			/* wilting point, mm */
		double		thermcap;			/* soil thermal capacity */
		double		total_n;			/* total soil N in 1m depth, kg N m-2 */
		double		total_c;			/* total soil C in 1m depth, kg C m-2 */
		double		plf_water;			/*  */
		double		slope_angle;		/* average slope angle, degree */
		double		slope_length;		/* slope length, m */
		
		double		fieldcap30;			/* soil water holding capacity of above 300 mm soil, mm */
		double		soildpth;			/* soil rooting depth, mm */
		double		hyd_cond;			/* hydraulic conductivity */
		double		bulkdens;			/* bulk density, kg litle-1 */
		double		soil_ph;			/* soil pH */
		double		sand_frac;			/* sand fraction */
		double		clay_frac;			/* caly fraction */

		/* calculated: Saxton (1986) */
		double		a_sw;				/* a parameter for water potential function */
		double		b_sw;				/* b parameter for water potential function */
		double		pore_cap;			/* pore capacity in 300mm depth, mm */
		double		field_cap;			/* field capacity in 300mm depth, mm */
};

/* grid conditions, derived from submodules ********************************************************/
struct Loct{ 
		/* time variables ***********************************/
		long		adyear;				/* historical year (AD), in year from the BEGIN */
		long		doy;				/* day of the year (DOY), from 0 to 364/365 */
		long		month;				/* month of the year, from Jan. to Dec., 0 to 11 */
		long		mday;				/* day of the month, from 1 to 28/30/31 */
		long		hour;				/* hour of the day */
		double		hangle;				/* hour angle, degree */

		/* time counter */
		long		time;				/* time to reach stabilization of carbon budget scheme */
		long		climy;
    
        long        veg_state;          /* state of vegetation */
                                        /* 0: natural */
                                        /* 1: cropland */
		/* stand age */
		double		age_stand;
				
		/* GHG scenario ***************************************/
		/* source: http://crga.atmos.uiuc.edu/research/post-sres.html
			M.E.Schlesinger and S.Malyshev			*/
		/* atmospheric CO2, ppmv */
		double		atm_co2_a1[553];	/* SRES A1 */
		double		atm_co2_a2[553];	/* SRES A2 */
		double		atm_co2_b1[553];	/* SRES B1 */
		double		atm_co2_b2[553];	/* SRES B2 */
		/* atmospheric CH4, pptv*/
		double		atm_ch4_a1[553];	/* SRES A1 */
		double		atm_ch4_a2[553];	/* SRES A2 */
		double		atm_ch4_b1[553];	/* SRES B1 */
		double		atm_ch4_b2[553];	/* SRES B2 */
		/* atmospheric N2O, pptv*/
		double		atm_n2o_a1[553];	/* SRES A1 */
		double		atm_n2o_a2[553];	/* SRES A2 */
		double		atm_n2o_b1[553];	/* SRES B1 */
		double		atm_n2o_b2[553];	/* SRES B2 */
		
		/*******************************************************/
		double		depo_nh4_model[24][YSTEP];
		double		depo_no3_model[24][YSTEP];
		double		depo_nh4_model_av[YSTEP];
		double		depo_no3_model_av[YSTEP];
		
		/* climate conditions ***********************************/
		/* annual */
		double		prate_ann;			/* annual precipitation */
		/* daily */
		double		tmp_sfc;			/* ground surface temperature, degree Celcius */
		double		tmp_2m;				/* 2m air temperature, degree Celcius */
		double		tmp10_soil;			/* soil temperature at 0-10 cm depth, degree Celcius */
		double		tmp200_soil;		/* soil temperature at 10-200 cm depth, degree Celcius */
		double		tmp300_soil;		/* soil temperature at 300 cm depth, degree Celcius */
		double		dswrf_sfc;			/* downward shortwave radiation at the surface, W m-2 */
		double		tcdc_clm;			/* total cloudiness, fraction */
		double		prate_sfc;			/* precipitation, mm day-1 */
		double		spfh_2m;			/* specific humidity, kg kg-1 */
		double		wnd_10m;			/* wind velocity, m s-1 */
		
		/* orbital parameters */
		double		soldec[YSTEP]; 		/* solar declination, degree */
		double		daylen[YSTEP]; 		/* day length, hour */
		double		solhgt_h[YSTEP][DSTEP]; 		/* solar hight at midday, degree */
		double		toprad_h[YSTEP][DSTEP];		/* downward solar radiation at the atmosphere-top, W m-2 */
		double		sinb_h[YSTEP][DSTEP];		/* solar angle parameter */
		double		day_srad[YSTEP];		/* daily solar radiation at the surface, MJ m-2 day-1 */

		/* hourly */
		double		sfcrad_h[DSTEP];	/* global radiation at the canopy-top, W m-2 */
		double		grdrad_h[DSTEP];	/* global radiation under the canopy, W m-2 */
		double		par_h[DSTEP];		/* hourly total photosynthetically active radiation (PAR), W m-2 */
		double		parb_h[DSTEP];		/* hourly beam PAR, W m-2*/
		double		pard_h[DSTEP];		/* hourly diffuse PAR, W m-2*/
		double		ppfd_h[DSTEP]; 		/* hourly total photon flux density (PPFD), micro mol photon m-2 s-1 */
		double		ppfdb_h[DSTEP];		/* hourly beam PPFD, micro mol photon m-2 s-1 */
		double		ppfdd_h[DSTEP];		/* hourly diffuse PPFD, micro mol photon m-2 s-1 */
		
		double		ppfd_c_h[DSTEP];	/* hourly PPFD, constant conversion factor, micro mol photon m-2 s-1 */

		/* clear sky */
		double		cssfcrad_h[DSTEP]; 		/* global radiation at the canopy-top, W m-2 */
		double		cspar_h[DSTEP]; 		/* hourly total photosynthetically active radiation (PAR), W m-2 */
		double		csparb_h[DSTEP];		/* hourly beam PAR, W m-2*/
		double		cspard_h[DSTEP];		/* hourly diffuse PAR, W m-2*/
		double		csppfd_h[DSTEP];		/* hourly total photon flux density (PPFD), micro mol photon m-2 s-1 */
		double		csppfdb_h[DSTEP];		/* hourly beam PPFD, micro mol photon m-2 s-1 */
		double		csppfdd_h[DSTEP];		/* hourly diffuse PPFD, micro mol photon m-2 s-1 */

		/* atmosphere */
		long		CO2y;				/* year of atmospheric CO2 scenario */
		double		aCO2;				/* ambient (canopy air) CO2 concentration, ppmv */
		double		bCO2;				/* background (free atmosphere) CO2 concentration, ppmv */
		double		d13C_aCO2;			/* stable carbon isotope composition of ambient CO2, permille */
		double		d13C_bCO2;			/* stable carbon isotope composition of background CO2, permille */
		
		double		air_prsr;			/* air pressure, hPa */
		double		air_dns;			/* density of air, kg m-3 */
		double		lhvp;				/* latent heat MJ kg-1 */

		double		vp;					/* vapour pressure, hPa */
		double		vps;				/* saturation vapour pressure, hPa */
		double		slope_vps;			/* slope of saturation vapour pressure related to tempertaure, hPa deg C-1 */
		double		r_aero;				/* aerodynamic resistance, s m-1 */
		double		vpd;				/* vapour pressure deficit, hPa */
		
		/* radiation: W m–2 *****************************/
		double		albedo_sfc;			/* surface albedo, fraction */
		/* longwave */
		double		rn_long_eco;		/* land surface */
		double		rn_long_tree;		/* trees */
		double		rn_long_c3;			/* C3 grass */
		double		rn_long_c4;			/* C4 grass */
		double		rn_long_ground;		/* soil surface */
		/* shortwave */
		double		rn_short_eco;		/* surface */
		double		rn_short_tree;		/* trees */
		double		rn_short_c3;		/* C3 grass */
		double		rn_short_c4;		/* C4 grass */
		double		rn_short_ground;	/* soil surface */
		/* total net */
		double		rn_eco;				/* surface */
		double		rn_tree;			/* tree layer */
		double		rn_c3;				/* lower layer C3 */
		double		rn_c4;				/* lower layer C4 */
		double		rn_ground;			/* soil ground */
		
		/* vegetation *************************************/
		double		comp_over1;			/* relative occupancy of species 1, fractional area coverage */
		double		comp_over2;			/* relative occupancy of species 2, fractional area coverage */
		
		double		fcover_veg;
		double		fcover_tree;
		double		fcover_c3;
		double		fcover_c4;
		double		fcover_ground;
		double		funder_c3;
		double		funder_c4;

		double		lai;				/* leaf area index (LAI), m2 m-2 */
		/* canopy conductance, mmol m-2 s-1 */
		double		gc_tree;			/* tree */
		double		gc_c3;				/* C3 herbaceous */
		double		gc_c4;				/* C4 herbaceous */
        
        /* development index for crops: 2013/03/17 by A.Ito */
        double      dvr;
        double      dvi;
    
        double      rooting_depth;      /* rooting depth, m */
		
		/* hydrology **************************************************/
		/* Penman-Monteith potential intercepted evaporation rate, mm day-1 */
		double		pm_incep_tree;		/* tree */
		double		pm_incep_c3; 		/* C3 herbaceous */
		double		pm_incep_c4; 		/* c3 herbaceous */
		/* actual intercepted evaporation rate, tree, mm day-1 */
		double		incep;
		double		incep_tree;			/* tree */
		double		incep_c3;			/* c3 herbaceous */
		double		incep_c4;			/* c4 herbaceous */
		
		/* Penman-Monteith potential soil evaporation rate, mm day-1 */
		double		pm_evpr;			/* potential evaporation from soil */
		/* actual soil evaporation rate, mm day-1 */		
		double		evpr;				/* actual evaporation */
		/* soil resistance to evaporation, s m-1 */
		double		rsst_soil;			/* resistance to evaporation */
		
		/* Penman-Monteith potential transpiration rate, tree, mm day-1 */
		double		pm_trnsp_tree; 		
		double		pm_trnsp_c3; 		
		double		pm_trnsp_c4; 		
		/* actual transpiration rate, tree, mm day-1 */
		double		trnsp;				/* total vegetation */
		double		trnsp_tree;			/* tree */
		double		trnsp_c3;			/* C3 grass */
		double		trnsp_c4;			/* C4 grass */
		
		double		thaw;				/* snow thaw water, mm day-1 */
		double		ro1; 				/* daily runoff of upper, mm */
		double		ro2; 				/* daily runoff of lower, mm */

		double		snp; 				/* snow fraction of precipitation */
		double		snow_acc;			/* snow accumulation, mm water */
		double		soilwtr_l;			/* soil water upper layer, mm water */
		double		soilwtr_h;			/* soil water lower layer, mm water */
		double		soilappr_l;			/* upper soil aperture, fraction */
		double		soilappr_w;			/* lower soil aperture, fraction */
		
		double		pet;				/* potential evapotranspiration */
		double		aet;				/* actual evapotranspiration */
		double		wfps;				/* water-filled pore space */
		double		vmc;
		double		vmc_pre;

		double		iw_casa;			/* CASA moisture index */
		double		rdr_casa;			/* CASA dryness index */
		double		e_casa;				/* CASA parameter */
		double		m_casa;				/* CASA parameter */
		double		m_casa_pre;			/* CASA parameter */
	
		/* soil water potential, MPa */
		/* litter */
		double		pot_total_l;		/* total */
		double		pot_grav_l;			/* gravity */
		double		pot_matric_l;		/* matric */
		/* humus */
		double		pot_total_h;		/* total */
		double		pot_grav_h;			/* gravity */
		double		pot_matric_h;		/* matric */
	
		/*****************************/
		double		d_c_soil;
		double		d_c_min;
		double		d_n_soil;
		
		double		ndepo[3];			/* N deposition */
		double		depo_no3;
		double		depo_nh4;
		
		double		nfert_no3;
		double		nfert_nh4;
		
		/* 30-min observation data **************************/
		long		metdat_frag;		/* flag of data availability */
		/* +1: temperature */
		/* +2: air humidity */
		/* +4: PAR */		/* i.e., 7, when all data are available */
		
		/* erosion */
		double		f_rain;
		double		f_slope;
		double		f_erodibility;
		double		f_cover;
		double		f_conservation;
		
		/* greenhouse gas */
		double		water_table_depth;		/* water table depth, m */
        double      water_table_depth_pre;  /* WTD at the last step */
		double		prof_ch4[SOIL_LAYER+2];		/* soil CH4 concentration, micro M */
		double		npp_max;	/* maximum daily NPP, g C m-2 day-1 */
		double		tmp_soil_mean;	/* annual mean temperature */

		/* surrogate variables for debugging intermediate variables */
		double		x1;
		double		x2;
		double		x3;
		double		x4[SOIL_LAYER+2];
		double		x5[SOIL_LAYER+2];
};

/* vegetation characteristics *******************************************/
struct Pchar{ 
		short		c_layer;		/* canopy layer ID */

		/** radiation **/
		double		albedo; 		/* reflectivity or albedo, fraction */
		double		ppfd_t; 		/* PPFD at the canopy-top, micro mol photon m-2 s-1 */
		double		ppfd_b;			/* PPFD at the canopy-bottom, micro mol photon m-2 s-1 */

		/*** allocation ***/
		double		opt_lai; 		/* optimum leaf area index, m2 m-2 */
		double		alloc_ass;		/* allocation ration for assimilation organ, fraction */
		double		alloc_abg; 		/* allocation ration for aboveground non-assimilation organ, fraction */
		double		malloc_f;		/* to foliage, ratio */
		double		malloc_c;		/* to stem, ratio */
		double		malloc_r; 		/* to root, ratio */
        double      malloc_g;       /* to grain */
		double		max_nsch_storage;	/* maximum capacity of carbohydrate storage */
		
		/*** phenology ***/
		long		phenoltype;
				/* 0: no vegetation, bareland */
				/* 1: evergreen */
				/* 2: cold deciduous */
		double		ctmp_lfdsp;			/* cumulative temperature for leaf display, degree-day */
		double		ctmp_lfshd;			/* cumulative temperature for leaf shedding, degree-day */
		double		gdd; 				/* growing degree days */
		double		cdd;				/* coldness degree days */
		double		grw_pd;				/* growing period, days */
		long		season;				/* phenological stage as a function of season */
				/* 0: dormancy */
				/* 1: vegetative growth */
				/* 2: emergence of new leaf */
				/* 3: abandon of old leaf */
		long		day_flush;			/* day of year */
		long		day_shed;			/* day of year */
    
        long        crop_stage;         /* growth stage for crops: 2013/03/17 by A.Ito */
	
		/* revision: 2012/01/06 (A.Ito) */
		double		ncd;				/* number of chilling (<5deg C) days */
		double		gddth;				/* threshold growing degree-days */
		double		gddjan;				/* cumulative temperature from January 1 */
		
		/*** photosynthesis ***/
		long		phototype;			/* photosynthetic metabolic pathway, 3=C3, 4=c3, 5=CAM */
		double		pmax;				/* potential maximum rate, micro mol CO2 m-2 s-1 */
		double		ptop;				/* canopy-top photosynthetic rate, micro mol CO2 m-2 s-1 */
		double		psat;				/* light-saturated rate, micro mol CO2 m-2 s-1 */
		double		sla;				/* specific leaf area, cm2 g-1 dm */
		double		eK0;				/* light attenuation coefficient, no dimension */
		double		eK;					/* light attenuation coefficient, no dimension */
		double		lue0;				/* control light dependence coefficient, mol CO2 mol photon-1 */
		double		lue;				/* daily quantum yield, mol CO2 mol photon-1 */ 
		double		topt;				/* optimum temperature, deg C */
		double		topt0;				/* optimum temperature (control), deg C */
		double		tmin;				/* minimum temperature, deg C */
		double		tmax;				/* maximum temperature, deg C */
		double		ci;					/* daily intercellular CO2 concentration, ppmv */
		double		kmci;				/* dependence of photosynthesis on intercellular CO2 concentration, ppmv */
		double		cmpcd0;				/* CO2 compensation point, ppmv */
		double		cmpcd;				/* CO2 compensation point, ppmv */
		double		gs;					/* daily stomatal conductance, mmol H2O m-2 s-1 */
		double		gc;					/* daily canopy conductance, mmol H2O m-2 s-1 */
		double		gs_b0;				/* parameter of stomata model, mmol H2O m-2 s-1 */
		double		gs_b1;				/* parameter of stomata model (VPD dependence) */
		double		gs_b2;				/* parameters of stomatal conductance, mmol H2O m-2 s-1 */
		double		km_nstl;			/* non-stomatal (soil moinsture) photosynthesis limitation factor */
		
		/*** respiration ***/
		/* specific growth respiration rate, g C g C-1 alloc */
		double		rgf;				/* leaf */
		double		rgc;				/* stem */
		double		rgr; 				/* root */
		/* specific maintenance respiration rate at 15 degC, mg C g C-1 day-1 */
		double		rmf;				/* leaf */
		double		rmc;				/* stem */
		double		rmr; 				/* root */
		double		rmf0; 				/* leaf, basic */
		double		rmc_s;				/* stem sapwood */
		double		rmr_s; 				/* fine root */
		double		rmc_h;				/* stem heartwood */
		double		rmr_h;				/* coarse root */
		/* temperature dependence, dimensionless */
		double		qTf0;				/* leaf at 15 deg C */
		double		qTc0;				/* stem at 15 deg C */
		double		qTr0;				/* root at 15 deg C */
		double		qTf;				/* leaf */
		double		qTc;				/* stem */
		double		qTr; 				/* root */
		/* size-dependence parameter */
		double		f_sz_s;				/* stem */
		double		f_sz_r;				/* root */
		
		/*** specific litter fall ***/
		double		lf0;				/* leaf, basic */
		double		lc0;				/* stem, basic */
		double		lr0; 				/* root, basic */
		double		lf;					/* leaf */
		double		lc;					/* stem */
		double		lr; 				/* root */
		double		dcd; 				/* leaf shedding parameter, fraction */
		
		/*** root stratification (Zeng 2001) ***/
		double		root_dpt_a; 		/* root profile parameter  */
		double		root_dpt_b; 		/* root profile parameter  */
		
		/*** photosynthetgic 13C fractionation ***/
		double		photo_13c_frac;		/* photosynthetic d13C fractionation */
		
		/* nitrogen */
		double		n_salvage;			/* N re-uptake from senescent leaves */
		double		cn0_strg;
		double		amax_nphoto;
		double		kn_nphoto;
		double		rd_n;
		double		n_conc_larea;
		double		n_conc_lmass;
		double		n_opt_leaf;
		double		n_max_mobile;
		
		/* for experiments of N-limitation and LAI increase */
		double		lai_contemp[YSTEP];	/* contemporary LAI */
		double		f_downreg;			/* down-regulation factor: experimental */
};

/* canopy characteristics **************************************/
struct Cchar{
		/*** De Pury & Farquhar canopy photosynthesis model ***/
		long		photo_type;			/* C3/C4 pathway */

		double		ppfdb_top;			/* beam PPFD at top */
		double		ppfdd_top;			/* diffuse PPFD at top */
		double		ppfdb_btm;			/* beam PPFD at bottom */
		double		ppfdd_btm;			/* diffuse PPFD at bottom */
		
		double		lma;				/* leaf mass per area */
		double		lai; 				/* daily leaf area index(LAI), m2 m-2 */
		double		lai_sn;				/* sunny LAI */
		double		lai_sd;				/* shade LAI */

		double		vpd_a;				/* ambient VPD, Pa */
		double		co2_a;				/* ambient CO2 concentration, Pa */
		double		ci_sn;				/* leaf-internal CO2 concentration, Pa */
		double		ci_sd;				/* leaf-internal CO2 concentration, Pa */
		double		gs_sn;				/* stomatal conductance of sun leaves */
		double		gs_sd;				/* stomatal conductance of shade leaves */
		double		o2_i;				/* leaf-internal O2 concentration, Pa */
		double		tmp;				/* canopy temperature, K */
		double		comp_co2;			/* CO2 compensation point, Pa */
		double		ke_d;				/* light attenuation coefficient, diffused */
		double		ke_b1;				/* light attenuation coefficient, beam */
		double		ke_b2;				/* light attenuation coefficient, beam and scattered beam */
		double		scttr;				/* scattering factor */
		double		spect;				/* spectral correction factor (PS-II absortption) */
		double		rfl_d;				/* reflectance, diffuse */
		double		rfl_b;				/* reflectance, beam */
		double		appfd;				/* absorbed PPFD */
		double		appfd_sn;			/* absorbed PPFD by sun leaves */
		double		appfd_sd;			/* absorbed PPFD by shade leaves */
		double		apar;				/* absorbed PAR */
		double		apar_sn;			/* absorbed PAR by sun leaves */
		double		apar_sd;			/* absorbed PAR by shade leaves */
		double		ke_n;				/* attenuation coefficient, nitrogen allocation */
		double		nc_ct;				/* canopy-top N concentration, mmol m-2 */
		double		nc_np;				/* non-photosynthetic N concentration, mmol m-2 */
		double		vcmx_b;				/* */
		double		vcmx_con;			/* */
		double		vcmx;				/* */
		double		vcmx_sn;			/* */
		double		vcmx_sd;			/* */
		double		vcmx_sn0;			/* */
		double		vcmx_sd0;			/* */
		double		photocap_n0;		/* */
		double		photocap_n;			/* */
		double		kc0;				/* kc at 25 deg C */
		double		kc;					/* Michaelis-Menten constant of carboxylation */
		double		ko0;				/* ko at 25 deg C */
		double		ko;					/* Michaelis-Menten constant of oxygenation */
		double		acen_kc;			/* activation energy of Kc */
		double		acen_ko;			/* activation energy of Ko */
		double		k_effc;				/* effective Michaelis-Menten constant */
		double		ft_et;				/* temperature factor of electron transport */
		double		acen_et;			/* activation energy of electron transport */
		double		hc_et;				/* */
		double		str_et;				/* */
		double		acen_vc;			/* activation energy of carboxylation */
		double		acen_rd;			/* activation energy of dark respiration */
		double		jvr;				/* Jmax/Vcmax ratio */
		double		cv_j;				/* convexity of J-APAR */
		double		jmx_sn;				/* Jmax of sunny */
		double		jmx_sd;				/* Jmax of shade */
		double		j_sn;				/* electron transport rate, sunny */
		double		j_sd;				/* electron transport rate, shade */
		double		cv_p;				/* convexity of assimulation rate */
		double		pc_sn;				/* */
		double		pc_sd;				/* */
		double		pj_sn;				/* */
		double		pj_sd;				/* */
		double		p_sn;				/* net assimilation rate of sun leaf */
		double		p_sd;				/* net assimilation rate of shade leaf */
		double		a_sn;				/* gross assimilation rate of sun leaf */
		double		a_sd;				/* gross assimilation rate of shade leaf */
		
		double		gb;					/* leaf boundary-layer conductance */
		double		gss_b0; 
		double		gss_b1; 
		double		gss_b2;				/* parameters of stomatal conductance, mmol H2O m-2 s-1 */
		
		double		rd0;
		double		rd; 				/* canopy dark respiration */
		double		rd_sn;				/* canopy dark respiration: sunny */
		double		rd_sd;				/* canopy dark respiration: shaded */
		double		ft_rd;
    
        double      g_msp;              /* mesophyll conductance */
};

/* soil characteristics *************************************/
struct Schar{ 
		double		albedo0;			/* bare soil reflectivity, or albedo */
		double		albedo;				/* whole soil (with snow) reflectivity, or albedo */
		double		tmp;				/* soil temperature */
		
		/* base decomposition coefficient, fraction */
		double		sr_lf;				/*  leaf litter */
		double		sr_lc;				/*  stem litter */
		double		sr_lr;				/*  root litter */
		double		sr_ha;				/*  active humus */
		double		sr_hi;				/*  intermediate humus */
		double		sr_hp;				/*  reluctant humus */
		
		/* combined environmental factor */
		double		f_tm_l;				/* litter */
		double		f_tm_h;				/* humus */
		
		/* moisture dependence, fraction of soil water */
		double		kml;				
		double		kmh; 
		/* moisture dependence, fraction of soil water */
		double		kmsl;				
		double		kmsh; 
	
		/* CO2 emission ratio to litter decomposition, fraction */
		double		f_co2_lf;			
		double		f_co2_lc;
		double		f_co2_lr;
		
		/* mineral soil formation ratio to litter decomposition, fraction */
		double		f_hm_a;				
		double		f_hm_i;
		double		f_hm_p;
		
		/* maximum nitrate concentration */
		double 		ntrcon_max;			/* g N / kg H2O */
};

/* ecosystem characteristics *************************************/
struct Echar{ 
		char		para_ver_id[15];		/* parameter vesion ID */
		char		para_date_id[15];		/* parameter update date */

		struct		Pchar tree;			/* for tree plants */
		struct		Pchar c3;			/* for c3 herbaceous plants */
		struct		Pchar c4;			/* for c4 herbaceous plants */
		struct		Schar soil;			/* for soil organic matter */
		
		/* canopies */
		struct 		Cchar overs1;		/* overstory 1 */
		struct 		Cchar overs2;		/* overstory 2 */
		struct 		Cchar unders;		/* floor plant */
};

/* plant biomass ************************************************************/
struct Pmas{ 
		double		lai;				/* leaf area index, m2 m-2 */
        double      rooting_depth;      /* rooting depth, m */
    
		/* carbon mass, Mg C ha-1 */ 
		double		fol; 				/* foliage */
		double		stm; 				/* stem and branch */
		double		rot; 				/* root */
		
		double		stm_sp;				/* sub-compartment: sap */
		double		stm_ht;				/* sub-compartment: heart */
		double		rot_fn;				/* sub-compartment: fine */
		double		rot_tp;				/* sub-compartment: tap/coarse */
	
		double		nsch_storage;		/* storage of non-structural carbohydrate */
        double      grain;
		
		double		plant;				/* total plant */
		
		/* stable carbon isotope composition, dimensionless */
		double		d13c_fol;			/* foliage */
		double		d13c_stm;			/* stem and branch */
		double		d13c_rot;			/* root */
		
		double		d13c_plant; 		/* total plant */

		/* nitrogen mass, g N ha-1 */ 
		double		n_canopy;			/* foliage */
		double		n_storage;			/* stem and branch storage (mobile) */
		double		n_struct;			/* stem and branch structual */

		double		n_stm_sp;			/* sub-compartment: sap */
		double		n_stm_ht;			/* sub-compartment: heart */
		double		n_rot_fn;			/* sub-compartment: fine */
		double		n_rot_tp;			/* sub-compartment: tap/coarse */

		double		n_plant;			/* total plant */
};

/* soil organic matter ********************************************************************/
struct Smas{ 
		/* carbon mass, Mg C ha-1 */ 
		double		ltr_tf;				/* foliage litter */
		double		ltr_tc;				/* stem litter */
		double		ltr_tr;				/* root litter */
		
		double		ltr_gf;				/* foliage litter */
		double		ltr_gc;				/* stem litter */
		double		ltr_gr;				/* root litter */
		
		double		ltr_tc_ag;			/* sub-compartment: aboveground */
		double		ltr_tc_bg;			/* sub-compartment: belowground */
		
		double		msl_a;				/* active humus */
		double		msl_i;				/* intermediate humus */
		double		msl_p;				/* passive humus */
		
		double		mcrb;				/* microbe */
		double		doc;				/* dissolved organic carbon, DOC */
		
		double		soil;				/* total soil */
		
		/* stable carbon isotope composition */
		double		d13c_ltr_f;			/* foliage */
		double		d13c_ltr_c;			/* stem */
		double		d13c_ltr_r;			/* root */
		
		double		d13c_msl_a;			/* active */
		double		d13c_msl_i;			/* intermediate */
		double		d13c_msl_p;			/* passive */

		double		d13c_soil;			/* total soil */

		/* nitrogen mass, g N ha-1 */
		/* inorganic */
		double		n_nh4;				/* ammonium */
		double		n_no3;				/* nitrate */

		/* organic */
		double		n_lttr;				/* foliage litter */
		double		n_hums;				/* active humus */
		double		n_mcrb;				/* microbe */
		
		double		n_soil;				/* total soil */
};

/* ecosystem carbon storage *************************************************/
struct Mass{ 
		/* evolving variables ( x(i) = x(i-1) + d(i) ) ************************************/
		double		snwa;				/* water equivalent snow depth, mm */
		double		sw30;				/* soil water content of above 30cm soil, mm */
		double		sww;				/* whole soil water content, mm */

		/* phenology */
		long		gd, bbm; 		/* vegetative growing period, days */
		double		gdd; 			/* cumulative temperature, growth degree days, degC days */
		
		double		gdd_tree; 		/* growing degree days */
		double		cdd_tree;		/* coldness degree days */

		double		gdd_c3; 		/* growing degree days */
		double		cdd_c3;			/* coldness degree days */

		double		gdd_c4; 		/* growing degree days */
		double		cdd_c4;			/* coldness degree days */
		
		/* pool */
		struct		Pmas tree;		/* tree plant mass */
		struct		Pmas c3;		/* floor C3 plant mass */
		struct		Pmas c4;		/* floor C4 plant mass */
		struct		Smas soil;		/* soil mass */

		/* carbon mass */
		double		total_c;		/* ecosystem total carbon */
		/* nitrogen mass */
		double		total_n;		/* ecosystem total nitrogen */

		/* stable carbon isotope composition */
		double		d13c_total;		/* ecosystem total carbon */
};

/* plant carbon fluxes ***************************************************************/
struct Pflx{ 
		/* carbon flux, Mg C ha-1 day-1 */ 
		double		gpp;				/* gross primary production */
		double		spp;				/* surplus production */
		double		epp;				/* effective production */
		double		npp;				/* net primary production */
		
		double		gpp_df97;
		
		double		tpf;				/* translocation of photosynthate to foliage */
		double		tpc;				/* translocation of photosynthate to stem */
		double		tpr;				/* translocation of photosynthate to root */
		double		tpp;				/* translocation of photosynthate, total */
        double      tpg;
		
		double		rtpc;               
		double		rtpr;
		
		double		rp;					/* plant respiration, =rpm+rpg */
		double		rpg;				/* plant growth respiration */
		double		rpm;				/* plant maintenance respiration */
		double		rfg;				/* foliage growth respiration */
		double		rfm;				/* foliage maintenance respiration */
		double		rcg;				/* stem and branch growth respiration */
		double		rcm;				/* stem and branch maintenance respiration */
		double		rrg;				/* root growth respiration */
		double		rrm;				/* root maintenance respiration */
		
		double		lf;					/* foliage litterfall */
		double		lc;					/* stem and branch litterfall */
		double		lr;					/* root litterfall */
		double		lL;					/* total litterfall */
		
		double		hvst;				/* harvest of biomass */
		
		/* hourly */
		double		ar_h;				/* hourly autotrophic respiration rate */
		
		/* CH4 emit */
		double		ch4_emit_kepp;
		double		emit_ch4_kirschbaum_mass;
		double		emit_ch4_kirschbaum_photo;
		
		/* nitrogen uptake */
		/* symbiotic N fixation */
		double		n_biofix;	
		/* root uptake */
		double		uptake_nh4;		
		double		uptake_no3;
		
		double		n_abdn_cnpy;
		double		n_abdn_strg;
		double		n_salvage;
		double		n_alloc_cnpy;
		double		n_alloc_strg;
		double		n_realloc;
		
		/* stable carbon isotope composition */
		double		d13c_gpp; 
		double		d13c_spp; 
		double		d13c_epp; 
		double		d13c_npp; 
		
		double		d13c_tpf; 
		double		d13c_tpc; 
		double		d13c_tpr; 
		double		d13c_tpp; 
		
		double		d13c_rp; 
		double		d13c_rpg; 
		double		d13c_rpm; 
		double		d13c_rfg; 
		double		d13c_rfm; 
		double		d13c_rcg; 
		double		d13c_rcm; 
		double		d13c_rrg; 
		double		d13c_rrm; 
		
		double		d13c_lf; 
		double		d13c_lc; 
		double		d13c_lr; 
		double		d13c_lL; 
		
		double		d13c_lf_c; 

		double		d13c_hvst; 
};

/* soil carbon fluxes, all in Mg C ha-1 day-1 *************************************/
struct Sflx{		 
		/* carbon flux ************************/ 
		/* litter input */
		double		li_tf;				/* foliage litter */
		double		li_tc;				/* stem litter */
		double		li_tr;				/* root litter */
		
		double		li_gf;				/* foliage litter */
		double		li_gc;				/* stem litter */
		double		li_gr;				/* root litter */
		
		/* humification */
		double		hf_tfa;		/* */
		double		hf_tfi;
		double		hf_tfp;
		double		hf_tca;
		double		hf_tci;
		double		hf_tcp;
		double		hf_tra;
		double		hf_tri;
		double		hf_trp;

		double		hf_gfa;		/* */
		double		hf_gfi;
		double		hf_gfp;
		double		hf_gca;
		double		hf_gci;
		double		hf_gcp;
		double		hf_gra;
		double		hf_gri;
		double		hf_grp;

		/* mineralization */
		double		mr_tf;			/* tree leaf litter */
		double		mr_tc;			/* tree stem litter */
		double		mr_tr;			/* tree root litter */
		double		mr_gf;			/* herb leaf litter */
		double		mr_gc;			/* herb stem litter */
		double		mr_gr;			/* herb root litter */
		double		mr_a;			/* active humus */
		double		mr_i;			/* intermediate humus */
		double		mr_p;			/* passive humus */

		double		hr;				/* total microbial heterotrophic respiration */
		double		hr_h;			/* short-term heterotrophic respiration */
		
		double		doc_boyer;		/* DOC export by Boyer et al. (1996), gC/m2/day */
		double		erosion_soil;
		double		erosion_orgmat;
		double		erosion_carbon;
		
		/* methane oxydation */
		double		ch4oxy_ridg;		/* Ridgwell et al. (1999) */
		double		ch4oxy_casa;		/* methane oxydation */
		double		ch4oxy_delgrosso;		/* methane oxydation */
		double		ch4oxy_curry;		/* Curry (2007) */
		double		ch4gen;				/* methane generation */
		
		double		ch4_wh;
		double		ch4_wh_plant;
		double		ch4_wh_ebull;
		double		ch4_wh_diff;

		/* nitrogen flux *********************************/
		/* N flow, g N ha-1 day-1 **/
		double	d_n2o_ntr_ngas;			/* N2O by nitrification, NGAS scheme */
		double	d_n2o_dnt_ngas;			/* N2O by denitrification, NGAS scheme */
		double	d_n2o_ngas;				/* total N2O, NGAS scheme */
		double	d_n2_ngas;				/* N2O by denitrification, NGAS scheme */
		
		double	d_n2o_ntr_ngas_2;			
		double	d_n2o_dnt_ngas_2;			
		double	d_n2o_ngas_2;			
		double	d_n2_ngas_2;			
		
		double	d_n2o_casa;			
		double	d_no_casa;			
		double	d_n2_casa;			
		
		double	n_nh3vlt;			
		double	n_leach;			
		double	n_minerlz_lttr;
		double	n_minerlz_hums;
		double	n_nitrif;			/* in NGAS */
		double	n_immbl;
		double	n_mcrb_abdn;
			
		/* stable carbon isotope composition, in per mil */
		double		d13c_li_f;		/* leaf litter input */
		double		d13c_li_c;		/* stem litter input */
		double		d13c_li_r;		/* root litter input */
		
		double		d13c_hf_f;
		double		d13c_hf_c;
		double		d13c_hf_r;

		double		d13c_mr_f;
		double		d13c_mr_c;
		double		d13c_mr_r;

		double		d13c_mr_a;
		double		d13c_mr_i;
		double		d13c_mr_p;

		double		d13c_hr; /* total microbial heterotrophic respiration */
};

/* ecosystem carbon fluxes *******************************************************/
struct Flux{ 
		struct		Pflx tree;		/* tree plant fluxes */
		struct		Pflx c3;		/* floor plant fluxes */
		struct		Pflx c4;		/* floor plant fluxes */
		struct		Sflx soil;		/* soil fluxes */

		/* ecosystem carbon budget */ 
		double		gpp;			/* ecosystem gross primary production */
		double		npp;			/* ecosystem net primary production */
		double		er;				/* ecosystem respiration */
		double		sr;				/* soil respiration */
		double		rr;				/* plant root respiration */
		double		nep;			/* net ecosystem production, Mg C ha-1 day-1 */
		double		necb;			/* net carbon balance of grid, Mg C ha-1 day-1 */
        double      hvst;

		/* VOC emission, g C m-2 day-1 */
		double		voc_isopr_g97;		/* Isoprene */
		double		voc_monotrp_g97;	/* Monoterpene */
		double		voc_methanl_g97;	/* Methanol */
		double		voc_acetone_g97;	/* Acetone */
		double		voc_actaldhd_g97;	/* Acetoaldehyde */
		double		voc_frmardhd_g97;	/* Formaldehyde */
		double		voc_formacd_g97;	/* Formacid */
		double		voc_acetacd_g97;	/* Acetoacid */
		double		voc_co_g97;			/* Carbon Monoxide */

		/* stable carbon isotope composition */
		double		d13c_nep; 
		double		d13c_ncb; 
		double		d13c_er;
};

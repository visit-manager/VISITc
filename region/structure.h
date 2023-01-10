/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
		create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
		modify: 2006 / 08 / 04 by Akihiko Ito: Version 0.2.0.
		modify: 2007 / 11 / 22 by Akihiko Ito: Version 0.2.0.
		modify: 2007 / 12 / 01 by Akihiko Ito: Version 0.2.0.
		modify: 2008 / 05 / 28 by Akihiko Ito: Version 0.2.0.
*/

/*********************************/
/* declaration of structures	*/
/*********************************/

/* calculation step */
#define DSTEP 24		/* diurnal 24: hourly */
/* number of soil layers */
/* #define SOILWET_LAYER 20 */
#define SOILWET_LAYER 6

#define GHG_TS 553

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
		char		file_ncep[128];		/* NCEP file name for site simulation */
		char		file_site[128];		/* site file name */
		char		file_para[128];		/* parameter file name */
		
		float	lat;				/* latitude, degree */
		float	lon;				/* longitude, degree */
		long		calc_flag;			/* column number in the area */
		
		short	veg_type;			/* vegetation type */
		long		veg_glc;			/* vegetation type by GLC2000 */
		long		veg_nlc;			/* vegetation type by new land cover */
		long		veg_mod12;			/* vegetation type by MOD12 */
		float		area;				/* grid area, ha */
		float		agrPt;				/* fraction of agricultural area */
		short		potveg;				/* potential vegetation */
		
		/* climatology: *_a[] means the average during 1948 to 2005 ***********/
		float		tmp2m_ann;			/* mean temperature, degree C */
		float		tmp2m_max;			/* maximum temperature, degree C */
		float		tmp2m_min;			/* minimum temperature, degree C */
		float		stmp10cm_av;		/* average soil temperature, deg C */
		float		gp_atem;			/* average air temperature during the growing period, degree C */
		float		gp_tem;				/* average temperature during the growing period, degree C */
		float		gp_pre;				/* precipitation during the prowing period, mm */
		float		prec_ann;		/* precipitation, mm year-1 */
		
		/* climate data for regional simulstions */
		float		prec_region;        /* precipitation */
		float		srad_region;        /* short-wave radiation */
		float		tmax_region;        /* max. temperature */
		float		tmin_region;        /* min. temperature */
		float		humd_region;        /* air humidity */
		float		wind_region;        /* wind velocity */
    
        float       tmp_clim;
        float       prec_clim;
        float       aet_clim;
        float       pet_clim;
        float       n_clim;
		
		/* soil ****************************************************************/
		long		stexture;			/* soil texture */
		long		soil_usda;			/* soil type by USDA map */
		long		soil_type;			/* soil type (simplified) */
		float		topo;				/* topography, orology, and altitude, m above MSL */
		float		fieldcap;			/* soil water holding capacity, mm */
		float		wiltpoint;			/* wilting point, mm */
		float		thermcap;			/* soil thermal capacity */
		float		total_n;			/* total soil N in 1m depth, kg N m-2 */
		float		total_c;			/* total soil C in 1m depth, kg C m-2 */
		float		plf_water;			/*  */
		float		slope_angle;		/* average slope angle, degree */
		float		slope_length;		/* slope length, m */
		
		float		fieldcap30;			/* soil water holding capacity of above 300 mm soil, mm */
		float		soildpth;			/* soil rooting depth, mm */
		float		hyd_cond;			/* hydraulic conductivity */
		float		bulkdens;			/* bulk density, kg litle-1 */
		float		soil_ph;			/* soil pH */
		float		sand_frac;			/* sand fraction */
		float		clay_frac;			/* caly fraction */

		/* calculated: Saxton (1986) */
		float		a_sw;				/* a parameter for water potential function */
		float		b_sw;				/* b parameter for water potential function */
		float		pore_cap;			/* pore capacity in 300mm depth, mm */
		float		field_cap;			/* field capacity in 300mm depth, mm */

		long		time;				/* time to reach stabilization of carbon budget scheme */
		/* stand age */
		float		age_stand;

        /* land-use change [Hurtt's matrix] *************************/
        float       gcrop;              /* fraction of each gridcell in cropland */
        float       gflcp;              /* fraction of each gridcell that transitioned from cropland to pasture */
        float       gflcs;              /* fraction of each gridcell that transitioned from cropland to secondary land */
        float       gflcu;              /* fraction of each gridcell that transitioned from cropland to urban land */
        float       gflpc;              /* fraction of each gridcell that transitioned from pasture to cropland */
        float       gflps;              /* fraction of each gridcell that transitioned from pasture to secondary land */
        float       gflpu;              /* fraction of each gridcell that transitioned from pasture to urban land */
        float       gflsc;              /* fraction of each gridcell that transitioned from secondary land to cropland */
        float       gflsp;              /* fraction of each gridcell that transitioned from secondary land to pasture */
        float       gflsu;              /* fraction of each gridcell that transitioned from secondary land to urban land */
        float       gfluc;              /* fraction of each gridcell that transitioned from urban land to cropland */
        float       gflup;              /* fraction of each gridcell that transitioned from urban land to pasture */
        float       gflus;              /* fraction of each gridcell that transitioned from urban land to secondary land */
        float       gflvc;              /* fraction of each gridcell that transitioned from primary land to cropland */
        float       gflvp;              /* fraction of each gridcell that transitioned from primary land to pasture */
        float       gflvu;              /* fraction of each gridcell that transitioned from primary land to urban land */
        float       gfsh1;              /* fraction of each gridcell that had wood harvested from mature secondary forested land */
        float       gfsh2;              /* fraction of each gridcell that had wood harvested from young secondary forested land */
        float       gfsh3;              /* fraction of each gridcell that had wood harvested from secondary non-forested land */
        float       gfvh1;              /* fraction of each gridcell that had wood harvested from primary forested land */
        float       gfvh2;              /* fraction of each gridcell that had wood harvested from primary non-forested land */
        float       gothr;              /* fraction of each gridcell in primary land */
        float       gpast;              /* fraction of each gridcell in pasture */
        float       gsbh1;              /* mature secondary forest biomass harvested from each gridcell (in kgC) */
        float       gsbh2;              /* young secondary forest biomass harvested from each gridcell (in kgC) */
        float       gsbh3;              /* secondary non-forest biomass harvested from each gridcell (in kgC) */
        float       gsecd;              /* fraction of each gridcell in secondary land */
        float       gssma;              /* mean age of secondary land in each gridcell */
        float       gssmb;              /* mean biomass density of secondary land in each gridcell (in kgC/m^2) */
        float       gsumm;              /* the sum of fractions of each gridcell occupied by cropland, pasture, primary, secondary, urban land, and ice/water (should add to 1 everywhere) */
        float       gurbn;              /* fraction of each gridcell in urban land */
        float       gvbh1;              /* primary forest biomass harvested from each gridcell (in kgC) */
        float       gvbh2;              /* primary non-forest biomass harvested from each gridcell (in kgC) */
        float       gzdis;              /*  */
       /* pre-industrial deforestation rate used for spin-up */
        float       deforest_pi;        
        float       gflvc_pi;        
        float       gflvp_pi;        
        float       gflsc_pi;        
        float       gflsp_pi;        
        float       gm_crop_con;        /* contemporary cropland fraction */
        float       gm_past_con;        /* contemporary pasture fraction */
        float       gm_crop_pre;        /* contemporary cropland fraction */
        float       gm_past_pre;        /* contemporary pasture fraction */
        float       gm_crop_base;       /* base cropland fraction */
        float       gm_crop_trend;      /* change in cropland fraction, yr-1 */
        float       gm_paddy;
        float       gm_paddy_base;
        float       gm_wetland;
    
        float       p_lai;              /* prescribed LAI by fix-CO2 experiments */
};

/* grid conditions, derived from submodules ********************************************************/
struct Loct{ 
		/* time variables ***********************************/
		long		adyear;				/* historical year, in year from the BEGIN */
		long		doy;				/* day of the year, from 0 to 364/365 */
		long		month;				/* month of the year, from Jan. to Dec., 0 to 11 */
		long		mday;				/* day of the month, from 1 to 28/30/31 */
		long		hour;				/* hour of the day */
		float		hangle;				/* hour angle, degree */

		/* time counter */
		long		climy;              /* year of the current climate dataset  */
        short       phase;
        
		/*******************************************************/
		float		depo_nh4_model;
		float		depo_no3_model;
		float		depo_nh4_model_av;
		float		depo_no3_model_av;  
		
		/* climate conditions ***********************************/
		/* annual */
		float		prate_ann;			/* annual precipitation */
        
		/* orbital parameters */
		float		soldec; 		/* solar declination, degree */
		float		daylen; 		/* day length, hour */
		float		solhgt_h; 		/* solar hight at midday, degree */
		float		toprad_h;		/* downward solar radiation at the atmosphere-top, W m-2 */
		float		sinb_h;			/* solar angle parameter */

		/* daily => hourly */
		float		tmp_sfc;			/* ground surface temperature, degree Celcius */
		float		tmp_2m;				/* 2m air temperature, degree Celcius */
		float		tmp10_soil;			/* soil temperature at 0-10 cm depth, degree Celcius */
		float		tmp200_soil;		/* soil temperature at 10-200 cm depth, degree Celcius */
		float		tmp300_soil;		/* soil temperature at 300 cm depth, degree Celcius */
		float		dswrf_sfc;			/* downward shortwave radiation at the surface, W m-2 */
		float		tcdc_clm;			/* total cloudiness, fraction */
		float		prate_sfc;			/* precipitation, mm day-1 */
		float		spfh_2m;			/* specific humidity, kg kg-1 */
		float		wnd_10m;			/* wind velocity, m s-1 */
    
        float       tair_dayav_c;
        float       tair_dayav;
		
		/* hourly calculated */
		float		sfcrad_h;		/* global radiation at the canopy-top, W m-2 */
		float		grdrad_h;		/* global radiation under the canopy, W m-2 */
		float		par_h;			/* hourly total photosynthetically active radiation (PAR), W m-2 */
		float		parb_h;			/* hourly beam PAR, W m-2*/
		float		pard_h;			/* hourly diffuse PAR, W m-2*/
		float		ppfd_h; 		/* hourly total photon flux density (PPFD), micro mol photon m-2 s-1 */
		float		ppfdb_h;		/* hourly beam PPFD, micro mol photon m-2 s-1 */
		float		ppfdd_h;		/* hourly diffuse PPFD, micro mol photon m-2 s-1 */
		float		ppfd_c_h;		/* hourly PPFD, constant conversion factor, micro mol photon m-2 s-1 */

		/* clear sky */
		float		cssfcrad_h; 	/* global radiation at the canopy-top, W m-2 */
		float		cspar_h; 		/* hourly total photosynthetically active radiation (PAR), W m-2 */
		float		csparb_h;		/* hourly beam PAR, W m-2*/
		float		cspard_h;		/* hourly diffuse PAR, W m-2*/
		float		csppfd_h;		/* hourly total photon flux density (PPFD), micro mol photon m-2 s-1 */
		float		csppfdb_h;		/* hourly beam PPFD, micro mol photon m-2 s-1 */
		float		csppfdd_h;		/* hourly diffuse PPFD, micro mol photon m-2 s-1 */
        
        float       fapar;

		/* atmosphere */
		long		CO2y;				/* year of atmospheric CO2 scenario */
		float		aCO2;				/* ambient (canopy air) CO2 concentration, ppmv */
		float		bCO2;				/* background (free atmosphere) CO2 concentration, ppmv */
		float		d13C_aCO2;			/* stable carbon isotope composition of ambient CO2, permille */
		float		d13C_bCO2;			/* stable carbon isotope composition of background CO2, permille */
		
		float		air_prsr;			/* air pressure, hPa */
		float		air_dns;			/* density of air, kg m-3 */
		float		lhvp;				/* latent heat MJ kg-1 */

		float		vp;					/* vapour pressure, hPa */
		float		vps;				/* saturation vapour pressure, hPa */
		float		slope_vps;			/* slope of saturation vapour pressure related to tempertaure, hPa deg C-1 */
		float		r_aero;				/* aerodynamic resistance, s m-1 */
		float		vpd;				/* vapour pressure deficit, hPa */
		
		/* radiation */
		float		albedo_sfc;			/* surface albedo, fraction */
		float		rn_long_eco;		/* net longwave radiation at the surface, W m-2 */
		float		rn_long_tree;		/* net longwave radiation for trees, W m-2 */
		float		rn_long_c3;			/* net longwave radiation for C3 grass, W m-2 */
		float		rn_long_c4;			/* net longwave radiation for C4 grass, W m-2 */
		float		rn_long_ground;		/* net longwave radiation for soil surface, W m-2 */

		float		rn_short_eco;		/* net shortwave radiation at the surface, W m-2 */
		float		rn_short_tree;		/* net shortwave radiation for trees, W m-2 */
		float		rn_short_c3;		/* net shortwave radiation for C3 grass, W m-2 */
		float		rn_short_c4;		/* net shortwave radiation for C4 grass, W m-2 */
		float		rn_short_ground;	/* net shortwave radiation for soil surface, W m-2 */

		float		rn_eco;				/* net shortwave radiation at the surface, W m-2 */
		float		rn_tree;			    /* net shortwave radiation of tree layer, W m-2 */
		float		rn_c3;				/* net shortwave radiation of lower layer C3, W m-2 */
		float		rn_c4;				/* net shortwave radiation of lower layer C4, W m-2 */
		float		rn_ground;			/* net shortwave radiation of soil ground, W m-2 */
		
		/* vegetation */		
		float		fcover_veg;			/* fractional vegetation cover */
		float		fcover_tree;		    /* fractional tree cover */
		float		fcover_c3;			/* fractional C3 grass cover */
		float		fcover_c4;			/* fractional C4 grass cover */
		float		fcover_ground;		/* fractional ground (soil surface) cover */
		float		funder_c3;			/* fractional C3 grass of understory */
		float		funder_c4;			/* fractional C4 grass of understory */

		float		lai;				/* leaf area index (LAI), m2 m-2 */
		/* canopy conductance, mmol m-2 s-1 */
		float		gc_tree;			/* tree */
		float		gc_c3;				/* C3 herbaceous */
		float		gc_c4;				/* C4 herbaceous */
    
        float       laimax_meta;            /* max LAI estimated by Iio-san's meta analysis */
		
		/* hydrology */
		/* Penman-Monteith potential intercepted evaporation rate, mm day-1 */
		float		pm_incep_tree;		/* tree */
		float		pm_incep_c3; 		    /* C3 herbaceous */
		float		pm_incep_c4; 		    /* C4 herbaceous */
		/* actual intercepted evaporation rate, tree, mm day-1 */
		float		incep;
		float		incep_tree;			/* tree */
		float		incep_c3;			/* C3 herbaceous */
		float		incep_c4;			/* C4 herbaceous */
		
		/* Penman-Monteith potential soil evaporation rate, mm day-1 */
		float		pm_evpr;			    /* potential evaporation from soil */
		/* actual soil evaporation rate, mm day-1 */		
		float		evpr;				/* actual evaporation */
		/* soil resistance to evaporation, s m-1 */
		float		rsst_soil;			/* resistance to evaporation */
		
		/* Penman-Monteith potential transpiration rate, tree, mm day-1 */
		float		pm_trnsp_tree; 		/* tree */
		float		pm_trnsp_c3; 		    /* C3 herbaceous */
		float		pm_trnsp_c4; 		    /* C4 herbaceous */
		/* actual transpiration rate, tree, mm day-1 */
		float		trnsp;              /* total */
		float		trnsp_tree; 		    /* tree */
		float		trnsp_c3;           /* C3 herbaceous */
		float		trnsp_c4;           /* C4 herbaceous */
		
		float		thaw;				/* snow thaw water, mm day-1 */
		float		ro1; 				/* monthly runoff of upper, mm */
		float		ro2; 				/* monthly runoff of lower, mm */

		float		snp; 				/* snow fraction of precipitation */
		float		snow_acc;			/* snow accumulation, mm water */
		float		soilwtr_l;			/* soil water upper layer, mm water */
		float		soilwtr_h;			/* soil water lower layer, mm water */
		float		soilappr_l;			/* upper soil aperture, fraction */
		float		soilappr_w;			/* lower soil aperture, fraction */
		
		float		pet;				    /* potential evapotranspiration */
		float		aet;				    /* actual evapotranspiration */
		float		wfps;				/* water-filled pore space */
		float		vmc;
		float		vmc_pre;

		float		iw_casa;			    /* CASA moisture index */
		float		rdr_casa;			/* CASA dryness index */
		float		e_casa;				/* CASA parameter */
		float		m_casa;				/* CASA parameter */
		float		m_casa_pre;			/* CASA parameter: previous time-step */
		
		/* soil water potential, MPa */
		/* litter */
		float		pot_total_l;		    /* total */
		float		pot_grav_l;			/* gravity */
		float		pot_matric_l;		    /* matric */
		/* humus */
		float		pot_total_h;		    /* total */
		float		pot_grav_h;			/* gravity */
		float		pot_matric_h;		    /* matric */
									 
		float		d_c_soil;
		float		d_c_min;
		float		d_n_soil;
		
		float		ndepo[3];			/* N deposition */
		float		depo_no3;
		float		depo_nh4;
		
		float		nfert_no3;
		float		nfert_nh4;
		
		/* 30-min observation data */
		long		metdat_frag;		/* flag of data availability */
		/* +1: temperature */
		/* +2: air humidity */
		/* +4: PAR */		/* i.e., 7, when all data are available */
		
		/* erosion */
		float		f_rain;				/* rain factor: dimensionless */
		float		f_slope;			    /* slope factor: dimensionless */
		float		f_erodibility;		/* soil factor: dimensionless */
		float		f_cover;			    /* vegetation cover factor: dimensionless */
		float		f_conservation;		/* conservation factor: dimensionless */
		
		/* greenhouse gas */
		float		water_table_depth;		/* water table depth, m */
		float		prof_ch4[SOILWET_LAYER+2];		/* soil CH4 concentration, micro M */
		float		npp_max;	/* maximum daily NPP, g C m-2 day-1 */
		float		tsoil_annav;	/* annual mean temperature */

        /* land-use change *******/
        float       deforest;
        float       deforest_v;
        float       deforest_s;
};

/* vegetation characteristics *************************************************************/
struct Pchar{ 
		/** radiation **/
		float		albedo; 		    /* reflectivity or albedo, fraction */
		float		ppfd_top; 		/* PPFD at the canopy-top, micro mol photon m-2 s-1 */
		float		ppfd_btm;		/* PPFD at the canopy-bottom, micro mol photon m-2 s-1 */

		float		ppfdb_top;			/* beam PPFD at top */
		float		ppfdd_top;			/* diffuse PPFD at top */
		float		ppfdb_btm;			/* beam PPFD at bottom */
		float		ppfdd_btm;			/* diffuse PPFD at bottom */
		
		/*** allocation ***/
		float		opt_lai; 		/* optimum leaf area index, m2 m-2 */
		float		alloc_ass;		/* allocation ration for assimilation organ, fraction */
		float		alloc_abg; 		/* allocation ration for aboveground non-assimilation organ, fraction */
		float		malloc_f;		/* to foliage, ratio */
		float		malloc_c;		/* to stem, ratio */
		float		malloc_r; 		/* to root, ratio */
		float		max_nsch_storage;	/* maximum capacity of carbohydrate storage */
		
		/*** phenology ***/
		long		phenoltype;
				/* 0: no vegetation, bareland */
				/* 1: evergreen */
				/* 2: cold deciduous */
		float		ctmp_lfdsp;			/* cumulative temperature for leaf display, degree-day */
		float		ctmp_lfshd;			/* cumulative temperature for leaf shedding, degree-day */
		/* phenology parameters moved to Mass structure ******************************************* 2009/04/20 A.Ito **/
	/*	float		gdd; 		*/		/* growing degree days */
	/*	float		cdd;		*/		/* coldness degree days */
	/*	float		grw_pd;		*/		/* growing period, days */
	/*	long		season;		*/		/* phenological stage as a function of season */
				/* 0: dormancy */
				/* 1: vegetative growth */
				/* 2: emergence of new leaf */
				/* 3: abandon of old leaf */
	/*	long		day_frash;	*/		/* day of year */
	/*	long		day_shed;	*/		/* day of year */
		
		/*** photosynthesis ***/
		long		phototype;			/* photosynthetic metabolic pathway, 3=C3, 4=C4, 5=CAM */
		float		pmax;				/* potential maximum rate, micro mol CO2 m-2 s-1 */
		float		ptop;				/* canopy-top photosynthetic rate, micro mol CO2 m-2 s-1 */
		float		psat;				/* light-saturated rate, micro mol CO2 m-2 s-1 */
		float		sla;				/* specific leaf area, cm2 g-1 dm */
		float		eK0;				/* light attenuation coefficient, no dimension */
		float		eK;					/* light attenuation coefficient, no dimension */
		float		lue0;				/* control light dependence coefficient, mol CO2 mol photon-1 */
		float		lue;				/* monthly quantum yield, mol CO2 mol photon-1 */ 
		float		topt;				/* optimum temperature, deg C */
		float		topt0;				/* optimum temperature (control), deg C */
		float		tmin;				/* minimum temperature, deg C */
		float		tmax;				/* maximum temperature, deg C */
		float		ci;					/* monthly intercellular CO2 concentration, ppmv */
		float		kmci;				/* dependence of photosynthesis on intercellular CO2 concentration, ppmv */
		float		cmpcd0;				/* CO2 compensation point, ppmv */
		float		cmpcd;				/* CO2 compensation point, ppmv */
		float		gs;					/* daily stomatal conductance, mmol H2O m-2 s-1 */
		float		gc;					/* daily canopy conductance, mmol H2O m-2 s-1 */
		float		gs_b0;				/* parameter of stomata model, mmol H2O m-2 s-1 */
		float		gs_b1;				/* parameter of stomata model (VPD dependence) */
		float		gs_b2;				/* parameters of stomatal conductance, mmol H2O m-2 s-1 */
		float		km_nstl;			/* non-stomatal (soil moinsture) photosynthesis limitation factor */
		
		/*** respiration ***/
		/* specific growth respiration rate, g C g C-1 alloc */
		float		rgf;				/* leaf */
		float		rgc;				/* stem */
		float		rgr; 				/* root */
		/* specific maintenance respiration rate at 15 degC, mg C g C-1 day-1 */
		float		rmf;				/* leaf */
		float		rmc;				/* stem */
		float		rmr; 				/* root */
		float		rmf0; 				/* leaf, basic */
		float		rmc_s;				/* stem sapwood */
		float		rmr_s; 				/* fine root */
		float		rmc_h;				/* stem heartwood */
		float		rmr_h;				/* coarse root */
		/* temperature dependence, dimensionless */
		float		qTf0;				/* leaf at 15 deg C */
		float		qTc0;				/* stem at 15 deg C */
		float		qTr0;				/* root at 15 deg C */
		float		qTf;				    /* leaf */
		float		qTc;				    /* stem */
		float		qTr; 				/* root */
		/* size-dependence parameter */
		float		f_sz_s;				/* stem */
		float		f_sz_r;				/* root */
		
		/*** specific litter fall ***/
		float		lf0;				/* leaf, basic */
		float		lc0;				/* stem, basic */
		float		lr0; 				/* root, basic */
		float		lf;					/* leaf */
		float		lc;					/* stem */
		float		lr; 				/* root */
		float		dcd; 				/* leaf shedding parameter, fraction */
		
		/*** root stratification ***/
		float		root_dpt_a; 		/* root profile parameter (Zeng 2001) */
		float		root_dpt_b; 		/* root profile parameter (Zeng 2001) */
		
		/*** photosynthetgic 13C fractionation ***/
		float		photo_13c_frac;		/* photosynthetic d13C fractionation */
		
		/* nitrogen */
		float		n_salvage;			/* N re-uptake from senescent leaves */
		float		cn0_strg;
		float		amax_nphoto;
		float		kn_nphoto;
		float		rd_n;
		float		n_conc_larea;
		float		n_conc_lmass;
		float		n_opt_leaf;

		/*** De Pury & Farquhar canopy photosynthesis model ***/
		long		photo_type;			/* 3=C3/4=C4 photosynthetic pathway */

		float		lma;				/* leaf mass per area */
		float		lai; 				/* monthly leaf area index(LAI), m2 m-2 */
		float		lai_sn;				/* sunny */
		float		lai_sd;				/* shade */

		float		vpd_a;				/* ambient VPD, Pa */
		float		co2_a;				/* ambient CO2 concentration, Pa */
		float		ci_sn;				/* leaf-internal CO2 concentration, Pa */
		float		ci_sd;				/* leaf-internal CO2 concentration, Pa */
		float		gs_sn;				/* stomatal conductance of sun leaves */
		float		gs_sd;				/* stomatal conductance of shade leaves */
		float		o2_i;				/* leaf-internal O2 concentration, Pa */
		float		tmp;				/* canopy temperature, K */
		float		comp_co2;			/* CO2 compensation point, Pa */
		float		ke_d;				/* light attenuation coefficient, diffused */
		float		ke_b1;				/* light attenuation coefficient, beam */
		float		ke_b2;				/* light attenuation coefficient, beam and scattered beam */
		float		scttr;				/* scattering factor */
		float		spect;				/* spectral correction factor (PS-II absortption) */
		float		rfl_d;				/* reflectance, diffuse */
		float		rfl_b;				/* reflectance, beam */
		float		appfd;				/* absorbed PPFD */
		float		appfd_sn;			/* absorbed PPFD by sun leaves */
		float		appfd_sd;			/* absorbed PPFD by shade leaves */
		float		apar;				/* absorbed PAR */
		float		apar_sn;			/* absorbed PAR by sun leaves */
		float		apar_sd;			/* absorbed PAR by shade leaves */
		float		ke_n;				/* attenuation coefficient, nitrogen allocation */
		float		nc_ct;				/* canopy-top N concentration, mmol m-2 */
		float		nc_np;				/* non-photosynthetic N concentration, mmol m-2 */
		float		vcmx_b;				/* max. carboxylation rate (base) */
		float		vcmx_con;			/* max. carboxylation rate (control) */
		float		vcmx;				/* max. carboxylation rate (actual) */
		float		vcmx_sn;			/* max. carboxylation rate (sunny leaves) */
		float		vcmx_sd;			/* max. carboxylation rate (shade leaves) */
		float		vcmx_sn0;			/* max. carboxylation rate (sunny leaves) at 25 degC */
		float		vcmx_sd0;			/* max. carboxylation rate (shade leaves) at 25 degC */
		float		photocap_n0;		/* non-photosynthetic nitrogen concentration */
		float		photocap_n;			/* N-based photosynthetic capacity */
		float		kc0;				/* kc at 25 deg C */
		float		kc;					/* Michaelis-Menten constant of carboxylation */
		float		ko0;				/* ko at 25 deg C */
		float		ko;					/* Michaelis-Menten constant of oxygenation */
		float		acen_kc;			/* activation energy of Kc */
		float		acen_ko;			/* activation energy of Ko */
		float		k_effc;				/* effective Michaelis-Menten constant */
		float		ft_et;				/* temperature factor of electron transport */
		float		acen_et;			/* activation energy of electron transport */
		float		hc_et;				/* */
		float		str_et;				/* */
		float		acen_vc;			/* activation energy of carboxylation */
		float		acen_rd;			/* activation energy of dark respiration */
		float		jvr;				/* Jmax/Vcmax ratio */
		float		cv_j;				/* convexity of J-APAR */
		float		jmx_sn;				/* Jmax of sunny */
		float		jmx_sd;				/* Jmax of shade */
		float		j_sn;				/* electron transport rate, sunny */
		float		j_sd;				/* electron transport rate, shade */
		float		cv_p;				/* convexity of assimulation rate */
		float		pc_sn;				/* rubisco-limited rate: sunny */
		float		pc_sd;				/* rubisco-limited rate: shade */
		float		pj_sn;				/* electron transport-limited rate: sunny */
		float		pj_sd;				/* electron transport-limited rate: shade */
		float		p_sn;				/* net assimilation rate of sun leaf */
		float		p_sd;				/* net assimilation rate of shade leaf */
		float		a_sn;				/* gross assimilation rate of sun leaf */
		float		a_sd;				/* gross assimilation rate of shade leaf */		
		float		gb;					/* leaf boundary-layer conductance */
        float       f_sw;               /* scalar of soil moisture */
    
		/* parameters of Ball-Berry stomatal conductance model */
		float		gss_b0; 
		float		gss_b1; 
		float		gss_b2;				/* parameters of stomatal conductance, mmol H2O m-2 s-1 */
		
		float		rd0;                /* canopy dark respiration at 25 degC */
		float		rd; 				/* canopy dark respiration */
		float		rd_sn;				/* canopy dark respiration: sunny */
		float		rd_sd;				/* canopy dark respiration: shade */
		float		ft_rd;              /* temperature scalar */
};

/* soil characteristics ************************************************************/
struct Schar{ 
		float		albedo0;			/* bare soil reflectivity, or albedo */
		float		albedo;				/* whole soil (with snow) reflectivity, or albedo */
		float		tmp;				/* soil temperature */
		
		/* decomposition coefficient, fraction */
		float		sr_lf;				/* foliage litter */
		float		sr_lc;				/* stem&branch litter */
		float		sr_lr;				/* root litter */
		float		sr_ha;				/* active humus */
		float		sr_hi;				/* intermediate humus */
		float		sr_hp;				/* passive humus */
		
		float		kml;				/* moisture dependence, fraction of soil water */
		float		kmh; 
		float		kmsl;				/* moisture dependence, fraction of soil water */
		float		kmsh; 

		float		f_co2_lf;			/* CO2 emission ratio to litter decomposition, fraction */
		float		f_co2_lc;
		float		f_co2_lr;
		
		float		f_hm_a;				/* mineral soil formation ratio to litter decomposition, fraction */
		float		f_hm_i;
		float		f_hm_p;
		
		/* combined environmental factor */
		float		f_tm_l;				/* litter */
		float		f_tm_h;				/* humus */
		
		/* maximum nitrate concentration */
		float 		ntrcon_max;			/* g N / kg H2O */
};

/* ecosystem characteristics *******************************************************/
struct Echar{ 
		char		para_ver_id[15];		/* parameter vesion ID */
		char		para_date_id[15];		/* parameter update date */

		struct		Pchar tree;			/* for tree plants */
		struct		Pchar c3;			/* for c3 herbaceous plants */
		struct		Pchar c4;			/* for c4 herbaceous plants */
		struct		Schar soil;			/* for soil organic matter */
};

/* plant biomass *********************************************************************/
struct Pmas{ 
		float		lai;				/* leaf area index, m2 m-2 */

		/* carbon mass, Mg C ha-1 */ 
		float		fol; 				/* foliage */
		float		stm; 				/* stem and branch */
		float		rot; 				/* root */
		
		float		stm_sp;				/* sub-compartment: sap */
		float		stm_ht;				/* sub-compartment: heart */
		float		rot_fn;				/* sub-compartment: fine */
		float		rot_tp;				/* sub-compartment: tap/coarse */
		
		float		nsch_storage;		/* storage of non-structural carbohydrate */
	
		float		plant;				/* total plant */
		
		float		gdd;				/* growing degree days */
		float		cdd;				/* coldness degree days */
		long		    season;             /* seasonal stage identifer */
		float		grwpd;              /* growing period, days */
		long		    day_frush;          /* day of leaf frush */
		long		    day_shed;           /* day of leaf shedding */
	

		/* stable carbon isotope composition, dimensionless */
		float		d13c_fol;			/* foliage */
		float		d13c_stm;			/* stem and branch */
		float		d13c_rot;			/* root */
		
		float		d13c_plant; 		/* total plant */

		/* nitrogen mass, g N ha-1 */ 
		float		n_cnpy;				/* foliage */
		float		n_strg;				/* stem and branch */

		float		n_stm_sp;			/* sub-compartment: sap */
		float		n_stm_ht;			/* sub-compartment: heart */
		float		n_rot_fn;			/* sub-compartment: fine */
		float		n_rot_tp;			/* sub-compartment: tap/coarse */

		float		n_plant;			/* total plant */
};

/* soil organic matter ***********************************************************/
struct Smas{ 
		/* carbon mass, Mg C ha-1 */ 
		float		ltr_tf;				/* foliage litter */
		float		ltr_tc;				/* stem litter */
		float		ltr_tr;				/* root litter */
		
		float		ltr_gf;				/* foliage litter */
		float		ltr_gc;				/* stem litter */
		float		ltr_gr;				/* root litter */
		
		float		ltr_tc_ag;			/* sub-compartment: aboveground */
		float		ltr_tc_bg;			/* sub-compartment: belowground */
		
		float		msl_a;				/* active humus */
		float		msl_i;				/* intermediate humus */
		float		msl_p;				/* passive humus */
		
		float		mcrb;				/* microbe */
		float		doc;				/* dissolved organic carbon, DOC */
		
		float		soil;				/* total soil */
		
        /* for SASU: 2013/12/03 by A.Ito */
        float       sasu_li_tf;
        float       sasu_li_tc;
        float       sasu_li_tr;
        float       sasu_li_gf;
        float       sasu_li_gc;
        float       sasu_li_gr;
		
        float       sasu_hf_ha;
        float       sasu_hf_hi;
        float       sasu_hf_hp;
		
        float       sasu_mr_tf;
        float       sasu_mr_tc;
        float       sasu_mr_tr;
        float       sasu_mr_gf;
        float       sasu_mr_gc;
        float       sasu_mr_gr;

        float       sasu_mr_ha;
        float       sasu_mr_hi;
        float       sasu_mr_hp;

		/* stable carbon isotope composition */
		float		d13c_ltr_f;		/* foliage */
		float		d13c_ltr_c;		/* stem */
		float		d13c_ltr_r;		/* root */
		
		float		d13c_msl_a;		/* active */
		float		d13c_msl_i;		/* intermediate */
		float		d13c_msl_p;		/* passive */

		float		d13c_soil;		/* total soil */

		/* nitrogen mass, g N ha-1 */
		/* inorganic */
		float		n_nh4;			/* ammonium */
		float		n_no3;			/* nitrate */

		/* organic */
		float		n_lttr;			/* foliage litter */
		float		n_hums;			/* active humus */
		float		n_mcrb;			/* microbe */
		
		float		n_soil;			/* total soil */
};

/* ecosystem carbon storage *************************************************/
struct Mass{ 
		/* evolving variables ( x(i) = x(i-1) + d(i) ) ************************************/
		float		snwa;				/* water equivalent snow depth, mm */
		float		sw30;				/* soil water content of above 30cm soil, mm */
		float		sww;				/* whole soil water content, mm */

        /* detritus in a grid produced by land-use change, Mg C ha-1 */
        float       c_ludetr_1;
        float       c_ludetr_10;
        float       c_ludetr_100;

		/* phenology */
        /*	long		gd, bbm; */		/* vegetative growing period, days */
        /*	float		gdd; */			/* cumulative temperature, growth degree days, degC days */
		
		/* pool */
		struct		Pmas tree;		/* tree plant mass */
		struct		Pmas c3;		/* floor C3 plant mass */
		struct		Pmas c4;		/* floor C4 plant mass */
		struct		Smas soil;		/* soil mass */

		/* carbon mass */
		float		total_c;		/* ecosystem total carbon */
		/* nitrogen mass */
		float		total_n;		/* ecosystem total nitrogen */

		/* stable carbon isotope composition */
		float		d13c_total;		/* ecosystem total carbon */
};

/* plant carbon fluxes ***************************************************************/
struct Pflx{ 
		/* carbon flux, Mg C ha-1 day-1 */ 
		float		gpp;				/* gross primary production */
		float		spp;				/* surplus production */
		float		epp;				/* effective production */
		float		npp;				/* net primary production */
		
		float		tpf;				/* translocation of photosynthate to foliage */
		float		tpc;				/* translocation of photosynthate to stem */
		float		tpr;				/* translocation of photosynthate to root */
		float		tpp;				/* translocation of photosynthate, total */
		
		float		rtpc;				/* re-allovcation from stem */
		float		rtpr;				/* re-allovcation from root */
		
		float		rp;					/* plant respiration, =rpm+rpg */
		float		rpg;				/* plant growth respiration */
		float		rpm;				/* plant maintenance respiration */
		float		rfg;				/* foliage growth respiration */
		float		rfm;				/* foliage maintenance respiration */
		float		rcg;				/* stem and branch growth respiration */
		float		rcm;				/* stem and branch maintenance respiration */
		float		rrg;				/* root growth respiration */
		float		rrm;				/* root maintenance respiration */
		
		float		lf;					/* foliage litterfall */
		float		lc;					/* stem and branch litterfall */
		float		lr;					/* root litterfall */
		float		lL;					/* total litterfall */
		
		float		hvst;				/* harvest of biomass */
		
		/* hourly */
		float		ar_h;				/* hourly autotrophic respiration rate */
		
		/* CH4 emit */
		float		emit_ch4_kirschbaum_mass;       /* leaf mass-based */
		float		emit_ch4_kirschbaum_photo;      /* photosynthesis-based */
		
		/* nitrogen uptake */
		/* symbiotic N fixation */
		float		n_biofix;           /* biological N fixation */
		/* root uptake */
		float		uptake_nh4;         /* nitrate, NO3- */
		float		uptake_no3;         /* ammonium, NH4+ */
		
		float		n_abdn_cnpy;
		float		n_abdn_strg;
		float		n_salvage;
		float		n_alloc_cnpy;
		float		n_alloc_strg;
		float		n_realloc;
		
		/* stable carbon isotope composition */
		float		d13c_gpp; 
		float		d13c_spp; 
		float		d13c_epp; 
		float		d13c_npp; 
		
		float		d13c_tpf; 
		float		d13c_tpc; 
		float		d13c_tpr; 
		float		d13c_tpp; 
		
		float		d13c_rp; 
		float		d13c_rpg; 
		float		d13c_rpm; 
		float		d13c_rfg; 
		float		d13c_rfm; 
		float		d13c_rcg; 
		float		d13c_rcm; 
		float		d13c_rrg; 
		float		d13c_rrm; 
		
		float		d13c_lf; 
		float		d13c_lc; 
		float		d13c_lr; 
		float		d13c_lL; 
		
		float		d13c_lf_c; 

		float		d13c_hvst; 
};

/* soil carbon fluxes, all in Mg C ha-1 day-1 *************************************/
struct Sflx{		 
		/* carbon flux ***************************************************/ 
		/* litter input */
		float		li_tf;				/* foliage litter */
		float		li_tc;				/* stem litter */
		float		li_tr;				/* root litter */
		
		float		li_gf;				/* foliage litter */
		float		li_gc;				/* stem litter */
		float		li_gr;				/* root litter */
		
		/* humification */
		float		hf_tfa;		/* */
		float		hf_tfi;
		float		hf_tfp;
		float		hf_tca;
		float		hf_tci;
		float		hf_tcp;
		float		hf_tra;
		float		hf_tri;
		float		hf_trp;

		float		hf_gfa;		/* */
		float		hf_gfi;
		float		hf_gfp;
		float		hf_gca;
		float		hf_gci;
		float		hf_gcp;
		float		hf_gra;
		float		hf_gri;
		float		hf_grp;

		/* mineralization */
		float		mr_tf;			/* tree leaf litter */
		float		mr_tc;			/* tree stem litter */
		float		mr_tr;			/* tree root litter */
		float		mr_gf;			/* herb leaf litter */
		float		mr_gc;			/* herb stem litter */
		float		mr_gr;			/* herb root litter */
		float		mr_a;			/* active humus */
		float		mr_i;			/* intermediate humus */
		float		mr_p;			/* passive humus */

		float		hr;				/* total microbial heterotrophic respiration */
		float		hr_h;			/* short-term heterotrophic respiration */
		
		float		doc_boyer;		/* DOC export by Boyer et al. (1996), gC/m2/day */
		float		erosion_soil;
		float		erosion_orgmat;
		float		erosion_carbon;
    
		/* methane oxydation */
		float		ch4oxy_ridg;		/* Ridgwell et al. (1999) */
		float		ch4oxy_casa;		/* methane oxydation */
		float		ch4oxy_delgrosso;		/* methane oxydation */
		float		ch4oxy_curry;		/* Curry (2007) */
		float		ch4gen;				/* methane generation */
		
		float		ch4_wh;
		float		ch4_wh_plant;
		float		ch4_wh_ebull;
		float		ch4_wh_diff;

		/* nitrogen flux **********************************************/
		/* N flow, g N ha-1 day-1 **/
		float	d_n2o_ntr_ngas;			/* N2O by nitrification, NGAS scheme */
		float	d_n2o_dnt_ngas;			/* N2O by denitrification, NGAS scheme */
		float	d_n2o_ngas;				/* total N2O, NGAS scheme */
		float	d_n2_ngas;				/* N2O by denitrification, NGAS scheme */
		
		float	d_n2o_ntr_ngas_2;			
		float	d_n2o_dnt_ngas_2;			
		float	d_n2o_ngas_2;			
		float	d_n2_ngas_2;			
		
		float	d_n2o_casa;			
		float	d_no_casa;			
		float	d_n2_casa;			
		
		float	n_nh3vlt;			
		float	n_leach;			
		float	n_minerlz_lttr;
		float	n_minerlz_hums;
		float	n_nitrif;			/* in NGAS */
		float	n_immbl;
		float	n_mcrb_abdn;
			
		/* stable carbon isotope composition, in per mil */
		float		d13c_li_f;		/* leaf litter input */
		float		d13c_li_c;		/* stem litter input */
		float		d13c_li_r;		/* root litter input */
		
		float		d13c_hf_f;
		float		d13c_hf_c;
		float		d13c_hf_r;

		float		d13c_mr_f;
		float		d13c_mr_c;
		float		d13c_mr_r;

		float		d13c_mr_a;
		float		d13c_mr_i;
		float		d13c_mr_p;

		float		d13c_hr; /* total microbial heterotrophic respiration */
};

/* ecosystem carbon fluxes *************************************/
struct Flux{ 
		struct		Pflx tree;		/* tree plant fluxes */
		struct		Pflx c3;		/* floor plant fluxes */
		struct		Pflx c4;		/* floor plant fluxes */
		struct		Sflx soil;		/* soil fluxes */

		/* ecosystem carbon budget */ 
		float		gpp;			/* ecosystem gross primary production */
		float		npp;			/* ecosystem net primary production */
		float		er;				/* ecosystem respiration */
		float		sr;				/* soil respiration */
		float		rr;				/* plant root respiration */
		float		nep;			/* net ecosystem production, Mg C ha-1 day-1 */
		float		nbp;			/* net biome production, Mg C ha-1 day-1 */
		float		necb;			/* net carbon balance of grid, Mg C ha-1 day-1 */

		/* VOC emission, g C m-2 day-1 */
		float		voc_isopr_g97;		/* Isoprene */
		float		voc_monotrp_g97;	/* Monoterpene */
		float		voc_methanl_g97;	/* Methanol */
		float		voc_acetone_g97;	/* Acetone */
		float		voc_actaldhd_g97;	/* Acetoaldehyde */
		float		voc_frmardhd_g97;	/* Formaldehyde */
		float		voc_formacd_g97;	/* Formacid */
		float		voc_acetacd_g97;	/* Acetoacid */
		float		voc_co_g97;			/* Carbon Monoxide */

        /* land-use change */
        float       lue_conv;           /* instantaneous at conversion */
        float       lue_detr_1;         /* within 1 year */
        float       lue_detr_10;        /* within 10 year */
        float       lue_detr_100;       /* within 100 year */

		/* stable carbon isotope composition */
		float		d13c_nep; 
		float		d13c_ncb; 
		float		d13c_er;
};

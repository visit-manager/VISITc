/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include"definition.h"

/* global variables ****************************/
/* configuration */
extern long WMODE;				/* working mode, site or region */
extern long WROW;				/* working row of 1st layer */
extern long WCOL;				/* working column of 1st layer */
extern long NROW;				/* number of row */
extern long NCOL;				/* number of column */
extern long YLDIST;			/* year of the last disturbance */
extern long SPUPT;				/* spin-up time (years) */

extern long BYR;				/* begin year */
extern long EYR;				/* end year */
extern long PERIOD;			/* period */

/* climatology: *_d[] means the time-series during 1948 to 2005 */	
extern long	month_day[12];
extern double	tmp_sfc_d[HCLIMD][YSTEP];		/* ground surface temperature, degree Celcius */
extern double	tmp_2m_d[HCLIMD][YSTEP];		/* 2m air temperature, degree Celcius */
extern double	tmp10_soil_d[HCLIMD][YSTEP];	/* soil temperature at 10 cm depth, degree Celcius */
extern double	tmp200_soil_d[HCLIMD][YSTEP];	/* soil temperature at 200 cm depth, degree Celcius */
extern double	dswrf_sfc_d[HCLIMD][YSTEP];	/* downward shortwave radiation at the surface, W m-2 */
extern double	tcdc_clm_d[HCLIMD][YSTEP];	/* total cloudiness, fraction */
extern double	prate_sfc_d[HCLIMD][YSTEP];	/* precipitation, mm mon-1 */
extern double	spfh_2m_d[HCLIMD][YSTEP];		/* specific humidity, kg kg-1 */
extern double	wind_10m_d[HCLIMD][YSTEP];	/* wind velocity, m s-1 */
extern double	vpd_d[HCLIMD][YSTEP]; 		/* VPD, hPa */

extern double	tmp_2m_ncep_dav[YSTEP];		/* 2m air temperature, degree Celcius */
extern double	prate_sfc_ncep_mav[12];		/* precipitation, mm mon-1 */

extern double	tmp_sfc_dav[YSTEP];	
extern double	tmp_2m_dav[YSTEP];	
extern double	tmp10_soil_dav[YSTEP];	
extern double	tmp200_soil_dav[YSTEP];	
extern double	prate_sfc_dav[YSTEP];	
extern double	dswrf_sfc_dav[YSTEP];	
extern double	spfh_2m_dav[YSTEP];	
extern double	tcdc_clm_dav[YSTEP];	

extern double	paddy_wtd[5][YSTEP];

extern double prec_altdata[HCLIMD][YSTEP];

/* GCM climate projection data */
extern double gcm_temp[121][YSTEP];
extern double gcm_prec[121][YSTEP];
extern double gcm_ahmd[121][YSTEP];
extern double gcm_dswrf[121][YSTEP];
extern double gcm_temp_av[YSTEP];
extern double gcm_prec_av[YSTEP];
extern double gcm_ahmd_av[YSTEP];
extern double gcm_dswrf_av[YSTEP];

/* LARS-generated dummy weather data */
extern double lars_tmin[63][YSTEP];
extern double lars_tmax[63][YSTEP];
extern double lars_prec[63][YSTEP];
extern double lars_srad[63][YSTEP];

/* sensitivity analysis */
extern short SENSANS_DIST;
/* 1: no disturbance */
extern short SENSANS_ACO2;
/* 1: no atmospheric CO2 rise */
extern short SENSANS_TEMP;
/* 1: no temperature change */
/* 2: use LARS-generated dammy data */
extern short SENSANS_PREC;
/* 1: no precipitation change */
/* 2: use LARS-generated dammy data */
extern short SENSANS_SRAD;
/* 1: no solar radiation change */
/* 2: use LARS-generated dammy data */
extern short N_LARS;

/* parameter sensitivity analysis */
extern short SA_PARA;
/* 1: photosyntehtic optimal temperature */
/* 2: photosyntehtic CO2 - Michaelis constant */
/* 3: photosynthetic soil water dependence */
/* 4: photosyntehtic light use efficiency */
/* 5: maintenance respiration temperature */
/* 6: decomposition temperature dependence */
/* 7: decomposition soil water dependence */
/* 8: disturbance intensity */

/* ensemble number */
extern short SA_PARA_EN;
extern double SA_PARA_VAR;

/* fix phenology ****************/
extern long fixlp_tree_season[YSTEP];
extern long fixlp_c3_season[YSTEP];
extern long fixlp_c4_season[YSTEP];

extern long fixlp_tree_dayflush[YSTEP];
extern long fixlp_c3_dayflush[YSTEP];
extern long fixlp_c4_dayflush[YSTEP];

extern long fixlp_tree_dayshed[YSTEP];
extern long fixlp_c3_dayshed[YSTEP];
extern long fixlp_c4_dayshed[YSTEP];

extern double fixlp_tree_lai[YSTEP];
extern double fixlp_c3_lai[YSTEP];
extern double fixlp_c4_lai[YSTEP];

/***** INITIALIZATION *****/
void f_setting(				struct Grid *grid);
void f_open_files(			struct Grid *grid, short rnum, FILE *fp_r[NFILE]);
void f_close_files(			FILE *fp_r[NFILE]);

void f_initialize(			struct Grid grid[DROW][DCOL], struct Loct *loct, struct Echar *echar, 
							struct Echar echar_type[17], struct Mass mass[DROW][DCOL], 
							struct Flux *flux, FILE *fp_r[NFILE]);
void f_init_site(				struct Grid *grid); 
void f_loct_proc(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void f_loct_init(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void atmco2_trend(				struct Grid *grid,	struct Loct* loct);
void set_parameter(			long mode, char [], struct Echar echar[]);
void f_init_region(			struct Grid grid[DROW][DCOL]);
void f_parameter_perturbation(short pnum);

/***** CLEARANCE *****/
void clear_a(				struct Grid *grid, struct Mass *mass);
void clear_b(				struct Loct *loct, struct Echar *echar, struct Flux *flux);

void pflux_zero(			struct Pflx *flux);
void sflux_zero(			struct Sflx *flux);
void delete_plant(			struct Pmas *mass);
void delete_soil(			struct Smas *mass);

/***** MASS & PARAMETERS INITIALIZATION *****/
void initTree(				struct Mass *mass);
void initFloor(				struct Mass *mass);
void initSoil(				struct Mass *mass);

/***** RADIATION *****/
double f_soldec(			struct Grid *grid, struct Loct *loct);
double f_solhgt(			struct Grid *grid, struct Loct *loct);
double f_daylen(			struct Grid *grid, struct Loct *loct);
double f_toprad(			struct Grid *grid, struct Loct *loct);
void f_ppfd(				struct Grid *grid, struct Loct *loct);
void f_net_rad(				struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass);
double albedo_soil(			struct Loct *loct, struct Schar *schar);

/***** HYDROLOGY *****/
void f_hydrology(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass);
double f_airdens(			struct Grid *grid, struct Loct *loct);
double f_vap_pre_sat(		struct Grid *grid, struct Loct *loct);
double f_slope_vps(			struct Grid *grid, struct Loct *loct);
double f_r_aero(			struct Grid *grid, struct Loct *loct);
double pm_evap(				struct Grid *grid, struct Loct *loct);
double pm_transp(			struct Grid *grid, struct Loct *loct, long canopy);
double pm_incep(			struct Grid *grid, struct Loct *loct, long canopy);
void f_soil_saxton(			struct Grid *grid);

/***** ECOPHYSIOLOGY *****/
void f_ecophysiology(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass);
void f_incelco2(			struct Loct *loct, struct Pchar *pchar);
void f_photo_qy(			struct Loct *loct, struct Pchar *pchar);
void f_opt_lai(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
double lai_mass(			struct Pmas *mass, struct Pchar *pchar);
double irr_attn(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void f_q10_ar(				struct Loct *loct, struct Pchar *pchar);
void f_spcfc_resp(			struct Pchar *pchar, struct Pmas *mass);
void f_mortality(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void f_stom_cond(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
double f_canopy_cond(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass);
void f_pc_sat(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar);

void f_leaf_aging(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass);
void f_leaf_prop(			struct Grid *grid, struct Loct *loct, struct Cchar *cchar);
void f_apar_sunshade(		struct Grid *grid,	struct Loct *loct, struct Cchar *cchar);
void f_gpp_sunshade(		long stage, struct Grid *grid, struct Loct *loct, struct Cchar *cchar);
double f_gas_sun(			struct Cchar *cchar, double ci);
double f_gas_shade(			struct Cchar *cchar, double ci);
void leaf_ansolv(			long frac, struct Cchar *cchar);
void f_sla_change(			struct  Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass);

void f_resp_h(				struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);

/***** PHENOLOGICAL CYCLES *****/
void f_growth_period(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology(					struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_evergreen(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_colddeciduous(	struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_ariddeciduous(	struct Grid *grid, struct Pchar *pchar);
void phenology_bareland(		struct Grid *grid, struct Pchar *pchar);
void phenology_grass(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_crop(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar);

/***** C FLUXES *****/
double f_gpp(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
double f_rfm(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
double f_rfg(		struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
double f_rcm(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
double f_rcg(		struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
double f_rrm(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
double f_rrg(		struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
double f_lf(		struct Grid *grid, struct Pchar *pchar, struct Pmas *mas);
double f_lc(		struct Grid *grid, struct Pchar *pchar, struct Pmas *mas);
double f_lr(		struct Grid *grid, struct Pchar *pchar, struct Pmas *mas);

double f_gpp_df97(	struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
					struct Pmas *pmas, struct Pflx *pflx);
void leaf_ansolv_df97(		double vpd, double aco2, double lai, double ft_rd, double vcmax, 
							double k_effc, double comp_co2, double j_et, double *assim,
							double *rdark, double *stocon, double *inco2, double *debug1, 
							double *debug2, double *debug3);

void f_allocation(			struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void reallocation_survival(	struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);

double frl(					struct Grid *grid, struct Loct *loct, struct Schar *schar);
double frh(					struct Grid *grid, struct Loct *loct, struct Schar *schar);

/***** CROP CYCLES *****/
void planting(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);
void harvesting(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);
void interval(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);

/***** SCHEMES *****/
void f_spinup(				struct Grid grid[DROW][DCOL], struct Loct *loct, struct Echar *echar, 
							struct Echar echar_type[17], struct Mass mass[DROW][DCOL], 
							struct Flux *flux, FILE *fp_r);
void f_experiment(			struct Grid grid[DROW][DCOL], struct Loct *loct, struct Echar *echar, 
							struct Echar echar_type[17], struct Mass mass[DROW][DCOL], 
							struct Flux *flux, FILE *fp_r[NFILE]);
void daily_scheme(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);

void f_flux_site(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux, FILE *fp_o, FILE *fp_r[NFILE]);

void disturbance_regime(	long year, struct Grid *grid, struct Loct *loct, struct Echar *echar,
							struct Mass *mass, struct Flux *flux);
void logging_event(			struct Grid *grid, struct Mass *mass);
void plant_process(			struct Grid *grid, struct Loct* loct, struct Pflx *flux,
							struct Pchar *pchar, struct Pmas *mass);

void f_doyTmody(			long year, long doy, long *month, long *day);
void f_modyTdoy(			long month, long day, long *doy);
void f_clim_correct_TKY(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_TMK(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_TSE(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_LSH(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_MKL(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_SKT(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_FJY(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_QHB(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_KBU(	struct Grid *grid, struct Loct* loct);
void f_clim_correct_Tongyu(	struct Grid *grid, struct Loct* loct);

/***** CARBON DYNAMICS *****/
void f_cycle_soil(			struct Grid *grid, struct Loct *loct, struct Schar *schar, 
							struct Smas *mass,	struct Sflx *flux);
void f_agriculture(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);

/*** BIOGEOCHEMISTRY ***/
void f_casa_moisture(		struct Grid *grid, struct Loct *loct);
void f_ch4oxy_ridgewell(	struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
void f_ch4oxy_curry(		struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
void f_ch4oxy_casa(			struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4oxy_delgrosso(	struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_n2o_emit_ngas(		struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
void f_n2o_emit_ngas_2(		struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
void f_n2o_emit_casa(		struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
void f_ch4emit_plant(		struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void f_ch4emit_walter(		struct Grid *grid, struct Loct *loct, struct Flux *flux);

void n_budget(				struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
void f_biolfix(				struct Grid *gsrid, struct Loct *loct, struct Mass *mass, 
							struct Flux *sflux);
void f_nh3_volatilization(	struct Grid *grid, struct Loct *loct, struct Smas *mass, 
							struct Sflx *flux);
void f_n_leaching(			struct Grid *grid, struct Loct *loct, struct Smas *mass, 
							struct Sflx *flux);

void f_n_uptake(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void f_n_abandon_salvage(	struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);
void f_n_immoblz(			struct Grid *grid, struct Loct *loct, struct Schar *schar, 
							struct Smas *mass, struct Sflx *flux);
void f_n_mcrb_abdn(			struct Grid *grid, struct Loct *loct, struct Schar *schar, 
							struct Smas *mass, struct Sflx *flux);
void f_n_alloc(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);
void f_n_realloc(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);
void f_n_conc(				struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);
void f_n_mineralz(			struct Grid *grid, struct Loct *loct, struct Smas *mass, 
							struct Sflx *flux);
void f_doc_boyer(			struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
void f_voc_emit_guenther97(	struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void f_erosion_rusle(		struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);

/*** OUTPUT ***/
void output_ansis_daily(	struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux, FILE *fp);
void f_ansis_ann(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux, double ansis_ann[N_ANSIS]);
void output_ansis_ann(		long year, double ansis_ann[N_ANSIS], FILE *fp);

/***** STABLE CARBON ISOTOPE *****/
void f_init_d13c(				struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void photo_13c_frac(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
double deltaTratio(			double delta);
double ratioTdelta(			double ratio);
double d13c_addition(		double d13c_a, double mass_a, double d13c_b, double mass_b);
double d13c_addition3(		double d13c_a, double mass_a, double d13c_b, double mass_b, 
							double d13c_c, double mass_c);
void co2_in_canopy(			struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);	
void d13c_efflux(			struct Grid *grid, struct Loct *loct, struct Flux *flux);

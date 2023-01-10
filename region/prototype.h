/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2012 / 11 / 30 by Akihiko Ito: Version 0.2.0.
*/

#include"definition.h"

/* **** CLEARANCE *****/
void clear_a(				struct Grid *grid, struct Mass *mass);
void clear_b(				struct Loct *loct, struct Echar *echar, struct Flux *flux);

void pflux_zero(			struct Pflx *flux);
void sflux_zero(			struct Sflx *flux);
void delete_plant(			struct Pmas *mass);
void delete_soil(			struct Smas *mass);

/* **** INITIALIZATION *****/
void f_setting(				struct Grid *grid);
void f_terminate(			FILE *fp_r[N_OFILE]);

void f_initialize(			struct Grid grid[], struct Loct *loct, struct Echar *echar, 
							struct Echar echar_type[MAX_BIOME], struct Mass mass[], 
							struct Flux *flux, FILE *fp_r[N_OFILE]);
void init_site(				struct Grid *grid); 
void f_init_cond(					struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void f_init_loct(					struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void atmco2_trend(			struct Grid *grid,	struct Loct* loct);
void set_parameter(			long mode, char [], struct Echar echar[]);
void f_init_global_run(		struct Grid grid[]);
void f_init_bamiyan_run(		struct Grid grid[]);
void f_init_easia_run(          struct Grid grid[]);
void f_init_japan_run(          struct Grid grid[]);
void f_init_bb_run(          struct Grid grid[]);
void f_init_japanc_run(          struct Grid grid[]);
void f_init_japanh_run(          struct Grid grid[]);
void f_init_japank_run(          struct Grid grid[]);
void f_open_global_clim(      short stage, short yr, short mon, FILE *fp_clim[6]);
void f_open_bamiyan_clim(      short stage, short yr, FILE *fp_clim[6]);
void f_open_easia_clim(      short stage, short yr, short mon, FILE *fp_clim[6]);
void f_open_japan_clim(      short stage, short yr, short mon, short day, short hour, FILE *fp_clim[6]);
void f_open_bb_clim(        short stage, short yr, short mon, short day, short hour, FILE *fp_clim[6]);
void f_open_japanc_clim(      short stage, short yr, short mon, short day, short hour, FILE *fp_clim[6]);
void f_open_japanh_clim(      short stage, short yr, short mon, short day, short hour, FILE *fp_clim[6]);
void f_open_japank_clim(      short stage, short yr, short mon, short day, short hour, FILE *fp_clim[6]);
void f_create_dfile_global(     short phase, short year, struct Grid *grid,
                            char filename[128], FILE *fp[N_OFILE]);
float grid_area(			    float lat1, float lat2, float lon1, float lon2);

/* **** MASS & PARAMETERS INITIALIZATION *****/
void initTree(				struct Mass *mass);
void initFloor(				struct Mass *mass);
void initSoil(				struct Mass *mass);

/***** RADIATION *****/
float f_soldec(				struct Grid *grid, struct Loct *loct);
float f_solhgt(				struct Grid *grid, struct Loct *loct);
float f_daylen(				struct Grid *grid, struct Loct *loct);
float f_toprad(				struct Grid *grid, struct Loct *loct);
void f_ppfd(				struct Grid *grid, struct Loct *loct);
void f_net_rad(				struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass);
float albedo_soil(			struct Loct *loct, struct Schar *schar);

/* **** HYDROLOGY *****/
void f_hydrology(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
                            struct Mass *mass);
float f_airdens(			struct Grid *grid, struct Loct *loct);
float f_vap_pre_sat(		struct Grid *grid, struct Loct *loct);
float f_slope_vps(			struct Grid *grid, struct Loct *loct);
float f_r_aero(				struct Grid *grid, struct Loct *loct);
float pm_evap(				struct Grid *grid, struct Loct *loct);
float pm_transp(			struct Grid *grid, struct Loct *loct, long canopy);
float pm_incep(				struct Grid *grid, struct Loct *loct, long canopy);
void f_soil_saxton(			struct Grid *grid);

/* **** ECOPHYSIOLOGY *****/
void f_ecophysiology(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
                            struct Pmas *mass);
void f_incelco2(			struct Loct *loct, struct Pchar *pchar);
void f_photo_qy(			struct Loct *loct, struct Pchar *pchar);
void f_opt_lai(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
float lai_mass(				struct Pmas *mass, struct Pchar *pchar);
float irr_attn(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void f_q10_ar(				struct Loct *loct, struct Pchar *pchar);
void f_spcfc_resp(			struct Pchar *pchar, struct Pmas *mass);
void f_mortality(			struct Grid *grid, struct Pchar *pchar);
void stom_cond(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
float canopy_cond(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
                            struct Pmas *mass);
void f_pc_sat(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar);

void f_leaf_aging(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
                            struct Mass *mass);
void f_leaf_prop(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void f_apar_sunshade(		struct Grid *grid,	struct Loct *loct, struct Pchar *pchar);
void f_gpp_sunshade(		long stage, struct Grid *grid, struct Loct *loct, 
                            struct Pchar *pchar);
float f_gas_sun(			struct Pchar *pchar, float ci);
float f_gas_shade(			struct Pchar *pchar, float ci);
void leaf_ansolv(			long frac, struct Pchar *pchar);
void f_sla_change(			struct Loct *loct, struct Echar *echar, struct Mass *mass);

void f_resp_h(				struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);

/* **** PHENOLOGICAL CYCLES *****/
void growthperiod(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
                                struct Pmas *mass);
void phenology_evergreen(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
                                struct Pmas *mass);
void phenology_colddeciduous(	struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
                                struct Pmas *mass);
void phenology_bareland(		struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);
void phenology_grass(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_agriculture(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_ariddeciduous(	struct Grid *grid, struct Pchar *pchar);

/* **** C FLUXES *****/
float f_gpp(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
float f_rfm(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
float f_rfg(		struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
float f_rcm(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
float f_rcg(		struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
float f_rrm(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
float f_rrg(		struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
float f_lf(			struct Grid *grid, struct Pchar *pchar, struct Pmas *mas);
float f_lc(			struct Grid *grid, struct Pchar *pchar, struct Pmas *mas);
float f_lr(			struct Grid *grid, struct Pchar *pchar, struct Pmas *mas);

void allocation(			struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void reallocation_survival(	struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);

float frl(					struct Grid *grid, struct Loct *loct, struct Schar *schar);
float frh(					struct Grid *grid, struct Loct *loct, struct Schar *schar);

/* **** CROP CYCLES *****/
void planting(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);
void harvesting(			struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);
void interval(				struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
							struct Pmas *mass, struct Pflx *flux);

/* **** SCHEMES *****/
void f_spinup(				struct Grid grid[], struct Loct *loct, struct Echar *echar, 
							struct Echar echar_type[MAX_BIOME], struct Mass mass[], 
                            struct Flux *flux, 
							FILE *fp_r);
void f_experiment(			struct Grid grid[], struct Loct *loct, struct Echar *echar, 
							struct Echar echar_type[MAX_BIOME], struct Mass mass[], 
                            struct Flux *flux, 
							FILE *fp_r[N_OFILE]);
void f_daily_scheme(		struct Grid *grid, struct Loct *loct, struct Echar *echar, 
                            struct Mass *mass, 
							struct Flux *flux);

void f_flux_TKY(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux, FILE *fp_o, 
                            FILE *fp_r[N_OFILE]);

void disturbance_regime(	long year, struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
void logging_event(			struct Grid *grid, struct Mass *mass);
void grazing_event(			struct Grid *grid, struct Loct* loct, struct Mass *mass);

void plant_process(			struct Grid *grid, struct Loct* loct, struct Pflx *flux,
							struct Pchar *pchar, struct Pmas *mass);
void f_cycle_soil(			struct Grid *grid, struct Loct *loct, struct Schar *schar, 
							struct Smas *mass,	struct Sflx *flux);
void f_agriculture(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);

void doyTmody(				long year, long doy, long *month, long *day);
void modyTdoy(				long month, long day, long *doy);

/* ** BIOGEOCHEMISTRY ***/
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
void f_biolfix(				struct Grid *gsrid, struct Loct *loct, struct Flux *sflux);
void f_nh3_volatilization(	struct Grid *grid, struct Loct *loct, struct Smas *mass, 
                            struct Sflx *flux);
void f_n_leaching(			struct Grid *grid, struct Loct *loct, struct Smas *mass, 
                            struct Sflx *flux);

void f_n_uptake(			struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);
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
void f_n_leaf_conc(			struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);
void f_n_mineralz(			struct Grid *grid, struct Loct *loct, struct Smas *mass, 
                            struct Sflx *flux);

void f_erosion_rusle(		struct Grid *grid, struct Loct *loct, struct Mass *mass, 
                            struct Flux *flux);
void f_doc_boyer(			struct Grid *grid, struct Loct *loct, struct Mass *mass, 
                            struct Flux *flux);

void f_voc_emit_guenther97(	struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);

/* ** OUTPUT ***/
void output_ansis_daily(	struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux, FILE *fp);
void f_ansis_ann(			struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux, float ansis_ann[MAX_ANSIS]);
void output_ansis_ann(		long year, float ansis_ann[MAX_ANSIS], FILE *fp);

/***** STABLE CARBON ISOTOPE *****/
void init_d13c(				struct Grid *grid, struct Loct *loct, struct Echar *echar, 
							struct Mass *mass, struct Flux *flux);
void photo_13c_frac(		struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
float deltaTratio(			float delta);
float ratioTdelta(			float ratio);
float d13c_addition(		float d13c_a, float mass_a, float d13c_b, float mass_b);
float d13c_addition3(		float d13c_a, float mass_a, float d13c_b, float mass_b, 
							float d13c_c, float mass_c);
void co2_in_canopy(			struct Grid *grid, struct Loct *loct, struct Mass *mass, 
							struct Flux *flux);	
void d13c_efflux(			struct Grid *grid, struct Loct *loct, struct Flux *flux);

/* land use change ************************************************************/
void f_open_landuse_cfsr(    FILE *fp_lu[N_LUDATA], long yr_offset, struct Grid grid[]);
void f_close_landuse_cfsr(   FILE *fp_lu[N_LUDATA]);
void f_read_landuse_cfsr(    short yr, struct Grid grid[],FILE *fp_lu[N_LUDATA]);
void f_conv_landuse(            struct Grid *grid, struct Loct *loct);
void f_emit_landuse(            struct Grid *grid, struct Loct *loct, struct Mass *mass, 
                                struct Flux *flux);



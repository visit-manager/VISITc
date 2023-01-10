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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* climatology: *_d[] means the time-series during 1948 to 2004- */	
extern float	tmp_sfc_d[PERIOD][366];		/* ground surface temperature, degree Celcius */
extern float	tmp_2m_d[PERIOD][366];		/* 2m air temperature, degree Celcius */
extern float	tmp10_soil_d[PERIOD][366];	/* soil temperature at 10 cm depth, degree Celcius */
extern float	tmp200_soil_d[PERIOD][366];	/* soil temperature at 200 cm depth, degree Celcius */
extern float	dswrf_sfc_d[PERIOD][366];	/* downward shortwave radiation at the surface, W m-2 */
extern float	tcdc_clm_d[PERIOD][366];		/* total cloudiness, fraction */
extern float	prate_sfc_d[PERIOD][366];	/* precipitation, mm mon-1 */
extern float	spfh_2m_d[PERIOD][366];		/* specific humidity, kg kg-1 */
extern float	wind_10m_d[PERIOD][366];	/* u-wind velocity, m s-1 */
extern float	vpd_d[PERIOD][366]; 		/* VPD, hPa */
extern float	tmp_sfc_dav[366];	
extern float	tmp_2m_dav[366];	
extern float	tmp10_soil_dav[366];	
extern float	tmp200_soil_dav[366];	
extern float	prate_sfc_dav[366];	
extern float	dswrf_sfc_dav[366];	
extern float	spfh_2m_dav[366];	
extern float	tcdc_clm_dav[366];

/* clear all variables and parameters ******************************/
void clear_a(
	struct Grid *grid, 
	struct Mass *mass
){
	mass->snwa = mass->sww = mass->sw30 = 0.0;
	delete_plant(&(mass->tree));
	delete_plant(&(mass->c3));
	delete_plant(&(mass->c4));
	delete_soil(&(mass->soil));
	mass->total_c = 0.0;
	mass->total_n = 0.0;
	
	grid->time = 0;
}

/****************************************************************/
void clear_b(
	struct Loct *loct, 
	struct Echar *echar, 
	struct Flux *flux
){
	long e, f;
	
	loct->climy = loct->doy = 0;
	
	for(f=BYR; f<=EYR; f++){
		for(e=0; e<366; e++){
			tmp_sfc_d[f-BYR][e] = 0.0; 
			tmp_2m_d[f-BYR][e] = 0.0; 
			tmp10_soil_d[f-BYR][e] = 0.0; 
			tmp200_soil_d[f-BYR][e] = 0.0; 
			dswrf_sfc_d[f-BYR][e] = 0.0; 
			tcdc_clm_d[f-BYR][e] = 0.0; 
			prate_sfc_d[f-BYR][e] = 0.0; 
			spfh_2m_d[f-BYR][e] = 0.0; 
			wind_10m_d[f-BYR][e] = 0.0;
			vpd_d[f-BYR][e] = 0.0; 
		}
	}
	
	for(e=0;e<366;e++){
		tmp_sfc_dav[e] = 0.0;
		tmp_2m_dav[e] = 0.0;
		tmp10_soil_dav[e] = 0.0;
		tmp200_soil_dav[e] = 0.0;
		prate_sfc_dav[e] = 0.0;
		dswrf_sfc_dav[e] = 0.0;
		spfh_2m_dav[e] = 0.0;
		tcdc_clm_dav[e] = 0.0;
	}
	
	loct->soldec = 0.0;
	loct->daylen = 0.0;
	loct->solhgt_h = 0.0;
	loct->toprad_h = 0.0;
	loct->sfcrad_h = 0.0;
	loct->grdrad_h = 0.0;
	loct->par_h = 0.0;
	loct->parb_h = 0.0;
	loct->pard_h = 0.0;
	loct->ppfd_h = 0.0;
	loct->ppfdb_h = 0.0;
	loct->ppfdd_h = 0.0;
		
	pflux_zero(&(flux->tree));
	pflux_zero(&(flux->c3));
	pflux_zero(&(flux->c4));
	sflux_zero(&(flux->soil));
	flux->nep = 0.0;
	flux->necb = 0.0;
	flux->voc_isopr_g97 = 0.0;
	flux->voc_monotrp_g97 = 0.0;
	flux->voc_methanl_g97 = 0.0;
	flux->voc_acetone_g97 = 0.0;
	flux->voc_actaldhd_g97 = 0.0;
	flux->voc_frmardhd_g97 = 0.0;
	flux->voc_formacd_g97 = 0.0;
	flux->voc_acetacd_g97 = 0.0;
	flux->voc_co_g97 = 0.0;
	flux->er = 0.0;
	flux->sr = 0.0;
	flux->rr = 0.0;
	
	loct->lai = 0.0;
	loct->aCO2 = 0.0;
	loct->bCO2 = 0.0;
	loct->d13C_aCO2 = 0.0;
	loct->d13C_bCO2 = 0.0;

	loct->d13C_aCO2 = 0.0;
	loct->rn_long_eco = 0.0;
	loct->rn_long_tree = 0.0;
	loct->rn_long_c3 = 0.0;
	loct->rn_long_c4 = 0.0;
	loct->rn_long_ground = 0.0;
	loct->rn_short_eco = 0.0;
	loct->rn_short_tree = 0.0;
	loct->rn_short_c3 = 0.0;
	loct->rn_short_c4 = 0.0;
	loct->rn_short_ground = 0.0;
	loct->rn_eco = 0.0;
	loct->rn_tree = 0.0;
	loct->rn_c3 = 0.0;
	loct->rn_c4 = 0.0;
	loct->rn_ground = 0.0;
	loct->air_prsr = 0.0;
	loct->air_dns = 0.0;
	loct->vp = 0.0;
	loct->vps = 0.0;
	loct->vpd = 0.0;
	loct->slope_vps = 0.0;
	loct->r_aero = 0.0;
	loct->pm_incep_tree = 0.0;
	loct->pm_incep_c3 = 0.0;
	loct->pm_incep_c4 = 0.0;
	loct->incep = 0.0;
	loct->incep_tree = 0.0;
	loct->incep_c3 = 0.0;
	loct->incep_c4 = 0.0;
	loct->pm_evpr = 0.0; 		
	loct->evpr = 0.0;			
	loct->pm_trnsp_tree = 0.0;	
	loct->pm_trnsp_c3 = 0.0; 	
	loct->pm_trnsp_c4 = 0.0; 	
	loct->trnsp = 0.0;
	loct->trnsp_tree = 0.0; 	
	loct->trnsp_c3 = 0.0; 	
	loct->trnsp_c4 = 0.0; 	
	loct->thaw = 0.0;
	loct->ro1 = 0.0;
	loct->ro2 = 0.0;
	loct->depo_no3 = 0.0;
	loct->depo_nh4 = 0.0;
	
	(echar->tree).psat = (echar->c3).psat = (echar->c4).psat = 0.0;
	(echar->tree).eK = (echar->c3).eK = (echar->c4).eK = 0.0;
	(echar->tree).lue = (echar->c3).lue = (echar->c4).lue = 0.0;
	(echar->tree).gs = (echar->c3).gs = (echar->c4).gs = 0.0;
	(echar->tree).gc = (echar->c3).gc = (echar->c4).gc = 0.0;
	(echar->tree).ci = (echar->c3).ci = (echar->c4).ci = 0.0;
	(echar->tree).opt_lai = (echar->c3).opt_lai = (echar->c4).opt_lai = 0.0;
	(echar->tree).qTf = (echar->c3).qTf = (echar->c4).qTf = 0.0;
	(echar->tree).qTc = (echar->c3).qTc = (echar->c4).qTc = 0.0;
	(echar->tree).qTr = (echar->c3).qTr = (echar->c4).qTr = 0.0;
	(echar->tree).lf = (echar->c3).lf = (echar->c4).lf = 0.0;
	(echar->tree).lc = (echar->c3).lc = (echar->c4).lc = 0.0;
	(echar->tree).lr = (echar->c3).lr = (echar->c4).lr = 0.0;
	(echar->tree).malloc_f = (echar->c3).malloc_f = (echar->c4).malloc_f = 0.0;
	(echar->tree).malloc_c = (echar->c3).malloc_c = (echar->c4).malloc_c = 0.0;
	(echar->tree).malloc_r = (echar->c3).malloc_r = (echar->c4).malloc_r = 0.0;
	
	(echar->soil).albedo = 0.0;
}

/* make the biome type zero *********************/
void delete_plant(
	struct Pmas *mass
){
	/* C */
	mass->lai = 0.0;
	mass->fol = 0.0;
	mass->stm = 0.0;
	mass->rot = 0.0;
	mass->nsch_storage = 0.0;
	mass->plant = 0.0;
	
	mass->grwpd = 0.0;
	mass->gdd = 0.0;
	mass->cdd = 0.0;
	mass->season = 0;
	mass->day_frush = 0;
	mass->day_shed = 0;
	
	/* N */
	mass->n_cnpy = 0.0;
	mass->n_strg = 0.0;
	mass->n_plant = 0.0;
}

/* make the biome type zero **********************/
void delete_soil(
	struct Smas *mass
){
	/* C */
	mass->ltr_tf = 0.0;
	mass->ltr_tc = 0.0;
	mass->ltr_tr = 0.0;
	mass->ltr_gf = 0.0;
	mass->ltr_gc = 0.0;
	mass->ltr_gr = 0.0;
	mass->msl_a = 0.0;
	mass->msl_i = 0.0;
	mass->msl_p = 0.0;
	mass->soil = 0.0;
	
    mass->sasu_li_tf = 0.0;
    mass->sasu_li_tc = 0.0;
    mass->sasu_li_tr = 0.0;
    mass->sasu_li_gf = 0.0;
    mass->sasu_li_gc = 0.0;
    mass->sasu_li_gr = 0.0;

    mass->sasu_hf_ha = 0.0;
    mass->sasu_hf_hi = 0.0;
    mass->sasu_hf_hp = 0.0;

    mass->sasu_mr_tf = 0.0;
    mass->sasu_mr_tc = 0.0;
    mass->sasu_mr_tr = 0.0;
    mass->sasu_mr_gf = 0.0;
    mass->sasu_mr_gc = 0.0;
    mass->sasu_mr_gr = 0.0;

    mass->sasu_mr_ha = 0.0;
    mass->sasu_mr_hi = 0.0;
    mass->sasu_mr_hp = 0.0;

	/* N */
	mass->n_lttr = 0.0;
	mass->n_hums = 0.0;
	mass->n_mcrb = 0.0;
	mass->n_nh4 = 0.0;
	mass->n_no3 = 0.0;
	mass->n_soil = 0.0;
}

/* make plant fluxes vacant **********************/
void pflux_zero(
	struct Pflx *flux
){
	flux->gpp = 0.0; 
	flux->epp = 0.0; 
	flux->spp = 0.0; 
	flux->npp = 0.0; 
	
	flux->tpp = 0.0; 
	flux->tpf = 0.0; 
	flux->tpc = 0.0; 		
	flux->tpr = 0.0; 	

	flux->rtpc = 0.0; 		
	flux->rtpr = 0.0; 	
	
	flux->rp = 0.0; 
	flux->rfm = 0.0; 
	flux->rcm = 0.0; 
	flux->rrm = 0.0; 
	flux->rpm = 0.0; 
	flux->rfg = 0.0; 
	flux->rcg = 0.0; 
	flux->rrg = 0.0; 	
	flux->rpg = 0.0; 	
			
	flux->lf = 0.0; 
	flux->lc = 0.0; 
	flux->lr = 0.0; 
	flux->lL = 0.0; 
	
	flux->hvst = 0.0; 

	flux->n_biofix = 0.0; 
	flux->uptake_nh4 = 0.0;		/* root uptake */
	flux->uptake_no3 = 0.0;
	flux->n_abdn_cnpy = 0.0;
	flux->n_abdn_strg = 0.0;
	flux->n_salvage = 0.0;
	flux->n_alloc_cnpy = 0.0;
	flux->n_alloc_strg = 0.0;
	flux->n_realloc = 0.0;
}

/* make soil fluxes vacant ******************/
void sflux_zero(
	struct Sflx *flux
){
	/* C */
	flux->li_tf = 0.0;
	flux->li_tc = 0.0;
	flux->li_tr = 0.0;
	flux->li_gf = 0.0;
	flux->li_gc = 0.0;
	flux->li_gr = 0.0;
		
	flux->hf_tfa = 0.0;
	flux->hf_tfi = 0.0;
	flux->hf_tfp = 0.0;
	flux->hf_tca = 0.0;
	flux->hf_tci = 0.0;
	flux->hf_tcp = 0.0;
	flux->hf_tra = 0.0;
	flux->hf_tri = 0.0;
	flux->hf_trp = 0.0;

	flux->hf_gfa = 0.0;
	flux->hf_gfi = 0.0;
	flux->hf_gfp = 0.0;
	flux->hf_gca = 0.0;
	flux->hf_gci = 0.0;
	flux->hf_gcp = 0.0;
	flux->hf_gra = 0.0;
	flux->hf_gri = 0.0;
	flux->hf_grp = 0.0;

	flux->mr_tf = 0.0;
	flux->mr_tc = 0.0;
	flux->mr_tr = 0.0;
	flux->mr_gf = 0.0;
	flux->mr_gc = 0.0;
	flux->mr_gr = 0.0;
	flux->mr_a = 0.0;
	flux->mr_i = 0.0;
	flux->mr_p = 0.0;
	
	flux->hr = 0.0;
	flux->ch4oxy_ridg = 0.0;
	flux->ch4oxy_casa = 0.0;
	flux->ch4oxy_delgrosso = 0.0;
	flux->ch4oxy_curry = 0.0;
	flux->ch4gen = 0.0;
	
	flux->doc_boyer = 0.0;

	/* N */
	flux->d_n2o_ntr_ngas = 0.0;
	flux->d_n2o_dnt_ngas = 0.0;
	flux->d_n2o_ngas = 0.0;
	flux->d_n2_ngas = 0.0;
	flux->d_n2o_ntr_ngas_2 = 0.0;
	flux->d_n2o_dnt_ngas_2 = 0.0;
	flux->d_n2o_ngas_2 = 0.0;
	flux->d_n2_ngas_2 = 0.0;
	flux->d_no_casa = 0.0;
	flux->d_n2_casa = 0.0;
	flux->d_n2o_casa = 0.0;
	flux->n_nh3vlt = 0.0;
	flux->n_leach = 0.0;
	flux->n_minerlz_lttr = 0.0;
	flux->n_minerlz_hums = 0.0;
	flux->n_nitrif = 0.0;
	flux->n_immbl = 0.0;
	flux->n_mcrb_abdn = 0.0;
}

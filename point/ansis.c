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

/***************************************************************************/
void output_ansis_daily(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp
){
	double rr_d, sr_d, voc_d;
	/* double dppfd, dppfdd, dppfdb; */

	rr_d = (flux->tree).rrm + (flux->tree).rrg + loct->funder_c3*(flux->c3).rrm + 
			loct->funder_c3*(flux->c3).rrg + loct->funder_c4*(flux->c4).rrm + 
			loct->funder_c4*(flux->c4).rrg;
	sr_d = rr_d + (flux->soil).hr;
	
	/* total VOC, g C ha-1 day-1 */
	voc_d = (flux->voc_isopr_g97 + flux->voc_monotrp_g97 + flux->voc_methanl_g97 + 
			flux->voc_acetone_g97 + flux->voc_actaldhd_g97 + flux->voc_frmardhd_g97 + 
			flux->voc_formacd_g97 + flux->voc_acetacd_g97 + flux->voc_co_g97);

	/* EDIT here to write results into "daily" file */
	fprintf(fp,"%4ld ", loct->adyear);
	fprintf(fp,"%4ld ", loct->adyear-BYR+1);
	fprintf(fp,"%3ld ", loct->doy);
    
 	fprintf(fp,"%le ", (mass->tree).lai+loct->funder_c3*(mass->c3).lai+loct->funder_c4*(mass->c4).lai);
	fprintf(fp,"%le ", flux->npp *0.1/24.0/3600.0);
	fprintf(fp,"%le ", (flux->soil).hr *0.1/24.0/3600.0);
 	fprintf(fp,"%le ", ((mass->tree).fol+loct->funder_c3*(mass->c3).fol+loct->funder_c4*(mass->c4).fol)*0.1);
 	fprintf(fp,"%le ", ((mass->tree).stm+loct->funder_c3*(mass->c3).stm+loct->funder_c4*(mass->c4).stm)*0.1);
 	fprintf(fp,"%le ", ((mass->tree).rot+loct->funder_c3*(mass->c3).rot+loct->funder_c4*(mass->c4).rot)*0.1);
 	fprintf(fp,"%le ", ((mass->soil).ltr_gf+(mass->soil).ltr_gc+(mass->soil).ltr_gr+
                            (mass->soil).ltr_tf+(mass->soil).ltr_tc+(mass->soil).ltr_tr)*0.1);
  	fprintf(fp,"%le ", ((mass->soil).msl_a+(mass->soil).msl_i+(mass->soil).msl_p)*0.1);
    
    
    /* fprintf(fp,"%9.4lf ", flux->gpp * 100.0);
    fprintf(fp,"%9.4lf ", flux->nep * 100.0);
    fprintf(fp,"%9.4lf ", flux->er * 100.0);
    fprintf(fp,"%9.4lf ", flux->sr * 100.0);
    
    fprintf(fp,"%lf ", loct->tmp_2m);
    fprintf(fp,"%lf ", loct->dswrf_sfc);
    fprintf(fp,"%lf ", loct->prate_sfc);
    fprintf(fp,"%lf ", loct->vpd);
    
    fprintf(fp,"%ld ", (echar->tree).season);
    fprintf(fp,"%lf ", (echar->tree).gdd);
    fprintf(fp,"%lf ", (echar->tree).cdd); */
    
    /* fprintf(fp,"%9.4lf ", (flux->tree).rfm * 100.0);
    fprintf(fp,"%9.4lf ", (flux->tree).rfg * 100.0); */
    
    /* fprintf(fp,"%lf ", flux->nep * 100.0);
    fprintf(fp,"%lf ", loct->rn_eco); */
    
    /* fprintf(fp,"%10.2lf ", (echar->tree).opt_lai);
    fprintf(fp,"%10.2lf ", (mass->tree).lai); */
    
    /* fprintf(fp,"%ld ", (echar->tree).season);
    fprintf(fp,"%ld ", (echar->tree).day_flush);
    fprintf(fp,"%ld ", (echar->tree).day_shed);
    fprintf(fp,"%lf ", (mass->tree).lai);
    
    fprintf(fp,"%ld ", (echar->c3).season);
    fprintf(fp,"%ld ", (echar->c3).day_flush);
    fprintf(fp,"%ld ", (echar->c3).day_shed);
    fprintf(fp,"%lf ", (mass->c3).lai);
    
    fprintf(fp,"%ld ", (echar->c4).season);
    fprintf(fp,"%ld ", (echar->c4).day_flush);
    fprintf(fp,"%ld ", (echar->c4).day_shed);
    fprintf(fp,"%lf ", (mass->c4).lai);
    
    fprintf(fp,"%lf ", flux->gpp * 100.0);
    fprintf(fp,"%lf ", flux->er * 100.0);
    fprintf(fp,"%lf ", flux->nep * 100.0); */
    
    /* fprintf(fp,"%lf ", loct->tmp_2m);
    fprintf(fp,"%ld ", (echar->tree).season);
    fprintf(fp,"%lf ", (echar->tree).gdd);
    fprintf(fp,"%lf ", (echar->tree).cdd); */
    
	/* fprintf(fp,"%7.2lf ", loct->prate_sfc);
	fprintf(fp,"%7.2lf ", loct->tmp_2m); */
    
    /* MSE output: 2013/03/30 by A.Ito */
    /* fprintf(fp,"%10.4lf ", loct->funder_c3*(mass->c3).lai);
    fprintf(fp,"%10.4lf ", flux->gpp*100.0);
	fprintf(fp,"%10.4lf ", flux->nep*100.0);
    fprintf(fp,"%10.4lf ", loct->rooting_depth);
    fprintf(fp,"%10.4lf ", (flux->soil).ch4_wh_diff);
    fprintf(fp,"%10.4lf ", (flux->soil).ch4_wh_ebull);
    fprintf(fp,"%10.4lf ", (flux->soil).ch4_wh_plant);
    fprintf(fp,"%10.4lf ", (flux->soil).ch4_wh_diff + (flux->soil).ch4_wh_ebull + (flux->soil).ch4_wh_plant); */

    /* CEAMIP */
    /* fprintf(fp,"%10.4lf ", (mass->tree).lai+loct->funder_c3*(mass->c3).lai+loct->funder_c4*(mass->c4).lai);
    fprintf(fp,"%9.4lf ", flux->gpp * 100.0);
    fprintf(fp,"%9.4lf ", flux->er * 100.0);
    fprintf(fp,"%9.4lf ", flux->nep * 100.0);
    fprintf(fp,"%9.4lf ", flux->sr * 100.0);
    fprintf(fp,"%9.4lf ", (flux->soil).hr * 100.0);  */

    
    /* fprintf(fp,"%10.2lf ", (echar->tree).gdd);
    fprintf(fp,"%4ld ", (echar->tree).season);
    fprintf(fp,"%10.2lf ", (echar->tree).sla); */
    
    /* fprintf(fp,"%7.2lf ", flux->voc_isopr_g97);
    fprintf(fp,"%7.2lf ", flux->voc_monotrp_g97);
    fprintf(fp,"%7.2lf ", flux->voc_methanl_g97);
    fprintf(fp,"%7.2lf ", flux->voc_acetone_g97);
    fprintf(fp,"%7.2lf ", flux->voc_acetacd_g97);
    fprintf(fp,"%7.2lf ", flux->voc_actaldhd_g97);
    fprintf(fp,"%7.2lf ", flux->voc_formacd_g97);
    fprintf(fp,"%7.2lf ", flux->voc_frmardhd_g97);
    fprintf(fp,"%7.2lf ", flux->voc_co_g97); */
	/* 
	
	fprintf(fp,"%9.4lf ", flux->gpp * 100.0);
	
	fprintf(fp,"%10.4lf ", (mass->tree).lai);
	fprintf(fp,"%10.4lf ", loct->funder_c3*(mass->c3).lai);
	fprintf(fp,"%10.4lf ", loct->funder_c4*(mass->c4).lai);
	
	fprintf(fp,"%10.4lf: ", (mass->tree).nsch_storage);	 */

	/* fprintf(fp,"%10.4lf ", mass->sw30);
	fprintf(fp,"%10.4lf ", loct->pot_grav_l);
	fprintf(fp,"%10.4lf ", loct->pot_matric_l);

	fprintf(fp,"%10.4lf ", mass->sww);
	fprintf(fp,"%10.4lf ", loct->pot_grav_h);
	fprintf(fp,"%10.4lf: ", loct->pot_matric_h); */
	
	/* fprintf(fp,"%4ld ", (echar->c3).season); 
	fprintf(fp,"%6.3lf ", loct->funder_c3*(mass->c3).lai);		
	fprintf(fp,"%10.4lf: ", (mass->c3).nsch_storage*1000.0/10000.0);
	
	fprintf(fp,"%4ld ", (echar->c4).season); 
	fprintf(fp,"%6.3lf ", loct->funder_c4*(mass->c4).lai);
	fprintf(fp,"%10.4lf: ", (mass->c4).nsch_storage*1000.0/10000.0); */
	
    /* 
    fprintf(fp,"%7.2lf ", loct->bCO2);
    fprintf(fp,"%7.2lf ", loct->vps);
    fprintf(fp,"%7.2lf ", loct->vp);
    fprintf(fp,"%7.2lf ", loct->vpd); */
    
	/* ADMIP output *******/
	/* fprintf(fp,"%le ", mass->sw30);
	fprintf(fp,"%le ", mass->sww);
	fprintf(fp,"%le ", loct->tmp10_soil+ZAT);
	fprintf(fp,"%le ", loct->tmp200_soil+ZAT);
	fprintf(fp,"%le ", loct->aet/86400.0);
	fprintf(fp,"%le ", loct->aet*2500000.0/86400.0);
	fprintf(fp,"%le ", -flux->gpp * 1000.0/10000.0/86400.0);
	fprintf(fp,"%le ", -flux->nep * 1000.0/10000.0/86400.0); 
	fprintf(fp,"%le ", loct->funder_c3*(mass->c3).lai);
	fprintf(fp,"%le ", loct->funder_c4*(mass->c4).lai);
	fprintf(fp,"%le ", loct->funder_c3*((mass->c3).fol+(mass->c3).stm+(mass->c3).rot)*1000.0/10000.0);
	fprintf(fp,"%le ", loct->funder_c4*((mass->c4).fol+(mass->c4).stm+(mass->c4).rot)*1000.0/10000.0);
	fprintf(fp,"%le ", loct->funder_c3*(flux->c3).rp * 1000.0/10000.0/86400.0);
	fprintf(fp,"%le ", loct->funder_c4*(flux->c4).rp * 1000.0/10000.0/86400.0);
	fprintf(fp,"%le ", (flux->soil).hr * 1000.0/10000.0/86400.0);
	fprintf(fp,"%le ", loct->albedo_sfc);
    fprintf(fp,"%le ", loct->snow_acc); */
    
    //fprintf(fp,"%4ld ", (echar->tree).season); 
    //fprintf(fp,"%4ld ", (echar->c3).season); 
    //fprintf(fp,"%4ld ", (echar->c4).season); 
	
    //fprintf(fp,"%10.4lf ", (echar->c3).gdd); 
    //fprintf(fp,"%10.4lf ", (echar->c3).cdd); 

    /*************************/
	
	/* for KBY preliminary analysis
	fprintf(fp,"%7.2lf ", loct->funder_c3*(mass->c3).lai);
	fprintf(fp,"%7.2lf ", loct->funder_c4*(mass->c4).lai);
	fprintf(fp,"%9.4lf ", flux->gpp * 100.0);
	fprintf(fp,"%9.4lf ", flux->er * 100.0);
	fprintf(fp,"%9.4lf ", flux->nep * 100.0);
	
	fprintf(fp,"%9.6lf ", (echar->c3).opt_lai);
	
	fprintf(fp,"%7.2lf ", loct->rn_long_eco);
	fprintf(fp,"%7.2lf ", loct->rn_short_eco);
	fprintf(fp,"%7.2lf ", loct->aet); */
	
	/*  (SR comparison at IARC)
	fprintf(fp,"%.4lf ", flux->sr*44.0/12.0*100.0*1000.0/86400.0);
	fprintf(fp,"%.4lf ", flux->rr*44.0/12.0*100.0*1000.0/86400.0);
	fprintf(fp,"%.4lf ", (flux->soil).hr*44.0/12.0*100.0*1000.0/86400.0);
	*/
	
	/*  
	fprintf(fp,"%7.2lf ", loct->dswrf_sfc);  */
	
	/*	fprintf(fp,"%7.2lf ", loct->bCO2);
	fprintf(fp,"%7.2lf ", (echar->tree).psat);
	fprintf(fp,"%7.2lf ", ((flux->tree).gpp+(flux->c3).gpp)*100.0);
	
	fprintf(fp,"%7.2lf ", (mass->tree).lai+(mass->c3).lai); 
	fprintf(fp,"%7.2lf ", ((mass->tree).fol+(mass->tree).stm+(mass->tree).rot+(mass->c3).fol+(mass->c3).stm+(mass->c3).rot)/10.0);
	fprintf(fp,"%7.2lf ", ((mass->soil).ltr_tf+(mass->soil).ltr_tc+(mass->soil).ltr_tr+(mass->soil).ltr_gf+(mass->soil).ltr_gc+
						   (mass->soil).ltr_gr+(mass->soil).msl_a+(mass->soil).msl_i+(mass->soil).msl_p)/10.0); */
	
/*	
	fprintf(fp,"%9.4lf ", flux->gpp * 0.1);
	fprintf(fp,"%7.2lf ", loct->dswrf_sfc);
	fprintf(fp,"%9.4lf ", (echar->c3).psat); 
	fprintf(fp,"%3ld ", loct->month+1);
	fprintf(fp,"%7.2lf ", (mass->c3).lai); 
	fprintf(fp,"%7.2lf ", (mass->c4).lai);
	fprintf(fp,"%9.4lf ", (flux->tree).gpp * 100.0);
	fprintf(fp,"%9.4lf ", (flux->c3).gpp * 100.0);
	fprintf(fp,"%9.4lf ", (flux->c4).gpp * 100.0);
	fprintf(fp,"%8.3lf ", (flux->soil).ch4oxy_curry);  */

/*	fprintf(fp,"%9.6lf ", (echar->tree).n_conc_larea);
	fprintf(fp,"%9.6lf ", (echar->c3).n_conc_larea);
	fprintf(fp,"%9.6lf ", (echar->c4).n_conc_larea);  */
	
/*	fprintf(fp,"%9.6lf ", loct->prate_sfc);
	fprintf(fp,"%7.2lf ", loct->ro2); 
	fprintf(fp,"%9.6lf ", (flux->soil).doc_boyer); 
	fprintf(fp,"%9.6lf ", (flux->soil).ch4oxy_ridg);
	fprintf(fp,"%8.3lf ", (flux->soil).ch4oxy_curry);
	fprintf(fp,"%8.3lf ", (flux->soil).ch4oxy_casa);
	fprintf(fp,"%8.3lf ", (flux->soil).ch4oxy_delgrosso);
	fprintf(fp,"%9.6lf ", (flux->tree).emit_ch4_kirschbaum_mass+loct->funder_c3*(flux->c3).emit_ch4_kirschbaum_mass+
							loct->funder_c4*(flux->c4).emit_ch4_kirschbaum_mass);
	fprintf(fp,"%9.6lf ", (flux->tree).emit_ch4_kirschbaum_photo+loct->funder_c3*(flux->c3).emit_ch4_kirschbaum_photo+
							loct->funder_c4*(flux->c4).emit_ch4_kirschbaum_photo);
	
	fprintf(fp,"%8.3lf ", (flux->soil).d_n2o_ntr_ngas);
	fprintf(fp,"%8.3lf ", (flux->soil).d_n2o_dnt_ngas);
	fprintf(fp,"%8.3lf ", (flux->soil).d_n2o_casa);
	fprintf(fp,"%7.3lf ", loct->wfps);
	fprintf(fp,"%7.3lf ", loct->e_casa);
	fprintf(fp,"%7.3lf ", loct->iw_casa);
	fprintf(fp,"%7.3lf ", loct->m_casa);  
	fprintf(fp,"%7.2lf ", loct->tmp10_soil);
	fprintf(fp,"%7.2lf ", loct->tmp200_soil);
	fprintf(fp,"%7.2lf ", loct->vpd);
	fprintf(fp,"%7.2lf ", loct->bCO2);
	fprintf(fp,"%.3lf ", flux->er*100.0);
	fprintf(fp,"%7.2lf ", mass->sw30);
	fprintf(fp,"%7.2lf ", mass->sww);
	fprintf(fp,"%9.6lf ", (flux->soil).d_n2o_ngas);
	fprintf(fp,"%9.6lf ", (flux->soil).doc_boyer);
	fprintf(fp,"%9.6lf ", voc_d);
	fprintf(fp,"%7.2lf ", (echar->tree).n_opt_leaf);
	fprintf(fp,"%7.2lf ", (echar->c3).n_opt_leaf);
	fprintf(fp,"%7.2lf ", (echar->c4).n_opt_leaf);
	
	fprintf(fp,"%7.2lf ", (echar->tree).n_conc_larea);
	fprintf(fp,"%7.2lf ", (echar->c3).n_conc_larea);
	fprintf(fp,"%7.2lf ", (echar->c4).n_conc_larea);  
	
	fprintf(fp,"%7.2lf ", (mass->tree).n_cnpy);
	fprintf(fp,"%7.2lf ", (mass->c3).n_cnpy);
	fprintf(fp,"%7.2lf ", (mass->c4).n_cnpy);  

	fprintf(fp,"%7.2lf ", loct->depo_no3);
	fprintf(fp,"%7.2lf ", loct->depo_nh4); */

	/************************************************/

/*	
	fprintf(fp,"%8.3lf ", (flux->soil).d_n2o_ngas_2);
	fprintf(fp,"%7.2lf ", loct->ro2);
	fprintf(fp,"%8.6lf ", (mass->soil).doc);  */

/*	fprintf(fp,"%.3lf ",flux->gpp*-100.0);	
	fprintf(fp,"%7.2lf ", loct->prate_sfc);
	fprintf(fp,"%7.2lf ", mass->sw30);
	fprintf(fp,"%7.2lf ", mass->snwa);
	fprintf(fp,"%7.2lf ", loct->rdr_casa);
	fprintf(fp,"%7.2lf ", loct->pet);
	fprintf(fp,"%7.2lf ", loct->vmc); */
	
/*	fprintf(fp,"%7.2lf ", loct->x1);
	fprintf(fp,"%7.2lf ", loct->x2);
	fprintf(fp,"%7.2lf ", loct->x3);
	fprintf(fp,"%7.2lf ", loct->x4);
	fprintf(fp,"%7.2lf ", loct->x5); */

/*	fprintf(fp,"%7.2lf ", (mass->soil).n_nh4);
	fprintf(fp,"%7.2lf ", (mass->soil).n_no3);
	fprintf(fp,"%7.2lf ", (mass->soil).n_lttr);
	fprintf(fp,"%7.2lf ", (mass->soil).n_hums);
	fprintf(fp,"%7.2lf ", (mass->soil).n_mcrb);
	
	fprintf(fp,"%7.2lf ", (mass->tree).n_cnpy);
	fprintf(fp,"%7.2lf ", (mass->tree).n_strg);
	fprintf(fp,"%7.2lf ", (mass->c3).n_cnpy);
	fprintf(fp,"%7.2lf ", (mass->c3).n_strg);
	fprintf(fp,"%7.2lf ", (mass->c4).n_cnpy);
	fprintf(fp,"%7.2lf ", (mass->c4).n_strg);
	
	fprintf(fp,"%7.4lf ", (flux->soil).n_nh3vlt*1000.0);
	fprintf(fp,"%7.4lf ", (flux->soil).n_leach);
	fprintf(fp,"%7.4lf ", (flux->soil).n_minerlz_lttr);
	fprintf(fp,"%7.4lf ", (flux->soil).n_minerlz_hums);
	fprintf(fp,"%7.4lf ", (flux->soil).n_nitrif);
	fprintf(fp,"%7.4lf ", (flux->soil).n_immbl);
	fprintf(fp,"%7.4lf ", (flux->soil).n_mcrb_abdn);
	
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2o_ntr_ngas);
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2o_dnt_ngas);
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2o_ngas);
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2_ngas);
	
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2o_ntr_ngas_2);
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2o_dnt_ngas_2);
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2o_ngas_2);
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2_ngas_2);
	
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2o_casa);
	fprintf(fp,"%7.4lf ", (flux->soil).d_no_casa);
	fprintf(fp,"%7.4lf ", (flux->soil).d_n2_casa);
	
	fprintf(fp,"%7.4lf ", (flux->tree).uptake_no3+(flux->c3).uptake_no3+(flux->c4).uptake_no3);
	fprintf(fp,"%7.4lf ", (flux->tree).uptake_nh4+(flux->c3).uptake_nh4+(flux->c4).uptake_nh4); */

	/*
	fprintf(fp,"%7.2lf ", (flux->tree).n_alloc_cnpy);
	fprintf(fp,"%7.2lf ", (flux->tree).n_realloc);
	fprintf(fp,"%7.2lf ", (flux->c3).n_alloc_cnpy);
	fprintf(fp,"%7.2lf ", (flux->c3).n_realloc);
	fprintf(fp,"%7.2lf ", (flux->c4).n_alloc_cnpy);
	fprintf(fp,"%7.2lf ", (flux->c4).n_realloc);
	*/
	/*
	fprintf(fp,"%7.2lf ", (flux->tree).uptake_no3);
	fprintf(fp,"%7.2lf ", (flux->tree).uptake_nh4);
	fprintf(fp,"%7.2lf ", (flux->c3).uptake_no3);
	fprintf(fp,"%7.2lf ", (flux->c3).uptake_nh4);
	fprintf(fp,"%7.2lf ", (flux->c4).uptake_no3);
	fprintf(fp,"%7.2lf ", (flux->c4).uptake_nh4);
	*/
	/*
	fprintf(fp,"%7.4lf ", (flux->soil).n_nh3vlt);
	fprintf(fp,"%7.4lf ", (flux->soil).n_leach);
	fprintf(fp,"%7.4lf ", (flux->soil).n_minerlz_lttr);
	fprintf(fp,"%7.4lf ", (flux->soil).n_minerlz_hums);
	fprintf(fp,"%7.4lf ", (flux->soil).n_nitrif);
	fprintf(fp,"%7.4lf ", (flux->soil).n_immbl);
	fprintf(fp,"%7.4lf ", (flux->soil).n_mcrb_abdn);
	*/
	/*
	fprintf(fp,"%7.2lf ", (flux->tree).n_biofix);
	fprintf(fp,"%7.2lf ", (flux->c3).n_biofix);
	fprintf(fp,"%7.2lf ", (flux->c4).n_biofix);
	*/

/*	fprintf(fp,"%.2lf ", loct->tmp10_soil);
	fprintf(fp,"%.3lf ",flux->rr*100.0);
	fprintf(fp,"%.3lf ",(flux->soil).hr*100.0);	
	fprintf(fp,"%.3lf ",flux->sr*100.0);	
		
	fprintf(fp,"%7.2lf ", loct->prate_sfc);
	fprintf(fp,"%7.2lf ", loct->incep);
	fprintf(fp,"%7.2lf ", loct->trnsp);
	fprintf(fp,"%7.2lf ", loct->evpr);
	fprintf(fp,"%7.2lf ", loct->ro1);
	fprintf(fp,"%7.2lf ", loct->ro2); 
	 */

/*	fprintf(fp,"%7.2lf ", loct->tmp_2m); 
	fprintf(fp,"%7.2lf ", (echar->c3).gdd); 
	fprintf(fp,"%7.2lf ", (echar->c3).cdd); 
	fprintf(fp,"%7.2lf ", (mass->c3).lai); 
	fprintf(fp,"%7.2lf ", (echar->c3).opt_lai); */
	
	/* for PAR analysis */
/*	dppfd = dppfdd = dppfdb = 0.0;
	for(f=0;f<48;f++){
		dppfd += loct->ppfd_h[f]*1800.0/1000000.0;
		dppfdd += loct->ppfdd_h[f]*1800.0/1000000.0;
		dppfdb += loct->ppfdb_h[f]*1800.0/1000000.0;
	}
	
	fprintf(fp,"%lf ",loct->tcdc_clm);
	fprintf(fp,"%lf ",loct->daylen[loct->doy]);
	fprintf(fp,"%lf ",loct->par_h[24]);
	fprintf(fp,"%lf ",loct->ppfd_h[24]);
	fprintf(fp,"%lf ",loct->csppfdb_h[24]+loct->csppfdd_h[24]);
	
	fprintf(fp,"%lf ",dppfd);
	fprintf(fp,"%lf ",dppfdd);
	fprintf(fp,"%lf ",dppfdb); */
	
	/* fprintf(fp,"%.3lf ",(mass->tree).lai);		
	fprintf(fp,"%.3lf ", loct->funder_c3*(mass->c3).lai);		
	fprintf(fp,"%.3lf ", loct->funder_c4*(mass->c4).lai);		*/
	
/* 	fprintf(fp,"%.2lf ", loct->ppfd_h[24]); */
/* 	fprintf(fp,"%.2lf ", loct->daylen[loct->doy]); */
/*	fprintf(fp,"%.2lf ", loct->ppfd_h[24]*(3600.0*loct->daylen[loct->doy])/2.0*(1.0-loct->tcdc_clm)/1000000.0); */
	
/*	fprintf(fp,"%.2lf ", (flux->tree).gpp*100.0 + loct->funder_c3*(flux->c3).gpp*100.0 + 
			loct->funder_c4*(flux->c4).gpp*100.0);
	fprintf(fp,"%.2lf ", (flux->tree).npp*100.0 + loct->funder_c3*(flux->c3).npp*100.0 + 
			loct->funder_c4*(flux->c4).npp*100.0);
	fprintf(fp,"%.2lf ", flux->nep*100.0);
	fprintf(fp,"%.2lf ", sr_d*100.0);
	fprintf(fp,"%.2lf %.2lf %.2lf ", (mass->tree).lai, loct->funder_c3*(mass->c3).lai, 
			loct->funder_c4*(mass->c4).lai); */
			
/*	fprintf(fp,"%.3lf ",(mass->tree).lai);		
	fprintf(fp,"%.3lf ", loct->funder_c3*(mass->c3).lai);		
	fprintf(fp,"%.3lf ", loct->funder_c4*(mass->c4).lai);		
	fprintf(fp,"%ld %.2lf %.2lf ", (echar->c3).season, (echar->c3).cdd);
	fprintf(fp,"	%.4lf ", (flux->tree).gpp*100.0);
	fprintf(fp,"%.4lf ", loct->funder_c3*(flux->c3).gpp*100.0);
	fprintf(fp,"%.4lf ", loct->funder_c4*(flux->c4).gpp*100.0); */

/*			
	fprintf(fp,"%.2lf ", (echar->c3).psat);
	fprintf(fp,"%.2lf ", (echar->c3).opt_lai); */
	
	fprintf(fp,"\n");
}

/*****************************************************************************************************/
void f_ansis_ann(
	struct Grid	*grid, 
	struct Loct	*loct, 
	struct Echar	*echar, 
	struct Mass	*mass, 
	struct Flux	*flux, 
	double ansis_ann[N_ANSIS]
){
	double rr_d, sr_d, er_d, fscale, voc_d, ynday;
	
	fscale = 1.0;
	ynday = 365.0;
	if(loct->climy%4==0){
		ynday = 366.0;
	}
	
	rr_d = (flux->tree).rrm + (flux->tree).rrg + (flux->c3).rrm + (flux->c3).rrg;
	sr_d = rr_d + (flux->soil).hr;
	er_d = (flux->tree).rp + (flux->c3).rp + (flux->soil).hr;
	/* total VOC, g C ha-1 day-1 */
	voc_d = (flux->voc_isopr_g97 + flux->voc_monotrp_g97 + flux->voc_methanl_g97 + 
			flux->voc_acetone_g97 + flux->voc_actaldhd_g97 + flux->voc_frmardhd_g97 + 
			flux->voc_formacd_g97 + flux->voc_acetacd_g97 + flux->voc_co_g97);

	ansis_ann[0] += ((flux->tree).gpp+loct->funder_c3*(flux->c3).gpp+loct->funder_c4*(flux->c4).gpp)*fscale;
	ansis_ann[1] += ((flux->tree).npp+loct->funder_c3*(flux->c3).npp+loct->funder_c4*(flux->c4).npp)*fscale;
	ansis_ann[2] += flux->nep*fscale;
	ansis_ann[3] += sr_d*fscale;
	ansis_ann[4] += er_d*fscale;
	
	ansis_ann[5] += (mass->tree).lai/365.0;
	ansis_ann[6] += (mass->tree).fol/365.0*fscale;
	ansis_ann[7] += (mass->tree).stm/365.0*fscale;
	ansis_ann[8] += (mass->tree).rot/365.0*fscale;
	ansis_ann[9] += loct->funder_c3*(mass->c3).lai/365.0;
	ansis_ann[10] += loct->funder_c3*(mass->c3).fol/365.0*fscale;
	ansis_ann[11] += loct->funder_c3*(mass->c3).stm/365.0*fscale;
	ansis_ann[12] += loct->funder_c3*(mass->c3).rot/365.0*fscale;
	ansis_ann[13] += loct->funder_c4*(mass->c4).lai/365.0;
	ansis_ann[14] += loct->funder_c4*(mass->c4).fol/365.0*fscale;
	ansis_ann[15] += loct->funder_c4*(mass->c4).stm/365.0*fscale;
	ansis_ann[16] += loct->funder_c4*(mass->c4).rot/365.0*fscale;
	
	ansis_ann[17] += (mass->soil).ltr_tf/365.0*fscale;
	ansis_ann[18] += (mass->soil).ltr_tc/365.0*fscale;
	ansis_ann[19] += (mass->soil).ltr_tr/365.0*fscale;
	ansis_ann[20] += (mass->soil).ltr_gf/365.0*fscale;
	ansis_ann[21] += (mass->soil).ltr_gc/365.0*fscale;
	ansis_ann[22] += (mass->soil).ltr_gr/365.0*fscale;
	ansis_ann[23] += (mass->soil).msl_a/365.0*fscale;
	ansis_ann[24] += (mass->soil).msl_i/365.0*fscale;
	ansis_ann[25] += (mass->soil).msl_p/365.0*fscale;
	
	ansis_ann[30] += (mass->tree).n_canopy/365.0;
	ansis_ann[31] += (mass->tree).n_storage/365.0;
	ansis_ann[32] += loct->funder_c3*(mass->c3).n_canopy/365.0;
	ansis_ann[33] += loct->funder_c3*(mass->c3).n_storage/365.0;
	ansis_ann[34] += loct->funder_c4*(mass->c4).n_canopy/365.0;
	ansis_ann[35] += loct->funder_c4*(mass->c4).n_storage/365.0;
	
	ansis_ann[36] += (mass->soil).n_nh4/365.0;
	ansis_ann[37] += (mass->soil).n_no3/365.0;
	ansis_ann[38] += (mass->soil).n_lttr/365.0;
	ansis_ann[39] += (mass->soil).n_hums/365.0;
	ansis_ann[40] += (mass->soil).n_mcrb/365.0;
	
	ansis_ann[48] += loct->depo_no3;
	ansis_ann[49] += loct->depo_nh4;

	ansis_ann[50] += (flux->tree).n_biofix;
	ansis_ann[51] += (flux->tree).uptake_no3;
	ansis_ann[52] += (flux->tree).uptake_nh4;
	ansis_ann[53] += (flux->tree).n_abdn_cnpy;
	ansis_ann[54] += (flux->tree).n_abdn_strg;
	ansis_ann[55] += (flux->tree).n_salvage;
	ansis_ann[56] += (flux->tree).n_alloc_cnpy;
	ansis_ann[57] += (flux->tree).n_alloc_strg;
	ansis_ann[58] += (flux->tree).n_realloc;
	
	ansis_ann[59] += loct->funder_c3*(flux->c3).n_biofix;
	ansis_ann[60] += loct->funder_c3*(flux->c3).uptake_no3;
	ansis_ann[61] += loct->funder_c3*(flux->c3).uptake_nh4;
	ansis_ann[62] += loct->funder_c3*(flux->c3).n_abdn_cnpy;
	ansis_ann[63] += loct->funder_c3*(flux->c3).n_abdn_strg;
	ansis_ann[64] += loct->funder_c3*(flux->c3).n_salvage;
	ansis_ann[65] += loct->funder_c3*(flux->c3).n_alloc_cnpy;
	ansis_ann[66] += loct->funder_c3*(flux->c3).n_alloc_strg;
	ansis_ann[67] += loct->funder_c3*(flux->c3).n_realloc;

	ansis_ann[68] += loct->funder_c4*(flux->c4).n_biofix;
	ansis_ann[69] += loct->funder_c4*(flux->c4).uptake_no3;
	ansis_ann[70] += loct->funder_c4*(flux->c4).uptake_nh4;
	ansis_ann[71] += loct->funder_c4*(flux->c4).n_abdn_cnpy;
	ansis_ann[72] += loct->funder_c4*(flux->c4).n_abdn_strg;
	ansis_ann[73] += loct->funder_c4*(flux->c4).n_salvage;
	ansis_ann[74] += loct->funder_c4*(flux->c4).n_alloc_cnpy;
	ansis_ann[75] += loct->funder_c4*(flux->c4).n_alloc_strg;
	ansis_ann[76] += loct->funder_c4*(flux->c4).n_realloc;

	ansis_ann[78] += (flux->soil).d_n2o_ntr_ngas;
	ansis_ann[79] += (flux->soil).d_n2o_dnt_ngas;
	ansis_ann[80] += (flux->soil).d_n2o_ngas;
	ansis_ann[81] += (flux->soil).d_n2_ngas;
	ansis_ann[82] += (flux->soil).d_n2o_ntr_ngas_2;
	ansis_ann[83] += (flux->soil).d_n2o_dnt_ngas_2;
	ansis_ann[84] += (flux->soil).d_n2o_ngas_2;
	ansis_ann[85] += (flux->soil).d_n2_ngas_2;
	ansis_ann[86] += (flux->soil).d_no_casa;
	ansis_ann[87] += (flux->soil).d_n2_casa;
	ansis_ann[88] += (flux->soil).d_n2o_casa;
	ansis_ann[89] += (flux->soil).n_nh3vlt;
	ansis_ann[90] += (flux->soil).n_leach;
	ansis_ann[91] += (flux->soil).n_minerlz_lttr;
	ansis_ann[92] += (flux->soil).n_minerlz_hums;
	ansis_ann[93] += (flux->soil).n_nitrif;
	ansis_ann[94] += (flux->soil).n_immbl;
	ansis_ann[95] += (flux->soil).n_mcrb_abdn;

	ansis_ann[96] += loct->depo_no3 + loct->depo_nh4;
	ansis_ann[97] += (flux->soil).doc_boyer;

	ansis_ann[98] += (flux->soil).ch4oxy_ridg;
	
	ansis_ann[99] += (flux->soil).doc_boyer;
	ansis_ann[100] += voc_d;

	ansis_ann[101] += (flux->tree).emit_ch4_kirschbaum_mass+loct->funder_c3*(flux->c3).emit_ch4_kirschbaum_mass+
							loct->funder_c4*(flux->c4).emit_ch4_kirschbaum_mass;
	ansis_ann[102] += (flux->tree).emit_ch4_kirschbaum_photo+loct->funder_c3*(flux->c3).emit_ch4_kirschbaum_photo+
							loct->funder_c4*(flux->c4).emit_ch4_kirschbaum_photo;
							
	ansis_ann[103] += (flux->soil).ch4oxy_curry;
	ansis_ann[104] += (flux->soil).ch4oxy_casa;
	ansis_ann[105] += (flux->soil).ch4oxy_delgrosso;
	
	ansis_ann[106] += (flux->tree).emit_ch4_kirschbaum_mass+loct->funder_c3*(flux->c3).emit_ch4_kirschbaum_mass+
							loct->funder_c4*(flux->c4).emit_ch4_kirschbaum_mass;
	ansis_ann[107] += (flux->tree).emit_ch4_kirschbaum_photo+loct->funder_c3*(flux->c3).emit_ch4_kirschbaum_photo+
							loct->funder_c4*(flux->c4).emit_ch4_kirschbaum_photo;
							
	ansis_ann[108] += (flux->soil).ch4_wh_plant;
	ansis_ann[109] += (flux->soil).ch4_wh_ebull;
	ansis_ann[110] += (flux->soil).ch4_wh_diff;
	ansis_ann[111] = ansis_ann[108] + ansis_ann[109] + ansis_ann[110];
	
	ansis_ann[112] += loct->prate_sfc;
	if(loct->tmp_2m > ansis_ann[113]){
		ansis_ann[113] = loct->tmp_2m;
	}
	if(loct->tmp_2m < ansis_ann[114]){
		ansis_ann[114] = loct->tmp_2m;
	}
	ansis_ann[115] += loct->tmp_2m / ynday;
	
	ansis_ann[116] += loct->dswrf_sfc * 24.0 * 3600.0 / 1000000.0;

	ansis_ann[117] += (flux->soil).erosion_carbon/365.0;
    
    ansis_ann[118] += loct->aet;
    ansis_ann[119] += loct->ro2;

    ansis_ann[120] += flux->voc_isopr_g97;
    ansis_ann[121] += flux->voc_monotrp_g97;
    ansis_ann[122] += (flux->soil).hr;
}

/*********************************************************************/
void output_ansis_ann(
	long year, 
	double ansis_ann[N_ANSIS], 
	FILE *fp
){
	double c_litter, c_humus;
	
	c_litter = ansis_ann[17]+ansis_ann[18]+ansis_ann[19]+ansis_ann[20]+ansis_ann[21]+ansis_ann[22];
	c_humus = ansis_ann[23]+ansis_ann[24]+ansis_ann[25];
	
	fprintf(fp, "%ld ", year);
	fprintf(fp, "%ld ", year-BYR+1);
	
	/* fprintf(fp, "%lf,", ansis_ann[6]+ansis_ann[7]+ansis_ann[8] + ansis_ann[10]+ansis_ann[11]+ansis_ann[12]);
	fprintf(fp, "%lf,", c_litter + c_humus); */

    fprintf(fp, "%le ", ansis_ann[5]+ansis_ann[9]+ansis_ann[13]);

	fprintf(fp, "%le ", ansis_ann[1]*0.1/365.0/24.0/3600.0);
	fprintf(fp, "%le ", ansis_ann[122]*0.1/365.0/24.0/3600.0);
    
	fprintf(fp, "%le ", (ansis_ann[6]+ansis_ann[10]+ansis_ann[14])*0.1);
	fprintf(fp, "%le ", (ansis_ann[7]+ansis_ann[11]+ansis_ann[15])*0.1);
	fprintf(fp, "%le ", (ansis_ann[8]+ansis_ann[12]+ansis_ann[16])*0.1);
	fprintf(fp, "%le ", c_litter*0.1);
	fprintf(fp, "%le ", c_humus*0.1);
    
	
	/* fprintf(fp, "%lf,", c_litter);
	fprintf(fp, "%lf,", c_humus); */
	
	/* fprintf(fp, "%lf,", ansis_ann[23]);
	fprintf(fp, "%lf,", ansis_ann[24]);
	fprintf(fp, "%lf,", ansis_ann[25]);
	fprintf(fp, "%lf,", ansis_ann[4]); 
	fprintf(fp, "%lf,", ansis_ann[5]+ansis_ann[9]);
	fprintf(fp, "%lf,", ansis_ann[6]+ansis_ann[10]);
    fprintf(fp, "%lf ", ansis_ann[117]);*/
	

	/* fprintf(fp, "%ld ", year);
	
	fprintf(fp, "%lf ", ansis_ann[112]);
	fprintf(fp, "%lf ", ansis_ann[113]);
	fprintf(fp, "%lf ", ansis_ann[114]);
	fprintf(fp, "%lf ", ansis_ann[115]);
	fprintf(fp, "%lf ", ansis_ann[116]);
	
	fprintf(fp, "%lf ", ansis_ann[0]);
	fprintf(fp, "%lf ", ansis_ann[1]);
	fprintf(fp, "%lf ", ansis_ann[2]);
	fprintf(fp, "%lf ", ansis_ann[3]);
	fprintf(fp, "%lf ", ansis_ann[4]);
    
    fprintf(fp, "%lf ", ansis_ann[2]);
    fprintf(fp, "%lf ", ansis_ann[103]);
    fprintf(fp, "%lf ", ansis_ann[80]);
    
    fprintf(fp, "%lf ", ansis_ann[2]*100.0*44.0/12.0 
                        + ansis_ann[103]/1000.0*25.0
                        - ansis_ann[80]/1000.0*298.0);

    fprintf(fp, "%lf ", ansis_ann[117]);
	
    fprintf(fp, "%lf ", ansis_ann[118]);
    fprintf(fp, "%lf ", ansis_ann[119]);
    fprintf(fp, "%lf ", ansis_ann[120]);
    fprintf(fp, "%lf ", ansis_ann[121]); */
		
	
	/*
	fprintf(fp, "%lf ", ansis_ann[23]);
	fprintf(fp, "%lf ", ansis_ann[24]);
	fprintf(fp, "%lf ", ansis_ann[25]); */
	
	
	
	/*	fprintf(fp, "%10.3lf ", ansis_ann[3]);
	fprintf(fp, "%10.3lf ", ansis_ann[4]);
	
	fprintf(fp, "%10.3lf ", ansis_ann[5]);
	fprintf(fp, "%10.3lf ", ansis_ann[9]);
	fprintf(fp, "%10.3lf ", ansis_ann[13]);

	fprintf(fp, "%10.3lf ", ansis_ann[6]);
	fprintf(fp, "%10.3lf ", ansis_ann[7]);
	fprintf(fp, "%10.3lf ", ansis_ann[8]);

	fprintf(fp, "%10.3lf ", ansis_ann[10]);
	fprintf(fp, "%10.3lf ", ansis_ann[11]);
	fprintf(fp, "%10.3lf ", ansis_ann[12]);

	fprintf(fp, "%10.3lf ", ansis_ann[14]);
	fprintf(fp, "%10.3lf ", ansis_ann[15]);
	fprintf(fp, "%10.3lf ", ansis_ann[16]);

	fprintf(fp, "%10.3lf ", ansis_ann[17]);
	fprintf(fp, "%10.3lf ", ansis_ann[18]);
	fprintf(fp, "%10.3lf ", ansis_ann[19]);
	fprintf(fp, "%10.3lf ", ansis_ann[20]);
	fprintf(fp, "%10.3lf ", ansis_ann[21]);
	fprintf(fp, "%10.3lf ", ansis_ann[22]);

	
	fprintf(fp, "%10.3lf ", ansis_ann[112]);  */

	
	//	fprintf(fp, "%12.3lf ", ansis_ann[30]); /* tree leaf N */
//	fprintf(fp, "%12.3lf ", ansis_ann[31]); /* tree storage N */
//	fprintf(fp, "%12.3lf ", ansis_ann[32]); /* c3grass leaf N */
//	fprintf(fp, "%12.3lf ", ansis_ann[33]); /* c3grass storage N */
//	fprintf(fp, "%12.3lf ", ansis_ann[34]); /* c4grass leaf N */
//	fprintf(fp, "%12.3lf ", ansis_ann[35]); /* c4grass storage N */
	
//	fprintf(fp, "%12.3lf ", ansis_ann[36]); /* soil NH4+ */
//	fprintf(fp, "%12.3lf ", ansis_ann[37]); /* soil NO3- */
//	fprintf(fp, "%12.3lf ", ansis_ann[38]); /* litter N */
//	fprintf(fp, "%12.3lf ", ansis_ann[39]); /* humus N */
//	fprintf(fp, "%12.3lf ", ansis_ann[40]); /* microbial N */
	
//	fprintf(fp, "%12.3lf ", ansis_ann[6]*1000000.0/ansis_ann[30]); 
//	fprintf(fp, "%12.3lf ", (ansis_ann[7]+ansis_ann[8])*1000000.0/ansis_ann[31]); 
//	fprintf(fp, "%12.3lf ", c_litter*1000000.0/ansis_ann[38]); 
//	fprintf(fp, "%12.3lf  ", c_humus*1000000.0/ansis_ann[39]);
		
	
//	fprintf(fp, "%12.3lf ", ansis_ann[49]); /* NH4 deposition */
//	fprintf(fp, "%12.3lf ", ansis_ann[48]); /* NO3 deposition */
//	fprintf(fp, "%12.3lf ", ansis_ann[50]+ansis_ann[59]+ansis_ann[68]); /* biological N2 fixation */
	
//	fprintf(fp, "%12.3lf ", ansis_ann[90]); /* NO3- leaching */
//	fprintf(fp, "%12.3lf ", ansis_ann[78]); /* nitrif N2O emission */
//	fprintf(fp, "%12.3lf ", ansis_ann[79]); /* denitrif N2O emission */
//	fprintf(fp, "%12.3lf ", ansis_ann[81]); /* denitrif N2 emission */
//	fprintf(fp, "%12.3lf ", ansis_ann[89]); /* NH3 volatilization */

//	fprintf(fp, "%12.3lf ", ansis_ann[52]); /* plant NH4+ uptake */
//	fprintf(fp, "%12.3lf ", ansis_ann[61]); /* plant NH4+ uptake */
//	fprintf(fp, "%12.3lf ", ansis_ann[70]); /* plant NH4+ uptake */
//	fprintf(fp, "%12.3lf ", ansis_ann[51]); /* plant NO3- uptake */
//	fprintf(fp, "%12.3lf ", ansis_ann[60]); /* plant NO3- uptake */
//	fprintf(fp, "%12.3lf ", ansis_ann[69]); /* plant NO3- uptake */

//	fprintf(fp, "%12.3lf ", ansis_ann[56]);
//	fprintf(fp, "%12.3lf ", ansis_ann[57]);
//	fprintf(fp, "%12.3lf ", ansis_ann[65]);
//	fprintf(fp, "%12.3lf ", ansis_ann[66]);
//	fprintf(fp, "%12.3lf ", ansis_ann[74]);
//	fprintf(fp, "%12.3lf ", ansis_ann[75]);
	
//	fprintf(fp, "%12.3lf ", ansis_ann[53]);
//	fprintf(fp, "%12.3lf ", ansis_ann[54]);
//	fprintf(fp, "%12.3lf ", ansis_ann[62]);
//	fprintf(fp, "%12.3lf ", ansis_ann[63]);
//	fprintf(fp, "%12.3lf ", ansis_ann[71]);
//	fprintf(fp, "%12.3lf ", ansis_ann[72]);
	
//	fprintf(fp, "%12.3lf ", ansis_ann[55]);
//	fprintf(fp, "%12.3lf ", ansis_ann[64]);
//	fprintf(fp, "%12.3lf ", ansis_ann[73]);
	
//	fprintf(fp, "%12.3lf ", ansis_ann[58]);
//	fprintf(fp, "%12.3lf ", ansis_ann[67]);
//	fprintf(fp, "%12.3lf ", ansis_ann[76]);

//	fprintf(fp, "%12.3lf ", ansis_ann[91]);
//	fprintf(fp, "%12.3lf ", ansis_ann[92]);
//	fprintf(fp, "%12.3lf ", ansis_ann[93]);
//	fprintf(fp, "%12.3lf ", ansis_ann[94]);
//	fprintf(fp, "%12.3lf ", ansis_ann[95]);

/*	fprintf(fp, "%12.3lf ", ansis_ann[96]); 
	
	fprintf(fp, "%12.3lf ", ansis_ann[6]);
	fprintf(fp, "%12.3lf ", ansis_ann[7]);
	fprintf(fp, "%12.3lf ", ansis_ann[8]);
	
	fprintf(fp, "%12.5lf ", ansis_ann[99]/1000000.0);
	
	fprintf(fp, "%12.3lf ", ansis_ann[98]);
	fprintf(fp, "%12.3lf ", ansis_ann[103]);
	fprintf(fp, "%12.3lf ", ansis_ann[104]);
	fprintf(fp, "%12.3lf ", ansis_ann[105]);
	
	fprintf(fp, "%12.3lf ", ansis_ann[106]*1000.0);
	fprintf(fp, "%12.3lf ", ansis_ann[107]*1000.0);

	fprintf(fp, "%12.3lf ", ansis_ann[88]);  */

/*	fprintf(fp, "%12.3lf ", ansis_ann[6]+ansis_ann[10]+ansis_ann[14]);
	fprintf(fp, "%12.3lf ", ansis_ann[7]+ansis_ann[11]+ansis_ann[15]);
	fprintf(fp, "%12.3lf ", ansis_ann[8]+ansis_ann[12]+ansis_ann[16]);
	
	fprintf(fp, "%12.3lf ", c_litter);
	fprintf(fp, "%12.3lf ", c_humus);
	
	fprintf(fp, "%12.3lf ", ansis_ann[99]/1000000.0);
	fprintf(fp, "%12.3lf ", ansis_ann[100]/1000000.0);
	fprintf(fp, "%12.3lf ", ansis_ann[101]/100.0*12.0/16.0);
	fprintf(fp, "%12.3lf ", ansis_ann[102]/100.0*12.0/16.0);
	  */

//	fprintf(fp, "%12.3lf ", ansis_ann[91]+ansis_ann[92]-ansis_ann[94]); /* microbe N */

/*	fprintf(fp, "%12.3lf ", ansis_ann[98]); */ /* CH4 oxy */
/*	fprintf(fp, "%12.3lf ", ansis_ann[6]+ansis_ann[10]+ansis_ann[14]);*/	/* leaf */
/*	fprintf(fp, "%12.3lf ", ansis_ann[7]+ansis_ann[11]+ansis_ann[15]);*/	/* stem */
/*	fprintf(fp, "%12.3lf ", ansis_ann[8]+ansis_ann[12]+ansis_ann[16]);*/	/* root */

/*	fprintf(fp, "%12.3lf ", (ansis_ann[17]/15.0 + ansis_ann[18]/80.0 + ansis_ann[19]/50.0 
		+(ansis_ann[20]+ansis_ann[21]+ansis_ann[22])/40.0)*1000000.0); */

	/* fprintf(fp, "\n"); */
	
	fprintf(fp, "\n");
}

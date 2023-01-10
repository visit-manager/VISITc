/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in November 02, 2007						*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 21 by Akihiko Ito: Version 0.2.0.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/***************************************************************************/
void output_ansis_daily(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp
){
	float rr_d, sr_d, voc_d;
	/* float dppfd, dppfdd, dppfdb; */

	rr_d = (flux->tree).rrm + (flux->tree).rrg + loct->funder_c3*(flux->c3).rrm + 
			loct->funder_c3*(flux->c3).rrg + loct->funder_c4*(flux->c4).rrm + 
			loct->funder_c4*(flux->c4).rrg;
	sr_d = rr_d + (flux->soil).hr;
	
	/* total VOC, g C ha-1 day-1 */
	voc_d = (flux->voc_isopr_g97 + flux->voc_monotrp_g97 + flux->voc_methanl_g97 + 
			flux->voc_acetone_g97 + flux->voc_actaldhd_g97 + flux->voc_frmardhd_g97 + 
			flux->voc_formacd_g97 + flux->voc_acetacd_g97 + flux->voc_co_g97);

	fprintf(fp,"%4ld ", loct->climy);
	fprintf(fp,"%3ld ", loct->doy);
	
	fprintf(fp,"%8.3lf ", (flux->soil).ch4oxy_curry);
	
	
	
	
	
/*	fprintf(fp,"%9.6lf ", loct->prate_sfc);
	fprintf(fp,"%7.2lf ", loct->ro2); 
	fprintf(fp,"%9.6lf ", (flux->soil).doc_boyer);  */
	
/*	fprintf(fp,"%8.3lf ", flux->nep*100.0);
	
	fprintf(fp,"%9.6lf ", (flux->soil).ch4oxy_ridg);
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
	fprintf(fp,"%7.3lf ", loct->m_casa);  */
	
/*	fprintf(fp,"%7.2lf ", loct->dswrf_sfc);
	fprintf(fp,"%7.2lf ", loct->prate_sfc);
	fprintf(fp,"%7.2lf ", loct->tmp_2m);
	fprintf(fp,"%7.2lf ", loct->tmp10_soil);
	fprintf(fp,"%7.2lf ", loct->tmp200_soil);
	fprintf(fp,"%7.2lf ", loct->vpd);
	fprintf(fp,"%7.2lf ", loct->bCO2);
	fprintf(fp,"%.3lf ", flux->er*100.0);
	
	fprintf(fp,"%7.2lf ", mass->sw30);
	fprintf(fp,"%7.2lf ", mass->sww);
	
	fprintf(fp,"%9.4lf ", flux->gpp * 100.0);
	fprintf(fp,"%9.4lf ", flux->npp * 100.0);
	fprintf(fp,"%9.4lf ", flux->nep* 100.0);
		
	fprintf(fp,"%9.6lf ", (flux->soil).d_n2o_ngas);
	
	fprintf(fp,"%9.6lf ", (flux->soil).doc_boyer);
	fprintf(fp,"%9.6lf ", voc_d);

	fprintf(fp,"%7.2lf ", (echar->tree).n_opt_leaf);
	fprintf(fp,"%7.2lf ", (echar->c3).n_opt_leaf);
	fprintf(fp,"%7.2lf ", (echar->c4).n_opt_leaf);
	
	fprintf(fp,"%7.2lf ", (echar->tree).n_conc_larea);
	fprintf(fp,"%7.2lf ", (echar->c3).n_conc_larea);
	fprintf(fp,"%7.2lf ", (echar->c4).n_conc_larea);  */
	
/*	fprintf(fp,"%7.2lf ", (mass->tree).n_cnpy);
	fprintf(fp,"%7.2lf ", (mass->c3).n_cnpy);
	fprintf(fp,"%7.2lf ", (mass->c4).n_cnpy);  */

/*	fprintf(fp,"%7.2lf ", loct->depo_no3);
	fprintf(fp,"%7.2lf ", loct->depo_nh4); */

	/************************************************/
	/* EDIT here to write results into "daily" file */

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
	
	fprintf(fp,"%7.2lf ", (mass->tree).lai); 
	fprintf(fp,"%7.2lf ", (mass->c3).lai); 
	fprintf(fp,"%7.2lf ", (mass->c4).lai); */

/*	fprintf(fp,"%7.2lf ", loct->tmp_2m); 
	fprintf(fp,"%7.2lf ", (echar->c3).gdd); 
	fprintf(fp,"%7.2lf ", (echar->c3).cdd); 
	fprintf(fp,"%4ld ", (echar->c3).season); 
	fprintf(fp,"%7.2lf ", (mass->c3).lai); 
	fprintf(fp,"%7.2lf ", (echar->c3).opt_lai); */
	
	/* for PAR analysis */
/*	dppfd = dppfdd = dppfdb = 0.0;
	for(f=0;f<48;f++){
		dppfd += loct->ppfd_h[f]*1800.0/1000000.0;
		dppfdd += loct->ppfdd_h[f]*1800.0/1000000.0;
		dppfdb += loct->ppfdb_h[f]*1800.0/1000000.0;
	}
	
	fprintf(fp,"%f ",loct->tcdc_clm);
	fprintf(fp,"%f ",loct->daylen[loct->doy]);
	fprintf(fp,"%f ",loct->par_h);
	fprintf(fp,"%f ",loct->ppfd_h);
	fprintf(fp,"%f ",loct->csppfdb_h+loct->csppfdd_h);
	
	fprintf(fp,"%f ",dppfd);
	fprintf(fp,"%f ",dppfdd);
	fprintf(fp,"%f ",dppfdb); */
	
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
	float ansis_ann[256]
){
	float rr_d, sr_d, er_d, fscale, voc_d;
	
	fscale = 1.0;
	
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
	
	ansis_ann[30] += (mass->tree).n_cnpy/365.0;
	ansis_ann[31] += (mass->tree).n_strg/365.0;
	ansis_ann[32] += (mass->c3).n_cnpy/365.0;
	ansis_ann[33] += (mass->c3).n_strg/365.0;
	ansis_ann[34] += (mass->c4).n_cnpy/365.0;
	ansis_ann[35] += (mass->c4).n_strg/365.0;
	
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
	
	ansis_ann[59] += (flux->c3).n_biofix;
	ansis_ann[60] += (flux->c3).uptake_no3;
	ansis_ann[61] += (flux->c3).uptake_nh4;
	ansis_ann[62] += (flux->c3).n_abdn_cnpy;
	ansis_ann[63] += (flux->c3).n_abdn_strg;
	ansis_ann[64] += (flux->c3).n_salvage;
	ansis_ann[65] += (flux->c3).n_alloc_cnpy;
	ansis_ann[66] += (flux->c3).n_alloc_strg;
	ansis_ann[67] += (flux->c3).n_realloc;

	ansis_ann[68] += (flux->c4).n_biofix;
	ansis_ann[69] += (flux->c4).uptake_no3;
	ansis_ann[70] += (flux->c4).uptake_nh4;
	ansis_ann[71] += (flux->c4).n_abdn_cnpy;
	ansis_ann[72] += (flux->c4).n_abdn_strg;
	ansis_ann[73] += (flux->c4).n_salvage;
	ansis_ann[74] += (flux->c4).n_alloc_cnpy;
	ansis_ann[75] += (flux->c4).n_alloc_strg;
	ansis_ann[76] += (flux->c4).n_realloc;

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
}

/*********************************************************************/
void output_ansis_ann(
	long year, 
	float ansis_ann[256], 
	FILE *fp
){
	long f;
	float c_litter, c_humus;
	
	c_litter = ansis_ann[17]+ansis_ann[18]+ansis_ann[19]+ansis_ann[20]+ansis_ann[21]+ansis_ann[22];
	c_humus = ansis_ann[23]+ansis_ann[24]+ansis_ann[25];
	
	fprintf(fp, "%3ld ", year);
	
	fprintf(fp, "%12.3f ", ansis_ann[0]);
	fprintf(fp, "%12.3f ", ansis_ann[1]);
	fprintf(fp, "%12.3f ", ansis_ann[2]);
	fprintf(fp, "%12.3f ", ansis_ann[3]);
	fprintf(fp, "%12.3f ", ansis_ann[4]);
	
	fprintf(fp, "%12.3f ", ansis_ann[6]+ansis_ann[7]+ansis_ann[8]);
	fprintf(fp, "%12.3f ", ansis_ann[10]+ansis_ann[11]+ansis_ann[12]+ansis_ann[14]+ansis_ann[15]+ansis_ann[16]);
	fprintf(fp, "%12.3f ", c_litter);
	fprintf(fp, "%12.3f ", c_humus);
	
	
/*	
	fprintf(fp, "%12.3f ", ansis_ann[96]); 
	
	
	fprintf(fp, "%12.3f ", ansis_ann[6]);
	fprintf(fp, "%12.3f ", ansis_ann[7]);
	fprintf(fp, "%12.3f ", ansis_ann[8]);
	
	fprintf(fp, "%12.3f ", ansis_ann[2]);
	fprintf(fp, "%12.5lf ", ansis_ann[99]/1000000.0);
	
	fprintf(fp, "%12.3f ", ansis_ann[98]);
	fprintf(fp, "%12.3f ", ansis_ann[103]);
	fprintf(fp, "%12.3f ", ansis_ann[104]);
	fprintf(fp, "%12.3f ", ansis_ann[105]);
	
	fprintf(fp, "%12.3f ", ansis_ann[106]*1000.0);
	fprintf(fp, "%12.3f ", ansis_ann[107]*1000.0);

	fprintf(fp, "%12.3f ", ansis_ann[78]);
	fprintf(fp, "%12.3f ", ansis_ann[79]);
	fprintf(fp, "%12.3f ", ansis_ann[88]);  */

/*	fprintf(fp, "%12.3f ", ansis_ann[6]+ansis_ann[10]+ansis_ann[14]);
	fprintf(fp, "%12.3f ", ansis_ann[7]+ansis_ann[11]+ansis_ann[15]);
	fprintf(fp, "%12.3f ", ansis_ann[8]+ansis_ann[12]+ansis_ann[16]);
	
	fprintf(fp, "%12.3f ", c_litter);
	fprintf(fp, "%12.3f ", c_humus);
	
	fprintf(fp, "%12.3f ", ansis_ann[99]/1000000.0);
	fprintf(fp, "%12.3f ", ansis_ann[100]/1000000.0);
	fprintf(fp, "%12.3f ", ansis_ann[101]/100.0*12.0/16.0);
	fprintf(fp, "%12.3f ", ansis_ann[102]/100.0*12.0/16.0);

	
	fprintf(fp, "%12.3f ", ansis_ann[36]); 
	fprintf(fp, "%12.3f ", ansis_ann[37]); 
	fprintf(fp, "%12.3f ", ansis_ann[38]); 
	fprintf(fp, "%12.3f ", ansis_ann[39]); 
	fprintf(fp, "%12.3f ", ansis_ann[40]); 
	fprintf(fp, "%12.3f ", c_litter*1000000.0/ansis_ann[38]); 
	fprintf(fp, "%12.3f ", c_humus*1000000.0/ansis_ann[39]);  */

/*	fprintf(fp, "%12.3f ", ansis_ann[49]); */ /* NH4 deposition */
/*	fprintf(fp, "%12.3f ", ansis_ann[48]); */ /* NO3 deposition */
/*	fprintf(fp, "%12.3f ", ansis_ann[78]); */ /* N2O nitrif */
/*	fprintf(fp, "%12.3f ", ansis_ann[79]); */ /* N2O denit */
	
//	fprintf(fp, "%12.3f ", ansis_ann[36]); /* NH4 */
//	fprintf(fp, "%12.3f ", ansis_ann[37]); /* NO3 */
//	fprintf(fp, "%12.3f ", ansis_ann[38]); /* litter N */
//	fprintf(fp, "%12.3f ", ansis_ann[39]); /* humus N */
//	fprintf(fp, "%12.3f ", ansis_ann[40]); /* microbe N */
//	fprintf(fp, "%12.3f ", ansis_ann[91]+ansis_ann[92]-ansis_ann[94]); /* microbe N */

/*	fprintf(fp, "%12.3f ", ansis_ann[98]); */ /* CH4 oxy */
	
/*	fprintf(fp, "%12.3f ", ansis_ann[6]+ansis_ann[10]+ansis_ann[14]);*/	/* leaf */
/*	fprintf(fp, "%12.3f ", ansis_ann[7]+ansis_ann[11]+ansis_ann[15]);*/	/* stem */
/*	fprintf(fp, "%12.3f ", ansis_ann[8]+ansis_ann[12]+ansis_ann[16]);*/	/* root */

/*	fprintf(fp, "%12.3f ", (ansis_ann[17]/15.0 + ansis_ann[18]/80.0 + ansis_ann[19]/50.0 
		+(ansis_ann[20]+ansis_ann[21]+ansis_ann[22])/40.0)*1000000.0); */

/*	fprintf(fp, "%12.3f ", ansis_ann[30]); 
	fprintf(fp, "%12.3f ", ansis_ann[31]); 
	fprintf(fp, "%12.3f ", ansis_ann[32]); 
	fprintf(fp, "%12.3f ", ansis_ann[33]); 
	fprintf(fp, "%12.3f ", ansis_ann[34]); 
	fprintf(fp, "%12.3f ", ansis_ann[35]); */
	
/*	for(f=0;f<=95;f++){
		fprintf(fp, "%7.2lf ", ansis_ann[f]); 
	} */

	/* fprintf(fp, "\n"); */
}

/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/*   initialize    environmental   locations  */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

#define TERM_HYD  0.1

/** initialization of climatic conditions (Primary data) ********************/
void f_loct_proc(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	long h;
	double  aaa, bbb, ccc, alt, tcK;
	double temp_factor, prec_factor;
	double dat_0, dat_1, dat_2, rr;
	double temp_diff, prec_diff, ahmd_diff, dswrf_diff;
	double mi[12] = {31.0, 29.5, 29.5, 30.5, 30.5, 30.5, 30.5, 31.0, 30.5, 30.5, 30.5, 30.5};
	double mm[12] = {16.0, 16.0, 13.0, 16.0, 15.0, 16.0, 15.0, 16.0, 16.0, 15.0, 16.0, 15.0};
	double temp_wclim_d, temp_corr, prec_corr;
	
	/* vegetation  ********************************************************/
	if(loct->time==1 && loct->doy==0){
		/* herbaceous C3/C4 composition (fixed) *************************/
		switch(grid->veg_type){
			case 9: case 10:
				/* loct->funder_c4 = 0.5; */  /* assumption */
				
				temp_factor = (grid->atmp2m_av + 5.0)*0.032;
				prec_factor = 1.0/(grid->prate_sfc_ann * 0.00025 + 0.65);
				loct->funder_c4 = temp_factor * prec_factor;
				loct->funder_c4 = (loct->funder_c4 > 0.0)?loct->funder_c4:0.0;
				loct->funder_c4 = (loct->funder_c4 < 0.95)?loct->funder_c4:0.95;
				break;
			default:
				loct->funder_c4 = 0.0;
		}
		
		if(strcmp(grid->site_id, "QHB")==0){
			loct->funder_c4 = 0.0;
		}
		
		/* prescribed C4 fraction at KBU site: 2011/09/29 by A.Ito **/
		if(strcmp(grid->site_id, "KBU")==0){
			loct->funder_c4 = 0.5; /* */
		}
		if(strcmp(grid->site_id, "Tongyu")==0){
			loct->funder_c4 = 0.5; /* */
		}
		
		loct->funder_c3 = (1.0 - loct->funder_c4);
		
		if(NOTICE == 1){
			printf("Floor C3/C4 = %lf/%lf\n", loct->funder_c3, loct->funder_c4);
		}
	}
	
	/* clear fluxes ********************************************/
	pflux_zero(&(flux->tree));
	pflux_zero(&(flux->c3));
	pflux_zero(&(flux->c4));
	sflux_zero(&(flux->soil));
	flux->gpp = 0.0;
	flux->npp = 0.0;
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

	/*  leaf area index (LAI), m2 m-2 ***/
	if(grid->veg_type == 4){
		/* deciduous broadleaved forest (e.g. Takayama) */
		f_sla_change(grid, loct, echar, mass);
	}
	
	(mass->tree).lai = lai_mass(&(mass->tree), &(echar->tree));
	(mass->c3).lai = lai_mass(&(mass->c3), &(echar->c3));
	(mass->c4).lai = lai_mass(&(mass->c4), &(echar->c4));
	loct->lai = (mass->tree).lai + (mass->c3).lai*loct->funder_c3 + 
			(mass->c4).lai*loct->funder_c4;
	
	/* fractional cover by vegetation and soil */
	aaa = 1.0 - exp(-(echar->tree).eK0*(mass->tree).lai);
	loct->fcover_tree = aaa;
	bbb = 1.0 - exp(-(echar->c3).eK0*(mass->c3).lai);
	ccc = 1.0 - exp(-(echar->c4).eK0*(mass->c4).lai);
	loct->fcover_c3 = (1.0 - aaa)*loct->funder_c3*bbb;
	loct->fcover_c4 = (1.0 - aaa)*loct->funder_c4*ccc;
	loct->fcover_veg = loct->fcover_tree + loct->fcover_c3 + loct->fcover_c4;
	if(loct->fcover_veg<0.0 || loct->fcover_veg>1.0){
		printf("!! Wrong vegetation cover: %lf\n", loct->fcover_veg);
	}
	loct->fcover_ground = 1.0 - loct->fcover_veg;
    
    loct->rooting_depth = 0.1;
    if((mass->tree).rooting_depth > loct->rooting_depth){
        loct->rooting_depth = (mass->tree).rooting_depth;
    }
    if((mass->c3).rooting_depth > loct->rooting_depth){
        loct->rooting_depth = (mass->c3).rooting_depth;
    }
    if((mass->c4).rooting_depth > loct->rooting_depth){
        loct->rooting_depth = (mass->c4).rooting_depth;
    }
	
	/* climate **************************************************************/
	if(loct->doy == 0){
		loct->prate_ann = 0.0;
	}
	/*  sequential climatic conditions  */
	if(loct->climy >= BYR && loct->climy <= EYR){
		loct->dswrf_sfc = dswrf_sfc_d[loct->climy-BYR][loct->doy];
		loct->tcdc_clm = tcdc_clm_d[loct->climy-BYR][loct->doy]/100.0;  
		loct->tmp_2m = tmp_2m_d[loct->climy-BYR][loct->doy]-ZAT;
		loct->tmp_sfc = tmp_sfc_d[loct->climy-BYR][loct->doy]-ZAT;
		loct->tmp10_soil = tmp10_soil_d[loct->climy-BYR][loct->doy]-ZAT;
		loct->tmp200_soil = tmp200_soil_d[loct->climy-BYR][loct->doy]-ZAT;
		loct->spfh_2m = spfh_2m_d[loct->climy-BYR][loct->doy];
		loct->prate_sfc = prate_sfc_d[loct->climy-BYR][loct->doy];		
		loct->wnd_10m = wind_10m_d[loct->climy-BYR][loct->doy];
		
		if((strcmp(grid->site_id, "PSO")==0) && (loct->climy<1955)){
			loct->prate_sfc = prate_sfc_dav[loct->doy];	
		}
	}
	
	/* future prediction */
	if(PREDICT >= 1 && loct->climy > EYR){
		temp_diff = gcm_temp[loct->climy-1980][loct->doy] - gcm_temp_av[loct->doy];
		prec_diff = gcm_prec[loct->climy-1980][loct->doy] - gcm_prec_av[loct->doy];
		dswrf_diff = gcm_dswrf[loct->climy-1980][loct->doy] - gcm_dswrf_av[loct->doy];
		
		/* temp_diff = prec_diff = dswrf_diff = 0.0; */
		
		loct->dswrf_sfc = dswrf_sfc_dav[loct->doy] + dswrf_diff;
		loct->tmp_2m = tmp_2m_dav[loct->doy]-ZAT + temp_diff;
		loct->tmp_sfc = tmp_sfc_dav[loct->doy]-ZAT + temp_diff;
		loct->tmp10_soil = tmp10_soil_dav[loct->doy]-ZAT + temp_diff;
		loct->tmp200_soil = tmp200_soil_dav[loct->doy]-ZAT + temp_diff;
		loct->prate_sfc = prate_sfc_dav[loct->doy] + prec_diff;
		if(loct->prate_sfc < 0.0){
			loct->prate_sfc = 0.0;
		}
		
		loct->wnd_10m = wind_10m_d[EYR-BYR][loct->doy];
		loct->tcdc_clm = tcdc_clm_d[EYR-BYR][loct->doy]/100.0;
	}
	
	/* annual precipitation */
	loct->prate_ann += loct->prate_sfc;
	
	/* sensitivity analysis ***********************************************/
	if(SENSANS_TEMP==1 && loct->adyear >= YLDIST){  // BYR
		loct->tmp_2m = tmp_2m_dav[loct->doy]-ZAT;
		loct->tmp_sfc = tmp_sfc_dav[loct->doy]-ZAT;
		loct->tmp10_soil = tmp10_soil_dav[loct->doy]-ZAT;
		loct->tmp200_soil = tmp200_soil_dav[loct->doy]-ZAT;
	}
	
	if(SENSANS_PREC==1 && loct->adyear >= YLDIST){  // BYR
		loct->prate_sfc = prate_sfc_dav[loct->doy];
	}

	if(SENSANS_SRAD==1 && loct->adyear >= YLDIST){  // BYR
		loct->dswrf_sfc = dswrf_sfc_dav[loct->doy];
		loct->tcdc_clm = tcdc_clm_dav[loct->doy]/100.0;
	}
	
	/* LARS-generated dummy weather data (60-year) */
	if(SENSANS_TEMP==2){
		loct->tmp_2m = (lars_tmin[loct->climy-BYR][loct->doy] + lars_tmax[loct->climy-BYR][loct->doy])/2.0;
		loct->tmp_sfc = loct->tmp_2m + (tmp_sfc_d[loct->climy-BYR][loct->doy] - tmp_2m_d[loct->climy-BYR][loct->doy]);
		loct->tmp10_soil = loct->tmp_2m + (tmp10_soil_d[loct->climy-BYR][loct->doy] - tmp_2m_d[loct->climy-BYR][loct->doy]);
		loct->tmp200_soil = loct->tmp_2m + (tmp200_soil_d[loct->climy-BYR][loct->doy] - tmp_2m_d[loct->climy-BYR][loct->doy]);
	}
	if(SENSANS_PREC==2){
		loct->prate_sfc = lars_prec[loct->climy-BYR][loct->doy];
	}
	if(SENSANS_SRAD==2){
		
		loct->dswrf_sfc = lars_srad[loct->climy-BYR][loct->doy]*1000000.0/3600.0/24.0;
	}
	
	/** climate corrections for FLUX site *******************/
	if(WMODE == 1 && (strcmp(grid->site_id, "TKY")==0)){ 
		f_clim_correct_TKY(grid, loct);		
	}
	if(WMODE == 1 && (strcmp(grid->site_id, "TMK")==0)){
		f_clim_correct_TMK(grid, loct);
	}
	if(WMODE == 1 && (strcmp(grid->site_id, "LSH")==0)){
		f_clim_correct_LSH(grid, loct);
	}
	if(WMODE == 1 && (strcmp(grid->site_id, "TSE")==0)){
		f_clim_correct_TSE(grid, loct);
	}
	if(WMODE == 1 && (strcmp(grid->site_id, "SKT")==0)){
		f_clim_correct_SKT(grid, loct);
	}
	if(WMODE == 1 && (strcmp(grid->site_id, "MKL")==0)){
		f_clim_correct_MKL(grid, loct);
	}
	if(WMODE == 1 && (strcmp(grid->site_id, "QHB")==0)){
		f_clim_correct_QHB(grid, loct);
	}
	
	
	if(WMODE == 1 && (strcmp(grid->site_id, "KBU")==0)){
		//f_clim_correct_KBU(grid, loct);
	}
	if(WMODE == 1 && (strcmp(grid->site_id, "Tongyu")==0)){
		//f_clim_correct_Tongyu(grid, loct);
	}
	
	/* local climate by WorldClim ***************************************/
	if(WMODE == 2){
		/* temperature (daily intepolation) */
		if(loct->month==0){
			dat_0 = grid->temp_wclim[11];
		}else{
			dat_0 = grid->temp_wclim[loct->month-1];
		}
								
		if(loct->month==11){
			dat_1 = grid->temp_wclim[11];
			dat_2 = grid->temp_wclim[0];
		}else{
			dat_1 = grid->temp_wclim[loct->month];
			dat_2 = grid->temp_wclim[loct->month+1];
		}
		rr = (double)(loct->mday + mm[loct->month])/mi[loct->month];
					
		temp_wclim_d = dat_0 + rr * (dat_1 - dat_0) + rr * (rr - 1.0) /2.0 * 
								((dat_2 - dat_1) - (dat_1 - dat_0));
				
		temp_corr = temp_wclim_d - tmp_2m_ncep_dav[loct->doy];

		/* precipitation */
		if(prate_sfc_ncep_mav[loct->month]>0.0){
			prec_corr = grid->prec_wclim[loct->month] / prate_sfc_ncep_mav[loct->month];
		}else{
			prec_corr = 0.0;
		}
		
		/* check */
		if(NOTICE==1 && grid->flag_datavl==1 && (temp_corr>20.0 || temp_corr<-20.0)){
			printf("! Suspicious temp_corr:%ld %.2lf %.2lf %.2lf\n", loct->doy, temp_wclim_d, 
						tmp_2m_ncep_dav[loct->doy], temp_corr);
		}
		if(NOTICE==1 && grid->flag_datavl==1 && (prec_corr<0.0 || prec_corr>10.0)){
			printf("! Suspicious prec_corr: %.2lf %.2lf %.2lf\n", grid->prec_wclim[loct->month],
					prate_sfc_ncep_mav[loct->month], prec_corr);
		}
		
		/* correction */
		loct->tmp_2m += temp_corr;
		loct->tmp_sfc += temp_corr;
		loct->tmp10_soil += temp_corr;
		loct->tmp200_soil += temp_corr;
		
		loct->prate_sfc *= prec_corr;
	}
	
	/*  snow-packing effect on soil temperature  */
	if(mass->snwa > 0.2){
		aaa = mass->snwa/(3.0 + mass->snwa);
		loct->tmp10_soil = 0.0*aaa + loct->tmp10_soil*(1.0-aaa);
		loct->tmp200_soil = 2.0*aaa + loct->tmp200_soil*(1.0-aaa);
	}
	
	/* atmosphere ****************************************************/
	/* CO2 ********************************/
	atmco2_trend(grid, loct);
	co2_in_canopy(grid, loct, mass, flux);
	
	/* air conditions  **************/
	/*  altitude  */
	alt = (grid->topo>=0.0)?grid->topo:0.0;  
	/*  air pressure */
	loct->air_prsr = 1013.25 * exp(-1.0*(28.964*0.001)*9.8*alt/(8.3144*(loct->tmp_2m + ZAT)));  
	
	/*  saturated vapour pressure, hPa  */
	loct->vps = f_vap_pre_sat(grid, loct);    
	/*  slope of saturated vapour pressure as a function of temperature  */
	loct->slope_vps = f_slope_vps(grid, loct);	
	 /*  aerodynamic resistance  */    
	loct->r_aero = f_r_aero(grid, loct);	   
	
	/* Future change in specific humidity by CCCma or ECHAM */
	if(loct->climy > EYR && (PREDICT == 10 || PREDICT == 20 || PREDICT == 21 || PREDICT == 22 ||
        PREDICT == 101 || PREDICT == 102 || PREDICT == 103 || PREDICT == 104)){
		ahmd_diff = gcm_ahmd[loct->climy-1980][loct->doy] - gcm_ahmd_av[loct->doy];
		
		loct->spfh_2m += ahmd_diff;
		if(loct->spfh_2m < 0.0){
			loct->spfh_2m = 0.0;
		}
	}
	
	/*  vapour pressure, hPa  */
	if(strcmp(grid->site_id, "CEAMIP_TKY")==0 || strcmp(grid->site_id, "CEAMIP_TUR")==0
            || strcmp(grid->site_id, "CEAMIP_TKC")==0|| strcmp(grid->site_id, "CEAMIP_SKR")==0
            || strcmp(grid->site_id, "CEAMIP_TMK")|| strcmp(grid->site_id, "CEAMIP_LSH")==0
            || strcmp(grid->site_id, "CEAMIP_MMF")==0|| strcmp(grid->site_id, "CEAMIP_MBF")==0
            || strcmp(grid->site_id, "CEAMIP_MKL")==0|| strcmp(grid->site_id, "CEAMIP_YLF")==0
            || strcmp(grid->site_id, "CEAMIP_YPF")==0|| strcmp(grid->site_id, "CEAMIP_SKT")==0
            || strcmp(grid->site_id, "CEAMIP_GDK")==0|| strcmp(grid->site_id, "CEAMIP_SMF")==0
            || strcmp(grid->site_id, "CEAMIP_TSE")==0|| strcmp(grid->site_id, "CEAMIP_CBS")==0
            || strcmp(grid->site_id, "CEAMIP_QYZ")==0|| strcmp(grid->site_id, "CEAMIP_PDF")==0
            || strcmp(grid->site_id, "CEAMIP_KBU")==0|| strcmp(grid->site_id, "CEAMIP_QHB")==0
            || strcmp(grid->site_id, "CEAMIP_HBG")==0|| strcmp(grid->site_id, "CEAMIP_YCS")==0
            || strcmp(grid->site_id, "CEAMIP_HFK")==0|| strcmp(grid->site_id, "CEAMIP_MSE")==0){
		loct->vp = loct->vps * spfh_2m_d[loct->climy-BYR][loct->doy];
	}else{
		loct->vp = loct->air_prsr*loct->spfh_2m/(0.622+0.378*loct->spfh_2m);    
	}
	
	/*  vapour pressure deficit, hPa  */
	loct->vpd = vpd_d[loct->climy-BYR][loct->doy];
	loct->vpd = (loct->vpd>=0.0)?loct->vpd:0.0;     /*  */ 
	
	/* Future change in specific humidity by MIROC-high */
	if((PREDICT == 1 || PREDICT == 100) && loct->climy > EYR){
		ahmd_diff = gcm_ahmd[loct->climy-1980][loct->doy] - gcm_ahmd_av[loct->doy];
		
		/* ahmd_diff = 0.0; */
		
		loct->vpd = loct->vps * (loct->vpd/loct->vps - ahmd_diff/100.0);
		if(loct->vpd < 0.0){
			loct->vpd = 0.0;
		}
		if(loct->vpd > loct->vps){
			loct->vpd = loct->vps;
		}
        loct->vp = loct->vps - loct->vpd;
	}
		
	/*  air density, kg m-3   */
	loct->air_dns = f_airdens(grid, loct);    
	
	/* latent heat MJ kg-1 */
	loct->lhvp = 1000000.0* (2.501 - 0.012/5.0*loct->tmp_2m);
	
	/* N deposition */
	/* f_n_deposit(grid, loct); */
	/* experiment for AGU2008F */
	if(loct->CO2y<1980){
		loct->depo_no3 = loct->depo_no3_model_av[loct->doy];
		loct->depo_nh4 = loct->depo_nh4_model_av[loct->doy];
	}else if(loct->CO2y>=1980 && loct->CO2y<=2003){
		loct->depo_no3 = loct->depo_no3_model[loct->climy - 1980][loct->doy];
		loct->depo_nh4 = loct->depo_nh4_model[loct->climy - 1980][loct->doy];
	}else if(loct->CO2y>2003){
		loct->depo_no3 = loct->depo_no3_model_av[loct->doy];
		loct->depo_nh4 = loct->depo_nh4_model_av[loct->doy];
	} /* */

	/* loct->depo_no3 = loct->depo_no3_model[0][loct->doy];
	loct->depo_nh4 = loct->depo_nh4_model[0][loct->doy]; */
	
	/* loct->depo_no3 = loct->depo_no3_model_av[loct->doy];
	loct->depo_nh4 = loct->depo_nh4_model_av[loct->doy]; */
	
	/* N fertilizer */
	loct->nfert_no3 = 0.0;
	loct->nfert_nh4 = 0.0;

	/* radiation *********************************************************************/
	/*  radiation fluxes and solar height */
	for(h=0; h<DSTEP; h++){
		loct->hour = h;
		loct->hangle = -180.0 + (double)h*7.5;
		/* canopy-top PPFD */
		f_ppfd(grid, loct); 
		
		if(loct->solhgt_h[loct->doy][h] > 0.0){
			tcK = ((echar->tree).eK0) / sin(loct->solhgt_h[loct->doy][h]*dTr);
			loct->grdrad_h[h] = exp(-tcK*(mass->tree).lai);
		}else{
			loct->grdrad_h[h] = 0.0;
		}
	}
	loct->hour = 24;  /* midday */

	/* soil surface albedo **********/
	loct->snow_acc = mass->snwa;
	(echar->soil).albedo = albedo_soil(loct, &(echar->soil));

	/* give irradiance attenuation coefficient  */
	(echar->tree).eK = irr_attn(grid, loct, &(echar->tree));
	(echar->c3).eK = irr_attn(grid, loct, &(echar->c3));
	(echar->c4).eK = irr_attn(grid, loct, &(echar->c4));
	
	/* net radiation */	
	f_net_rad(grid, loct, echar, mass);
	
	/* hydrological water budget  ************************************************/
    /* soil degradation */
    if(strcmp(grid->site_id, "PSO")==0 && DGSOIL==1){
    
        if(loct->veg_state==1 && loct->age_stand>=1 && loct->doy==0){
            /* loss of soil clay */
            aaa = 0.01 * grid->clay_frac;
            
            grid->clay_frac -= aaa;
            grid->sand_frac += aaa;
            
            if(grid->clay_frac < 0.05){
                grid->clay_frac = 0.05;
            }
            if(grid->sand_frac > 0.75){
                grid->sand_frac = 0.75;
            }
            
            grid->bulkdens += 0.0075;
           
            f_soil_saxton(grid);
            grid->fieldcap30 = grid->field_cap;
            //grid->fieldcap = grid->fieldcap - grid->fieldcap30;
            grid->fieldcap = grid->field_cap * 1.0/0.3;
        }
    }
    
	f_hydrology(grid, loct, echar, mass);
	
	loct->snow_acc = mass->snwa;
	loct->soilwtr_l = mass->sw30;
	loct->soilwtr_h = mass->sww;
	
	/* soil apperture **/
	loct->soilappr_l = (grid->fieldcap30 - mass->sw30)/grid->fieldcap30;
	loct->soilappr_l = (loct->soilappr_l < 1.0)?loct->soilappr_l:1.0;
	loct->soilappr_l = (loct->soilappr_l > 0.0)?loct->soilappr_l:0.0;
	loct->soilappr_w = (grid->fieldcap - mass->sww)/grid->fieldcap;
	loct->soilappr_w = (loct->soilappr_w < 1.0)?loct->soilappr_w:1.0;
	loct->soilappr_w = (loct->soilappr_w > 0.0)?loct->soilappr_w:0.0;
	
	/* water-filled pore space ***********/
	loct->wfps = ((((mass->sw30/10.0)/30.0)*100.0) / ((1.0 - grid->bulkdens/2.65)*100.0));  /**/
/*	loct->wfps = ((((loct->sw30/10.0)/30.0)*100.0) / ((1.0 - grid->bulkdens/1.325)*100.0)); */
	
	if(loct->wfps < 0.05){
		loct->wfps = 0.05;
	}
	
	/* water potential: 2011/10/25 by A.Ito **********************/
	/* litter ***********/
	loct->pot_grav_l = -0.05;
	if(mass->sw30>=0.2){
		if(grid->stexture==0){
			loct->pot_matric_l = -0.121 * pow(mass->sw30/grid->fieldcap30,-4.05);
		}else if(grid->stexture==1){
			loct->pot_matric_l = -0.478 * pow(mass->sw30/grid->fieldcap30,-5.39);
		}else if(grid->stexture==2){
			loct->pot_matric_l = -0.405 * pow(mass->sw30/grid->fieldcap30,-11.4);
		}
	}else{
		if(grid->stexture==0){
			loct->pot_matric_l = -0.121 * pow(0.2/grid->fieldcap30,-4.05);
		}else if(grid->stexture==1){
			loct->pot_matric_l = -0.478 * pow(0.2/grid->fieldcap30,-5.39);
		}else if(grid->stexture==2){
			loct->pot_matric_l = -0.405 * pow(0.2/grid->fieldcap30,-11.4);
		}
	}
	loct->pot_total_l = loct->pot_grav_l + loct->pot_matric_l;
	/* humus ***********/
	loct->pot_grav_h = -1.00;
	if(mass->sww>=0.2){
		if(grid->stexture==0){
			loct->pot_matric_h = -0.121 * pow(mass->sww/grid->fieldcap,-4.05);
		}else if(grid->stexture==1){
			loct->pot_matric_h = -0.478 * pow(mass->sww/grid->fieldcap,-5.39);
		}else if(grid->stexture==2){
			loct->pot_matric_h = -0.405 * pow(mass->sww/grid->fieldcap,-11.4);
		}
	}else{
		if(grid->stexture==0){
			loct->pot_matric_h = -0.121 * pow(0.2/grid->fieldcap,-4.05);
		}else if(grid->stexture==1){
			loct->pot_matric_h = -0.478 * pow(0.2/grid->fieldcap,-5.39);
		}else if(grid->stexture==2){
			loct->pot_matric_h = -0.405 * pow(0.2/grid->fieldcap,-11.4);
		}
	}
	loct->pot_total_h = loct->pot_grav_h + loct->pot_matric_h;

	/* ecophysiology: changing  ***************************************************/
	f_ecophysiology(grid, loct, &(echar->tree), &(mass->tree));
	f_ecophysiology(grid, loct, &(echar->c3), &(mass->c3));
	f_ecophysiology(grid, loct, &(echar->c4), &(mass->c4));
	
	if(WMODE == 1 && (strcmp(grid->site_id, "TMK")==0)){
		(echar->tree).lc = (echar->tree).lc0 * (0.3 + 0.9*exp(-0.03*loct->age_stand));
		(echar->tree).lr = (echar->tree).lr0 * (0.3 + 0.9*exp(-0.03*loct->age_stand));
	}
    
    /* developmental process in rice***********************************************/
	if(loct->doy==0){
		loct->dvi = 0.0;
	}
	
	if(loct->dvi<0.20){
		loct->dvr = (1/51.37) * (1/(1+exp(-0.3504*(loct->tmp_2m-18.81))));
	}else if(loct->dvi>=0.20&&loct->daylen[loct->doy]<15.64){
		loct->dvr = (1/51.37) * ((1-exp(0.6736*(loct->daylen[loct->doy]-15.64)))/(1+exp(-0.3504*(loct->tmp_2m-18.81))));
	}else if(loct->dvi>=0.20&&loct->daylen[loct->doy]>=15.64){
		loct->dvr = 0.0;
	}
	
	if(loct->doy>=100){
		loct->dvi += loct->dvr;
	}
}

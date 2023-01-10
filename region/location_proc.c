/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in November 02, 2007						*/

/*  [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
*/

/*   initialize    environmental   locations  */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

#define TERM_HYD  0.1

/* climatology: *_d[] means the time-series during 1948 to 2004- */	
extern float tmp_sfc_d[PERIOD][366];		/* ground surface temperature, degree Celcius */
extern float tmp_2m_d[PERIOD][366];		/* 2m air temperature, degree Celcius */
extern float tmp10_soil_d[PERIOD][366];	/* soil temperature at 10 cm depth, degree Celcius */
extern float tmp200_soil_d[PERIOD][366];	/* soil temperature at 200 cm depth, degree Celcius */
extern float dswrf_sfc_d[PERIOD][366];		/* downward shortwave radiation at the surface, W m-2 */
extern float tcdc_clm_d[PERIOD][366];		/* total cloudiness, fraction */
extern float prate_sfc_d[PERIOD][366];		/* precipitation, mm mon-1 */
extern float spfh_2m_d[PERIOD][366];		/* specific humidity, kg kg-1 */
extern float wind_10m_d[PERIOD][366];		/* u-wind velocity, m s-1 */
extern float vpd_d[PERIOD][366];   		/* VPD, hPa */
extern float tmp_2m_ncep_dav[366];			/* 2m air temperature, degree Celcius */
extern float prate_sfc_ncep_mav[12];		/* precipitation, mm mon-1 */
extern float tmp_sfc_dav[366];	
extern float tmp_2m_dav[366];	
extern float tmp10_soil_dav[366];	
extern float tmp200_soil_dav[366];	
extern float prate_sfc_dav[366];	
extern float dswrf_sfc_dav[366];	
extern float spfh_2m_dav[366];	
extern float tcdc_clm_dav[366];
extern long	WMODE;

/** initialization of climatic conditions (Primary data) *************************/
void f_init_cond(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	float  aaa, bbb, ccc, alt;
	float temp_factor, prec_factor;
	
	/* clear fluxes ****************************************************/
	pflux_zero(&(flux->tree));
	pflux_zero(&(flux->c3));
	pflux_zero(&(flux->c4));
	sflux_zero(&(flux->soil));
	flux->gpp = 0.0;
	flux->npp = 0.0;
	flux->nep = 0.0;
	flux->nbp = 0.0;
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
    flux->lue_conv = 0.0; 
    flux->lue_detr_1 = 0.0; 
    flux->lue_detr_10 = 0.0;
    flux->lue_detr_100 = 0.0;

    /* orbital radiation components *********************/
    /* hour angle: local time */
	/* loct->hangle = -180.0 + (float)(loct->hour)*15.0 + (float)(grid->lon) +90.0; */ /* 2022/01/13 */
	loct->hangle = -180.0 + (float)(loct->hour)*15.0 + (float)(grid->lon); /* 2022/08/01 */
    if(loct->hangle > 180.0){
        loct->hangle -= 360.0;
    }
    if(loct->hangle < -180.0){
        loct->hangle += 360.0;
    }
    
    if(loct->hour == 0){
        /* solar declination */
        loct->soldec = f_soldec(grid, loct); 
        /* day length */
        loct->daylen = f_daylen(grid, loct);
    }
		
    /* solar height */
    loct->solhgt_h = f_solhgt(grid, loct);
    
    /* atmosphere-top solar radiation */
    loct->toprad_h = f_toprad(grid, loct); 
	
	/* vegetation  *************************************************/
	/* if(grid->time==1 && loct->doy==0){ */
		/* C3/C4 composition (fixed) *************************/
		switch(grid->veg_type){
			case 9: case 10:
				/* loct->funder_c4 = 0.5; */  /* assumption */
				
				temp_factor = (grid->tmp2m_ann + 5.0)*0.032;
				prec_factor = 1.0/(grid->prec_ann * 0.00025 + 0.65);
				loct->funder_c4 = temp_factor * prec_factor;
				loct->funder_c4 = (loct->funder_c4 > 0.0)?loct->funder_c4:0.0;
				loct->funder_c4 = (loct->funder_c4 < 0.95)?loct->funder_c4:0.95;
				break;
			default:
				;
				loct->funder_c4 = 0.0;
		}
    
        if(strcmp(grid->site_id, "BAMIYAN")==0){
            loct->funder_c4 = 0.0;
        }
        
		loct->funder_c3 = (1.0 - loct->funder_c4);
		
	/* } */
	
	/*  leaf area index (LAI), m2 m-2 ***/
	if(grid->veg_type == 4){
		/* deciduous broadleaved forest (e.g. Takayama) */
		/* f_sla_change(loct, echar, mass); */
	}
	
	(mass->tree).lai = lai_mass(&(mass->tree), &(echar->tree));
	(mass->c3).lai = lai_mass(&(mass->c3), &(echar->c3));
	(mass->c4).lai = lai_mass(&(mass->c4), &(echar->c4));
    
    /* constrain max LAI of woody plants using meta-analysis */
    if(CONST_LAIMAX == 1 && loct->phase==1){
        if((mass->tree).lai > loct->laimax_meta){
            (mass->tree).lai = loct->laimax_meta;
        }
    }
    
	/* DEBUG *********************************************************** 2009/04/18 **** !!! for debug !!! */
	/*  (mass->tree).lai = 3.0;
	(mass->c3).lai = 1.0;
	(mass->c4).lai = 0.0;  */
	
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
		printf("!! Wrong vegetation cover: %f\n", loct->fcover_veg);
	}
	loct->fcover_ground = 1.0 - loct->fcover_veg;
	
	/* climate **************************************************************/
    /* 	if(loct->doy == 0){
		loct->prate_ann = 0.0;
	} */
	/*  sequential climatic conditions  */
	/* deleted 2009/05/25 by A.Ito **/
	
	/* regional simulation data */
	if(WMODE == 2){
		/* USE spatial data -**********************/
		loct->dswrf_sfc = grid->srad_region;
        if(loct->dswrf_sfc < 0.0){
            loct->dswrf_sfc = 0.0;
        }
        
		loct->tmp_2m = (grid->tmax_region + grid->tmin_region)/2.0;
        
		loct->prate_sfc = grid->prec_region;
        if(loct->prate_sfc < 0.0){
            loct->prate_sfc = 0.0;
        }
        
		loct->wnd_10m = grid->wind_region;
		
		/*!!!!!!!!! assumed */
		loct->tmp_sfc = loct->tmp_2m;
		//loct->tmp_sfc = loct->tmp_2m + 0.01*loct->dswrf_sfc;
		loct->tmp10_soil = loct->tmp_2m;
		//loct->tmp10_soil = loct->tmp_2m + 0.005*loct->dswrf_sfc;
		loct->tmp200_soil = loct->tmp_2m;
        
		/*!!!!!!!!! assumed */
        if(loct->toprad_h > 0.0){
            aaa = loct->dswrf_sfc / loct->toprad_h;
            aaa = (aaa>1.0)?1.0:aaa;
            aaa = (aaa<0.0)?0.0:aaa;

            bbb = (0.8964 - aaa) / 0.5392;
            bbb = (bbb>1.0)?1.0:bbb;
            bbb = (bbb<0.0)?0.0:bbb;
        }else{
            aaa = 0.0;
            bbb = 0.5;
        }
        
		/* loct->tcdc_clm = 0.5; */  /* for DEBUG */
		loct->tcdc_clm = bbb;
        
        /* sensitivity analysis for cloudiness: 2017/05/30 by A.Ito */
        /* revised 2017/09/14 by A.Ito */
        if(EX_CLOUD == 1 || EX_CLOUD == 2){
            
            if(EX_CLOUD == 1){
                loct->tcdc_clm *= 1.1;
            }else if(EX_CLOUD == 2){
                loct->tcdc_clm *= 0.9;
            }
            if(loct->tcdc_clm >= 1.0){
                loct->tcdc_clm = 1.0;
            }
            if(loct->tcdc_clm <= 0.0){
                loct->tcdc_clm = 0.0;
            }
            
            aaa = loct->toprad_h * (0.8964 - 0.5392 * loct->tcdc_clm);
        
            if(aaa < 0.60 * loct->dswrf_sfc){
                loct->dswrf_sfc = 0.6 * loct->dswrf_sfc;
            }
            if(aaa > 1.4 * loct->dswrf_sfc){
                loct->dswrf_sfc = 1.4 * loct->dswrf_sfc;
            }
        }
	}
    
    if(loct->hour == 0){
        loct->tair_dayav = loct->tair_dayav_c;
        loct->tair_dayav_c = 0.0;
    }
    loct->tair_dayav_c += loct->tmp_2m / (float)DSTEP;
	
	/*  snow-packing effect on soil temperature  */
	if(mass->snwa > 0.2){
		aaa = mass->snwa/(3.0 + mass->snwa);
		loct->tmp10_soil = 0.0*aaa + loct->tmp10_soil*(1.0 - aaa);
		loct->tmp200_soil = 2.0*aaa + loct->tmp200_soil*(1.0 - aaa);
	}
	
	/* atmosphere *********************************************************/
	/* CO2 **************************************/
    atmco2_trend(grid, loct);
	co2_in_canopy(grid, loct, mass, flux);
	
    /* hour angle: local time */
	loct->hangle = -180.0 + (float)(loct->hour)*15.0 + (float)(grid->lon);
    if(loct->hangle > 180.0){
        loct->hangle -= 360.0;
    }
    if(loct->hangle < -180.0){
        loct->hangle += 360.0;
    }
    
	/* air conditions  **********************************/
	/*  altitude  */
	alt = (grid->topo >= 0.0)?grid->topo:0.0;  
	/*  air pressure */
	loct->air_prsr = 1013.25 * exp(-1.0*(28.964 * 0.001) * 9.8 * alt / (8.3144*(loct->tmp_2m + ZAT)));  
	
	/*  saturated vapour pressure, hPa  */
	loct->vps = f_vap_pre_sat(grid, loct);    
	/*  slope of saturated vapour pressure as a function of temperature  */
	loct->slope_vps = f_slope_vps(grid, loct);	
	 /*  aerodynamic resistance  */    
	loct->r_aero = f_r_aero(grid, loct);	   
	
	/* Future change in specific humidity by CCCma or ECHAM */
	/* deleted 2009/05/25 by A.Ito **/
	
	/*  vapour pressure, hPa  */
	// loct->vp = loct->air_prsr*loct->spfh_2m/(0.622+0.378*loct->spfh_2m);    
	/*  vapour pressure deficit, hPa  */
	// loct->vpd = vpd_d[loct->climy-BYR][loct->doy];
	// loct->vpd = (loct->vpd>=0.0)?loct->vpd:0.0;     /*  */ 
	
	/* Future change in specific humidity by MIROC-high */
	/* deleted 2009/05/25 by A.Ito **/
	
	/*******************************************************/
    /* data: specific humidity */
    if(grid->humd_region < 0.0){
        grid->humd_region = 0.0;
    }
    if(strcmp(grid[0].site_id, "GLOBAL")==0){
        loct->vp = loct->air_prsr * grid->humd_region/(0.622 + 0.378 * grid->humd_region);
    }else if(strcmp(grid[0].site_id, "BAMIYAN")==0){
        loct->vp = loct->slope_vps * grid->humd_region;
    }else if(strcmp(grid[0].site_id, "JAPAN")==0 || strcmp(grid[0].site_id, "BB")==0
            || strcmp(grid[0].site_id, "JAPANc")==0|| strcmp(grid[0].site_id, "JAPANh")==0
            || strcmp(grid[0].site_id, "JAPANk")==0){
        loct->vp = grid->humd_region;
    }
    loct->vpd = loct->vps - loct->vp;
    loct->vpd = (loct->vpd >= 0.0)?loct->vpd:0.0;
	
	/*  air density, kg m-3 */
	loct->air_dns = f_airdens(grid, loct);    
	
	/* latent heat MJ kg-1 */
	loct->lhvp = 1000000.0 * (2.501 - 0.012/5.0*loct->tmp_2m);
	
	/* N deposition */
	/* f_n_deposit(grid, loct); */
	/* experiment for AGU2008F */
	if(loct->CO2y < 1980){
		loct->depo_no3 = loct->depo_no3_model_av;
		loct->depo_nh4 = loct->depo_nh4_model_av;
	}else if(loct->CO2y >= 1980 && loct->CO2y <= 2003){
		loct->depo_no3 = loct->depo_no3_model;
		loct->depo_nh4 = loct->depo_nh4_model;
	}else if(loct->CO2y > 2003){
		loct->depo_no3 = loct->depo_no3_model_av;
		loct->depo_nh4 = loct->depo_nh4_model_av;
	}

	/* loct->depo_no3 = loct->depo_no3_model[0][loct->doy];
	loct->depo_nh4 = loct->depo_nh4_model[0][loct->doy]; */
	
	/* loct->depo_no3 = loct->depo_no3_model_av[loct->doy];
	loct->depo_nh4 = loct->depo_nh4_model_av[loct->doy]; */
	
	/* N fertilizer */
	loct->nfert_no3 = 0.0;
	loct->nfert_nh4 = 0.0;

	/* radiation *********************************************************************/
	/* canopy-top PPFD */
	f_ppfd(grid, loct); 

	/* soil surface albedo **********/
	loct->snow_acc = mass->snwa;
	(echar->soil).albedo = albedo_soil(loct, &(echar->soil));

	/* give irradiance attenuation coefficient  */
	(echar->tree).eK = irr_attn(grid, loct, &(echar->tree));
	(echar->c3).eK = irr_attn(grid, loct, &(echar->c3));
	(echar->c4).eK = irr_attn(grid, loct, &(echar->c4));
	
	/*  net radiation  ***************/	
	f_net_rad(grid, loct, echar, mass);
	
	/* CURRENT soil moisture based on the previous soil water content:  2009/05/01 by A.Ito **/
	loct->snow_acc = mass->snwa;
	loct->soilwtr_l = mass->sw30;
	loct->soilwtr_h = mass->sww;
	loct->soilappr_l = (grid->fieldcap30 - mass->sw30)/grid->fieldcap30;
	loct->soilappr_l = (loct->soilappr_l < 1.0)?loct->soilappr_l:1.0;
	loct->soilappr_l = (loct->soilappr_l > 0.0)?loct->soilappr_l:0.0;
	loct->soilappr_w = (grid->fieldcap - mass->sww)/grid->fieldcap;
	loct->soilappr_w = (loct->soilappr_w < 1.0)?loct->soilappr_w:1.0;
	loct->soilappr_w = (loct->soilappr_w > 0.0)?loct->soilappr_w:0.0;

	/* ecophysiology: changing ***************************************************/
	f_ecophysiology(grid, loct, &(echar->tree), &(mass->tree));
	f_ecophysiology(grid, loct, &(echar->c3), &(mass->c3));
	f_ecophysiology(grid, loct, &(echar->c4), &(mass->c4));
	
	/**  hydrological water budget ************************************************/
	f_hydrology(grid, loct, echar, mass);
	
	/** UPDATED soil moisture : 2009/05/01 by A.Ito **/
	loct->snow_acc = mass->snwa;
	loct->soilwtr_l = mass->sw30;
	loct->soilwtr_h = mass->sww;
	loct->soilappr_l = (grid->fieldcap30 - mass->sw30)/grid->fieldcap30;
	loct->soilappr_l = (loct->soilappr_l < 1.0)?loct->soilappr_l:1.0;
	loct->soilappr_l = (loct->soilappr_l > 0.0)?loct->soilappr_l:0.0;
	loct->soilappr_w = (grid->fieldcap - mass->sww)/grid->fieldcap;
	loct->soilappr_w = (loct->soilappr_w < 1.0)?loct->soilappr_w:1.0;
	loct->soilappr_w = (loct->soilappr_w > 0.0)?loct->soilappr_w:0.0;
	
	loct->wfps = ((((mass->sw30/10.0)/30.0)*100.0) / ((1.0 - grid->bulkdens/2.65)*100.0));  /**/
	/* loct->wfps = ((((loct->sw30/10.0)/30.0)*100.0) / ((1.0 - grid->bulkdens/1.325)*100.0)); */
	
	/* to avoid invalid calculation */
	if(loct->wfps < 0.05){
		loct->wfps = 0.05;
	} 
	
	/* water potential: 2011/10/25 by A.Ito **********************/
	/* litter ***********/
	loct->pot_grav_l = -0.05;
	if(mass->sw30>=0.2){
		if(grid->stexture == 0){
			loct->pot_matric_l = -0.121 * pow(mass->sw30/grid->fieldcap30,-4.05);
		}else if(grid->stexture == 1){
			loct->pot_matric_l = -0.478 * pow(mass->sw30/grid->fieldcap30,-5.39);
		}else if(grid->stexture == 2){
			loct->pot_matric_l = -0.405 * pow(mass->sw30/grid->fieldcap30,-11.4);
		}
	}else{
		if(grid->stexture == 0){
			loct->pot_matric_l = -0.121 * pow(0.2/grid->fieldcap30,-4.05);
		}else if(grid->stexture == 1){
			loct->pot_matric_l = -0.478 * pow(0.2/grid->fieldcap30,-5.39);
		}else if(grid->stexture == 2){
			loct->pot_matric_l = -0.405 * pow(0.2/grid->fieldcap30,-11.4);
		}
	}
	loct->pot_total_l = loct->pot_grav_l + loct->pot_matric_l;
	/* humus ***********/
	loct->pot_grav_h = -1.00;
	if(mass->sww>=0.2){
		if(grid->stexture == 0){
			loct->pot_matric_h = -0.121 * pow(mass->sww/grid->fieldcap,-4.05);
		}else if(grid->stexture == 1){
			loct->pot_matric_h = -0.478 * pow(mass->sww/grid->fieldcap,-5.39);
		}else if(grid->stexture == 2){
			loct->pot_matric_h = -0.405 * pow(mass->sww/grid->fieldcap,-11.4);
		}
	}else{
		if(grid->stexture == 0){
			loct->pot_matric_h = -0.121 * pow(0.2/grid->fieldcap,-4.05);
		}else if(grid->stexture == 1){
			loct->pot_matric_h = -0.478 * pow(0.2/grid->fieldcap,-5.39);
		}else if(grid->stexture == 2){
			loct->pot_matric_h = -0.405 * pow(0.2/grid->fieldcap,-11.4);
		}
	}
	loct->pot_total_h = loct->pot_grav_h + loct->pot_matric_h;

	/* deleted 2009/05/25 by A.Ito **/
}

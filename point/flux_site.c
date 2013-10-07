/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/*
Ito, A., Inatomi, M., Mo, W., Lee, M., Koizumi, H., Saigusa, N., Murayama, S., and Yamamoto, S.: 
  Examination of model-estimated ecosystem respiration by use of flux measurement data from a 
  cool-temperate deciduous broad-leaved forest in central Japan, Tellus, 59B, 616–624, 2007.
*/

/* NEE flux simulation using the 30-min scheme *****************************************/
void f_flux_site(
	struct  Grid *grid, 
	struct  Loct  *loct,  
	struct  Echar *echar, 
	struct  Mass  *mass,  
	struct  Flux  *flux, 
	FILE *fp_fxd, 
	FILE *fp_r[NFILE]
){
	long h;
	long yr, day,need_flag, needt_flag, flux_byr, flux_eyr;
	double hr, assim_h, ecoresp_h, apar_h, scale;
	double vp, vps, vpd;
	double fgppd, ferd, fneed, fneed_obs;
	double fgppdt, ferdt, fneedt, fneedt_obs, fapardt;
	double obs_neeh, obs_uw, ustar;
	double obs_atmp, obs_vpd, obs_wind;
	double obs_stmp[6], obs_swtr[4];
	double obs_par_in;
	double obs_rhmd;
	double fatmp, fppfd;
	
	double apar_bn, apar_qn, apar_un;
	double apar_bd, apar_qd, apar_ud;
	double gpp_bn, gpp_qn, gpp_un;
	double gpp_bd, gpp_qd, gpp_ud;
	double pc_bn, pc_qn, pc_un;
	double pc_bd, pc_qd, pc_ud;
	double pj_bn, pj_qn, pj_un;
	double pj_bd, pj_qd, pj_ud;
	double rd_bn, rd_qn, rd_un;
	double rd_bd, rd_qd, rd_ud;
	
	double ipar_d, apar_d, gpp_d, tmp_d, er_d;
	double vcmax_av, jmax_av, Iabs_av, j_av, ci_av, gs_av, ac_av, aj_av, rd_av;
	double limit_ac[4], limit_aj[4];
	
	scale = 1800.0/1000000.0;
	/* 1800.0 seconds: 30 minute */
	/* 1000000.0 micro mol: mol */
    
    flux_byr = -9999;
    flux_eyr = -9999;
    
    if(strcmp(grid->site_id, "TKY")==0){
        flux_byr = 1998;
        flux_eyr = 2005;
    }
	
    if(strcmp(grid->site_id, "CEAMIP_TSE")==0){
        flux_byr = 2002;
        flux_eyr = 2005;
    }
	
	if(loct->climy>=flux_byr && loct->climy<=flux_eyr){
		fgppd=ferd=fneed=fneed_obs=0.0;
		fgppdt=ferdt=fneedt=fapardt=fneedt_obs=0.0;
		need_flag=needt_flag=0;
		
		apar_bn = apar_qn = apar_un = 0.0;
		apar_bd = apar_qd = apar_ud = 0.0;
		gpp_bn = gpp_qn = gpp_un = 0.0;
		gpp_bd = gpp_qd = gpp_ud = 0.0;
		pc_bn = pc_qn = pc_un = 0.0;
		pc_bd = pc_qd = pc_ud = 0.0;
		pj_bn = pj_qn = pj_un = 0.0;
		pj_bd = pj_qd = pj_ud = 0.0;
		rd_bn = rd_qn = rd_un = 0.0;
		rd_bd = rd_qd = rd_ud = 0.0;
		
		ipar_d = apar_d = gpp_d = tmp_d = er_d = 0.0;
		vcmax_av = jmax_av = Iabs_av = j_av = ci_av = gs_av = ac_av = aj_av = rd_av = 0.0;
		limit_ac[0] = limit_aj[0] = 0.0;
		limit_ac[1] = limit_aj[1] = 0.0;
		limit_ac[2] = limit_aj[2] = 0.0;
		limit_ac[3] = limit_aj[3] = 0.0;
		
		/*  seasonal change in leaf properties  ********************************/
		if(strcmp(grid->site_id, "TKY")==0 || strcmp(grid->site_id, "CEAMIP_TSE")==0){
			f_leaf_aging(grid, loct, echar, mass);
		}
		if(AG_VMX==1){
			(echar->overs1).photocap_n = (echar->overs1).photocap_n0* (echar->overs1).vcmx_b;
			(echar->overs2).photocap_n = (echar->overs2).photocap_n0* (echar->overs2).vcmx_b;
		}else{
			(echar->overs1).photocap_n = (echar->overs1).photocap_n0;
			(echar->overs2).photocap_n = (echar->overs2).photocap_n0;
		}
		(echar->unders).photocap_n = (echar->unders).photocap_n0;
		
		/* 30-min step *************************************************************/
		for(h=0;h<DSTEP;h++){
			loct->hour = h;
			loct->hangle = -180.0+(double)h*7.5;
			loct->metdat_frag = 0;
			
			/*** environmental condition ************************************/
            if(strcmp(grid->site_id, "TKY")==0){
                fscanf(fp_fxd,"%ld %ld %lf", &yr, &day, &hr);
                fscanf(fp_fxd,"%lf %lf", &obs_atmp, &obs_rhmd);
                fscanf(fp_fxd,"%lf %lf %lf %lf %lf", &obs_stmp[0], &obs_stmp[1], 
                        &obs_stmp[2], &obs_stmp[3], &obs_stmp[4]);
                fscanf(fp_fxd,"%lf %lf %lf %lf", &obs_swtr[0], &obs_swtr[1],
                        &obs_swtr[2], &obs_swtr[3]);
                fscanf(fp_fxd,"%lf %lf %lf", &obs_par_in, &obs_uw, &obs_neeh);
                                
                if(yr!=loct->climy){
                    printf("TKY: WRONG HOURLY DATA <%ld - %ld:%ld>\n", loct->climy, yr, day);
                }
                
                /* u* (u-star): friction velocity related to atmospheric stability, m s-1 */
                if(obs_uw<0.0 && obs_uw>-50.0){
                    ustar = sqrt(-1.0*obs_uw);
                }else{
                    ustar = -999.999;
                }
			
            }
            
            if(strcmp(grid->site_id, "CEAMIP_TSE")==0){
                fscanf(fp_fxd,"%lf", &obs_par_in);
                fscanf(fp_fxd,"%lf", &obs_atmp);
                fscanf(fp_fxd,"%lf", &obs_vpd);
                fscanf(fp_fxd,"%lf", &obs_wind);
                fscanf(fp_fxd,"%lf", &obs_stmp[0]);
                fscanf(fp_fxd,"%lf", &ustar);
                fscanf(fp_fxd,"%lf", &obs_neeh);
                
                /* printf("%lf\n", obs_par_in); */
                
                obs_stmp[1] = obs_stmp[2] = obs_stmp[3] = obs_stmp[0];
            }
			
			/* temperature data ***************************************************************/
			if(obs_atmp>-30.0 && obs_atmp<50.0){
				/*  normal: tower-observed  */
				(echar->overs1).tmp = (echar->overs2).tmp = obs_atmp;
				(echar->unders).tmp = obs_atmp;
				loct->metdat_frag += 1;
				fatmp = obs_atmp;
			}else{
				/*  abnormal: NCEP/NCAR  */
				(echar->overs1).tmp = (echar->overs2).tmp = loct->tmp_2m;
				(echar->unders).tmp = loct->tmp_2m;		
				fatmp = loct->tmp_2m;
			}
			if(((obs_stmp[0]+obs_stmp[1])*0.5)>-30.0 && ((obs_stmp[0]+obs_stmp[1])*0.5)<50.0){
				(echar->soil).tmp = (obs_stmp[0]+obs_stmp[1])*0.5;
			}else{
				(echar->soil).tmp = loct->tmp10_soil;
			}
			
			/* vapor data ***************************************************************/
            if(strcmp(grid->site_id, "TKY")==0){
                if((fatmp>-30.0&&fatmp<50.0) && (obs_rhmd>=0.0&&obs_rhmd<120.0)){
                    if(fatmp>0.0){  /*  at water surface  */
                        vps = 6.1078*pow(10.0, (7.5*fatmp)/(237.3+fatmp));
                    }else if(fatmp<=0.0){  /*  at ice surface  */
                        vps = 6.1078*pow(10.0, (9.5*fatmp)/(265.3+fatmp));
                    }
                    if(obs_rhmd>=100.0) obs_rhmd = 100.0;	
                    vp = vps * obs_rhmd/100.0;
                    vpd = (vps > vp)?(vps-vp):0.0;
                    (echar->overs1).vpd_a = (echar->overs2).vpd_a = vpd;
                    (echar->unders).vpd_a = vpd;
                    loct->metdat_frag += 2;
                }else{
                    (echar->overs1).vpd_a = (echar->overs2).vpd_a = loct->vpd;
                    (echar->unders).vpd_a = loct->vpd;				
                }
            }
            if(strcmp(grid->site_id, "CEAMIP_TSE")==0){
                loct->metdat_frag += 2;
                (echar->overs1).vpd_a = (echar->overs2).vpd_a = obs_vpd/10.0;
                (echar->unders).vpd_a = obs_vpd/10.0;				
            }
			
			/* par data ***************************************************************/
			if(obs_par_in>-50.0 && obs_par_in<3000.0){
				obs_par_in = (obs_par_in>0.0)?obs_par_in:0.0;
				
				if(loct->ppfd_h[h]>0.0 && obs_par_in>=0.0){
					(echar->overs1).ppfdb_top = (echar->overs2).ppfdb_top = 
												obs_par_in * loct->ppfdb_h[h]/loct->ppfd_h[h];
					(echar->overs1).ppfdd_top = (echar->overs2).ppfdd_top = 
												obs_par_in * loct->ppfdd_h[h]/loct->ppfd_h[h];
				}else{
					(echar->overs1).ppfdb_top = (echar->overs2).ppfdb_top = 0.0;
					(echar->overs1).ppfdd_top = (echar->overs2).ppfdd_top = 0.0;
				}
				loct->metdat_frag += 4;
				fppfd = obs_par_in;
			}else{
				(echar->overs1).ppfdb_top = (echar->overs2).ppfdb_top = loct->ppfdb_h[h];
				(echar->overs1).ppfdd_top = (echar->overs2).ppfdd_top = loct->ppfdd_h[h];
				fppfd = loct->ppfdb_h[h] + loct->ppfdd_h[h];
			}
						
			/* leaf biochemistry **********************************************/
			/* initial Ci */
			(echar->overs1).ci_sn = (echar->tree).ci/10.0;
			(echar->overs1).ci_sd = (echar->tree).ci/10.0;
			(echar->overs2).ci_sn = (echar->tree).ci/10.0;
			(echar->overs2).ci_sd = (echar->tree).ci/10.0;
			(echar->unders).ci_sn = (echar->c3).ci/10.0;
			(echar->unders).ci_sd = (echar->c3).ci/10.0;
			
			/* single-leaf properties */
			f_leaf_prop(grid, loct, &(echar->overs1));
			f_leaf_prop(grid, loct, &(echar->overs2));
			f_leaf_prop(grid, loct, &(echar->unders));
						
			/* sun/shade PPFD absorption */
			f_apar_sunshade(grid, loct, &(echar->overs1));
			f_apar_sunshade(grid, loct, &(echar->overs2));
			(echar->unders).ppfdb_top = (echar->overs1).ppfdb_btm*loct->comp_over1 
										+ (echar->overs2).ppfdb_btm*loct->comp_over2;
			(echar->unders).ppfdd_top = (echar->overs1).ppfdd_btm*loct->comp_over1 
										+ (echar->overs2).ppfdd_btm*loct->comp_over2;
			f_apar_sunshade(grid, loct, &(echar->unders));
			
			/* canopy CO2 assimilation */
			f_gpp_sunshade((echar->tree).season, grid, loct, &(echar->overs1));
			f_gpp_sunshade((echar->tree).season, grid, loct, &(echar->overs2));
			if(mass->snwa<20.0){
				/* snow-free */
				f_gpp_sunshade((echar->c3).season, grid, loct, &(echar->unders));
				//f_gpp_sunshade((echar->c3).season, grid, loct, &(echar->unders));
			}else{
				/* snow-covered */
				(echar->unders).p_sn = 0.0;
				(echar->unders).p_sd = 0.0;
			}
			
			/**  half-hourly respiration  **********************************/
			f_resp_h(grid, loct, echar, mass, flux);
			
			/****************************************************************************************/
			apar_h = ((echar->overs1).apar_sn + (echar->overs1).apar_sd)*loct->comp_over1 
					+ ((echar->overs2).apar_sn + (echar->overs2).apar_sd)*loct->comp_over2 
					+ (echar->unders).apar_sn + (echar->unders).apar_sd;
			
			assim_h = ((echar->overs1).p_sn + (echar->overs1).p_sd)*loct->comp_over1 
					+ ((echar->overs2).p_sn + (echar->overs2).p_sd)*loct->comp_over2 
					+ (echar->unders).p_sn + (echar->unders).p_sd;
						
			ecoresp_h = (flux->tree).ar_h + (flux->c3).ar_h + (flux->soil).hr_h; 
			/*  ecoresp_h = (flux->soil).hr_h;  */
			
			if(ecoresp_h>=0.0&&ecoresp_h<=100.0){
				;
			}else{
				printf("%lf %lf %lf\n", (flux->tree).ar_h, (flux->c3).ar_h, (flux->soil).hr_h);
			}
			
			ipar_d += ((echar->overs1).ppfdb_top + (echar->overs1).ppfdd_top)*scale;
			//ipar_d += loct->ppfdb_h[h]*scale;
			
			apar_d += (((echar->overs1).appfd_sn + (echar->overs1).appfd_sd)*loct->comp_over1 + 
					((echar->overs2).appfd_sn + (echar->overs2).appfd_sd)*loct->comp_over2 + 
					((echar->unders).appfd_sn + (echar->unders).appfd_sd))*scale;			
			gpp_d += assim_h*scale;
			tmp_d += (echar->overs1).tmp/(double)DSTEP;
			er_d += ecoresp_h*scale;
			if(h>=18 && h<=23){
				vcmax_av += ((echar->overs1).vcmx_sn*loct->comp_over1 + (echar->overs2).vcmx_sn*loct->comp_over2 + 
						(echar->overs2).vcmx_sn*loct->comp_over1 + (echar->overs2).vcmx_sd*loct->comp_over2)/6.0;
				jmax_av += ((echar->overs1).jmx_sn*loct->comp_over1 + (echar->overs2).jmx_sd*loct->comp_over2 + 
						(echar->overs2).jmx_sn*loct->comp_over1 + (echar->overs2).jmx_sd*loct->comp_over2)/6.0;
				Iabs_av += (((echar->overs1).apar_sn + (echar->overs1).apar_sd)*loct->comp_over1 + 
					 	((echar->overs2).apar_sn + (echar->overs2).apar_sd)*loct->comp_over2)/6.0;
				j_av += (((echar->overs1).j_sn + (echar->overs1).j_sd)*loct->comp_over1 + 
					 	((echar->overs2).j_sn + (echar->overs2).j_sd)*loct->comp_over2)/6.0;
				ac_av += (((echar->overs1).pc_sn + (echar->overs1).pc_sd)*loct->comp_over1 + 
					 	((echar->overs2).pc_sn + (echar->overs2).pc_sd)*loct->comp_over2)/6.0;
				aj_av += (((echar->overs1).pj_sn + (echar->overs1).pj_sd)*loct->comp_over1 + 
					 	((echar->overs2).pj_sn + (echar->overs2).pj_sd)*loct->comp_over2)/6.0;
				rd_av += (((echar->overs1).rd_sn + (echar->overs1).rd_sd)*loct->comp_over1 + 
					 	((echar->overs2).rd_sn + (echar->overs2).rd_sd)*loct->comp_over2)/6.0;
				if((echar->overs1).pc_sn < (echar->overs1).pj_sn){
					limit_ac[0] += 1.0;
				}else{
					limit_aj[0] += 1.0;
				}
				if((echar->overs1).pc_sd < (echar->overs1).pj_sd){
					limit_ac[1] += 1.0;
				}else{
					limit_aj[1] += 1.0;
				}
				if((echar->overs2).pc_sn < (echar->overs2).pj_sn){
					limit_ac[2] += 1.0;
				}else{
					limit_aj[2] += 1.0;
				}
				if((echar->overs2).pc_sd < (echar->overs2).pj_sd){
					limit_ac[3] += 1.0;
				}else{
					limit_aj[3] += 1.0;
				}
			}
			
			/************************/
			/*  daily  */
			if(loct->metdat_frag==7 && obs_neeh>-100.0 && obs_neeh<100.0){ /*  g C m-2 day-1  */
				need_flag++;
			/* 	fgppd += assim_h*scale*12.0;
				ferd += ecoresp_h*scale*12.0;
				fneed += ecoresp_h*scale*12.0 - assim_h*scale*12.0;  */
				fneed_obs += obs_neeh*scale*12.0;
			}
			
			fgppd += assim_h*scale*12.0;
			ferd += ecoresp_h*scale*12.0;
			fneed += ecoresp_h*scale*12.0 - assim_h*scale*12.0;  /*  */
						
			/*  daytime:  g C m-2 day-1  */
			if(fppfd>50.0 && 
			   loct->metdat_frag==7 && 
			   ustar>0.3 && 
			   obs_neeh>-100.0&&obs_neeh<100.0){ 
				/*  &&ustar>0.3   &&obs_neeh>-100.0&&obs_neeh<100.0 */
				needt_flag++;  
				fgppdt += assim_h*scale*12.0;
				ferdt += ecoresp_h*scale*12.0;
				fneedt += ecoresp_h*scale*12.0 - assim_h*scale*12.0;
				fneedt_obs += obs_neeh*scale*12.0;
				
				fapardt += ((echar->overs1).apar_sn*loct->comp_over1+(echar->overs2).apar_sn*loct->comp_over2)*scale;
				fapardt += (echar->unders).apar_sn*scale;
				fapardt += ((echar->overs1).apar_sd*loct->comp_over1+(echar->overs2).apar_sd*loct->comp_over2)*scale;
				fapardt += (echar->unders).apar_sd*scale;
			}
			
			/*  daytime  canopy properties  */
		/*	if(ippfd>50.0){  */
			if(fppfd>=0.0){
				apar_bn += (echar->overs1).apar_sn*loct->comp_over1*scale;
				apar_qn += (echar->overs2).apar_sn*loct->comp_over2*scale;
				apar_un += (echar->unders).apar_sn*scale;
				
				apar_bd += (echar->overs1).apar_sd*loct->comp_over1*scale;
				apar_qd += (echar->overs2).apar_sd*loct->comp_over2*scale;
				apar_ud += (echar->unders).apar_sd*scale;
				
				gpp_bn += (echar->overs1).p_sn*loct->comp_over1*scale*12.0;
				gpp_qn += (echar->overs2).p_sn*loct->comp_over2*scale*12.0;
				gpp_un += (echar->unders).p_sn*scale*12.0;
				
				gpp_bd += (echar->overs1).p_sd*loct->comp_over1*scale*12.0;
				gpp_qd += (echar->overs2).p_sd*loct->comp_over2*scale*12.0;
				gpp_ud += (echar->unders).p_sd*scale*12.0;
				
				pc_bn += (echar->overs1).pc_sn*loct->comp_over1*scale;
				pc_qn += (echar->overs2).pc_sn*loct->comp_over2*scale;
				pc_un += (echar->unders).pc_sn*scale;
				
				pc_bd += (echar->overs1).pc_sd*loct->comp_over1*scale;
				pc_qd += (echar->overs2).pc_sd*loct->comp_over2*scale;
				pc_ud += (echar->unders).pc_sd*scale;

				pj_bn += (echar->overs1).pj_sn*loct->comp_over1*scale;
				pj_qn += (echar->overs2).pj_sn*loct->comp_over2*scale;
				pj_un += (echar->unders).pj_sn*scale;
				
				pj_bd += (echar->overs1).pj_sd*loct->comp_over1*scale;
				pj_qd += (echar->overs2).pj_sd*loct->comp_over2*scale;
				pj_ud += (echar->unders).pj_sd*scale;

				rd_bn += (echar->overs1).rd*loct->comp_over1*scale;
				rd_qn += (echar->overs2).rd*loct->comp_over2*scale;
				rd_un += (echar->unders).rd*scale;

				rd_bd += (echar->overs1).rd*loct->comp_over1*scale;
				rd_qd += (echar->overs2).rd*loct->comp_over2*scale;
				rd_ud += (echar->unders).rd*scale;
			}
			
			/**  hourly property  *****************************/
			if(loct->month==5 && loct->mday<10){
				/* fprintf(fp_r[loct->climy-1998+9],"%ld %ld %ld %ld %ld %ld ",
					loct->climy, loct->doy, loct->month+1, loct->mday+1, h, loct->metdat_frag);  
				
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf %lf %lf ",(echar->overs1).tmp, 
					(echar->overs1).ppfdb_top, (echar->overs1).ppfdd_top, (echar->overs1).vpd_a);
				
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).lai_sn,(echar->overs1).lai_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).lai_sn,(echar->overs2).lai_sd);
				
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).vcmx_sn,(echar->overs1).vcmx_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).vcmx_sn,(echar->overs2).vcmx_sd);
				
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).apar_sn,(echar->overs1).apar_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).apar_sn,(echar->overs2).apar_sd); */
				
				
				/*	fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).jmx_sn,(echar->overs1).jmx_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).jmx_sn,(echar->overs2).jmx_sd);
				
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).ci_sn,(echar->overs1).ci_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).ci_sn,(echar->overs2).ci_sd);
				
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).gs_sn,(echar->overs1).gs_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).gs_sn,(echar->overs2).gs_sd);
				
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).j_sn,(echar->overs1).j_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).j_sn,(echar->overs2).j_sd);
				
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).pc_sn,(echar->overs1).pc_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).pc_sn,(echar->overs2).pc_sd);

				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).pj_sn,(echar->overs1).pj_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).pj_sn,(echar->overs2).pj_sd);

				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).p_sn,(echar->overs1).p_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).p_sn,(echar->overs2).p_sd);

				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).rd_sn,(echar->overs1).rd_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).rd_sn,(echar->overs2).rd_sd);

				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs1).a_sn,(echar->overs1).a_sd);
				fprintf(fp_r[loct->climy-1998+9],"%lf %lf ",(echar->overs2).a_sn,(echar->overs2).a_sd); */

				/*  fprintf(fp_r[loct->climy-1998+9],"\n"); */
			}		
			
		/*	fprintf(fp_r[loct->climy-1998+9],"%4ld %3ld %2ld %2ld ", 
					loct->climy, loct->doy, h, loct->metdat_frag);
			fprintf(fp_r[loct->climy-1998+9],"%12.5lf ", ecoresp_h*scale*12.0);
			fprintf(fp_r[loct->climy-1998+9],"\n"); */
			
		/*	fprintf(fp_r[grid->year-1998+3],"%lf %lf %lf ",(echar->overs1).a_sd,(echar->overs1).p_sd, (echar->overs1).rd_sd);
			fprintf(fp_r[grid->year-1998+3],"%lf %lf ",(echar->overs1).pc_sd,(echar->overs1).pj_sd);
			fprintf(fp_r[grid->year-1998+3],"%lf %lf ",(echar->overs1).gs_sd,(echar->overs1).gb);
			fprintf(fp_r[grid->year-1998+3],"%lf %lf ",(echar->overs1).co2_a,(echar->overs1).ci_sd);
			
			fprintf(fp_r[grid->year-1998+3],"%lf %lf ",(echar->overs1).vcmx_sd,(echar->overs1).comp_co2);
			fprintf(fp_r[grid->year-1998+3],"%lf %lf ",(echar->overs1).k_effc,(echar->overs1).j_sd);   */
			
		/*	fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf %.2lf ",(echar->overs1).a_sn,(echar->overs1).p_sn, (echar->overs1).rd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs1).pc_sn,(echar->overs1).pj_sn);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs1).gs_sn,(echar->overs1).gb);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs1).co2_a,(echar->overs1).ci_sn);  */
			
		/*	fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs1).lai_sn, (echar->overs1).lai_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs2).lai_sn, (echar->overs2).lai_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->unders).lai_sn, (echar->unders).lai_sd);

			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs1).apar_sn, (echar->overs1).apar_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs2).apar_sn, (echar->overs2).apar_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->unders).apar_sn, (echar->unders).apar_sd);
			
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs1).j_sn,(echar->overs1).j_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs2).j_sn,(echar->overs2).j_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->unders).j_sn,(echar->unders).j_sd);
			
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf %.2lf ",(echar->overs1).pc_sn,(echar->overs1).pj_sn, (echar->overs1).p_sn);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf %.2lf ",(echar->overs1).pc_sd,(echar->overs1).pj_sd, (echar->overs1).p_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf %.2lf ",(echar->overs2).pc_sn,(echar->overs2).pj_sn, (echar->overs2).p_sn);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf %.2lf ",(echar->overs2).pc_sd,(echar->overs2).pj_sd, (echar->overs2).p_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf %.2lf ",(echar->unders).pc_sn,(echar->unders).pj_sn, (echar->unders).p_sn);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf %.2lf ",(echar->unders).pc_sd,(echar->unders).pj_sd, (echar->unders).p_sd);

			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs1).gs_sn,(echar->overs1).gs_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs2).gs_sn,(echar->overs2).gs_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->unders).gs_sn,(echar->unders).gs_sd);
			
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs1).ci_sn,(echar->overs1).ci_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->overs2).ci_sn,(echar->overs2).ci_sd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf %.2lf ",(echar->unders).ci_sn,(echar->unders).ci_sd);
			
			fprintf(fp_r[grid->year-1998+3],"%.2lf ",(echar->overs1).rd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf ",(echar->overs2).rd);
			fprintf(fp_r[grid->year-1998+3],"%.2lf ",(echar->unders).rd);  */

			/**  hourly NEE  *****************************/
			fprintf(fp_r[3],"%4ld ", loct->climy);
			fprintf(fp_r[3],"%4ld ", loct->doy);
			fprintf(fp_r[3],"%4ld ", loct->hour);
			fprintf(fp_r[3],"%4ld ", loct->metdat_frag);
			
			fprintf(fp_r[3],"%12.4lf ",obs_atmp);
			fprintf(fp_r[3],"%12.4lf ",obs_par_in);
			fprintf(fp_r[3],"%12.4lf ",(echar->overs1).ppfdb_top);
			fprintf(fp_r[3],"%12.4lf ",(echar->overs1).ppfdd_top);
			
            fprintf(fp_r[3],"%12.4lf ", apar_h);
			fprintf(fp_r[3],"%12.4lf ", assim_h);
            fprintf(fp_r[3],"%12.4lf ", ecoresp_h);
            fprintf(fp_r[3],"%12.4lf ", ecoresp_h-assim_h);
            fprintf(fp_r[3],"%12.4lf ", obs_neeh);
            
            
			//fprintf(fp_r[3],"%12.2lf %12.2lf ",(echar->overs1).ci_sn,(echar->overs1).ci_sd);
			//fprintf(fp_r[3],"%12.4lf ", ustar);
			//fprintf(fp_r[3],"%12.4lf ", apar_h);
			//fprintf(fp_r[3],"%12.4lf ",loct->ppfdb_h[h]+loct->ppfdd_h[h]);
			//fprintf(fp_r[3],"%12.4lf ", ecoresp_h*scale*12.0);
			//fprintf(fp_r[3],"%12.4lf ", ecoresp_h*scale*12.0 - assim_h*scale*12.0);
			//fprintf(fp_r[3],"%12.4lf ", obs_neeh*scale*12.0);
			//fprintf(fp_r[3],"%12.4lf ", apar_h*scale);

			//fprintf(fp_r[3],"%12.4lf ", obs_rhmd);
			//fprintf(fp_r[3],"%12.4lf ", vpd);

			/* if(loct->metdat_frag==7&&ustar>0.3&&obs_neeh>-100.0&&obs_neeh<100.0){
				fprintf(fp_r[3],"%12.4lf ", obs_neeh*scale*12.0);
				fprintf(fp_r[3],"%12.4lf ", ecoresp_h*scale*12.0 - assim_h*scale*12.0);
			} */
			fprintf(fp_r[3],"\n");		
		} /**  end of hourly roop  ******************************************************/
		
		/**  daily NEE  *****************************/
		fprintf(fp_r[4],"%4ld %4ld %2ld ", loct->climy, loct->doy, need_flag);
		fprintf(fp_r[4],"%4ld ", (echar->tree).season);
		fprintf(fp_r[4],"%lf ", ipar_d);
		fprintf(fp_r[4],"%lf ", apar_d);
	/* 	fprintf(fp_r[4],"%12.6lf ", (echar->overs1).vcmx_b*(echar->overs1).vcmx_con); 
		fprintf(fp_r[4],"%12.6lf ", (echar->overs1).jvr); 
		fprintf(fp_r[4],"%12.6lf ", (echar->overs1).lma); 
		fprintf(fp_r[4],"%12.6lf ", (echar->overs1).rd0);    */
	/*	fprintf(fp_r[4],"%12.6lf ", (echar->overs2).vcmx_b*(echar->overs1).vcmx_con); 
		fprintf(fp_r[4],"%12.6lf ", (echar->overs2).jvr); 
		fprintf(fp_r[4],"%12.6lf ", (echar->overs2).lma); 
		fprintf(fp_r[4],"%12.6lf ", (echar->overs2).rd0);     */
		/*  continuous 1 */
		fprintf(fp_r[4],"%12.6lf %12.6lf %12.6lf ", fgppd, ferd, fneed);   /*  */
		/*  continuous 2 */
		if(need_flag>=45){  /*  need_flag==48  */
			/*  good day:  output by 30-min model  */
			fprintf(fp_r[4],"%12.6lf %12.6lf %12.6lf ", fgppd, ferd, fneed);   /*  */
		}else{
			/*  bad day:  surrogate by daily-model  */
			fprintf(fp_r[4],"%12.6lf %12.6lf %12.6lf ", 
					((flux->tree).gpp+(flux->c3).gpp)*100.0, 
					((flux->tree).rp+(flux->c3).rp+(flux->soil).hr)*100.0, 
					(((flux->tree).rp+(flux->c3).rp+(flux->soil).hr)-((flux->tree).gpp+(flux->c3).gpp))*100.0);  /* */
		}
		/*  fprintf(fp_r[4],"%12.6lf %12.6lf %12.6lf ", (echar->overs1).lai, (echar->overs2).lai, (echar->unders).lai);   */ 
		/*  only good day  */
		if(need_flag>=45){  /*  need_flag==48  */
			fprintf(fp_r[4],"%12.6lf %12.6lf %12.6lf %12.6lf ", fgppd, ferd, fneed, fneed_obs);   /*  */ 
		}else{
			;
		}
		fprintf(fp_r[4],"\n");   /*  */
		
		/* daily daytime NEE *****************************/
		/* if(needt_flag>=5){  */
		if(needt_flag>=1){
			fprintf(fp_r[5],"%4ld %4ld %4ld %4ld ", loct->climy, loct->month, loct->doy, needt_flag);
			fprintf(fp_r[5],"%lf ", fapardt);
			fprintf(fp_r[5],"%lf %lf ", fgppdt, ferdt);
			fprintf(fp_r[5],"%lf %lf ", fneedt, fneedt_obs);
			fprintf(fp_r[5],"\n");   /*  */
		}else{
			fprintf(fp_r[5],"%4ld %4ld %4ld %4ld ", loct->climy, loct->month, loct->doy, needt_flag);
			fprintf(fp_r[5],"\n");   /*  */
		}
		
		/**  daily canopy  *****************************/
		/* fprintf(fp_r[6],"%4ld %4ld %4ld %4ld ", loct->climy, loct->doy, (echar->tree).day_flush, needt_flag);  
		fprintf(fp_r[6],"%lf %lf %lf ", apar_bn, apar_qn, apar_un);
		fprintf(fp_r[6],"%lf %lf %lf ", apar_bd, apar_qd, apar_ud);
		fprintf(fp_r[6],"%lf %lf %lf ", gpp_bn, gpp_qn, gpp_un);
		fprintf(fp_r[6],"%lf %lf %lf ", gpp_bd, gpp_qd, gpp_ud);
		fprintf(fp_r[6],"%lf %lf %lf ", pc_bn, pc_qn, pc_un);
		fprintf(fp_r[6],"%lf %lf %lf ", pc_bd, pc_qd, pc_ud);
		fprintf(fp_r[6],"%lf %lf %lf ", pj_bn, pj_qn, pj_un);
		fprintf(fp_r[6],"%lf %lf %lf ", pj_bd, pj_qd, pj_ud);
		fprintf(fp_r[6],"%lf %lf %lf ", rd_bn, rd_qn, rd_un);
		fprintf(fp_r[6],"%lf %lf %lf ", rd_bd, rd_qd, rd_ud);
		
		fprintf(fp_r[6],"%lf %lf %lf %lf %lf ", ipar_d, apar_d, gpp_d, tmp_d, er_d);
		fprintf(fp_r[6],"%lf %lf %lf %lf %lf %lf %lf ", vcmax_av, jmax_av, Iabs_av, j_av, ac_av, aj_av, rd_av); */

	/*	fprintf(fp_r[6],"%lf %lf ", (gpp_bn+gpp_bd)*12.0, (gpp_qn+gpp_qd)*12.0);
		fprintf(fp_r[6],"%lf %lf ", (gpp_bn+gpp_bd-rd_bn-rd_bd)*12.0, (gpp_qn+gpp_qd-rd_qn-rd_qd)*12.0);  */
		
	/*	fprintf(fp_r[6],"%lf %lf %lf %lf %lf %lf %lf %lf ", 
				limit_ac[0], limit_aj[0], limit_ac[1], limit_aj[1], limit_ac[2], limit_aj[2], limit_ac[3], limit_aj[3]);  */
		/*  fprintf(fp_r[6],"\n");   */
	}else{
		for(h=0;h<DSTEP;h++){
			loct->hour = h;
			loct->hangle = -180.0+(double)h*7.5;
			loct->metdat_frag = 0;
			
			(echar->overs1).ppfdb_top = (echar->overs2).ppfdb_top = loct->ppfdb_h[h];
			(echar->overs1).ppfdd_top = (echar->overs2).ppfdd_top = loct->ppfdd_h[h];
			(echar->overs1).tmp = (echar->overs2).tmp = loct->tmp_2m;
			(echar->unders).tmp = loct->tmp_2m;
			(echar->overs1).vpd_a = (echar->overs2).vpd_a = loct->vpd;
			(echar->unders).vpd_a = loct->vpd;
		}
	}
}

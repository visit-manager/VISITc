/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* SIMULATION ********************************************************************/
/* area */
void f_experiment(
	struct Grid grid[DROW][DCOL], 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Echar echar_type[NFILE], 
	struct Mass mass[DROW][DCOL], 
	struct Flux *flux, 
	FILE *fp_r[NFILE]
){
	long e, f, g, h, i;
	char filename[100];
	long ndy, end_year, ddummy;
	double ansis_ann[256],fdummy;
	double aet_a, flux_mon[6][12],dmon[12],aaa;
	FILE *fp_fxd, *fp_lp;
	FILE *fp_restart;
	
	if(NOTICE==1){
		printf("Start experimental phase\n");
	}
	
/*	if(WMODE==1){
		for(f=1;f<=16;f++){
			fprintf(fp_r[f],"%s %s\n", echar->para_ver_id, echar->para_date_id);
		}
	}else if(WMODE==2){
		for(f=1;f<=16;f++){
			fprintf(fp_r[f],"%s %s\n", echar_type[0].para_ver_id, echar_type[0].para_date_id);
		}
	} */
	
	if(USE_RESTART == 1){
		strcpy(filename, grid[0][0].site_id);
		strcat(filename, "_restart.txt");
		if( (fp_restart = fopen(filename,"rt")) == NULL){
			printf("! No restart data !\n");
			exit (1);
		}
	}
    
    if(FIX_PHENOLOGY == 1){
        fp_lp = fopen("TKY_131024_phenology.txt","rt");
        
        for(f=0; f<366; f++){
            fscanf(fp_lp,"%ld", &ddummy);
            fscanf(fp_lp,"%ld", &ddummy);
            
            fscanf(fp_lp,"%ld", &fixlp_tree_season[f]);
            fscanf(fp_lp,"%ld", &fixlp_tree_dayflush[f]);
            fscanf(fp_lp,"%ld", &fixlp_tree_dayshed[f]);
            fscanf(fp_lp,"%lf", &fixlp_tree_lai[f]);
            
            fscanf(fp_lp,"%ld", &fixlp_c3_season[f]);
            fscanf(fp_lp,"%ld", &fixlp_c3_dayflush[f]);
            fscanf(fp_lp,"%ld", &fixlp_c3_dayshed[f]);
            fscanf(fp_lp,"%lf", &fixlp_c3_lai[f]);
            
            fscanf(fp_lp,"%ld", &fixlp_c4_season[f]);
            fscanf(fp_lp,"%ld", &fixlp_c4_dayflush[f]);
            fscanf(fp_lp,"%ld", &fixlp_c4_dayshed[f]);
            fscanf(fp_lp,"%lf", &fixlp_c4_lai[f]);
            
            fscanf(fp_lp,"%lf", &fdummy);
            fscanf(fp_lp,"%lf", &fdummy);
            fscanf(fp_lp,"%lf", &fdummy);
        }
    }
	
	/***************************************************************************/
	for(f=0; f<NROW; f++){
		for(g=0; g<NCOL; g++){
			if(WMODE==2){
				echar = &(echar_type[grid[f][g].veg_type]);
			}
			
		/*	for(h=1;h<NFILE;h++){
				fprintf(fp_r[h],"%ld %ld\n", f, g);
			} */
			
			if(USE_RESTART == 1){
				fscanf(fp_restart,"%lf", &(mass[f][g].snwa));
				fscanf(fp_restart,"%lf", &(mass[f][g].sw30));
				fscanf(fp_restart,"%lf", &(mass[f][g].sww));
				fscanf(fp_restart,"%lf", &((mass[f][g].tree).fol));
				fscanf(fp_restart,"%lf", &((mass[f][g].tree).stm));
				fscanf(fp_restart,"%lf", &((mass[f][g].tree).rot));
				fscanf(fp_restart,"%lf", &((mass[f][g].c3).fol));
				fscanf(fp_restart,"%lf", &((mass[f][g].c3).stm));
				fscanf(fp_restart,"%lf", &((mass[f][g].c3).rot));
				fscanf(fp_restart,"%lf", &((mass[f][g].c4).fol));
				fscanf(fp_restart,"%lf", &((mass[f][g].c4).stm));
				fscanf(fp_restart,"%lf", &((mass[f][g].c4).rot));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).ltr_tf));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).ltr_tc));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).ltr_tr));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).ltr_gf));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).ltr_gc));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).ltr_gr));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).msl_a));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).msl_i));
				fscanf(fp_restart,"%lf", &((mass[f][g].soil).msl_p));
			}
			
			if(grid[f][g].flag_datavl == 1){
				/* initialize N cycle */
				/* init_nitrogen(mass, flux); */
				
				/* 30-min microclimate ***************/
				if(FLUX_SCHEME == 1){
                    if(strcmp(grid[f][g].site_id, "TKY")==0){
                        if( (fp_fxd = fopen("tky_30min_9805.dat","rt")) == NULL){
                        //if( (fp_fxd = fopen("tky_30min_9805_ex2.txt","rt")) == NULL){
                            printf("No input data (TKY FLUX) !\n");
                            exit (1);
                        }
                    }
				}
				
				if(PREDICT == 0){
					end_year = EYR;
				}else if(PREDICT >= 1){
					end_year = 2050;
                    
                    if(PREDICT==100){
                        end_year = 2100;
                    }
				}
				
				/* roop for experimental stage ************************************************/
				for(e=BYR; e<=end_year; e++){	
					loct->climy = e;
					if(SENSANS_ACO2 == 1){
						loct->CO2y = BYR;
					}else{
						loct->CO2y = e;
					}
					loct->time++;
					loct->adyear = e;
					loct->age_stand += 1.0;
					
					/* analysis of annual values */
					aet_a = 0.0;
					for(i=0;i<N_ANSIS;i++){
						ansis_ann[i] = 0.0;
					}
                    ansis_ann[114] = 50.0;
                    
                    if(FLUX_SCHEME == 1){
                        if(strcmp(grid[f][g].site_id, "CEAMIP_TSE")==0){
                            if(e==2002){
                                if((fp_fxd = fopen("tse_30min_2002b.txt","rt"))==NULL){
                                    printf("No tse_30min_2002b.txt.txt\n");
                                    exit(1);
                                }
                                
                                /* fscanf(fp_fxd,"%lf",&aaa); printf("%lf\n",aaa);
                                fscanf(fp_fxd,"%lf",&aaa); printf("%lf\n",aaa);
                                fscanf(fp_fxd,"%lf",&aaa); printf("%lf\n",aaa);
                                fscanf(fp_fxd,"%lf",&aaa); printf("%lf\n",aaa);
                                fscanf(fp_fxd,"%lf",&aaa); printf("%lf\n",aaa);
                                fscanf(fp_fxd,"%lf",&aaa); printf("%lf\n",aaa);
                                fscanf(fp_fxd,"%lf",&aaa); printf("%lf\n",aaa); */
                            }
                            if(e==2003){
                                if((fp_fxd = fopen("tse_30min_2003b.txt","rt"))==NULL){
                                    printf("No tse_30min_2003b.txt\n");
                                    exit(1);
                                }
                            }
                            if(e==2004){
                                if((fp_fxd = fopen("tse_30min_2004b.txt","rt"))==NULL){
                                    printf("No tse_30min_2004b.txt\n");
                                    exit(1);
                                }
                            }
                            if(e==2005){
                                if((fp_fxd = fopen("tse_30min_2005b.txt","rt"))==NULL){
                                    printf("No tse_30min_2005b.txt\n");
                                    exit(1);
                                }
                            }
                        }
                    }
					
					/* seasonal roop */
					ndy = (e%4==0)?366:365;
					for(h=0;h<ndy;h++){
						loct->doy = h;
						/* simulation suing fixed climate data: 2010/09/06 by A.Ito */
						/* no diurnal and no  */
						if(FIX_CLIM==1){
							loct->climy = 2001;
							loct->CO2y = BYR;
							loct->doy = 120;
							
							if(e>=(BYR+10)){
								loct->CO2y = 2050;
							}
						}
						f_doyTmody(e, h, &(loct->month), &(loct->mday));
						loct->hour = 24;
						
						if(loct->mday == 0){
                            dmon[loct->month]= 0.0;
							flux_mon[0][loct->month] = flux_mon[1][loct->month] = flux_mon[2][loct->month] = 0.0;
							flux_mon[3][loct->month] = flux_mon[4][loct->month] = flux_mon[5][loct->month] = 0.0;
						}

						/*** local condition ************************************************/
						f_loct_proc(&(grid[f][g]), loct, echar, &(mass[f][g]), flux);
									
						/* disturbance *****/
						/* loct->climy = e; */
						if(SENSANS_DIST == 1){
							;
						}else{
							disturbance_regime(e, &(grid[f][g]), loct, echar, &(mass[f][g]), flux); /* */
						}
						/* loct->climy = BYR; */
						
						/* basic scheme *****/
						daily_scheme(&(grid[f][g]), loct, echar, &(mass[f][g]), flux);
                        
                        /**/
                        if(loct->veg_state == 1){
                            f_agriculture(&(grid[f][g]), loct, echar, &(mass[f][g]), flux);
                        }

						/* 30min flux ************************************************************/
						if(FLUX_SCHEME == 1){
							f_flux_site(&(grid[f][g]), loct, echar, &(mass[f][g]), flux, fp_fxd, fp_r);
						}
						loct->hour = 24;

						/* monthly NEE */
                        dmon[loct->month]+=1.0;
						flux_mon[0][loct->month] += (mass[f][g].tree).lai
                                    +loct->funder_c3*(mass[f][g].c3).lai+loct->funder_c4*(mass[f][g].c4).lai;
						flux_mon[1][loct->month] += flux->gpp*100.0;
						flux_mon[2][loct->month] += flux->er*100.0;
						flux_mon[3][loct->month] += flux->nep*100.0;
						flux_mon[4][loct->month] += flux->sr*100.0;
						flux_mon[5][loct->month] += (flux->soil).hr*100.0;

						aet_a += loct->aet;
						
						/* GHG schemes *******************************************************************/
						if(GHG_CALC==1 || N_CYCLE ==1){				
							/* CH4 */
							if(GHG_OUT == 1){
								fprintf(fp_r[7], "%4ld %2ld %2ld ", e, loct->month, h);
								fprintf(fp_r[7], "%10.3lf ", loct->vmc);
								fprintf(fp_r[7], "%10.3lf ", loct->wfps);
								fprintf(fp_r[7], "%10.3lf ", loct->water_table_depth);
								fprintf(fp_r[7], "%10.3lf   ", loct->tmp10_soil);
								
								fprintf(fp_r[7], "%10.3lf ", (flux->soil).ch4oxy_ridg);
								fprintf(fp_r[7], "%10.3lf ", (flux->soil).ch4oxy_casa);
								fprintf(fp_r[7], "%10.3lf ", (flux->soil).ch4oxy_delgrosso);
								fprintf(fp_r[7], "%10.3lf   ", (flux->soil).ch4oxy_curry);
								
								fprintf(fp_r[7], "%10.3lf ", (flux->soil).ch4_wh_diff);
								fprintf(fp_r[7], "%10.3lf ", (flux->soil).ch4_wh_ebull);
								fprintf(fp_r[7], "%10.3lf   ", (flux->soil).ch4_wh_plant);

							/*	fprintf(fp_r[7], "%10.3lf ", loct->prof_ch4[0]);
								fprintf(fp_r[7], "%10.3lf ", loct->prof_ch4[1]);
								fprintf(fp_r[7], "%10.3lf ", loct->prof_ch4[10]);
								
								fprintf(fp_r[7], "%10.3lf ", loct->x1);
								fprintf(fp_r[7], "%10.3lf ", loct->x2);
								fprintf(fp_r[7], "%10.3lf ", loct->x3);	*/
								
								for(i=0;i<=SOIL_LAYER;i++){
									fprintf(fp_r[7], "%10.2lf ", loct->prof_ch4[i]);
								}
								for(i=0;i<=SOIL_LAYER;i++){
									fprintf(fp_r[7], "%10.3lf ", loct->x4[i]);
								}
								
								fprintf(fp_r[7], "\n");
							}
									
							/* N2O */
							if(GHG_OUT == 1){
								fprintf(fp_r[8], "%4ld %2ld %2ld ", e, loct->month, loct->mday);
								
								fprintf(fp_r[8], "%7.2lf ", loct->prate_sfc);
								fprintf(fp_r[8], "%7.2lf ", (mass[f][g]).sw30);
								fprintf(fp_r[8], "%7.3lf ", loct->rdr_casa);
								fprintf(fp_r[8], "%lf ", loct->m_casa);
								fprintf(fp_r[8], "%lf ", loct->e_casa);
								fprintf(fp_r[8], "%lf ", loct->iw_casa);
								fprintf(fp_r[8], "%lf ", loct->wfps);
				
								fprintf(fp_r[8], "%lf ", (flux->soil).d_n2o_ngas);
								fprintf(fp_r[8], "%lf ", (flux->soil).d_n2o_ntr_ngas);
								fprintf(fp_r[8], "%lf ", (flux->soil).d_n2o_dnt_ngas);
								
								fprintf(fp_r[8], "%lf ", (flux->soil).d_n2o_ngas_2);
								fprintf(fp_r[8], "%lf ", (flux->soil).d_n2o_ntr_ngas_2);
								fprintf(fp_r[8], "%lf ", (flux->soil).d_n2o_dnt_ngas_2);
								
								fprintf(fp_r[8], "%lf ", (flux->soil).d_n2o_casa);
								fprintf(fp_r[8], "%lf ", (flux->soil).d_n2_casa);
								fprintf(fp_r[8], "%lf ", (flux->soil).d_no_casa);
							
								fprintf(fp_r[8], "\n"); /* */
							}
						}
						
						/******************************************************************************/
						f_ansis_ann(&(grid[f][g]), loct, echar, &(mass[f][g]), flux, ansis_ann);									
						
						/* daily result iutput */
						if(WMODE==1 && (loct->climy>=1990 && loct->climy<=2012)){
							output_ansis_daily(&(grid[f][g]), loct, echar, &(mass[f][g]), flux, fp_r[2]);
						}
						if((WMODE==2) && (loct->climy>=1990 && loct->climy<=2012)){
							output_ansis_daily(&(grid[f][g]), loct, echar, &(mass[f][g]), flux, fp_r[loct->climy - 1998 + 9]);
						}
					}
					f_erosion_rusle(&grid[f][g], loct, &mass[f][g], flux);
					
					(mass[f][g].soil).msl_a -= (flux->soil).erosion_carbon;
					
					/* annual result iutput */
					output_ansis_ann(loct->climy, ansis_ann, fp_r[1]);
				
                /*	fprintf(fp_r[1],"%lf ", loct->f_rain);
					fprintf(fp_r[1],"%lf ", loct->f_slope);
					fprintf(fp_r[1],"%lf ", loct->f_erodibility);
					fprintf(fp_r[1],"%lf ", loct->f_cover);
					fprintf(fp_r[1],"%lf ", loct->f_conservation);
					fprintf(fp_r[1],"%lf ", (flux->soil).erosion_soil);
					fprintf(fp_r[1],"%lf ", (flux->soil).erosion_orgmat);
					fprintf(fp_r[1],"%lf ", (flux->soil).erosion_carbon); */
                    
                   /*  for(h=0;h<12;h++){
                        fprintf(fp_r[2],"%ld %ld %lf %lf %lf %lf %lf %lf\n",e,h+1,flux_mon[0][h]/dmon[h],
                            flux_mon[1][h],flux_mon[2][h],flux_mon[3][h],flux_mon[4][h],flux_mon[5][h]);
                    } */
					
					/* monitoring */
					if(WMODE==2){
						printf("%4d %4d ", grid[f][g].dg_row, grid[f][g].dg_col);
					}
					printf("%4ld GPP:%5.1lf NPP:%5.1lf NEP:%6.3lf  LAI:%4.1lf TREE:%5.1lf C3:%5.1lf C4:%5.1lf  LITTER:%5.1lf SOIL:%6.1lf\n", 
						e, ansis_ann[0], ansis_ann[1], ansis_ann[2], 
						ansis_ann[5]+ansis_ann[9]+ansis_ann[13], 
						ansis_ann[6]+ansis_ann[7]+ansis_ann[8], 
						ansis_ann[10]+ansis_ann[11]+ansis_ann[12], 
						ansis_ann[14]+ansis_ann[15]+ansis_ann[16], 
						ansis_ann[17]+ansis_ann[18]+ansis_ann[19]+ansis_ann[20]+ansis_ann[21]+ansis_ann[22], 
						ansis_ann[23]+ansis_ann[24]+ansis_ann[25]); /* */
                    
                    if(FLUX_SCHEME == 1){
                        if(strcmp(grid[f][g].site_id, "CEAMIP_TSE")==0 && (e==2002||e==2003||e==2004||e==2005)){
                            fclose(fp_fxd);
                        }
                    }
				}
				/****** end of experimental roop *******/
				
				if(FLUX_SCHEME == 1){
                    if(strcmp(grid[f][g].site_id, "TKY")==0){
                        fclose(fp_fxd);
                    }
				}
			}
		}
	}
	
	if(USE_RESTART == 1){
		fclose(fp_restart);
	}

	if(NOTICE==1){
		printf("done\n");
	}
    
    if(FIX_PHENOLOGY == 1){
        fclose(fp_lp);
    }
}

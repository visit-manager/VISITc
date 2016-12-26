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
	struct Echar echar_type[NBIOME], 
	struct Mass mass[DROW][DCOL], 
	struct Flux *flux, 
	FILE *fp_r[NFILE]
){
	long e, f, g, h, i;
	char filename[128];
	long ndy, end_year, ddummy;
	double ansis_ann[N_ANSIS],ansis_mon[N_ANSIS],fdummy;
	double aet_a, dmon[12];
    long monday[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	FILE *fp_fxd, *fp_lp;
	FILE *fp_restart;
	
	if(NOTICE==1){
		printf("Start experimental phase\n");
	}
	
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
                
                /* LUC-MIP 2014/11/11 by A.Ito */
                if(EX_LUCMIP==11){
                    end_year = 2940;
                }
                if(EX_LUCMIP==12){
                    end_year = 2990;
                }
                if(EX_LUCMIP==13){
                    end_year = 3130;
                }
                if(EX_LUCMIP==2){
                    end_year = 2930;
                }
                if(EX_LUCMIP==3){
                    end_year = 2930;
                }
                if(EX_LUCMIP==4){
                    end_year = 1930+3000;
                }
                if(EX_LUCMIP==5){
                    end_year = 1930+3000;
                }
                if(EX_LUCMIP==6){
                    end_year = 2930;
                }
                if(EX_LUCMIP==7){
                    end_year = 2930;
                }
				
				/* roop for experimental stage ************************************************/
				for(e=BYR; e<=end_year; e++){	
					
					if(SENSANS_ACO2 == 1){
						loct->CO2y = BYR;
					}else if(SENSANS_ACO2 == 2){
						loct->CO2y = 2005;
					}else{
						loct->CO2y = e;
					}
                    
                    loct->climy = e;
                    /* LUC-MIP: random 1983-2012: 2014/11/06 by A.Ito */
                    if((strcmp(grid[f][g].site_id, "LUCMIP0")==0) || (strcmp(grid[f][g].site_id, "LUCMIP0")==1)
                        || (strcmp(grid[f][g].site_id, "LUCMIP0")==2) || (strcmp(grid[f][g].site_id, "LUCMIP0")==3)){
                        
                        //loct->climy = 1983 + (long)(30.0 * (double)rand() / (double)RAND_MAX);
                        loct->climy = 1983 + (e - BYR)%30;
                        
                        loct->CO2y = 2005;
                    }
                    
                    /* ASIAMIP: 2015/12/08 by A.Ito */
                    /* ASIAMIP: 2016/01/20 by A.Ito */
                    if(EX_ASIAMIP == 1){
                        if(e<=1900){
                            loct->CO2y = 1901;
                        }else{
                            loct->CO2y = e;
                        }
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
					
					/* seasonal roop ****************************/
                    if(e%4!=0 || LEAPYEAR == 0 || EX_ASIAMIP == 1){
                        ndy = 365;
                        monday[1] = 28;
                    }else{
                        ndy = 366;
                        monday[1] = 29;
                    }
                    
					for(h=0;h<ndy;h++){
                    
						loct->doy = h;
						/* simulation suing fixed climate data: 2010/09/06 by A.Ito */
						/* no diurnal and no  */
						if(FIX_CLIM == 1){
							loct->climy = 2001;
							loct->CO2y = BYR;
							loct->doy = 120;
							
							if(e >= (BYR+10)){
								loct->CO2y = 2050;
							}
						}
                        
                        if(EX_ASIAMIP == 1){
                            f_doyTmody(2001, h, &(loct->month), &(loct->mday));
                        }else{
                            f_doyTmody(e, h, &(loct->month), &(loct->mday));
                        }
						loct->hour = 24;
                        
						if(loct->mday == 0){
                            dmon[loct->month]= 0.0;
                            
                            for(i=0;i<N_ANSIS;i++){
                                ansis_mon[i] = 0.0;
                            }
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
						f_ansis_mon(&(grid[f][g]), loct, echar, &(mass[f][g]), flux, ansis_mon);

						f_ansis_ann(&(grid[f][g]), loct, echar, &(mass[f][g]), flux, ansis_ann);
						
						/* daily result iutput */
						//if(WMODE==1 && (loct->climy>=1995 && loct->climy<=2013)){
							output_ansis_daily(&(grid[f][g]), loct, echar, &(mass[f][g]), flux, fp_r[1]);
						//}

                        if(loct->mday == (monday[loct->month]-1)){
                            output_ansis_mon(1, loct->adyear, loct->month, ansis_mon, fp_r[2]);
                            output_ansis_mon(2, loct->adyear, loct->month, ansis_mon, fp_r[3]);
                        }
					}
                    
					f_erosion_rusle(&grid[f][g], loct, &mass[f][g], flux);
                    if(NECB_POC==1){
                        (mass[f][g].soil).msl_a -= (flux->soil).erosion_carbon;
                    }
                    
					/* annual result output */
					output_ansis_ann(loct->adyear, ansis_ann, fp_r[4]);
				
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

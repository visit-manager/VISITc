/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/* header files */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/*********************************************************************/
void f_initialize(
	struct Grid grid[DROW][DCOL], 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Echar echar_type[17], 
	struct Mass mass[DROW][DCOL], 
	struct Flux *flux, 
	FILE *fp_r[17]
){
	char filename[100], num[4];
	long h,i,yr,ndy,year,day,nn;
	double dd[4], dummy, sunshine, ss;
	FILE *fp_ghg, *fp_ndepo, *fp_gwp, *fp_dmw, *fp_wtd, *fp_altp;
	
	/* Number of days for each month */
	month_day[0] = 31;	/* January */
	month_day[1] = 28;	/* February */
	month_day[2] = 31;	/* March */
	month_day[3] = 30;	/* April */
	month_day[4] = 31;	/* May */
	month_day[5] = 30;	/* June */
	month_day[6] = 31;	/* July */
	month_day[7] = 31;	/* August */
	month_day[8] = 30;	/* September */
	month_day[9] = 31;	/* October */
	month_day[10] = 30;	/* November */
	month_day[11] = 31;	/* December */

	/******************************************************************/	
	/** clear all parameters *********************/
	if(NOTICE==1){
		printf("Clearing...");
	}
	for(h=0; h<DROW; h++){
		for(i=0; i<DCOL; i++){
			clear_a(&(grid[h][i]), &(mass[h][i]));

			grid[h][i].dg_row = h;
			grid[h][i].dg_col = i;
		}
	}
	clear_b(loct, echar, flux);	
	if(NOTICE==1){
		printf("done\n");
	}

	/******************************************************************/	
	/** initialize grid condition **/
	if(NOTICE==1){
		printf("Initializing site...");
	}
    
    f_init_site(&grid[0][0]);

    if(NOTICE==1){
		printf("done\n");
	}
		
	/******************************************************************/	
	if(NOTICE==1){
		printf("Reading GHG scenario...");
	}
	if((fp_ghg=fopen("AtmGHG_timeseries.txt","rt"))==NULL){
		printf("No atmospheric GHG scenario!\n");
		exit(1);
	}
	for(h=0;h<553;h++){
		fscanf(fp_ghg,"%ld", &yr);
		fscanf(fp_ghg,"%lf", &loct->atm_co2_a1[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_co2_b1[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_co2_a2[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_co2_b2[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_ch4_a1[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_ch4_b1[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_ch4_a2[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_ch4_b2[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_n2o_a1[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_n2o_b1[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_n2o_a2[h]);
		fscanf(fp_ghg,"%lf", &loct->atm_n2o_b2[h]);
	}
	fclose(fp_ghg);
	if(NOTICE==1){
		printf("done\n");
	}	

	/***************************************************************/
	if(NOTICE==1){
		printf("Reading N deposition scenario...");
	}
	
	/* initialization */
	for(i=0;i<366;i++){
		loct->depo_nh4_model_av[i] = 0.0;
		loct->depo_no3_model_av[i] = 0.0;
	}
	
	if(strcmp(grid[0][0].site_id, "TKY")==0){
		if((fp_ndepo = fopen("ndepo_1980-2003_tky.txt","rt"))==NULL){
			printf("No atmospheric N deposition data!\n");
			exit(1);
		}
	}else if(strcmp(grid[0][0].site_id, "FJY")==0){
		if((fp_ndepo = fopen("ndepo_1980-2003_fjy.txt","rt"))==NULL){
			printf("No atmospheric N deposition data!\n");
			exit(1);
		}
	}else if(strcmp(grid[0][0].site_id, "TMK")==0){
		if((fp_ndepo = fopen("ndepo_1980-2003_tmk.txt","rt"))==NULL){
			printf("No atmospheric N deposition data!\n");
			exit(1);
		}
	}else if(strcmp(grid[0][0].site_id, "QHB")==0){
		if((fp_ndepo = fopen("ndepo_1980-2003_qhb.txt","rt"))==NULL){
			printf("No atmospheric N deposition data!\n");
			exit(1);
		}
	}
	
	if((strcmp(grid[0][0].site_id, "TKY")==0) || (strcmp(grid[0][0].site_id, "FJY")==0)
		 || (strcmp(grid[0][0].site_id, "TMK")==0)){
		for(h=1980;h<=2003;h++){
			ndy = (h%4==0)?366:365;
			for(i=0;i<ndy;i++){
				fscanf(fp_ndepo,"%ld %ld", &year, &day);
				fscanf(fp_ndepo,"%lf %lf %lf %lf", &dd[0], &dd[1], &dd[2], &dd[3]);
				
				/* g N/ha/day */
				loct->depo_nh4_model[h-1980][i] = (dd[0] + dd[2])*1000.0;
				loct->depo_no3_model[h-1980][i] = (dd[1] + dd[3])*1000.0;
				
				/* 24-year average */
				loct->depo_nh4_model_av[i] += (dd[0] + dd[2])*1000.0/24.0;
				loct->depo_no3_model_av[i] += (dd[1] + dd[3])*1000.0/24.0;
			}
		}
		fclose(fp_ndepo);
	}

	if(NOTICE==1){
		printf("done\n");
	}
	
	/* water-table depth of paddy field *********************/
	if((strcmp(grid[0][0].site_id, "MSE")==0) && WH_CH4==1){
		if((fp_wtd = fopen("MSE_wtd_2001-2005.txt","rt"))==NULL){
			printf("No water-table depth data!\n");
			exit(1);
		}
		for(i=0;i<366;i++){
			fscanf(fp_wtd,"%ld %lf %lf %lf %lf %lf", &day,
            &paddy_wtd[0][i], &paddy_wtd[1][i], &paddy_wtd[2][i], &paddy_wtd[3][i], &paddy_wtd[4][i]);
		}
	}

	/*******************************************************/
	if(PREDICT >= 1){
		if(NOTICE==1){
			printf("Reading future climate scenario...");
		}
		
		for(h=1980;h<=2050;h++){
			for(i=0;i<366;i++){
				gcm_temp[h-1980][i] = 0.0;
				gcm_prec[h-1980][i] = 0.0;
				gcm_ahmd[h-1980][i] = 0.0;
				gcm_dswrf[h-1980][i] = 0.0;
			}
		}
		for(i=0;i<366;i++){
			gcm_temp_av[i] = 0.0;
			gcm_prec_av[i] = 0.0;
			gcm_ahmd_av[i] = 0.0;
			gcm_dswrf_av[i] = 0.0;
		}
		
		/* climate scanario 1980-2050 */
		if(PREDICT == 1){ /* MIROC-high A1B ************************/
			if(strcmp(grid[0][0].site_id, "TKY")==0){
				if((fp_gwp = fopen("ext_MIROCh_tky.dat","rt"))==NULL){
					printf("No MIROC climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "FJY")==0){
				if((fp_gwp = fopen("ext_MIROCh_fjy.dat","rt"))==NULL){
					printf("No MIROC climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "TMK")==0){
				if((fp_gwp = fopen("ext_MIROCh_tmk.dat","rt"))==NULL){
					printf("No MIROC climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "QHB")==0){
				if((fp_gwp = fopen("ext_MIROCh_qhb.dat","rt"))==NULL){
					printf("No MIROC climate data!\n");
					exit(1);
				}
			}
		}else if(PREDICT == 10){ /* MPI-ECHAM A1B ************************/
			if(strcmp(grid[0][0].site_id, "TKY")==0){
				if((fp_gwp = fopen("int_ext_MPI_A1B_TKY.dat","rt"))==NULL){
					printf("No MPI-ECHAM climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "FJY")==0){
				if((fp_gwp = fopen("int_ext_MPI_A1B_FJY.dat","rt"))==NULL){
					printf("No MPI-ECHAM climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "TMK")==0){
				if((fp_gwp = fopen("int_ext_MPI_A1B_TMK.dat","rt"))==NULL){
					printf("No MPI-ECHAM climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "QHB")==0){
				if((fp_gwp = fopen("int_ext_MPI_A1B_QHB.dat","rt"))==NULL){
					printf("No MPI-ECHAM climate data!\n");
					exit(1);
				}
			}
		}else if(PREDICT == 20){ /* CCCma A1B ************************/
			if(strcmp(grid[0][0].site_id, "TKY")==0){
				if((fp_gwp = fopen("int_ext_CCC_A1B_TKY.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "FJY")==0){
				if((fp_gwp = fopen("int_ext_CCC_A1B_FJY.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "TMK")==0){
				if((fp_gwp = fopen("int_ext_CCC_A1B_TMK.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "QHB")==0){
				if((fp_gwp = fopen("int_ext_CCC_A1B_QHB.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}
		}else if(PREDICT == 21){ /* CCCma A2 ************************/
			if(strcmp(grid[0][0].site_id, "TKY")==0){
				if((fp_gwp = fopen("int_ext_CCC_A2_TKY.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "FJY")==0){
				if((fp_gwp = fopen("int_ext_CCC_A2_FJY.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "TMK")==0){
				if((fp_gwp = fopen("int_ext_CCC_A2_TMK.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "QHB")==0){
				if((fp_gwp = fopen("int_ext_CCC_A2_QHB.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}
		}else if(PREDICT == 22){ /* CCCma B1 ************************/
			if(strcmp(grid[0][0].site_id, "TKY")==0){
				if((fp_gwp = fopen("int_ext_CCC_B1_TKY.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "FJY")==0){
				if((fp_gwp = fopen("int_ext_CCC_B1_FJY.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "TMK")==0){
				if((fp_gwp = fopen("int_ext_CCC_B1_TMK.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}else if(strcmp(grid[0][0].site_id, "QHB")==0){
				if((fp_gwp = fopen("int_ext_CCC_B1_QHB.dat","rt"))==NULL){
					printf("No CCCma climate data!\n");
					exit(1);
				}
			}
		}else if(PREDICT == 100){ /* MIROC-high A1B ************************/
			if(strcmp(grid[0][0].site_id, "PSO")==0){
				if((fp_gwp = fopen("ext_MIROCh_pso.dat","rt"))==NULL){
					printf("No MIROC-PSO climate data!\n");
					exit(1);
				}
            }
        }else if(PREDICT == 101){ /* MPI_ECHAM A1B ************************/
			if(strcmp(grid[0][0].site_id, "PSO")==0){
				if((fp_gwp = fopen("mpi_a1b_pso.txt","rt"))==NULL){
					printf("No MPI climate data!\n");
					exit(1);
				}
            }
        }else if(PREDICT == 102){ /* CCCma A1B ************************/
			if(strcmp(grid[0][0].site_id, "PSO")==0){
				if((fp_gwp = fopen("cccma_a1b_pso.txt","rt"))==NULL){
					printf("No CCC A1B climate data!\n");
					exit(1);
				}
            }
        }else if(PREDICT == 103){ /* CCCma A2 ************************/
			if(strcmp(grid[0][0].site_id, "PSO")==0){
				if((fp_gwp = fopen("cccma_a2_pso.txt","rt"))==NULL){
					printf("No CCC A2 climate data!\n");
					exit(1);
				}
            }
        }else if(PREDICT == 104){ /* CCCma B1 ************************/
			if(strcmp(grid[0][0].site_id, "PSO")==0){
				if((fp_gwp = fopen("cccma_b1_pso.txt","rt"))==NULL){
					printf("No CCC B1 climate data!\n");
					exit(1);
				}
            }
        }else{
			printf("Wrong GCM scenario number !!!!!");
			exit(1);
		}

		if((strcmp(grid[0][0].site_id, "TKY")==0) || (strcmp(grid[0][0].site_id, "FJY")==0)
			 || (strcmp(grid[0][0].site_id, "TMK")==0) || (strcmp(grid[0][0].site_id, "QHB")==0)){
			for(h=1980;h<=2050;h++){
				ndy = (h%4==0)?366:365;
				for(i=0;i<ndy;i++){
					fscanf(fp_gwp,"%ld %ld", &year, &day);
					fscanf(fp_gwp,"%lf %lf %lf %lf", &gcm_temp[h-1980][i], &gcm_prec[h-1980][i], 
							&gcm_ahmd[h-1980][i], &gcm_dswrf[h-1980][i]);
					
					gcm_prec[h-1980][i] *= 24.0*3600.0;
					
					if(h>=1980 && h<2000){
						gcm_temp_av[i] += gcm_temp[h-1980][i] / 20.0;
						gcm_prec_av[i] += gcm_prec[h-1980][i] / 20.0;
						gcm_ahmd_av[i] += gcm_ahmd[h-1980][i] / 20.0;
						gcm_dswrf_av[i] += gcm_dswrf[h-1980][i] / 20.0;
					}
				}
			}
            gcm_temp_av[365] = gcm_temp_av[364];
            gcm_prec_av[365] = gcm_prec_av[364];
            gcm_ahmd_av[365] = gcm_ahmd_av[364];
            gcm_dswrf_av[365] = gcm_dswrf_av[364];
			fclose(fp_gwp);
		}
        
        if((strcmp(grid[0][0].site_id, "PSO")==0)){
        
            if(PREDICT==100){
                for(h=1980;h<=2100;h++){
                    ndy = (h%4==0)?366:365;
                    for(i=0;i<ndy;i++){
                        fscanf(fp_gwp,"%ld %ld", &year, &day);
                        fscanf(fp_gwp,"%lf %lf %lf %lf", &gcm_temp[h-1980][i], &gcm_prec[h-1980][i], 
                                &gcm_ahmd[h-1980][i], &gcm_dswrf[h-1980][i]);
                        
                        gcm_prec[h-1980][i] *= 24.0*3600.0;
                        
                        if(h>=1980 && h<2000){
                            gcm_temp_av[i] += gcm_temp[h-1980][i] / 20.0;
                            gcm_prec_av[i] += gcm_prec[h-1980][i] / 20.0;
                            gcm_ahmd_av[i] += gcm_ahmd[h-1980][i] / 20.0;
                            gcm_dswrf_av[i] += gcm_dswrf[h-1980][i] / 20.0;
                        }
                    }
                }
                gcm_temp_av[365] = gcm_temp_av[364];
                gcm_prec_av[365] = gcm_prec_av[364];
                gcm_ahmd_av[365] = gcm_ahmd_av[364];
                gcm_dswrf_av[365] = gcm_dswrf_av[364];
            }else{
                for(h=1980;h<=2050;h++){
                    ndy = (h%4==0)?366:365;
                    for(i=0;i<ndy;i++){
                        fscanf(fp_gwp,"%ld %ld", &year, &day);
                        fscanf(fp_gwp,"%lf %lf %lf %lf", &gcm_temp[h-1980][i], &gcm_prec[h-1980][i], 
                                &gcm_ahmd[h-1980][i], &gcm_dswrf[h-1980][i]);
                        
                        gcm_prec[h-1980][i] *= 24.0*3600.0;
                        
                        if(h>=1980 && h<2000){
                            gcm_temp_av[i] += gcm_temp[h-1980][i] / 20.0;
                            gcm_prec_av[i] += gcm_prec[h-1980][i] / 20.0;
                            gcm_ahmd_av[i] += gcm_ahmd[h-1980][i] / 20.0;
                            gcm_dswrf_av[i] += gcm_dswrf[h-1980][i] / 20.0;
                        }
                    }
                }
                gcm_temp_av[365] = gcm_temp_av[364];
                gcm_prec_av[365] = gcm_prec_av[364];
                gcm_ahmd_av[365] = gcm_ahmd_av[364];
                gcm_dswrf_av[365] = gcm_dswrf_av[364];
            }
			fclose(fp_gwp);
		}
		
		if(NOTICE==1){
			printf("done\n");
		}
	}
	
	/* alternative precipitation data by APHRODITE *********************/
	/* if(strcmp(grid[0][0].site_id, "KBU")==0){
		if((fp_altp = fopen("ext_aprodite_kbu.txt","rt"))==NULL){
			printf("No APHRODITE-KBU data!\n");
			exit(1);
		}
		
		for(h=1951;h<=2007;h++){
			ndy = (h%4==0)?366:365;
			
			for(i=0;i<ndy;i++){
				fscanf(fp_altp,"%ld %ld %lf %lf", &year, &day, 
									&prec_altdata[h-1951][i],&dummy);
			}
		}
		fclose(fp_altp);
	} */
	/* if(strcmp(grid[0][0].site_id, "Tongyu")==0){
		if((fp_altp = fopen("ext_aprodite_tongyu.txt","rt"))==NULL){
			printf("No APHRODITE-Tongyu data!\n");
			exit(1);
		}
		
		for(h=1951;h<=2007;h++){
			ndy = (h%4==0)?366:365;
			
			for(i=0;i<ndy;i++){
				fscanf(fp_altp,"%ld %ld %lf %lf", &year, &day, 
					   &prec_altdata[h-1951][i],&dummy);
			}
		}
		fclose(fp_altp);
	} */
	
	/*****************************************************************/
	if(SENSANS_TEMP==2 || SENSANS_PREC==2 || SENSANS_SRAD==2){
		if(NOTICE==1){
			printf("Reading LARS-generated dummy weather scenario...");
		}

		if(strcmp(grid[0][0].site_id, "TKY")==0){
			strcpy(filename, "./tky_lars/TakayamaWG_");
			nn = (long)N_LARS;
			snprintf(num, 4, "%03ld", nn);
			strcat(filename, num);
			strcat(filename, ".dat");
			if((fp_dmw = fopen(filename,"rt"))==NULL){
				printf("No TakayamaWG_XXX.dat!\n");
				exit(1);
			}
		}else if(strcmp(grid[0][0].site_id, "FJY")==0){
			strcpy(filename, "fjy_BaseWG");
			nn = (long)N_LARS;
			snprintf(num, 4, "%02ld", nn);
			strcat(filename, num);
			strcat(filename, ".dat");
			if((fp_dmw = fopen(filename,"rt"))==NULL){
				printf("No fjy_BaseWGXX.dat!\n");
				exit(1);
			}
		}else if(strcmp(grid[0][0].site_id, "TMK")==0){
			strcpy(filename, "tmk_BaseWG");
			nn = (long)N_LARS;
			snprintf(num, 4, "%02ld", nn);
			strcat(filename, num);
			strcat(filename, ".dat");
			if((fp_dmw = fopen(filename,"rt"))==NULL){
				printf("No tmk_BaseWGXX.dat!\n");
				exit(1);
			}
		}else if(strcmp(grid[0][0].site_id, "QHB")==0){
			strcpy(filename, "qhb_BaseWG");
			nn = (long)N_LARS;
			snprintf(num, 4, "%02ld", nn);
			strcat(filename, num);
			strcat(filename, ".dat");
			if((fp_dmw = fopen(filename,"rt"))==NULL){
				printf("No qhb_BaseWGXX.dat!\n");
				exit(1);
			}
		}else{
			printf("!! Wrong site ID for LARS dummy data\n");
			exit(1);
		}
		
		for(h=0;h<HCLIMD;h++){
			for(i=0;i<365;i++){
				fscanf(fp_dmw,"%ld %ld", &year, &day);
				fscanf(fp_dmw,"%lf", &lars_tmin[h][i]);
				fscanf(fp_dmw,"%lf", &lars_tmax[h][i]);
				fscanf(fp_dmw,"%lf", &lars_prec[h][i]);
				fscanf(fp_dmw,"%lf", &lars_srad[h][i]);
				fscanf(fp_dmw,"%lf", &sunshine);
			}
		}
		
		fclose(fp_dmw);
		
		if(NOTICE==1){
			printf("done\n");
		}
	}
    
    if(SPINUP==4){
        fp_dmw = fopen("TakayamaWG_10000y.txt","rt");
        
        for(h=0;h<N_SU_LARS;h++){
			for(i=0;i<365;i++){
				fscanf(fp_dmw,"%ld %ld", &year, &day);
				fscanf(fp_dmw,"%lf", &su_lars_tmin[h][i]);
				fscanf(fp_dmw,"%lf", &su_lars_tmax[h][i]);
				fscanf(fp_dmw,"%lf", &su_lars_prec[h][i]);
				fscanf(fp_dmw,"%lf", &su_lars_srad[h][i]);
				fscanf(fp_dmw,"%lf", &sunshine);
			}
            
            su_lars_tmin[h][365] = su_lars_tmin[h][364];
            su_lars_tmax[h][365] = su_lars_tmax[h][364];
            su_lars_prec[h][365] = su_lars_prec[h][364];
            su_lars_srad[h][365] = su_lars_srad[h][364];
		}
        
        fclose(fp_dmw);
    }

	/******************************************************************/	
	/** initialize vegatation and soil conditions **/ 
	if(NOTICE==1){
		printf("Initializing vegetation and soil...");
	}
	/*** parameters ***/
	strcpy(filename, "parameter_S1b.txt");
    
    if(EX_LUCMIP==3 || EX_LUCMIP==5 || EX_LUCMIP==7){
        
        if(strcmp(grid[0][0].site_id, "LUCMIP0")==0){
            strcpy(grid[0][0].file_para,"parameter_LUCMIP0_crop.txt");
        }
        if(strcmp(grid[0][0].site_id, "LUCMIP1")==0){
            strcpy(grid[0][0].file_para,"parameter_LUCMIP1_crop.txt");
        }
        if(strcmp(grid[0][0].site_id, "LUCMIP2")==0){
            strcpy(grid[0][0].file_para,"parameter_LUCMIP2_crop.txt");
        }
        if(strcmp(grid[0][0].site_id, "LUCMIP3")==0){
            strcpy(grid[0][0].file_para,"parameter_LUCMIP3_crop.txt");
        }
    }
    
    set_parameter(1, grid[0][0].file_para, echar);
    set_parameter(2, filename, echar_type);
	
	/* using default parameter set ************************ (2008/12/07 A.Ito) */
	if(PARA_CHOICE == 1){
		*echar = echar_type[grid[0][0].veg_type];
	}
	
    /* parameter perturbation *****************************/
	if(SA_PARA_EN==1){
		(echar->tree).topt0 *= SA_PARA_VAR;
		(echar->c3).topt0 *= SA_PARA_VAR;
		(echar->c4).topt0 *= SA_PARA_VAR;
	}
	
	if(SA_PARA_EN==2){
		(echar->tree).kmci *= SA_PARA_VAR;
		(echar->c3).kmci *= SA_PARA_VAR;
		(echar->c4).kmci *= SA_PARA_VAR;
	}
	
	if(SA_PARA_EN==3){
		(echar->tree).km_nstl *= SA_PARA_VAR;
		(echar->c3).km_nstl *= SA_PARA_VAR;
		(echar->c4).km_nstl *= SA_PARA_VAR;
	}
	
	if(SA_PARA_EN==4){
		(echar->tree).lue0 *= SA_PARA_VAR;
		(echar->c3).lue0 *= SA_PARA_VAR;
		(echar->c4).lue0 *= SA_PARA_VAR;
	}
	
	if(SA_PARA_EN==5){
		(echar->tree).qTf0 *= SA_PARA_VAR;
		(echar->c3).qTf0 *= SA_PARA_VAR;
		(echar->c4).qTf0 *= SA_PARA_VAR;

		(echar->tree).qTc0 *= SA_PARA_VAR;
		(echar->c3).qTc0 *= SA_PARA_VAR;
		(echar->c4).qTc0 *= SA_PARA_VAR;

		(echar->tree).qTr0 *= SA_PARA_VAR;
		(echar->c3).qTr0 *= SA_PARA_VAR;
		(echar->c4).qTr0 *= SA_PARA_VAR;
	}
	
    /* multiple perturbation: 2012/02/20 by A.Ito */
    if(PERTURB_MPARA == 1){
        /* pmax */
        ss = 0.0;
        for(h=0;h<12;h++){
            ss += (double)rand() / RAND_MAX;
        }
        ss -= 6.0;
        ss = ss*0.05 + 1.0;
        if(ss<0.85){
            ss = 0.85;
        }
        if(ss>1.15){
            ss = 1.15;
        }
        (echar->tree).pmax *= ss;
        (echar->c3).pmax *= ss;
        (echar->c4).pmax *= ss;
        
        /* lue0 */
        ss = 0.0;
        for(h=0;h<12;h++){
            ss += (double)rand() / RAND_MAX;
        }
        ss -= 6.0;
        ss = ss*0.05 + 1.0;
        if(ss<0.85){
            ss = 0.85;
        }
        if(ss>1.15){
            ss = 1.15;
        }
         (echar->tree).lue0 *= ss;
        (echar->c3).lue0 *= ss;
        (echar->c4).lue0 *= ss; /**/
        
        /* Q10 */
        ss = 0.0;
        for(h=0;h<12;h++){
            ss += (double)rand() / RAND_MAX;
        }
        ss -= 6.0;
        ss = ss*0.1 + 1.0;
        if(ss<0.8){
            ss = 0.8;
        }
        if(ss>1.2){
            ss = 1.2;
        }
        (echar->tree).qTf0 *= ss;
        (echar->c3).qTf0 *= ss;
        (echar->c4).qTf0 *= ss;
        (echar->tree).qTc0 *= ss;
        (echar->c3).qTc0 *= ss;
        (echar->c4).qTc0 *= ss;
        (echar->tree).qTr0 *= ss;
        (echar->c3).qTr0 *= ss;
        (echar->c4).qTr0 *= ss; 
        
        /* soil decomposition */
        ss = 0.0;
        for(h=0;h<12;h++){
            ss += (double)rand() / RAND_MAX;
        }
        ss -= 6.0;
        ss = ss*0.1 + 1.0;
        if(ss<0.7){
            ss = 0.7;
        }
        if(ss>1.3){
            ss = 1.3;
        }
        (echar->soil).sr_lf *= ss;
        (echar->soil).sr_lc *= ss;
        (echar->soil).sr_lr *= ss;
        (echar->soil).sr_ha *= ss;
        (echar->soil).sr_hi *= ss;
        (echar->soil).sr_hp *= ss;
    }
    
    /* initialize stand age */
	loct->age_stand = 0.0;
	
	/*** mass ***/
	/* initialize */
	for(h=0;h<NROW;h++){
		for(i=0;i<NCOL;i++){

			switch(WMODE){
				case 1:
					break;
				case 2: 
					echar = &(echar_type[grid[h][i].veg_type]);
					break;
			}
			if((echar->tree).phenoltype != 0){
				/* forests */
				initTree(&(mass[h][i]));
			}
			initFloor(&(mass[h][i]));
			initSoil(&(mass[h][i]));
		}
	}
	if(NOTICE==1){
		printf("done\n");
	}

	/** initialize location conditions **/
	if(NOTICE==1){
		printf("Initializing local environment...\n");
	}
	for(h=0; h<NROW; h++){
		for(i=0; i<NCOL; i++){
			if(NOTICE==1){
				printf("%ld - %ld\n", h,i);
			}
			
			if(grid[h][i].flag_datavl == 1){
				f_loct_init(&grid[h][i], loct, echar, &mass[h][i], flux);	
			}
		}
	}
	if(NOTICE==1){
		printf("done\n");
	}
	
	/** initialize stable carbon isotope **/ 
	if(NOTICE==1){
		printf("Initializing stable carbon isotope...");
	}
    
    f_init_d13c(&grid[0][0], loct, echar, &mass[0][0], flux);

	if(NOTICE==1){
		printf("done\n");
	}

}

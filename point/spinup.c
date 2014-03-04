/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */
/*	SASU-implemented version in October 29, 2013                        */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* spin-up *****************************************************************/
void f_spinup(
	struct Grid grid[DROW][DCOL], 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Echar echar_type[NFILE], 
	struct Mass mass[DROW][DCOL], 
	struct Flux *flux, 
	FILE *fp_spinup
){
	long e, f, g, h, i, j, k;
	char filename[128];
	long ndy, nroop;
	double gpp_a, npp_a, nep_a, lai_a, plant_a, soil_a;
	double aaa[20], xx;
	FILE *fp_restart;
    
    /* variables for SASU: 2013/10/29 */
    double ss_v[9], ss_l[6], ss_h[3];
    double nsasu_v[10], nsasu[10];
    double nppi[10][9],lf[10][9], tpc[10][9];
    double li[10][6][9],ld[10][6];
    double hf[10][3][6],hd[10][3];
	
	/* roop to stable stage *****************************/	
	strcpy(filename, grid[0][0].site_id);
	strcat(filename, "_restart.txt");
	fp_restart = fopen(filename,"wt");

	if(NOTICE==1){
		printf("Start spin-up phase\n");
	}
	
	if(WMODE==1){
		fprintf(fp_spinup,"%s %s\n", echar->para_ver_id, echar->para_date_id);
	}else if(WMODE==2){
		fprintf(fp_spinup,"%s %s\n", echar_type[0].para_ver_id, echar_type[0].para_date_id);
	}
    
    srand(113);

	for(h=0; h<NROW; h++){
		for(i=0; i<NCOL; i++){
			if(WMODE==2){
				echar = &(echar_type[grid[h][i].veg_type]);
			}
			
			loct->time = 0;
			loct->age_stand = 0.0;
			if(grid[h][i].flag_datavl==1){
				loct->time = 0;
				nroop = (long)(SPUPT/PERIOD);
				
				for(g=0; g<=nroop; g++){
					
                    /* yearly roop **********************************/
					for(e=BYR ; e<=EYR ; e++){	
						
                        /**/
                        if(RANDCLIM_SU == 1){
                            loct->climy = BYR + (long)(PERIOD * (double)rand() / (double)RAND_MAX);
                        }else{
                            loct->climy = e;
                        }
                        
						loct->CO2y = BYR;
						loct->time++;
						
						/* stand age, year */
						loct->age_stand += 1.0;
				
						loct->adyear = ((BYR-1) - PERIOD*(nroop - g) - (EYR - e));
						
						/* ndy = (e%4==0)?366:365; */
						ndy = 365;
						/** initialize climate conditions **/
						gpp_a = npp_a = nep_a = lai_a = plant_a = soil_a = 0.0;
						for(f=0;f<20;f++){
							aaa[f] = 0.0;
						}
						
                        /* SASU: 2013/10/29 by A.Ito */
                        if(USE_SASU ==2 && (loct->time>=20 && loct->time<=50)){
                            nsasu_v[loct->time%10] = 0.0;
                            /* initialization: vegetation */
                            for(j=0;j<9;j++){
                                ss_v[j] = 0.0;
                                nppi[loct->time%10][j] = tpc[loct->time%10][j] = 0.0;
                                lf[loct->time%10][j] = 0.0;
                            }
                        }

                        if(USE_SASU >=1 && ((loct->time>=30 && loct->time<=60)
                                || loct->time==80 || loct->time==110)){
                            
                            nsasu[loct->time%10] = 0.0;
                            /* initialization: litter */
                            for(j=0;j<6;j++){
                                ss_l[j] = 0.0;
                                for(k=0;k<9;k++){
                                    li[loct->time%10][j][k] = 0.0;
                                }
                                ld[loct->time%10][j] = 0.0;
                            }

                            /* initialization: humus */
                            for(j=0;j<3;j++){
                                ss_h[j] = 0.0;
                                for(k=0;k<6;k++){
                                    hf[loct->time%10][j][k] = 0.0;
                                }
                                hd[loct->time%10][j] = 0.0;
                            }
                        }
                        
                        //printf("YYY %ld\n", loct->time%10);
                        
                        /* daily roop **********************************/
                        xx = 0;
                        for(f=0;f<ndy;f++){
							loct->doy = f;
							if(FIX_CLIM==1){
								loct->climy = 2001;
								loct->CO2y = BYR;
								loct->doy = 120;
							}
							f_doyTmody(e, f, &(loct->month), &(loct->mday));
															
							/*** set environmental condition ***/
							f_loct_proc(&grid[h][i], loct, echar, &mass[h][i], flux);	
							
							/***** disturbance *****/
							/* loct->climy = e; */
							if(SENSANS_DIST==1){
								;
							}else{
								disturbance_regime(loct->adyear, &grid[h][i], loct, echar, &mass[h][i], flux); /* */
							}
							/* loct->climy = BYR; */
							
							/* basic scheme ******************************************/
							daily_scheme(&grid[h][i], loct, echar, &mass[h][i], flux);
							
							if(flux->npp > loct->npp_max){
								loct->npp_max = flux->npp;
							}
							if(loct->npp_max < 1.0){
								loct->npp_max = 1.0;
							}

							/* fprintf(fp_spinup, "%4ld %4ld ", e, f); 
							fprintf(fp_spinup, "%6.2lf ", (flux->tree).gpp*100.0); 
							fprintf(fp_spinup, "%6.2lf ", (echar->tree).psat); 
							fprintf(fp_spinup, "%6.2lf ", (echar->tree).opt_lai); 
							fprintf(fp_spinup, "%6.3lf ", (mass->tree).lai); 
							fprintf(fp_spinup, "%6.2lf ", (mass->tree).plant); 
							fprintf(fp_spinup, "\n"); */
							
							gpp_a += (flux->tree).gpp + loct->funder_c3 * (flux->c3).gpp + loct->funder_c4 * (flux->c4).gpp;
							npp_a += (flux->tree).npp + loct->funder_c3 * (flux->c3).npp + loct->funder_c4 * (flux->c4).npp;
							nep_a += flux->nep;
							lai_a += ((mass[h][i].tree).lai + loct->funder_c3 * (mass[h][i].c3).lai + 
								loct->funder_c4 * (mass[h][i].c4).lai)/(double)ndy;
							plant_a += ((mass[h][i].tree).plant + loct->funder_c3 * (mass[h][i].c3).plant + 
								loct->funder_c4 * (mass[h][i].c4).plant)/(double)ndy;
							soil_a += (mass[h][i].soil).soil/(double)ndy;
                            
                            xx += (flux->tree).tpr - (flux->tree).rrg;
							
							loct->m_casa_pre = loct->m_casa;
							loct->vmc_pre = loct->vmc;
							
							aaa[0] += (mass[h][i].soil).n_lttr/365.0;
							aaa[1] += (flux->tree).n_abdn_cnpy + (flux->tree).n_abdn_strg + 
									(flux->c3).n_abdn_cnpy + (flux->c3).n_abdn_strg + 
									(flux->c4).n_abdn_cnpy + (flux->c4).n_abdn_strg;
							aaa[2] += (flux->soil).n_minerlz_lttr;
							aaa[3] += (mass[h][i].tree).n_canopy;
							aaa[4] += (mass[h][i].tree).n_storage/365.0;
							aaa[5] += (mass[h][i].c3).n_canopy/365.0;
							aaa[6] += (mass[h][i].c3).n_storage/365.0;
							aaa[7] += (flux->tree).uptake_no3 + (flux->tree).uptake_nh4;
							aaa[8] += (flux->c3).uptake_no3 + (flux->c3).uptake_nh4;
							aaa[9] += (flux->tree).n_alloc_cnpy;
							aaa[10] += (flux->c3).n_alloc_cnpy;
							aaa[11] += (flux->tree).n_alloc_strg;
							aaa[12] += (flux->c3).n_alloc_strg;
							aaa[13] += (mass[h][i].soil).n_hums/365.0;
							aaa[14] += (mass[h][i].soil).n_mcrb/365.0;
							
							if(g==nroop){
								(echar->tree).lai_contemp[f] += (mass[h][i].tree).lai / (double)(EYR-BYR+1);
								(echar->c3).lai_contemp[f] += (mass[h][i].c3).lai / (double)(EYR-BYR+1);
								(echar->c4).lai_contemp[f] += (mass[h][i].c4).lai / (double)(EYR-BYR+1);
							}
                            
                            /* SASU: 2013/10/29 by A.Ito */
                            if(USE_SASU ==2 && (loct->time>=20 && loct->time<=59)){
                                
                                /* vegetation */
                                for(j=0;j<10;j++){
                                    nsasu_v[j] += 1.0;
                                
                                    nppi[j][0] += (flux->tree).tpf + (flux->tree).rtpc + (flux->tree).rtpr - (flux->tree).rfg;
                                    nppi[j][1] += (flux->tree).tpc - (flux->tree).rcg;
                                    nppi[j][2] += (flux->tree).tpr - (flux->tree).rrg;
                                    nppi[j][3] += (flux->c3).tpf - (flux->c3).rfg;
                                    nppi[j][4] += (flux->c3).tpc - (flux->c3).rcg;
                                    nppi[j][5] += (flux->c3).tpr - (flux->c3).rrg;
                                    nppi[j][6] += (flux->c4).tpf - (flux->c4).rfg;
                                    nppi[j][7] += (flux->c4).tpc - (flux->c4).rcg;
                                    nppi[j][8] += (flux->c4).tpr - (flux->c4).rrg;

                                    lf[j][0] += (echar->tree).lf;
                                    lf[j][1] += (echar->tree).lc + (flux->tree).rtpc/(mass[h][i].tree).stm;
                                    lf[j][2] += (echar->tree).lr + (flux->tree).rtpr/(mass[h][i].tree).rot;
                                    lf[j][3] += (echar->c3).lf;
                                    lf[j][4] += (echar->c3).lc;
                                    lf[j][5] += (echar->c3).lr;
                                    lf[j][6] += (echar->c4).lf;
                                    lf[j][7] += (echar->c4).lc;
                                    lf[j][8] += (echar->c4).lr;
                                    
                                    tpc[j][0] += (flux->tree).tpf;
                                    tpc[j][1] += (flux->tree).tpc;
                                    tpc[j][2] += (flux->tree).tpr;
                                    tpc[j][3] += (flux->c3).tpf;
                                    tpc[j][4] += (flux->c3).tpc;
                                    tpc[j][5] += (flux->c3).tpr;
                                    tpc[j][6] += (flux->c4).tpf;
                                    tpc[j][7] += (flux->c4).tpc;
                                    tpc[j][8] += (flux->c4).tpr;
                                }
                               
                                /* if((flux->tree).tpc < 0.0){
                                    printf("%ld %ld %lf %lf\n", loct->adyear, f, (flux->tree).tpc);
                                } */
                            }
                            
                            if(USE_SASU >=1 && (loct->time>=30 && loct->time<=69)
                                || (loct->time>=80 && loct->time<=99)|| (loct->time>=110 && loct->time<=139)){
                            /* if(USE_SASU ==1 && ((loct->time>=30 && loct->time<=50) ||
                                                (loct->time>=60 && loct->time<=80) ||
                                                (loct->time>=90 && loct->time<=110) ||
                                                (loct->time>=120 && loct->time<=140) ||
                                                (loct->time>=150 && loct->time<=170)) ){ */
                                for(j=0;j<10;j++){
                                    nsasu[j] += 1.0;
                                    /* litter */
                                    li[j][0][3] += loct->fcover_c3 * (flux->c3).lf;
                                    li[j][1][4] += loct->fcover_c3 * (flux->c3).lc;
                                    li[j][2][5] += loct->fcover_c3 * (flux->c3).lr;
                                    li[j][0][6] += loct->fcover_c4 * (flux->c4).lf;
                                    li[j][1][7] += loct->fcover_c4 * (flux->c4).lc;
                                    li[j][2][8] += loct->fcover_c4 * (flux->c4).lr;
                                    li[j][3][0] += (flux->tree).lf;
                                    li[j][4][1] += (flux->tree).lc;
                                    li[j][5][2] += (flux->tree).lr;
                                    ld[j][0] += (echar->soil).sr_lf/1000.0 * (echar->soil).f_tm_l;
                                    ld[j][1] += (echar->soil).sr_lc/1000.0 * (echar->soil).f_tm_l;
                                    ld[j][2] += (echar->soil).sr_lr/1000.0 * (echar->soil).f_tm_l;
                                    ld[j][3] += (echar->soil).sr_lf/1000.0 * (echar->soil).f_tm_l;
                                    ld[j][4] += (echar->soil).sr_lc/1000.0 * (echar->soil).f_tm_l;
                                    ld[j][5] += (echar->soil).sr_lr/1000.0 * (echar->soil).f_tm_l;
                                    
                                    /* humus */
                                    hf[j][0][0] += (flux->soil).hf_gfa;
                                    hf[j][0][1] += (flux->soil).hf_gca;
                                    hf[j][0][2] += (flux->soil).hf_gra;
                                    hf[j][0][3] += (flux->soil).hf_tfa;
                                    hf[j][0][4] += (flux->soil).hf_tca;
                                    hf[j][0][5] += (flux->soil).hf_tra;
                                    hf[j][1][0] += (flux->soil).hf_gfi;
                                    hf[j][1][1] += (flux->soil).hf_gci;
                                    hf[j][1][2] += (flux->soil).hf_gri;
                                    hf[j][1][3] += (flux->soil).hf_tfi;
                                    hf[j][1][4] += (flux->soil).hf_tci;
                                    hf[j][1][5] += (flux->soil).hf_tri;
                                    hf[j][2][0] += (flux->soil).hf_gfp;
                                    hf[j][2][1] += (flux->soil).hf_gcp;
                                    hf[j][2][2] += (flux->soil).hf_grp;
                                    hf[j][2][3] += (flux->soil).hf_tfp;
                                    hf[j][2][4] += (flux->soil).hf_tcp;
                                    hf[j][2][5] += (flux->soil).hf_trp;
                                    hd[j][0] += (echar->soil).sr_ha/1000.0 * (echar->soil).f_tm_h;
                                    hd[j][1] += (echar->soil).sr_hi/1000.0 * (echar->soil).f_tm_h;
                                    hd[j][2] += (echar->soil).sr_hp/1000.0 * (echar->soil).f_tm_h;
                                }
                            }
						}
						f_erosion_rusle(&grid[h][i], loct, &mass[h][i], flux);
						if(loct->time > 100){
							(mass[h][i].soil).msl_a -= (flux->soil).erosion_carbon;
						}
						
						/* fprintf(fp_spinup, "%4d %4d %4ld %4ld ", grid[h][i].dg_row, grid[h][i].dg_col, g, e);
						
						fprintf(fp_spinup, "%8.4lf ", aaa[0]);
						fprintf(fp_spinup, "%8.4lf ", aaa[1]);
						fprintf(fp_spinup, "%8.4lf ", aaa[2]);
						fprintf(fp_spinup, "%8.4lf ", aaa[3]);
						fprintf(fp_spinup, "%8.4lf ", aaa[4]);
						fprintf(fp_spinup, "%8.4lf ", aaa[5]);
						fprintf(fp_spinup, "%8.4lf ", aaa[6]);
						fprintf(fp_spinup, "%8.4lf ", aaa[7]);
						fprintf(fp_spinup, "%8.4lf ", aaa[8]);
						fprintf(fp_spinup, "%8.4lf ", aaa[9]);
						fprintf(fp_spinup, "%8.4lf ", aaa[10]);
						fprintf(fp_spinup, "%8.4lf ", aaa[11]);
						fprintf(fp_spinup, "%8.4lf ", aaa[12]);
						fprintf(fp_spinup, "%8.4lf ", aaa[13]);
						fprintf(fp_spinup, "%8.4lf ", aaa[14]); */
						
						fprintf(fp_spinup,"%ld ", loct->adyear); 
						if(WMODE==2){
							fprintf(fp_spinup,"%4d %4d ", grid[h][i].dg_row, grid[h][i].dg_col); 
						}
						fprintf(fp_spinup,"%4ld %4ld ", g, e);
						fprintf(fp_spinup,"%6.2lf ", gpp_a);
						fprintf(fp_spinup,"%6.2lf ", npp_a);
						fprintf(fp_spinup,"%6.2lf ", nep_a);
						fprintf(fp_spinup,"%6.2lf ", lai_a);
						fprintf(fp_spinup,"%6.2lf ", plant_a);
						fprintf(fp_spinup,"%6.2lf ", soil_a);
                        
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].tree).fol);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].tree).stm);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].tree).rot);

                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].c3).fol);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].c3).stm);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].c3).rot);

                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].c4).fol);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].c4).stm);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].c4).rot);
                        
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).ltr_gf);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).ltr_gc);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).ltr_gr);

                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).ltr_tf);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).ltr_tc);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).ltr_tr);

                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).msl_a);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).msl_i);
                        fprintf(fp_spinup,"%7.3lf ", (mass[h][i].soil).msl_p);
                        
                        //fprintf(fp_spinup,"%7.5lf ", xx);
                        
						printf("AD:%4ld ", loct->adyear); 
						if(WMODE==2){
							printf("%4d %4d ", grid[h][i].dg_row, grid[h][i].dg_col); 
						}
						printf("R:%4ld Y:%4ld ", g, loct->climy);
						printf("GPP:%6.2lf ", gpp_a);
						printf("NPP:%6.2lf ", npp_a);
						printf("NEP:%6.2lf ", nep_a);
						printf("LAI:%6.2lf ", lai_a);
						printf("PLANT:%6.2lf ", plant_a);
						printf("SOIL:%6.2lf ", soil_a);
						printf("\n");
                        
                        /* SASU: 2013/10/29 by A.Ito */
                        if(USE_SASU ==2 && (loct->time>=29 && loct->time<=59) ){
                            /* vegetation */
                            /* average paramepers *****************/
                            for(j=0;j<9;j++){
                                ss_v[j] = 0.0;
                                //nppi[loct->time%10][j] *= 1.0 / nsasu_v[loct->time%10];
                                //lf[loct->time%10][j] *= 1.0 / nsasu_v[loct->time%10];
                            }
 
                            /* estimate steady-state carbon stock *****************/
                            /* vegetation */
                            for(j=0;j<9;j++){
                                ss_v[j] += nppi[(loct->time+1)%10][j];
                                if(lf[loct->time%10][j] > 0.0){
                                    ss_v[j] *= 1.0 / lf[(loct->time+1)%10][j];
                                }else{
                                    ss_v[j] = 0.0;
                                }

                                fprintf(fp_spinup,"%lf %lf %lf ", ss_v[j],nppi[(loct->time+1)%10][j],lf[(loct->time+1)%10][j]);
                            }
                            
                            /* update carbon stock by using SASU-estimated one */    
                            /* vegetation */
                            if(nppi[(loct->time+1)%10][0] > 0.0){
                                (mass[h][i].tree).fol = ss_v[0];
                            }
                            if(nppi[(loct->time+1)%10][1] > 0.0){
                                (mass[h][i].tree).stm = ss_v[1];
                            }
                            if(nppi[(loct->time+1)%10][2] > 0.0){
                                (mass[h][i].tree).rot = ss_v[2];
                            }
                            if(nppi[(loct->time+1)%10][3] > 0.0){
                                (mass[h][i].c3).fol = ss_v[3];
                            }
                            if(nppi[(loct->time+1)%10][4] > 0.0){
                                (mass[h][i].c3).stm = ss_v[4];
                            }
                            if(nppi[(loct->time+1)%10][5] > 0.0){
                                (mass[h][i].c3).rot = ss_v[5];
                            }
                            if(nppi[(loct->time+1)%10][6] > 0.0){
                                (mass[h][i].c4).fol = ss_v[6];
                            }
                            if(nppi[(loct->time+1)%10][7] > 0.0){
                                (mass[h][i].c4).stm = ss_v[7];
                            }
                            if(nppi[(loct->time+1)%10][8] > 0.0){
                                (mass[h][i].c4).rot = ss_v[8];
                            }
                        }
                        
                        if(USE_SASU >=1 && ((loct->time>=39 && loct->time<=69) || loct->time == 99 || loct->time == 139) ){
                        /* if(USE_SASU ==1 && (loct->time==50 || loct->time==80 || loct->time==110 ||
                                            loct->time==140 || loct->time==170)){ */
                            
                            /* average paramepers *****************/
                           
                            /* litter */
                            for(j=0;j<6;j++){
                                ss_l[j] = 0.0;
                                for(k=0;k<9;k++){
                                    //li[loct->time%10][j][k] *= 1.0 / nsasu[loct->time%10];
                                }
                                //ld[loct->time%10][j] *= 1.0 / nsasu[loct->time%10];
                            }
                            /* humus */
                            for(j=0;j<3;j++){
                                ss_h[j] = 0.0;
                                for(k=0;k<6;k++){
                                    //hf[loct->time%10][j][k] *= 1.0 / nsasu[loct->time%10];
                                }
                                //hd[loct->time%10][j] *= 1.0 / nsasu[loct->time%10];
                            }
                            
                            /* estimate steady-state carbon stock *****************/
                            /* litter */
                            for(j=0;j<6;j++){
                                for(k=0;k<9;k++){
                                    ss_l[j] += li[(loct->time+1)%10][j][k];
                                }
                                if(ld[(loct->time+1)%10][j] > 0.0){
                                    ss_l[j] *= 1.0 / ld[(loct->time+1)%10][j];
                                }else{
                                    ss_l[j] = 0.0;
                                }

                                fprintf(fp_spinup,"%7.3lf ", ss_l[j]);
                            }
                            
                            /* humus */
                            for(j=0;j<3;j++){
                                for(k=0;k<6;k++){
                                    ss_h[j] += hf[(loct->time+1)%10][j][k];
                                }
                                if(hd[(loct->time+1)%10][j] > 0.0){
                                    ss_h[j] *= 1.0 / hd[(loct->time+1)%10][j];
                                }else{
                                    ss_h[j] = 0.0;
                                }

                                fprintf(fp_spinup,"%7.3lf ", ss_h[j]);
                            }
                            
                            /* update carbon stock by using SASU-estimated one */                            
                            /* litter */
                            (mass[h][i].soil).ltr_gf = ss_l[0];
                            (mass[h][i].soil).ltr_gc = ss_l[1];
                            (mass[h][i].soil).ltr_gr = ss_l[2];
                            (mass[h][i].soil).ltr_tf = ss_l[3];
                            (mass[h][i].soil).ltr_tc = ss_l[4];
                            (mass[h][i].soil).ltr_tr = ss_l[5];
                            /* humus */
                            (mass[h][i].soil).msl_a = ss_h[0];
                            (mass[h][i].soil).msl_i = ss_h[1];
                            (mass[h][i].soil).msl_p = ss_h[2];
                        }
                        
						fprintf(fp_spinup, "\n"); /* */
					}
				}
			}
			
			/* output for restart file */
			fprintf(fp_restart,"%lf ", mass[h][i].snwa);
			fprintf(fp_restart,"%lf ", mass[h][i].sw30);
			fprintf(fp_restart,"%lf ", mass[h][i].sww);

			fprintf(fp_restart,"%lf ", (mass[h][i].tree).fol);
			fprintf(fp_restart,"%lf ", (mass[h][i].tree).stm);
			fprintf(fp_restart,"%lf ", (mass[h][i].tree).rot);
			fprintf(fp_restart,"%lf ", (mass[h][i].c3).fol);
			fprintf(fp_restart,"%lf ", (mass[h][i].c3).stm);
			fprintf(fp_restart,"%lf ", (mass[h][i].c3).rot);
			fprintf(fp_restart,"%lf ", (mass[h][i].c4).fol);
			fprintf(fp_restart,"%lf ", (mass[h][i].c4).stm);
			fprintf(fp_restart,"%lf ", (mass[h][i].c4).rot);

			fprintf(fp_restart,"%lf ", (mass[h][i].soil).ltr_tf);
			fprintf(fp_restart,"%lf ", (mass[h][i].soil).ltr_tc);
			fprintf(fp_restart,"%lf ", (mass[h][i].soil).ltr_tr);
			fprintf(fp_restart,"%lf ", (mass[h][i].soil).ltr_gf);
			fprintf(fp_restart,"%lf ", (mass[h][i].soil).ltr_gc);
			fprintf(fp_restart,"%lf ", (mass[h][i].soil).ltr_gr);
			fprintf(fp_restart,"%lf ", (mass[h][i].soil).msl_a);
			fprintf(fp_restart,"%lf ", (mass[h][i].soil).msl_i);
			fprintf(fp_restart,"%lf\n", (mass[h][i].soil).msl_p);
		}
	}
	
	fclose(fp_restart);	

	if(NOTICE==1){
		printf("done\n");
	}
}

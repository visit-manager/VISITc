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
	long e, f, g, h, i;
	char filename[100];
	long ndy, nroop;
	double gpp_a, npp_a, nep_a, lai_a, plant_a, soil_a;
	double aaa[20];
	FILE *fp_restart;
	
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
					
					for(e=BYR ; e<=EYR ; e++){	
						loct->climy = e;
						loct->CO2y = BYR;
						loct->time++;
						
						/* stand age, year */
						loct->age_stand += 1.0;
				
						loct->adyear = ((BYR-1) - PERIOD*(nroop - g) - (EYR - e));
						
						/* ndy = (e%4==0)?366:365; */
						ndy = 365;
						/** initialize climate conditions **/
						gpp_a=npp_a=nep_a=lai_a=plant_a=soil_a = 0.0;
						for(f=0;f<20;f++){
							aaa[f] = 0.0;
						}
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
						
						fprintf(fp_spinup,"%4ld ", loct->adyear); 
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
						fprintf(fp_spinup, "\n"); /* */
						
						printf("AD:%4ld ", loct->adyear); 
						if(WMODE==2){
							printf("%4d %4d ", grid[h][i].dg_row, grid[h][i].dg_col); 
						}
						printf("R:%4ld Y:%4ld ", g, e);
						printf("GPP:%6.2lf ", gpp_a);
						printf("NPP:%6.2lf ", npp_a);
						printf("NEP:%6.2lf ", nep_a);
						printf("LAI:%6.2lf ", lai_a);
						printf("PLANT:%6.2lf ", plant_a);
						printf("SOIL:%6.2lf ", soil_a);
						printf("\n");
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

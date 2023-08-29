/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
*/

/* header files */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* climatology: *_d[] means the time-series during 1948 to 2007- */	
extern float tmp_sfc_d[PERIOD][366];		/* ground surface temperature, degree Celcius */
extern float tmp_2m_d[PERIOD][366];			/* 2m air temperature, degree Celcius */
extern float tmp10_soil_d[PERIOD][366];		/* soil temperature at 10 cm depth, degree Celcius */
extern float tmp200_soil_d[PERIOD][366];	/* soil temperature at 200 cm depth, degree Celcius */
extern float dswrf_sfc_d[PERIOD][366];		/* downward shortwave radiation at the surface, W m-2 */
extern float tcdc_clm_d[PERIOD][366];		/* total cloudiness, fraction */
extern float prate_sfc_d[PERIOD][366];		/* precipitation, mm mon-1 */
extern float spfh_2m_d[PERIOD][366];		/* specific humidity, kg kg-1 */
extern float wind_10m_d[PERIOD][366];		/* u-wind velocity, m s-1 */
extern float vpd_d[PERIOD][366];			/* VPD, hPa */
extern float paddy_wtd[366];

extern long	month_day[12];
extern long	WMODE;
extern long	WGRIDS;

/* GHG scenario ***************************************/
/* source: http://crga.atmos.uiuc.edu/research/post-sres.html
 M.E.Schlesinger and S.Malyshev			*/
/* atmospheric CO2, ppmv */
extern float		atm_co2_a1[553];	/* SRES A1 */
extern float		atm_co2_a2[553];	/* SRES A2 */
extern float		atm_co2_b1[553];	/* SRES B1 */
extern float		atm_co2_b2[553];	/* SRES B2 */
/* atmospheric CH4, pptv*/
extern float		atm_ch4_a1[553];	/* SRES A1 */
extern float		atm_ch4_a2[553];	/* SRES A2 */
extern float		atm_ch4_b1[553];	/* SRES B1 */
extern float		atm_ch4_b2[553];	/* SRES B2 */
/* atmospheric N2O, pptv*/
extern float		atm_n2o_a1[553];	/* SRES A1 */
extern float		atm_n2o_a2[553];	/* SRES A2 */
extern float		atm_n2o_b1[553];	/* SRES B1 */
extern float		atm_n2o_b2[553];	/* SRES B2 */

/******************************************************************************/
void f_initialize(
	struct Grid grid[], 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Echar echar_type[MAX_BIOME], 
	struct Mass mass[], 
	struct Flux *flux, 
	FILE *fp_r[N_OFILE]
){
	char filename[100];
	long h, i, yr, day;
	FILE *fp_ghg, *fp_wtd;
	
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
	for(i=0; i<WGRIDS; i++){
		grid[i] = grid[0];
	}
	
	/** clear all parameters **********/
	if(NOTICE == 1){
		printf("Clearing variables...");
	}
	for(i=0; i<WGRIDS; i++){
		clear_a(&(grid[i]), &(mass[i]));
	}
	clear_b(loct, echar, flux);	
	if(NOTICE == 1){
		printf("done\n");
	}

	/******************************************************************/	
	/** initialize grid condition **/
	if(NOTICE == 1){
		printf("Initializing site...\n");
	}
	switch(WMODE){
		case 1:
			/* point */
			init_site(&grid[0]); 
			break;
		case 2: case 3:
			/* region */
            if(strcmp(grid->site_id, "GLOBAL")==0){
                f_init_global_run(grid);
            }else if(strcmp(grid->site_id, "BAMIYAN")==0){
                f_init_bamiyan_run(grid);
            }else if(strcmp(grid->site_id, "EASIA")==0){
                f_init_easia_run(grid);
            }else if(strcmp(grid->site_id, "JAPAN")==0){
                f_init_japan_run(grid);
            }else if(strcmp(grid->site_id, "BB")==0){
                f_init_bb_run(grid);
            }else if(strcmp(grid->site_id, "JAPANc")==0){
                f_init_japanc_run(grid);
            }else if(strcmp(grid->site_id, "JAPANh")==0){
                f_init_japanh_run(grid);
            }else if(strcmp(grid->site_id, "JAPANk")==0){
                f_init_japank_run(grid);
            }else if(strcmp(grid->site_id, "PAWCs")==0){
                f_init_pawcs_run(grid);
            }
			break;
		default:
			break;
	}
	if(NOTICE == 1){
		printf("done\n");
	}
		
	/******************************************************************/	
	if(NOTICE == 1){
		printf("Reading GHG scenario...");
	}
	if((fp_ghg=fopen("./data/AtmGHG_timeseries.txt","rt"))==NULL){
		printf("No atmospheric GHG scenario!\n");
		exit(1);
	}
	for(h=0;h<553;h++){
		fscanf(fp_ghg,"%ld", &yr);
		fscanf(fp_ghg,"%f", &atm_co2_a1[h]);
		fscanf(fp_ghg,"%f", &atm_co2_a2[h]);
		fscanf(fp_ghg,"%f", &atm_co2_b1[h]);
		fscanf(fp_ghg,"%f", &atm_co2_b2[h]);
		fscanf(fp_ghg,"%f", &atm_ch4_a1[h]);
		fscanf(fp_ghg,"%f", &atm_ch4_a2[h]);
		fscanf(fp_ghg,"%f", &atm_ch4_b1[h]);
		fscanf(fp_ghg,"%f", &atm_ch4_b2[h]);
		fscanf(fp_ghg,"%f", &atm_n2o_a1[h]);
		fscanf(fp_ghg,"%f", &atm_n2o_a2[h]);
		fscanf(fp_ghg,"%f", &atm_n2o_b1[h]);
		fscanf(fp_ghg,"%f", &atm_n2o_b2[h]);
	}
	fclose(fp_ghg);
	if(NOTICE == 1){
		printf("done\n");
	}	

	/***************************************************************/
	if(NOTICE == 1){
		printf("Reading N deposition scenario...");
	}
	
	/* initialization */
	for(i=0;i<366;i++){
		loct->depo_nh4_model_av = 0.0;
		loct->depo_no3_model_av = 0.0;
	}
	
    /* site simulation */
    
	if(NOTICE == 1){
		printf("done\n");
	}
	
	/* water-table depth of paddy field: MASE ***************************/
	if(WH_CH4 == 1 && strcmp(grid->site_id, "MSE")==0){
		if((fp_wtd = fopen("paddy_wtd.dat","rt"))==NULL){
			printf("No water-table depth data!\n");
			exit(1);
		}
		for(i=0;i<366;i++){
			fscanf(fp_wtd,"%ld %f", &day, &paddy_wtd[i]);
		}
	}

	/* ******************************************************/
	
	/* *****************************************************************/
	/** initialize vegatation and soil conditions **/ 
	if(NOTICE == 1){
		printf("Reading parameters...");
	}
	/*** parameters ***/
	/* strcpy(filename, "./data/parameter_S1b.dat"); */ /* S1 */
	/* strcpy(filename, "./data/parameter_CEAMIP.txt"); */ /* 2011/11/03 by A.Ito */
    
    /* 2012/02/07 by A.Ito */
    strcpy(filename, grid->file_para);
    
	switch(WMODE){
		case 1:
			set_parameter(1, grid[0].file_para, echar);
			set_parameter(2, filename, echar_type);
			break;
		case 2: case 3:
			set_parameter(2, filename, echar_type);
			break;
	}
	
	if(NOTICE == 1){
		printf("done\n");
	}
	
	/* mass **************************************************************/
	if(NOTICE == 1){
		printf("Initializing vegetation and soil mass...");
	}
	/*** initialize ***/
	for(i=0;i<WGRIDS;i++){
		grid[i].age_stand = 0.0;

		switch(WMODE){
			case 1:
				break;
			case 2: 
				echar = &(echar_type[grid[i].veg_type]);
				break;
		}
        
        /* initialize mass *************/
        if(grid[i].flag_datavl == 1){
            if((echar->tree).phenoltype != 0){
                /* forests */
                initTree(&(mass[i]));
            }
            initFloor(&(mass[i]));
            initSoil(&(mass[i]));
        }
	}
	if(NOTICE == 1){
		printf("done\n");
	}

	/* initialize location conditions **********************/
	if(NOTICE == 1){
		printf("Initializing local environment...\n");
	}
	for(i=0; i<WGRIDS; i++){
		
		/* initialize location *************/
		if(grid[i].flag_datavl == 1){
			f_init_loct(&grid[i], loct, echar, &mass[i], flux);	
			if(NOTICE == 1){
				if(i%1000==0) printf("%ld: %f %f %f\n", i, mass[i].sw30, mass[i].sww, mass[i].snwa);
			}
		}
	}
	if(NOTICE == 1){
		printf("done\n");
	}
	
	/* initialize stable carbon isotope **********************/ 
	if(NOTICE == 1){
		printf("Initializing stable carbon isotope...");
	}
	switch(WMODE){
		case 1:
			init_d13c(&grid[0], loct, echar, &mass[0], flux);
			break;
		case 2: case 3: 
			for(i=0;i<WGRIDS;i++){
				init_d13c(&grid[i], loct, echar, &mass[i], flux);
			}
			break;
		default:
			break;
	}
	if(NOTICE == 1){
		printf("done\n");
	}

	/******************************************************************/	
	if(NOTICE == 1){
		printf("Creating result files...");
	}
	switch(WMODE){
		case 1: /* point */
			for(h=0;h<N_OFILE;h++){
				strcpy(filename, grid[0].site_id);
				strcat(filename, grid[0].file_name);
				
				switch(h){
					case 0:		strcat(filename, "_spinup");		break;
					case 1:		strcat(filename, "_hourly");		break;
					case 2:		strcat(filename, "_daily");			break;
					case 3:		strcat(filename, "_monthly");		break;
					case 4:		strcat(filename, "_annual");		break;
					case 5:		strcat(filename, "_ansis1");	break;
					case 6:		strcat(filename, "_ansis2");		break;
					case 7:		strcat(filename, "_ansis3");			break;
					case 8:		strcat(filename, "_ansis4");			break;
					case 9:		strcat(filename, "_ansis5");		break;
					case 10:	strcat(filename, "_99-30min");		break;
					case 11:	strcat(filename, "_00-30min");		break;
					case 12:	strcat(filename, "_01-30min");		break;
					case 13:	strcat(filename, "_02-30min");		break;
					case 14:	strcat(filename, "_03-30min");		break;
					case 15:	strcat(filename, "_04-30min");		break;
					case 16:	strcat(filename, "_05-30min");		break;
				}
				strcat(filename, ".txt");
				fp_r[h] = fopen(filename,"wt");
			}
			break;
		case 2: case 3:  /* region */
			
            for(h=0;h<N_OFILE;h++){
                strcpy(filename, grid[0].site_id);
                strcat(filename, grid[0].file_name);
                switch(h){
                    case 0:		strcat(filename, "_spinup");		break;
                    case 1:		strcat(filename, "_hourly");		break;
                    case 2:		strcat(filename, "_daily");		break;
                    case 3:		strcat(filename, "_monthly");		break;
                    case 4:		strcat(filename, "_annual");		break;
                        
                    case 5:		strcat(filename, "_ansis1");		break;
                    case 6:		strcat(filename, "_ansis2");		break;
                    case 7:		strcat(filename, "_ansis3");		break;
                    case 8:		strcat(filename, "_ansis4");		break;
                    case 9:		strcat(filename, "_ansis5");		break;
                    case 10:	strcat(filename, "_04-daily");		break;
                    case 11:	strcat(filename, "_05-daily");		break;
                    case 12:	strcat(filename, "_06-daily");		break;
                    case 13:	strcat(filename, "_07-daily");		break;
                    case 14:	strcat(filename, "_dat1");			break;
                    case 15:	strcat(filename, "_dat2");			break;
                    case 16:	strcat(filename, "_dat3");			break;
                }
                strcat(filename, ".txt");
                fp_r[h] = fopen(filename,"wt");
            }
			
			break;
	}
	if(NOTICE == 1){
		printf("done\n");
	}
}

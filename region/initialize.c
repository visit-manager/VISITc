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
extern float paddy_wtd[366];

extern long	WMODE;
extern long	WGRIDS;
extern struct Echar echar_type[MAX_BIOME];

/* GHG scenario ***************************************/
/* source: http://crga.atmos.uiuc.edu/research/post-sres.html
 M.E.Schlesinger and S.Malyshev			*/
/* atmospheric CO2, ppmv */
extern float		atm_co2[4][553];	/* CO2 scenarios */
/* atmospheric CH4, pptv*/
extern float		atm_ch4[4][553];	/* CH4 scenarios */
/* atmospheric N2O, pptv*/
extern float		atm_n2o[4][553];	/* n2O scenarios */

/* *****************************************************************************/
void f_initialize(
	struct Grid grid[], 
	struct Loct loct[],
	struct Echar echar[],
	struct Mass mass[],
	struct Flux flux[]
){
	char filename[128];
	long i, yr, day;
	FILE *fp_ghg, *fp_wtd;
	
	/* Number of days for each month */
	/* month_day[0] = 31; */	/* January */
    /* month_day[1] = 28; */	/* February */
    /* month_day[2] = 31; */	/* March */
    /* month_day[3] = 30; */	/* April */
    /* month_day[4] = 31; */	/* May */
    /* month_day[5] = 30; */	/* June */
    /* month_day[6] = 31; */	/* July */
    /* month_day[7] = 31; */	/* August */
    /* month_day[8] = 30; */	/* September */
    /* month_day[9] = 31; */	/* October */
    /* month_day[10] = 30; */	/* November */
    /* month_day[11] = 31; */	/* December */

	/* *****************************************************************/	
	for(i=0; i<WGRIDS; i++){
		grid[i] = grid[0];
	}
	
	/* clear all parameters **********/
	if(NOTICE == 1){
		printf("Clearing variables...");
	}
	for(i=0; i<WGRIDS; i++){
		f_clear_a(&(grid[i]), &(mass[i]));
        f_clear_b(&(loct[i]), &(echar[i]), &(flux[i]));
	}
	
	if(NOTICE == 1){
		printf("done\n");
	}

	/* *****************************************************************/
	/* initialize grid condition **/
	if(NOTICE == 1){
		printf("Initializing site...\n");
	}
	switch(WMODE){
		case 1:
			/* point */
			/* init_site(&grid[0]);  */
            printf("BAD WMODE (site)\n");
            exit(1);
			break;
		case 2: case 3:
			/* region */
            if(strcmp(grid[0].area_id, "GLOBAL")==0){
                f_init_global_run(grid);
            }else if(strcmp(grid[0].area_id, "BAMIYAN")==0){
                f_init_bamiyan_run(grid);
            }else if(strcmp(grid[0].area_id, "EASIA")==0){
                f_init_easia_run(grid);
            }else if(strcmp(grid[0].area_id, "JAPAN")==0){
                f_init_japan_run(grid);
            }else if(strcmp(grid[0].area_id, "BB")==0){
                f_init_bb_run(grid);
            }else if(strcmp(grid[0].area_id, "JAPANc")==0){
                f_init_japanc_run(grid);
            }else if(strcmp(grid[0].area_id, "JAPANh")==0){
                f_init_japanh_run(grid);
            }else if(strcmp(grid[0].area_id, "JAPANk")==0){
                f_init_japank_run(grid);
            }else if(strcmp(grid[0].area_id, "PAWCs")==0){
                f_init_pawcs_run(grid);
            }
			break;
		default:
			break;
	}
	if(NOTICE == 1){
		printf("done\n");
	}
		
	/* *****************************************************************/
	if(NOTICE == 1){
		printf("Reading GHG scenario...");
	}
	if((fp_ghg=fopen("./data/AtmGHG_timeseries.txt","rt"))==NULL){
		printf("No atmospheric GHG scenario!\n");
		exit(1);
	}
	for(i=0;i<553;i++){
		fscanf(fp_ghg,"%ld", &yr);
		fscanf(fp_ghg,"%f", &atm_co2[0][i]);
		fscanf(fp_ghg,"%f", &atm_co2[1][i]);
		fscanf(fp_ghg,"%f", &atm_co2[2][i]);
		fscanf(fp_ghg,"%f", &atm_co2[3][i]);
		fscanf(fp_ghg,"%f", &atm_ch4[0][i]);
		fscanf(fp_ghg,"%f", &atm_ch4[1][i]);
		fscanf(fp_ghg,"%f", &atm_ch4[2][i]);
		fscanf(fp_ghg,"%f", &atm_ch4[3][i]);
		fscanf(fp_ghg,"%f", &atm_n2o[0][i]);
		fscanf(fp_ghg,"%f", &atm_n2o[1][i]);
		fscanf(fp_ghg,"%f", &atm_n2o[2][i]);
		fscanf(fp_ghg,"%f", &atm_n2o[3][i]);
	}
	fclose(fp_ghg);
	if(NOTICE == 1){
		printf("done\n");
	}	

	/* **************************************************************/
	if(NOTICE == 1){
		printf("Reading N deposition scenario...");
	}
	
	/* initialization */
    for(i=0; i<WGRIDS; i++){
        loct[i].depo_nh4_model_av = 0.0;
        loct[i].depo_no3_model_av = 0.0;
    }
	
    /* site simulation */
    
	if(NOTICE == 1){
		printf("done\n");
	}
	
	/* water-table depth of paddy field: MASE ***************************/
	if(WH_CH4 == 1 && strcmp(grid[0].area_id, "MSE")==0){
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
	/* initialize vegatation and soil conditions **/
	if(NOTICE == 1){
		printf("Reading parameters...");
	}
	/*** parameters ***/
	/* strcpy(filename, "./data/parameter_S1b.dat"); */ /* S1 */
	/* strcpy(filename, "./data/parameter_CEAMIP.txt"); */ /* 2011/11/03 by A.Ito */
    
    /* 2012/02/07 by A.Ito */
    strcpy(filename, grid[0].file_para);
    
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
				echar[i] = echar_type[grid[i].veg_type];
				break;
		}
        
        /* initialize mass *************/
        if(grid[i].flag_datavl == 1){
            if((echar[i].tree).phenoltype != 0){
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
			f_init_loct(&grid[i], &loct[i], &echar[i], &mass[i], &flux[i]);
			if(NOTICE == 1){
				if(i%1000 == 0) printf("%ld: %f %f %f\n", i, mass[i].sw30, mass[i].sww, mass[i].snwa);
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
			init_d13c(&grid[0], &loct[0], &echar[0], &mass[0], &flux[0]);
			break;
		case 2: case 3: 
			for(i=0;i<WGRIDS;i++){
				init_d13c(&grid[i], &loct[i], &echar[i], &mass[i], &flux[i]);
			}
			break;
		default:
			break;
	}
	if(NOTICE == 1){
		printf("done\n");
	}
}

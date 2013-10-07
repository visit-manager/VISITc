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

void f_setting(
	struct Grid *grid
){
	char varread[16];
	char smode[] = "MODE";
	char sfsfx[] = "FSFX";
	char ssite[] = "SITE";
	char sreg[] = "REG";
	char syrdist[] = "LDYEAR";
	char sbyear[] = "BYEAR";
	char seyear[] = "EYEAR";
	char sspupt[] = "SPUPT";
	char sclimf[] = "CLIMF";
	char ssitef[] = "SITEF";
	char sparaf[] = "PARAF";
	double dist_min;
	FILE *fp_set;
	
	WMODE = 0;
	
	if(NOTICE==1){
		printf("Setting ");
	}

	if((fp_set = fopen("setting.txt","rt")) == NULL){
		printf("No setting data !\n");
		exit (1);
	}
	
	/* MODE:  working mode */
	fscanf(fp_set,"%s", &varread);
	if(strcmp(varread, smode)==0){
		fscanf(fp_set,"%ld", &WMODE);
	}
	
	/* file suffix string */
	fscanf(fp_set,"%s", &varread);
	if(strcmp(varread, sfsfx)==0){
		fscanf(fp_set,"%s", &grid->file_name);
	}
	
	dist_min = 999.999;
	YLDIST = 0;
	/*  disturbance year */
	/* Takayama */
	/* #define YR_CUT 1964 */
	/* S1 region */
	/* #define YR_CUT 1965  */
	/* #define YR_CUT 1955  */
	/* #define YR_CUT 1965  */
	/* #define YR_CUT 1975  */
	switch(WMODE){
		case 1: /* site ********************************/
			/* site ID */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, ssite)==0){
				fscanf(fp_set,"%s", &grid->site_id);
			}
			/* spin-up time */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sspupt)==0){
				fscanf(fp_set,"%ld", &SPUPT);
			}
			/* begin year */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sbyear)==0){
				fscanf(fp_set,"%ld", &BYR);
			}
			/* end year */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, seyear)==0){
				fscanf(fp_set,"%ld", &EYR);
			}
			/* last disturbance year */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, syrdist)==0){
				fscanf(fp_set,"%ld", &YLDIST);
				/* printf("DEBUG: %ld\n", EYR); */
			}
			
			PERIOD = EYR - BYR + 1;

			/* climate file name */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sclimf)==0){
				fscanf(fp_set,"%s", &grid->file_clim);
			}
			/* site description filename */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, ssitef)==0){
				fscanf(fp_set,"%s", &grid->file_site);
			}
			/* parameter file name */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sparaf)==0){
				fscanf(fp_set,"%s", &grid->file_para);
			}

			NCOL = 1;
			NROW = 1;
			WROW = 0;
			WCOL = 0;
			break;
		case 2: /* region **********************************/
			/* study area */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sreg)==0){
				fscanf(fp_set,"%ld %ld", &WROW, &WCOL);
				fscanf(fp_set,"%ld %ld", &NROW, &NCOL);
			}
			/* last disturbance year */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, syrdist)==0){
				fscanf(fp_set,"%ld", &YLDIST);
			}
			/* spin-up time */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sspupt)==0){
				fscanf(fp_set,"%ld", &SPUPT);
			}
			break;
	}
	
	if(NROW>DROW){
		printf("! Wrong NROW: %ld\n", NROW);
		exit (1);
	}
	if(NCOL>DROW){
		printf("! Wrong NCOL: %ld\n", NCOL);
		exit (1);
	}
		
	if(NOTICE==1){
		printf("Mode:%ld Grid:%ld-%ld Site:%s", WMODE, 
			WROW, WCOL, grid->site_id);
	}
	
	if(NOTICE==1){
		printf("...done\n");
	}
	
	/* close setting.dat */
	fclose(fp_set);
}

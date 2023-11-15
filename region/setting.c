/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/*  [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 07 by Akihiko Ito: Version 0.2.0.
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

extern long WMODE;
extern long WGRIDS;
extern long YLDIST;
extern long SPUPT;

/* ************************************************/
void f_setting(
	struct Grid *grid
){
	char varread[10];
	char smode[] = "MODE";
	char sfsfx[] = "FSFX";
	char ssite[] = "SITE";
	char sregion[] = "REGION";
	char sreg[] = "GRIDS";
	char syrdist[] = "YLDIST";
	char sspupt[] = "SPUPT";
	char sncepf[] = "NCEPF";
	char ssitef[] = "SITEF";
	char sparaf[] = "PARAF";
	float dist_min;
	FILE *fp_set;
	
	WMODE = 0;
	
	if(NOTICE == 1){
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
	
	/*  disturbance year */
	dist_min = 999.999;
	YLDIST = 0;
	switch(WMODE){
		case 1: /* site ********************************/
			/* site ID */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, ssite)==0){
				fscanf(fp_set,"%s", &grid->area_id);
			}
			printf("Site ID: %s\n", grid->area_id);
			
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

			/* NCEP climate file name */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sncepf)==0){
				fscanf(fp_set,"%s", &grid->file_ncep);
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

			WGRIDS = 1;
			break;
		case 2: case 3: /* region **********************************/
			/* region ID */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sregion)==0){
				fscanf(fp_set,"%s", &grid->area_id);
			}
			printf("Region ID: %s\n", grid->area_id);
			/* strcpy(grid->area_id, "AsiaMIP"); */

			/* study area: number of grids */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sreg)==0){
				fscanf(fp_set,"%ld", &WGRIDS);
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
            
            /* parameter file name: added by A.Ito (2012/02/07) */
			fscanf(fp_set,"%s", &varread);
			if(strcmp(varread, sparaf)==0){
				fscanf(fp_set,"%s", &grid->file_para);
			}
			
			break;
	}
	
	if(WGRIDS > (NROW*NCOL)){
		printf("! Too much grids: %ld\n", WGRIDS);
		exit (1);
	}
		
	if(NOTICE==1){
		printf("Mode:%ld Grids:%ld ", WMODE, WGRIDS);
	}
	
	if(NOTICE==1){
		printf("...done\n");
	}
	
	/* close setting.dat */
	fclose(fp_set);
}

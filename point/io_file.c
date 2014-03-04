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

/* IO output files *****************************************************/
void f_open_files(
	struct Grid *grid, 
    short rnum, 
	FILE *fp_r[NFILE]
){
	long h;
	char filename[100], num[4];
	
	/******************************************************************/	
	if(NOTICE==1){
		printf("Creating result files...");
	}
	switch(WMODE){
		case 1:
			for(h=0;h<NFILE;h++){
				strcpy(filename, grid->site_id);
				strcat(filename, grid->file_name);
				
				if(SENSANS_DIST==1){
					strcat(filename, "_no-distrb");
				}
				if(SENSANS_ACO2==1){
					strcat(filename, "_fix-co2");
				}
				if(SENSANS_TEMP==1){
					strcat(filename, "_fix-temp");
				}else if(SENSANS_TEMP==2){
					strcat(filename, "_lars-temp");
					strcat(filename, num);
				}
				if(SENSANS_PREC==1){
					strcat(filename, "_fix-prec");
				}else if(SENSANS_PREC==2){
					strcat(filename, "_lars-prec");
					strcat(filename, num);
				}
				if(SENSANS_SRAD==1){
					strcat(filename, "_fix-srad");
				}else if(SENSANS_SRAD==2){
					strcat(filename, "_lars-srad");
					strcat(filename, num);
				}
				
				if(SA_PARA==1){
					strcat(filename, "_sa-ptopt");
				}else if(SA_PARA==2){
					strcat(filename, "_sa-pco2");
				}else if(SA_PARA==3){
					strcat(filename, "_sa-psw");
				}else if(SA_PARA==4){
					strcat(filename, "_sa-plue");
				}else if(SA_PARA==5){
					strcat(filename, "_sa-rq10");
				}else if(SA_PARA==6){
					strcat(filename, "_sa-dtemp");
				}else if(SA_PARA==7){
					strcat(filename, "_sa-dsw");
				}else if(SA_PARA==8){
					strcat(filename, "_sa-dist");
				}
				
				if(SA_PARA!=0){
					snprintf(num, 4, "%03d", rnum);
					strcat(filename, num);
				}
				
				if(PREDICT==1 || PREDICT==100){
					strcat(filename, "_pred-miroc-a1b");
				}
				if(PREDICT==10 || PREDICT==101){
					strcat(filename, "_pred-mpi-a1b");
				}
				if(PREDICT==20 || PREDICT==102){
					strcat(filename, "_pred-cccma-a1b");
				}
				if(PREDICT==21 ||  PREDICT==103){
					strcat(filename, "_pred-cccma-a2");
				}
				if(PREDICT==22 || PREDICT==104){
					strcat(filename, "_pred-cccma-b1");
				}
				
				if(F_DOWNREG==1){
					strcat(filename, "_dreg1");
				}
				
				if(FIX_CLIM==1){
					strcat(filename, "_fclim");
				}
                
                if(DGSOIL==1){
					strcat(filename, "_dgsoil");
				}
            
                if(LCCONV==1){
					strcat(filename, "_lcconv");
				}

                if(EX_N_ECOPHYS==1){
					strcat(filename, "_nex1");
				}
                if(EX_N_ECOPHYS==2){
					strcat(filename, "_nex2");
				}
            
				switch(h){
					case 0:		strcat(filename, "_spinup");		break;
					case 1:		strcat(filename, "_annual");		break;
					case 2:		strcat(filename, "_daily");			break;
					case 3:		strcat(filename, "_nee30min");		break;
					case 4:		strcat(filename, "_needaily");		break;
					case 5:		strcat(filename, "_needaytime");	break;
					case 6:		strcat(filename, "_canopy");		break;
					case 7:		strcat(filename, "_ch4");			break;
					case 8:		strcat(filename, "_n2o");			break;
					case 9:		strcat(filename, "_98-30min");		break;
					case 10:	strcat(filename, "_99-30min");		break;
					case 11:	strcat(filename, "_00-30min");		break;
					case 12:	strcat(filename, "_01-30min");		break;
					case 13:	strcat(filename, "_02-30min");		break;
					case 14:	strcat(filename, "_03-30min");		break;
					case 15:	strcat(filename, "_04-30min");		break;
					case 16:	strcat(filename, "_05-30min");		break;
				}
				
				strcat(filename, ".txt");
				/* if(h==1){
				 strcat(filename, ".csv");
				 }else{
				 strcat(filename, ".txt");
				 } */
				
				fp_r[h] = fopen(filename,"wt");
			}
			break;
		case 2: 
			for(h=0;h<NFILE;h++){
				strcpy(filename, grid->site_id);
				if(WMODE==3){
					strcat(filename, "-s1");
				}
				strcat(filename, grid->file_name);
				switch(h){
					case 0:		strcat(filename, "_spinup");		break;
					case 1:		strcat(filename, "_annual");		break;
					case 2:		strcat(filename, "_daily");			break;
					case 3:		strcat(filename, "_dat1");			break;
					case 4:		strcat(filename, "_dat2");			break;
					case 5:		strcat(filename, "_dat3");			break;
					case 6:		strcat(filename, "_dat4");			break;
					case 7:		strcat(filename, "_ch4");			break;
					case 8:		strcat(filename, "_n2o");			break;
					case 9:		strcat(filename, "_98-daily");		break;
					case 10:	strcat(filename, "_99-daily");		break;
					case 11:	strcat(filename, "_00-daily");		break;
					case 12:	strcat(filename, "_01-daily");		break;
					case 13:	strcat(filename, "_02-daily");		break;
					case 14:	strcat(filename, "_03-daily");		break;
					case 15:	strcat(filename, "_04-daily");		break;
					case 16:	strcat(filename, "_05-daily");		break;
				}
				strcat(filename, ".txt");
				fp_r[h] = fopen(filename,"wt");
			}
			break;
	}
	if(NOTICE==1){
		printf("done\n");
	}
	
}

/******************************************************/
void f_close_files(
	FILE *fp_r[NFILE]
){
	long h;
	
	if(NOTICE==1){
		printf("Close output files...");
	}
	
	for(h=0;h<NFILE;h++){
		fclose(fp_r[h]); 
	}

	if(NOTICE==1){
		printf("done\n");
	}
	
}

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
            strcat(filename, "_sa-dst");
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
        
        if(EX_DGSOIL==1){
            strcat(filename, "_dgsoil");
        }
    
        if(EX_LCCONV==1){
            strcat(filename, "_EX_LCCONV");
        }

        if(EX_N_ECOPHYS==1){
            strcat(filename, "_nex1");
        }
        if(EX_N_ECOPHYS==2){
            strcat(filename, "_nex2");
        }
        
        if(EX_SPINUP==1){
            strcat(filename, "_su-1");
        }
        if(EX_SPINUP==2){
            strcat(filename, "_su-2");
        }
        if(EX_SPINUP==3){
            strcat(filename, "_su-3");
        }
        if(EX_SPINUP==4){
            strcat(filename, "_su-4");
        }
        if(EX_SPINUP==5){
            strcat(filename, "_su-5");
        }
        
        /* LUC-MIP: 2014/11/11 */
        if(EX_LUCMIP==11){
           strcat(filename, "_lucmip-1a");
        }
        if(EX_LUCMIP==12){
           strcat(filename, "_lucmip-1b");
        }
        if(EX_LUCMIP==13){
           strcat(filename, "_lucmip-1c");
        }
        if(EX_LUCMIP==2){
           strcat(filename, "_lucmip-2");
        }
        if(EX_LUCMIP==3){
           strcat(filename, "_lucmip-3");
        }
        if(EX_LUCMIP==4){
           strcat(filename, "_lucmip-4");
        }
        if(EX_LUCMIP==5){
           strcat(filename, "_lucmip-5");
        }
        if(EX_LUCMIP==6){
           strcat(filename, "_lucmip-6");
        }
        if(EX_LUCMIP==7){
           strcat(filename, "_lucmip-7");
        }
    
        switch(h){
            case 0:		strcat(filename, "_E01");		break;
            case 1:		strcat(filename, "_E02");		break;
            case 2:		strcat(filename, "_E03");		break;
            case 3:		strcat(filename, "_E04");		break;
            case 4:		strcat(filename, "_E05");		break;
            case 5:		strcat(filename, "_E06");       break;
            case 6:		strcat(filename, "_E07");		break;
            case 7:		strcat(filename, "_E08");		break;
            case 8:		strcat(filename, "_E09");		break;
            case 9:		strcat(filename, "_E10");		break;
            case 10:	strcat(filename, "_E11");		break;
            case 11:	strcat(filename, "_E12");		break;
            case 12:	strcat(filename, "_E13");		break;
            case 13:	strcat(filename, "_E14");		break;
            case 14:	strcat(filename, "_E15");		break;
            case 15:	strcat(filename, "_E16");		break;
            case 16:	strcat(filename, "_E17");		break;
        }
        
        strcat(filename, ".csv");
        /* strcat(filename, ".txt"); */
        
        /* if(h==1){
        strcat(filename, ".txt");
         }else{
         strcat(filename, ".txt");
         } */
        
        fp_r[h] = fopen(filename,"wt");
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

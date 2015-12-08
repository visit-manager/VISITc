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

/* parameter perturbation ************************************/
void f_parameter_perturbation(
	short rnum
){
	/* sensitivity analysis */
	extern short SENSANS_DIST;
	/* 1: no disturbance */
	extern short SENSANS_ACO2;
	/* 1: no atmospheric CO2 rise */
	extern short SENSANS_TEMP;
	/* 1: no temperature change */
	/* 2: use LARS-generated dammy data */
	extern short SENSANS_PREC;
	/* 1: no precipitation change */
	/* 2: use LARS-generated dammy data */
	extern short SENSANS_SRAD;
	/* 1: no solar radiation change */
	/* 2: use LARS-generated dammy data */
	extern short N_LARS;
	
	/* parameter sensitivity analysis */
	extern short SA_PARA;
	/* 1: photosyntehtic optimal temperature */
	/* 2: photosyntehtic CO2 - Michaelis constant */
	/* 3: photosynthetic soil water dependence */
	/* 4: photosyntehtic light use efficiency */
	/* 5: maintenance respiration temperature */
	/* 6: decomposition temperature dependence */
	/* 7: decomposition soil water dependence */
	/* 8: disturbance intensity */
	
	/* ensemble number */
	extern short SA_PARA_EN;
	extern double SA_PARA_VAR;
		
	long f;
	double ss;
	
	/* attaribution analysis *****************/
	
	if(rnum == 0){
		SENSANS_DIST = 0;
		SENSANS_ACO2 = 0;
		SENSANS_TEMP = 0;
		SENSANS_PREC = 0;
		SENSANS_SRAD = 0;
		N_LARS = 0;
	}else{
		/*	if(argc<7){
		printf("Insufficient augments !!!\n");
		exit(1);
		}
		 
		SENSANS_DIST = (short)atol(argv[1]);
		SENSANS_ACO2 = (short)atol(argv[2]);
		SENSANS_TEMP = (short)atol(argv[3]);
		SENSANS_PREC = (short)atol(argv[4]);
		SENSANS_SRAD = (short)atol(argv[5]);
		N_LARS = (short)atol(argv[6]); */
		
		/* SENSANS_DIST = 0;
		SENSANS_ACO2 = 1;
		SENSANS_TEMP = 2;
		SENSANS_PREC = 2;
		SENSANS_SRAD = 2; */
		/* N_LARS = 1; */
		/* N_LARS = rnum; */

        SENSANS_DIST = 0;
		SENSANS_ACO2 = 0;
		SENSANS_TEMP = 0;
		SENSANS_PREC = 0;
		SENSANS_SRAD = 0;
		N_LARS = 0;
    }
    
    if(PERTURB_MPARA==1){
        srand(rnum);
        SA_PARA = 10;
        SA_PARA_EN = 0;
    }
	
	/* parameter sensitivity analysis *****************/
	/*	if(argc<3){
	 printf("Insufficient augments !!!\n");
	 exit(1);
	 }
	 
	 SA_PARA = (long)atol(argv[1]);
	 SA_PARA_EN = (long)atol(argv[2]);
	 
	 srand(SA_PARA_EN);
	 ss = 0.0;
	 for(f=0;f<12;f++){
	 ss += (double)rand() / RAND_MAX;
	 }
	 ss -= 6.0;
	 SA_PARA_VAR = ss*0.1 + 1.0;  */
	
	SA_PARA_VAR = 0.0; /**/
	
	printf("PARASA: %d %d %lf\n", SA_PARA, SA_PARA_EN, SA_PARA_VAR);
}


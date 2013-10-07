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

/* initialize vegetation conditions ********************************************/
void set_parameter(
	long mode, 
	char filename[], 
	struct Echar echar[]
){
	long f, NREC;
	double in_dat;
	FILE *fp_parameter;
	
	/* open parameter file ************/
	if((fp_parameter = fopen(filename,"rt")) == NULL){
		printf("No parameter data !: %s\n", filename);
		exit (1);
	}
	
	/* scan version and date ID */
	fscanf(fp_parameter,"%s", &(echar[0].para_ver_id));
	fscanf(fp_parameter,"%s", &(echar[0].para_date_id));
	
	/* number of biome types */
	switch(mode){
		case 1:
			NREC = 1;
			break;
		case 2: 
			NREC = 16;  /* for new land cover data */
			break;
		default:
			printf("!! Wrong mode\n");
			exit(1);
			break;
	}
	
	/* check biome number */
	if(mode==2){
		for(f=0;f<NREC;f++){
			fscanf(fp_parameter, "%lf", &in_dat);
			if((long)in_dat!=f){
				printf("!! ERROR: WRONG PARAMETER FILE !!\n");
				exit(1);
			}
		}
	}
	
	/* check divider */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}
	
	for(f=0;f<NREC;f++){
		(echar[f].tree).c_layer = 0;
		(echar[f].c3).c_layer = 1;
		(echar[f].c4).c_layer = 1;
	}
	
	/* trees **********************************************/
	/*  daily model  parameter */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).albedo = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).alloc_ass = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).alloc_abg = in_dat;
	}
	
	/* phenology */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).phenoltype = (long)in_dat;
		/* (echar[f].tree).phenoltype = 0;   */ /* test ******************/
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).ctmp_lfdsp = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).ctmp_lfshd = in_dat;
	}
	
	/*   photosynthesis parameters  */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).phototype = (long)in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).sla = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).eK0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).lue0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).pmax = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).topt0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).tmin = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).tmax = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).gs_b0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).gs_b1 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).gs_b2 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).km_nstl = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).kmci = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).cmpcd0 = in_dat;
	}
	
	/*   respiration parameters  */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).rgf = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).rgc = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).rgr = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).rmf0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).rmc_s = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).rmr_s = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).rmc_h = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).rmr_h = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).qTf0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).qTc0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).qTr0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).f_sz_s = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).f_sz_r = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).lf0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).lc0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).lr0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).dcd = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).root_dpt_a = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].tree).root_dpt_b = in_dat;
	}
	
	/* check */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}

	/* herb c3 **********************************************/
	/*  daily model  parameter */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).albedo = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).alloc_ass = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).alloc_abg = in_dat;
	}
	
	/* phenology */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).phenoltype = (long)in_dat;
		/* (echar[f].c3).phenoltype = 0;   */ /* test ******************/
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).ctmp_lfdsp = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).ctmp_lfshd = in_dat;
	}
	
	/*   photosynthesis parameters  */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).phototype = (long)in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).sla = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).eK0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).lue0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).pmax = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).topt0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).tmin = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).tmax = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).gs_b0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).gs_b1 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).gs_b2 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).km_nstl = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).kmci = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).cmpcd0 = in_dat;
	}
	
	/*   respiration parameters  */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).rgf = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).rgc = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).rgr = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).rmf0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).rmc_s = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).rmr_s = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).rmc_h = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).rmr_h = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).qTf0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).qTc0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).qTr0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).f_sz_s = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).f_sz_r = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).lf0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).lc0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).lr0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).dcd = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).root_dpt_a = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c3).root_dpt_b = in_dat;
	}
	
	/* check */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}

	/* herb c3 **********************************************/
	/*  daily model  parameter */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).albedo = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).alloc_ass = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).alloc_abg = in_dat;
	}
	
	/* phenology */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).phenoltype = (long)in_dat;
		/* (echar[f].c4).phenoltype = 0;   */ /* test ******************/
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).ctmp_lfdsp = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).ctmp_lfshd = in_dat;
	}
	
	/*   photosynthesis parameters  */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).phototype = (long)in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).sla = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).eK0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).lue0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).pmax = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).topt0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).tmin = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).tmax = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).gs_b0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).gs_b1 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).gs_b2 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).km_nstl = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).kmci = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).cmpcd0 = in_dat;
	}
	
	/*   respiration parameters  */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).rgf = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).rgc = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).rgr = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).rmf0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).rmc_s = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).rmr_s = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).rmc_h = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).rmr_h = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).qTf0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).qTc0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).qTr0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).f_sz_s = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).f_sz_r = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).lf0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).lc0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).lr0 = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).dcd = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).root_dpt_a = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].c4).root_dpt_b = in_dat;
	}
	
	/* check */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}

	/* soil **********************************************/
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).albedo0 = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).sr_lf = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).sr_lc = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).sr_lr = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).sr_ha = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).sr_hi = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).sr_hp = in_dat;
		
		/*  printf("DEBUD:  %lf\n", in_dat);  */
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).kml = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).kmh = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).kmsl = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).kmsh = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).f_co2_lf = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).f_co2_lc = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).f_co2_lr = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).f_hm_a = in_dat;
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).f_hm_i = in_dat;
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].soil).f_hm_p = in_dat; 
	}
	
	for(f=0;f<NREC;f++){
		if( ((echar[f].soil).f_hm_a+(echar[f].soil).f_hm_i+(echar[f].soil).f_hm_p)<0.95 ||
		 ((echar[f].soil).f_hm_a+(echar[f].soil).f_hm_i+(echar[f].soil).f_hm_p)>1.05){
			printf("Warning: soil humification imbalance: %lf\n", 
				(echar[f].soil).f_hm_a+(echar[f].soil).f_hm_i+(echar[f].soil).f_hm_p);
		}
	
		(echar[f].soil).ntrcon_max=		0.0001;
	}
	
	/* check */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}

	/* 30-min De Pury & Farquhar model's parameters *****************/
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).ke_d = (echar[f].overs2).ke_d = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).ke_n = (echar[f].overs2).ke_n = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).scttr = (echar[f].overs2).scttr = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).spect = (echar[f].overs2).spect = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).rfl_d = (echar[f].overs2).rfl_d = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).photocap_n0 = (echar[f].overs2).photocap_n0 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).nc_ct = (echar[f].overs2).nc_ct = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).nc_np = (echar[f].overs2).nc_np = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).acen_et = (echar[f].overs2).acen_et = in_dat; 
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).hc_et = (echar[f].overs2).hc_et = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).str_et = (echar[f].overs2).str_et = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).kc0 = (echar[f].overs2).kc0 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).ko0 = (echar[f].overs2).ko0 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).acen_kc = (echar[f].overs2).acen_kc = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).acen_ko = (echar[f].overs2).acen_ko = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).acen_vc = (echar[f].overs2).acen_vc = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).acen_rd = (echar[f].overs2).acen_rd = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).cv_j = (echar[f].overs2).cv_j = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).cv_p = (echar[f].overs2).cv_p = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).jvr = (echar[f].overs2).jvr = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).gss_b0 = (echar[f].overs2).gss_b0 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).gss_b1 = (echar[f].overs2).gss_b1 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].overs1).gss_b2 = (echar[f].overs2).gss_b2 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		if((echar[f].overs1).ke_n > 0.0){
			(echar[f].overs1).vcmx_con = (echar[f].overs1).photocap_n0*((echar[f].overs1).nc_ct-
				(echar[f].overs1).nc_np)*(1.0-exp(-(echar[f].overs1).ke_n))/(echar[f].overs1).ke_n;
		}else{
			(echar[f].overs1).vcmx_con = 0.0;
		}
		
		if((echar[f].overs2).ke_n > 0.0){
			(echar[f].overs2).vcmx_con = (echar[f].overs2).photocap_n0*((echar[f].overs2).nc_ct-
				(echar[f].overs2).nc_np)*(1.0-exp(-(echar[f].overs2).ke_n))/(echar[f].overs2).ke_n;
		}else{
			(echar[f].overs2).vcmx_con = 0.0;
		}
	}

	/* check */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).ke_d = in_dat; 
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).ke_n = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).scttr = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).spect = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).rfl_d = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).photocap_n0 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).nc_ct = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).nc_np = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).acen_et = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).hc_et = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).str_et = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).kc0 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).ko0 = in_dat; 
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).acen_kc = in_dat; 
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).acen_ko = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).acen_vc = in_dat; 
	}
	
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).acen_rd = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).cv_j = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).cv_p = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).jvr = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).gss_b0 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).gss_b1 = in_dat; 
	}

	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		(echar[f].unders).gss_b2 = in_dat; 
	}
	
	for(f=0;f<NREC;f++){
		if((echar[f].unders).ke_n > 0.0){
			(echar[f].unders).vcmx_con = (echar[f].unders).photocap_n0*((echar[f].unders).nc_ct-
				(echar[f].unders).nc_np)*(1.0-exp(-(echar[f].unders).ke_n))/(echar[f].unders).ke_n;
		}else{
			(echar[f].unders).vcmx_con = 0.0;
		}
	}
	
	/* check */
	for(f=0;f<NREC;f++){
		fscanf(fp_parameter, "%lf", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}
	
	fclose(fp_parameter);
}

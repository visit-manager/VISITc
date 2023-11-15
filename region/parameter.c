/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in November 02, 2007						*/

/*  [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

extern long NBIOME;

/* initialize vegetation conditions ********************************************/
void set_parameter(
	long mode, 
	char filename[], 
	struct Echar echar[]
){
	long f;
	float in_dat;
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
			NBIOME = 1;
			break;
		case 2: 
            /* variable biome numbers: 2012/02/23 by A.Ito */
            fscanf(fp_parameter,"%ld", &NBIOME);
            /* NBIOME = 16; //old setting */
			break;
		default:
			printf("!! Wrong mode\n");
			exit(1);
			break;
	}
	
	/* check biome number */
	if(mode == 2){
		for(f=0;f<NBIOME;f++){
			fscanf(fp_parameter, "%f", &in_dat);
			if((long)in_dat!=f){
				printf("!! ERROR: WRONG PARAMETER FILE !!\n");
				exit(1);
			}
		}
	}
	
	/* check divider */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}
	
	/* trees **********************************************/
	/*  daily model  parameter */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).albedo = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).alloc_ass = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).alloc_abg = in_dat;
	}
	
	/* phenology */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).phenoltype = (long)in_dat;
		/* (echar[f].tree).phenoltype = 0;   */ /* test ******************/
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).ctmp_lfdsp = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).ctmp_lfshd = in_dat;
	}
	
	/*   photosynthesis parameters  */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).phototype = (long)in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).sla = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).eK0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).lue0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).pmax = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).topt0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).tmin = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).tmax = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).gs_b0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).gs_b1 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).gs_b2 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).km_nstl = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).kmci = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).cmpcd0 = in_dat;
	}
	
	/*   respiration parameters  */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rgf = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rgc = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rgr = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rmf0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rmc_s = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rmr_s = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rmc_h = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rmr_h = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).qTf0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).qTc0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).qTr0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).f_sz_s = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).f_sz_r = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).lf0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).lc0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).lr0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).dcd = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).root_dpt_a = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).root_dpt_b = in_dat;
	}
	
	/* check */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}

	/* herb c3 **********************************************/
	/*  daily model  parameter */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).albedo = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).alloc_ass = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).alloc_abg = in_dat;
	}
	
	/* phenology */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).phenoltype = (long)in_dat;
		/* (echar[f].c3).phenoltype = 0;   */ /* test ******************/
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).ctmp_lfdsp = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).ctmp_lfshd = in_dat;
	}
	
	/*   photosynthesis parameters  */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).phototype = (long)in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).sla = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).eK0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).lue0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).pmax = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).topt0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).tmin = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).tmax = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).gs_b0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).gs_b1 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).gs_b2 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).km_nstl = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).kmci = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).cmpcd0 = in_dat;
	}
	
	/*   respiration parameters  */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rgf = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rgc = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rgr = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rmf0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rmc_s = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rmr_s = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rmc_h = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rmr_h = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).qTf0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).qTc0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).qTr0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).f_sz_s = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).f_sz_r = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).lf0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).lc0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).lr0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).dcd = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).root_dpt_a = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).root_dpt_b = in_dat;
	}
	
	/* check */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}

	/* herb c3 **********************************************/
	/*  daily model  parameter */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).albedo = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).alloc_ass = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).alloc_abg = in_dat;
	}
	
	/* phenology */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).phenoltype = (long)in_dat;
		/* (echar[f].c4).phenoltype = 0;   */ /* test ******************/
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).ctmp_lfdsp = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).ctmp_lfshd = in_dat;
	}
	
	/*   photosynthesis parameters  */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).phototype = (long)in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).sla = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).eK0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).lue0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).pmax = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).topt0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).tmin = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).tmax = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).gs_b0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).gs_b1 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).gs_b2 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).km_nstl = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).kmci = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).cmpcd0 = in_dat;
	}
	
	/*   respiration parameters  */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).rgf = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).rgc = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).rgr = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).rmf0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).rmc_s = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).rmr_s = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).rmc_h = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).rmr_h = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).qTf0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).qTc0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).qTr0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).f_sz_s = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).f_sz_r = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).lf0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).lc0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).lr0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).dcd = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).root_dpt_a = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c4).root_dpt_b = in_dat;
	}
	
	/* check */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}

	/* soil **********************************************/
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).albedo0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).sr_lf = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).sr_lc = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).sr_lr = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).sr_ha = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).sr_hi = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).sr_hp = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).kml = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).kmh = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).kmsl = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).kmsh = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).f_co2_lf = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).f_co2_lc = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).f_co2_lr = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).f_hm_a = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).f_hm_i = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].soil).f_hm_p = in_dat; 
	}
	
	for(f=0;f<NBIOME;f++){
		if( ((echar[f].soil).f_hm_a+(echar[f].soil).f_hm_i+(echar[f].soil).f_hm_p)<0.95 ||
		 ((echar[f].soil).f_hm_a+(echar[f].soil).f_hm_i+(echar[f].soil).f_hm_p)>1.05){
			printf("Warning: soil humification imbalance: %ld %f\n", f,
				(echar[f].soil).f_hm_a+(echar[f].soil).f_hm_i+(echar[f].soil).f_hm_p);
		}
	
		(echar[f].soil).ntrcon_max = 0.0001;
	}
	
	/* check */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		if((long)in_dat != 9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}

	/* 30-min De Pury & Farquhar model's parameters *****************/
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).ke_d = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).ke_n = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).scttr = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).spect = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).rfl_d = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).photocap_n0 = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).nc_ct = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).nc_np = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).acen_et = in_dat; 
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).hc_et = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).str_et = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).kc0 = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).ko0 = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).acen_kc = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).acen_ko = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).acen_vc = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).acen_rd = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).cv_j = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).cv_p = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).jvr = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).gss_b0 = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).gss_b1 = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].tree).gss_b2 = in_dat; 
	}

	for(f=0;f<NBIOME;f++){
		if((echar[f].tree).ke_n > 0.0){
			(echar[f].tree).vcmx_con = (echar[f].tree).photocap_n0*((echar[f].tree).nc_ct-
				(echar[f].tree).nc_np)*(1.0-exp(-(echar[f].tree).ke_n))/(echar[f].tree).ke_n;
		}else{
			(echar[f].tree).vcmx_con = 0.0;
		}
	}

	/* check */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).ke_d = (echar[f].c4).ke_d = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).ke_n = (echar[f].c4).ke_n = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).scttr = (echar[f].c4).scttr = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).spect = (echar[f].c4).spect = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).rfl_d = (echar[f].c4).rfl_d = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).photocap_n0 = (echar[f].c4).photocap_n0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).nc_ct = (echar[f].c4).nc_ct = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).nc_np = (echar[f].c4).nc_np = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).acen_et = (echar[f].c4).acen_et = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).hc_et = (echar[f].c4).hc_et = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).str_et = (echar[f].c4).str_et = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).kc0 = (echar[f].c4).kc0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).ko0 = (echar[f].c4).ko0 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).acen_kc = (echar[f].c4).acen_kc = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).acen_ko = (echar[f].c4).acen_ko = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).acen_vc = (echar[f].c4).acen_vc = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).acen_rd = (echar[f].c4).acen_rd = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).cv_j = (echar[f].c4).cv_j = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).cv_p = (echar[f].c4).cv_p = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).jvr = (echar[f].c4).jvr = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).gss_b0 = (echar[f].c4).gss_b0 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).gss_b1 = (echar[f].c4).gss_b1 = in_dat;
	}

	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		(echar[f].c3).gss_b2 = (echar[f].c4).gss_b2 = in_dat;
	}
	
	for(f=0;f<NBIOME;f++){
		if((echar[f].c3).ke_n > 0.0){
			(echar[f].c3).vcmx_con = (echar[f].c3).photocap_n0*((echar[f].c3).nc_ct-
				(echar[f].c3).nc_np)*(1.0-exp(-(echar[f].c3).ke_n))/(echar[f].c3).ke_n;
		}else{
			(echar[f].c3).vcmx_con = 0.0;
		}
        
        if((echar[f].c4).ke_n > 0.0){
			(echar[f].c4).vcmx_con = (echar[f].c4).photocap_n0*((echar[f].c4).nc_ct-
				(echar[f].c4).nc_np)*(1.0-exp(-(echar[f].c4).ke_n))/(echar[f].c4).ke_n;
		}else{
			(echar[f].c4).vcmx_con = 0.0;
		}
	}
	
	/* check */
	for(f=0;f<NBIOME;f++){
		fscanf(fp_parameter, "%f", &in_dat);
		if((long)in_dat!=9999999){
			printf("!! ERROR: WRONG PARAMETER FILE !!\n");
			exit(1);
		}
	}
	
	fclose(fp_parameter);
}

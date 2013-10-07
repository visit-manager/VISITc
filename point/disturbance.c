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

/*  disturbance  *************************************************/
void disturbance_regime(
	long dyear, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar,
	struct Mass *mass, 
	struct Flux *flux
){
	double defol, abn_leaves, abn_branches;
    char filename[128];

	abn_leaves = abn_branches = 0.0;
	
	if(strcmp(grid->site_id, "TKY")==0){
		/* logging *************************/
		if(dyear==YLDIST && loct->doy==354){
			logging_event(grid, mass);  	   
		}
		
		/* TYPHOON *************************/
        /*
        Ito, A.: Evaluation of defoliation impacts of tropical cyclones on the forest 
        carbon budget using flux data and a process-based model, 
        Journal of Geophysical Research, 115, doi:10.1029/2010JG001314, 
        doi:10.1029/2010JG001314, 2010. */

		if(dyear==2004 && loct->doy==162){	/* CONSON ***/
			defol = 0.024429;
		}else if(dyear==2004 && loct->doy==173){	/* DIANMU *****/
			defol = 0.024285;
		}else if(dyear==2004 && loct->doy==212){	/* NAMTHEUM ***/
			defol = 0.065341;
		}else if(dyear==2004 && loct->doy==216){	/* MALOU ***/
			defol = 0.066581;
		}else if(dyear==2004 && loct->doy==231){	/* MEGI ***/
			defol = 0.069883;
		}else if(dyear==2004 && loct->doy==243){	/* CHABA *****/
			defol = 0.133493;
		}else if(dyear==2004 && loct->doy==251){	/* SONGDA *****/
			defol = 0.143948;
		}else if(dyear==2004 && loct->doy==272){	/* MEARI ***/
			defol = 0.081585;
		}else if(dyear==2004 && loct->doy==282){	/* MA-ON ***/
			defol = 0.028193;
		}else if(dyear==2004 && loct->doy==293){	/* TOKAGE ***/
			defol = 0.028201;
		}else{
			defol = 0.0;
		}
		
		/* light */
		abn_leaves = (mass->tree).fol * defol;
		abn_branches = (mass->tree).fol * defol * 5.296;
		/* coefficient 5.296 is obtained from Ohtsuka et al.(2005)  
		 leaves: 3.38 Mg/ha, stems: 17.9 Mg/ha
		 */
		if(abn_branches > (mass->tree).stm*0.9){
			abn_branches = (mass->tree).stm*0.9;
		}
		
		(mass->tree).fol -= abn_leaves;
		(mass->tree).stm -= abn_branches;
		(mass->soil).ltr_tf += abn_leaves;
		(mass->soil).ltr_tc += abn_branches;
	}
	
	if((strcmp(grid->site_id, "TMK")==0) || (strcmp(grid->site_id, "FJY")==0)){
		if(dyear==YLDIST && loct->doy==354){
            logging_event(grid, mass);  	 /*   */
			
			loct->age_stand = 0.0;
		}
	}else if((strcmp(grid->site_id, "TSE")==0) || (strcmp(grid->site_id, "CEAMIP_TSE")==0)){
		if(dyear==YLDIST && loct->doy==1){
			logging_event(grid, mass);  	   /*   */
			loct->age_stand = 0.0;
		}
        if(dyear==(YLDIST+1) && loct->doy==270){
            (mass->tree).fol += 1.0;
            (mass->tree).stm += 1.0;
            (mass->tree).rot += 1.0;
            (mass->tree).nsch_storage += 1.0;
		}
	}else if(strcmp(grid->site_id, "PSO")==0){
		//if((dyear==1980 || dyear==2000 || dyear==2020 || dyear==2040) && loct->doy==364){
		if((dyear==YLDIST) && loct->doy==364){
			logging_event(grid, mass);  	   /*   */
			loct->age_stand = 0.0;
			loct->veg_state = 1;
            
            if(LCCONV==1){
                strcpy(filename,"parameter_PSO_op.txt");		
                set_parameter(1, filename, echar);
            }
		}
        
	}else{
		if(dyear==YLDIST && loct->doy==354){
			logging_event(grid, mass);  	   /* */
			loct->age_stand = 0.0;
		}
	}
	
	if(WMODE==2||WMODE==3){
		if(dyear==YLDIST && loct->doy==354){
			logging_event(grid, mass); 
		}
	}
}

/* logging ***************************************************************/
void logging_event(
	struct Grid *grid, 
	struct Mass *mass
){
    long f;
	double ems_fol, ems_stm, ems_rot, ems_ltr_f, ems_ltr_c, ems_ltr_r;
	double brn_fol, brn_stm, brn_rot;
	double nloss_fol, nloss_storage;
    double ss;
	
/*	ems_fol=(mass->tree).fol*0.45;
	ems_stm=(mass->tree).stm*0.35;
	ems_rot=(mass->tree).rot*0.35;
	ems_ltr=(mass->soil).ltr*0.65;

	brn_fol=(mass->tree).fol*0.50;
	brn_stm=(mass->tree).stm*0.60;
	brn_rot=(mass->tree).rot*0.60;  */
    
    /* perturbation in disturbance intensity: 2012/02/20 by A.Ito */
    ss = 0.0;
    for(f=0;f<12;f++){
        ss += (double)rand() / RAND_MAX;
    }
    ss -= 6.0;
    ss = ss*0.1 + 1.0;
    if(ss<0.7){
        ss = 0.7;
    }
    if(ss>1.3){
        ss = 1.3;
    }
	
    ems_ltr_f = 0.0;
    ems_ltr_c = 0.0;
    ems_ltr_r = 0.0;

	/* emission */
	if(strcmp(grid->site_id, "TKY")==0){
		ems_fol = (mass->tree).fol*0.45;
		ems_stm = (mass->tree).stm*0.65;
		ems_rot = (mass->tree).rot*0.65;
		ems_ltr_f = (mass->soil).ltr_tf*0.40;
		ems_ltr_c = (mass->soil).ltr_tc*0.40;
		ems_ltr_r = (mass->soil).ltr_tr*0.40;
		brn_fol = (mass->tree).fol*0.35;
		brn_stm = (mass->tree).stm*0.20;
		brn_rot = (mass->tree).rot*0.25; /* */
		
        /* perturbed simulation: 2012/02/20 by A.Ito */
        if(PERTURB_MPARA==1){
            ems_fol = (mass->tree).fol*0.15 * ss;
            ems_stm = (mass->tree).stm*0.50 * ss;
            ems_rot = (mass->tree).rot*0.05 * ss;
            ems_ltr_f = (mass->soil).ltr_tf*0.00 * ss;
            ems_ltr_c = (mass->soil).ltr_tc*0.10 * ss;
            ems_ltr_r = (mass->soil).ltr_tr*0.00 * ss;
            brn_fol = (mass->tree).fol*0.60 * ss;
            brn_stm = (mass->tree).stm*0.25 * ss;
            brn_rot = (mass->tree).rot*0.60 * ss;
        }
        
        /* burnt residue */
		/* brn_fol = (mass->tree).fol*0.25;
		brn_stm = (mass->tree).stm*0.10;
		brn_rot = (mass->tree).rot*0.15; */

	/*	ems_fol = (mass->tree).fol*0.55;
		ems_stm = (mass->tree).stm*0.85;
		ems_rot = (mass->tree).rot*0.75;
		ems_ltr_f = (mass->soil).ltr_tf*0.65;
		ems_ltr_c = (mass->soil).ltr_tc*0.65;
		ems_ltr_r = (mass->soil).ltr_tr*0.65;  */
		/* burnt residue */
	/*	brn_fol = (mass->tree).fol*0.30;
		brn_stm = (mass->tree).stm*0.10;
		brn_rot = (mass->tree).rot*0.20;  */
	}else if((strcmp(grid->site_id, "TMK")==0)){
		ems_fol = (mass->tree).fol*0.45;
		ems_stm = (mass->tree).stm*0.55;
		ems_rot = (mass->tree).rot*0.45;
		ems_ltr_f = (mass->soil).ltr_tf*0.70;
		ems_ltr_c = (mass->soil).ltr_tc*0.70;
		ems_ltr_r = (mass->soil).ltr_tr*0.70;
		/* burnt residue */
		brn_fol = (mass->tree).fol*0.20;
		brn_stm = (mass->tree).stm*0.10;
		brn_rot = (mass->tree).rot*0.10;
		
		(mass->soil).msl_a *= 0.5;
		(mass->soil).msl_i *= 0.5;
		(mass->soil).msl_p *= 0.5;
	}else if((strcmp(grid->site_id, "FJY")==0)){
		ems_fol = (mass->tree).fol*0.9;
		ems_stm = (mass->tree).stm*0.9;
		ems_rot = (mass->tree).rot*0.9;
		ems_ltr_f = (mass->soil).ltr_tf*0.85;
		ems_ltr_c = (mass->soil).ltr_tc*0.85;
		ems_ltr_r = (mass->soil).ltr_tr*0.85;
		/* burnt residue */
		brn_fol = (mass->tree).fol*0.05;
		brn_stm = (mass->tree).stm*0.05;
		brn_rot = (mass->tree).rot*0.05;
		
		(mass->soil).msl_a *= 0.1;
		(mass->soil).msl_i *= 0.1;
		(mass->soil).msl_p *= 0.1; /* */
	}else if((strcmp(grid->site_id, "TSE")==0) || (strcmp(grid->site_id, "CEAMIP_TSE")==0)){
		ems_fol = (mass->tree).fol*0.49;
		ems_stm = (mass->tree).stm*0.65;
		ems_rot = (mass->tree).rot*0.65;
		ems_ltr_f = (mass->soil).ltr_tf*0.0;
		ems_ltr_c = (mass->soil).ltr_tc*0.0;
		ems_ltr_r = (mass->soil).ltr_tr*0.0;
		/* burnt residue */
		brn_fol = (mass->tree).fol*0.20;
		brn_stm = (mass->tree).stm*0.10;
		brn_rot = (mass->tree).rot*0.10;
		
		(mass->tree).fol = 0.1;
		(mass->tree).stm = 0.1;
		(mass->tree).rot = 0.1;
		(mass->tree).nsch_storage = 0.05;
        
		(mass->c3).fol = 0.1;
		(mass->c3).stm = 0.1;
		(mass->c3).rot = 0.1;
		
	/*	(mass->soil).msl_a *= 0.95;
		(mass->soil).msl_i *= 0.95;
		(mass->soil).msl_p *= 0.95;  */
	}else if((strcmp(grid->site_id, "PSO")==0) && LCCONV==0){
	/*	ems_fol = (mass->tree).fol*0.2;
		ems_stm = (mass->tree).stm*0.2;
		ems_rot = (mass->tree).rot*0.2;
		brn_fol = (mass->tree).fol*0.2;
		brn_stm = (mass->tree).stm*0.2;
		brn_rot = (mass->tree).rot*0.2; */
		
		ems_fol = (mass->tree).fol*0.4;
		ems_stm = (mass->tree).stm*0.4;
		ems_rot = (mass->tree).rot*0.4;
		brn_fol = (mass->tree).fol*0.4;
		brn_stm = (mass->tree).stm*0.4;
		brn_rot = (mass->tree).rot*0.4;
		
	}else if((strcmp(grid->site_id, "PSO")==0) && LCCONV==1){
	/*	ems_fol = (mass->tree).fol*0.2;
		ems_stm = (mass->tree).stm*0.2;
		ems_rot = (mass->tree).rot*0.2;
		brn_fol = (mass->tree).fol*0.2;
		brn_stm = (mass->tree).stm*0.2;
		brn_rot = (mass->tree).rot*0.2; */
		
		ems_fol = (mass->tree).fol*0.475;
		ems_stm = (mass->tree).stm*0.475;
		ems_rot = (mass->tree).rot*0.475;
		brn_fol = (mass->tree).fol*0.475;
		brn_stm = (mass->tree).stm*0.475;
		brn_rot = (mass->tree).rot*0.475;
		
	}else if((strcmp(grid->site_id, "CEAMIP_MKL")==0)){
		ems_fol = (mass->tree).fol*0.25;
		ems_stm = (mass->tree).stm*0.25;
		ems_rot = (mass->tree).rot*0.25;
		brn_fol = (mass->tree).fol*0.25;
		brn_stm = (mass->tree).stm*0.25;
		brn_rot = (mass->tree).rot*0.25;
	}else{
		ems_fol = (mass->tree).fol*0.55;
		ems_stm = (mass->tree).stm*0.85;
		ems_rot = (mass->tree).rot*0.75;
		ems_ltr_f = (mass->soil).ltr_tf*0.65;
		ems_ltr_c = (mass->soil).ltr_tc*0.65;
		ems_ltr_r = (mass->soil).ltr_tr*0.65;
		/* burnt residue */
		brn_fol = (mass->tree).fol*0.30;
		brn_stm = (mass->tree).stm*0.10;
		brn_rot = (mass->tree).rot*0.20;
	}
	
	/********************************************/
	if(SA_PARA==8){
		ems_fol *= SA_PARA_VAR;
		ems_stm *= SA_PARA_VAR;
		ems_rot *= SA_PARA_VAR;
		ems_ltr_f *= SA_PARA_VAR;
		ems_ltr_c *= SA_PARA_VAR;
		ems_ltr_r *= SA_PARA_VAR;
	}
	
	/* mass balance */
	(mass->tree).fol -= ems_fol+brn_fol;
	if((mass->tree).fol < 0.0){
		(mass->tree).fol = INIT_C;
	}
	(mass->tree).stm -= ems_stm+brn_stm;
	if((mass->tree).stm < 0.0){
		(mass->tree).stm = INIT_C;
	}
	(mass->tree).rot -= ems_rot+brn_rot;
	if((mass->tree).rot < 0.0){
		(mass->tree).rot = INIT_C;
	}

	(mass->soil).ltr_tf -= ems_ltr_f;
	if((mass->soil).ltr_tf < 0.0){
		(mass->soil).ltr_tf = INIT_C;
	}
	(mass->soil).ltr_tc -= ems_ltr_c;
	if((mass->soil).ltr_tc < 0.0){
		(mass->soil).ltr_tc = INIT_C;
	}
	(mass->soil).ltr_tr -= ems_ltr_r;
	if((mass->soil).ltr_tr < 0.0){
		(mass->soil).ltr_tr = INIT_C;
	}
	
	(mass->soil).ltr_tf += brn_fol;
	(mass->soil).ltr_tc += brn_stm;
	(mass->soil).ltr_tr += brn_rot;
	
	/* nitrogen flows *****/
	nloss_fol = 0.3 * (mass->tree).n_canopy;
	nloss_storage = 0.3 * (mass->tree).n_storage;
	(mass->tree).n_canopy -= nloss_fol;
	(mass->tree).n_storage -= nloss_storage;
	(mass->soil).n_lttr += nloss_fol + nloss_storage;

	(mass->c3).fol += 0.01;
	(mass->c3).stm += 0.01;
	(mass->c3).rot += 0.01;
}

/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in Nvember 02, 2007							*/

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

extern long	WMODE;
extern long YLDIST;

void disturbance_regime(
	long dyear, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	if(strcmp(grid->site_id, "TKY")==0){
		/*******  disturbance  *******/
		if(dyear==YLDIST && loct->doy==354){
			logging_event(grid, mass);  	   /*   */
		}
		
		/*******  TYPHOON  *******/
	/*	if(e==1998&&f==264){   */
		if(dyear&&loct->doy==171){
		/* 	typhoon_lf = (mass->tree).fol*0.3;
			(mass->tree).fol -= typhoon_lf;
			(flux->tree).lf += typhoon_lf;	*/
		/*  	(mass->tree).fol-=1.2;
			(flux->tree).lf+=1.2;  */
		/*	(mass->tree).fol-=0.5;
			(flux->tree).lf+=0.5;    */
		}  
		if(dyear==2004 && loct->doy==240){
		
		}
		if(dyear==2004 && loct->doy==249){
		
		}
	}
	
	if((strcmp(grid->site_id, "TMK")==0) || (strcmp(grid->site_id, "FJY")==0)){
		if(dyear==YLDIST && loct->doy==354){
			logging_event(grid, mass);  	   /*   */
			
			grid->age_stand = 0.0;
		}
	}else{
		if(dyear==YLDIST && loct->doy==354){
			logging_event(grid, mass);  	   /*   */
			grid->age_stand = 0.0;
		}
	}
	
	if(WMODE == 2 || WMODE == 3){
		if(dyear==YLDIST && loct->doy==354){
			logging_event(grid, mass); 
		}
	}
}

/* logging **************************************************************************/
void logging_event(
	struct Grid *grid, 
	struct Mass *mass
){
	float ems_fol, ems_stm, ems_rot, ems_ltr_f, ems_ltr_c, ems_ltr_r;
	float brn_fol, brn_stm, brn_rot;
	
/*	ems_fol=(mass->tree).fol*0.45;
	ems_stm=(mass->tree).stm*0.35;
	ems_rot=(mass->tree).rot*0.35;
	ems_ltr=(mass->soil).ltr*0.65;

	brn_fol=(mass->tree).fol*0.50;
	brn_stm=(mass->tree).stm*0.60;
	brn_rot=(mass->tree).rot*0.60;  */
	
	/* emission */
	if(strcmp(grid->site_id, "TKY")==0){
		ems_fol = (mass->tree).fol*0.45;
		ems_stm = (mass->tree).stm*0.55;
		ems_rot = (mass->tree).rot*0.45;
		ems_ltr_f = (mass->soil).ltr_tf*0.40;
		ems_ltr_c = (mass->soil).ltr_tc*0.40;
		ems_ltr_r = (mass->soil).ltr_tr*0.40;
		/* burnt residue */
		brn_fol = (mass->tree).fol*0.25;
		brn_stm = (mass->tree).stm*0.10;
		brn_rot = (mass->tree).rot*0.15;

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
	}else if((strcmp(grid->site_id, "TMK")==0) || (strcmp(grid->site_id, "FJY")==0)){
		ems_fol = (mass->tree).fol*0.65;
		ems_stm = (mass->tree).stm*0.85;
		ems_rot = (mass->tree).rot*0.85;
		ems_ltr_f = (mass->soil).ltr_tf*0.75;
		ems_ltr_c = (mass->soil).ltr_tc*0.75;
		ems_ltr_r = (mass->soil).ltr_tr*0.75;
		/* burnt residue */
		brn_fol = (mass->tree).fol*0.20;
		brn_stm = (mass->tree).stm*0.10;
		brn_rot = (mass->tree).rot*0.10;

		(mass->soil).msl_a *= 0.1;
		(mass->soil).msl_i *= 0.1;
		(mass->soil).msl_p *= 0.1;
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
	
	/* mass balance */
	(mass->tree).fol -= ems_fol+brn_fol;
	(mass->tree).stm -= ems_stm+brn_stm;
	(mass->tree).rot -= ems_rot+brn_rot;
	(mass->soil).ltr_tf -= ems_ltr_f;
	(mass->soil).ltr_tc -= ems_ltr_c;
	(mass->soil).ltr_tr -= ems_ltr_r;
	
	(mass->soil).ltr_tf += brn_fol;
	(mass->soil).ltr_tc += brn_stm;
	(mass->soil).ltr_tr += brn_rot;
}

/* Bamiyan grazing: 2017/09/03 by A.Ito ********************/
void grazing_event(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass
){
    float abg_mass, f_grazing;
    float intensity_grazing;
    float scenario_grazing_2[36]={1.0, 1.017413023, 1.041780486,
            1.065151572, 1.088113039, 1.111098532, 1.134164424,
            1.157162386, 1.180086962, 1.20289565, 1.225557045,
            1.248083924, 1.270502262, 1.292802974, 1.314970266,
            1.336985785, 1.358841328, 1.380518512, 1.401974132,
            1.423154621, 1.444014976, 1.464528947, 1.484681746,
            1.504454013, 1.52382965, 1.542792315, 1.561326153,
            1.579412414, 1.597029848, 1.614156257, 1.630774021,
            1.646870119, 1.662439338, 1.67748302, 1.692007477,
            1.70601725};
    float scenario_grazing_3[36]={1.0, 0.997064762, 1.000525979,
            1.002512138, 1.003640822, 1.00434507, 1.004690913,
            1.004562214, 1.003974351, 1.002911556, 1.001369335,
            0.999379926, 0.996984373, 0.994194432, 0.991016744,
            0.987456391, 0.983526256, 0.979231829, 0.974561752,
            0.969499357, 0.964035909, 0.958176576, 0.951934428,
            0.94531959, 0.938344335, 0.931020747, 0.92336116,
            0.915376146, 0.907074906, 0.898466253, 0.889561657,
            0.880374989, 0.870923941, 0.861228971, 0.851312188,
            0.841193808};
    float scenario_grazing_4[36]={1.0, 0.986890632, 0.980211259,
            0.972135081, 0.96329865, 0.954138115, 0.944727214,
            0.934967357, 0.92488532, 0.914478633, 0.903755334,
            0.892756182, 0.881528289, 0.870091424, 0.858460294,
            0.846647838, 0.83467325, 0.822548883, 0.81027272,
            0.79783859, 0.785247175, 0.77251048, 0.759646468,
            0.746670177, 0.733597837, 0.720444986, 0.707226826,
            0.693956699, 0.680646464, 0.667307288, 0.65395191,
            0.640594362, 0.627250915, 0.6139392, 0.600677332,
            0.587481395};
    float scenario_grazing_5[36]={1.0, 0.976716502, 0.960104896,
            0.942377941, 0.924185067, 0.905959407, 0.887776037,
            0.869546757, 0.851302437, 0.833046131, 0.814790321,
            0.796576248, 0.778449133, 0.760428477, 0.742528616,
            0.724761784, 0.707145001, 0.689688836, 0.67239147,
            0.655247715, 0.638258134, 0.621432345, 0.604784282,
            0.588324964, 0.572065836, 0.556017285, 0.540188947,
            0.524588589, 0.509222469, 0.494095993, 0.47921541,
            0.464587557, 0.450220498, 0.436122832, 0.422303034,
            0.408767703};
    
    intensity_grazing = 1.0;
    if(EX_BAMIYAN_S == 0){
        /* no grazing */ ;
        intensity_grazing = 0.0;
    }else if(EX_BAMIYAN_S == 1){
        /* constant grazing */ ;
        intensity_grazing = 1.0;
    }else if(EX_BAMIYAN_S == 2){
        /* population-proportional grazing */ ;
        intensity_grazing = 1.0;
        if(loct->adyear <=2014){
            intensity_grazing = 1.0;
        }else if(loct->adyear >=2015 && loct->adyear <=2050){
            intensity_grazing = scenario_grazing_2[loct->adyear - 2015];
        }else if(loct->adyear >=2051){
            intensity_grazing = scenario_grazing_2[2050 - 2015];
        }
    }else if(EX_BAMIYAN_S == 3){
        /* -2% grazing */ ;
        intensity_grazing = 1.0;
        if(loct->adyear <=2014){
            intensity_grazing = 1.0;
        }else if(loct->adyear >=2015 && loct->adyear <=2050){
            intensity_grazing = scenario_grazing_3[loct->adyear - 2015];
        }else if(loct->adyear >=2051){
            intensity_grazing = scenario_grazing_3[2050 - 2015];
        }
    }else if(EX_BAMIYAN_S == 4){
        /* -3% grazing */ ;
        intensity_grazing = 1.0;
        if(loct->adyear <=2014){
            intensity_grazing = 1.0;
        }else if(loct->adyear >=2015 && loct->adyear <=2050){
            intensity_grazing = scenario_grazing_4[loct->adyear - 2015];
        }else if(loct->adyear >=2051){
            intensity_grazing = scenario_grazing_4[2050 - 2015];
        }
    }else if(EX_BAMIYAN_S == 5){
        /* -4% grazing */ ;
        intensity_grazing = 1.0;
        if(loct->adyear <=2014){
            intensity_grazing = 1.0;
        }else if(loct->adyear >=2015 && loct->adyear <=2050){
            intensity_grazing = scenario_grazing_5[loct->adyear - 2015];
        }else if(loct->adyear >=2051){
            intensity_grazing = scenario_grazing_5[2050 - 2015];
        }
    }
    
    /**********************************************/
    if(loct->funder_c4 > 0.01){
        abg_mass = loct->funder_c3 * ((mass->c3).fol + (mass->c3).stm)
                    + loct->funder_c4 * ((mass->c4).fol + (mass->c4).stm);
        /* 0.01: assumed grazing demand MgC/ha/day */
    
        if(loct->hour == 0 && abg_mass > (1.50/dmTc) && loct->tmp_2m>5.0){
            /* grazing intensity assumed */
            f_grazing = intensity_grazing * 0.001 * (mass->c3).fol;
            (mass->c3).fol -= f_grazing;

            f_grazing = intensity_grazing * 0.001 * (mass->c3).stm;
            (mass->c3).stm -= f_grazing;

            f_grazing = intensity_grazing * 0.001 * (mass->c4).fol;
            (mass->c4).fol -= f_grazing;

            f_grazing = intensity_grazing * 0.001 * (mass->c4).stm;
            (mass->c4).stm -= f_grazing;
        }
    }else{
        abg_mass = (mass->c3).fol + (mass->c3).stm;
        
        if(loct->hour == 0 && abg_mass > (1.50/dmTc) && loct->tmp_2m>5.0){
            f_grazing = intensity_grazing * 0.001 * (mass->c3).fol;
            (mass->c3).fol -= f_grazing;

            f_grazing = intensity_grazing * 0.001 * (mass->c3).stm;
            (mass->c3).stm -= f_grazing;
        }
    }
    
}

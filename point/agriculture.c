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

void f_agriculture(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Echar *echar, 
	struct Mass *mass,  
	struct Flux *flux
){
    /* Pasoh */
    if(strcmp(grid->site_id, "PSO")==0){
        
        /* harvest */
        if(loct->age_stand > 5.0 && loct->doy==180 && (mass->tree).stm > 5.0){
            (mass->tree).stm -= 3.3;
        }
        
        /* fertilizer input */
        if(loct->age_stand>=1 && ((long)(loct->age_stand)%5)==0 && loct->doy==180){
            loct->nfert_no3 = 10000.0;
            loct->nfert_nh4 = 10000.0;
        }else{
            loct->nfert_no3 = 0.0;
            loct->nfert_nh4 = 0.0;
        }
    }
}

/* CROP **************************************/
void phenology_crop(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	double crit_gdd_b, crit_cdd_b;
	double crit_gdd_e, crit_cdd_e;
	
	/* critical GDD for the beginning of leaf flush */
	crit_gdd_b = 260.0; 
	/* critical GDD for the end of leaf flush */
	crit_gdd_e = 580.0; 
	/* critical GDD for the beginning of leaf shedding */
	crit_cdd_b = -20.0; 
	/* critical GDD for the end of leaf shedding */
	crit_cdd_e = -190.0;
	
	if(pchar->gdd<crit_gdd_b){
		/* dormancy */
		pchar->season = 0;
		pchar->day_flush = 0;
		pchar->day_shed ++;
	}else if( (pchar->gdd >= crit_gdd_b)&&(pchar->gdd<crit_gdd_e) ){
		/* leaf display */
		pchar->season = 2;
		pchar->cdd = 0.0;
		pchar->day_flush ++;
		pchar->day_shed = 0;
	}else if(pchar->gdd >= crit_gdd_e){
		/* growing period */
		pchar->season = 1;
		pchar->day_flush ++;
		
		if( (pchar->cdd<crit_cdd_b)&&(pchar->cdd >= crit_cdd_e) ){
			/* leaf shedding */
			pchar->season = 3;
			pchar->day_shed ++;
		}else if(pchar->cdd<crit_cdd_e){
			/* dormancy */
			pchar->season = 0;
			pchar->gdd = 0.0;
			pchar->day_shed ++;
		}
	}
	
	/**************************/
    if(grid->veg_type==12){
        if(loct->doy<120){
            pchar->crop_stage = 0;
        }else if(loct->doy>=120 && loct->doy<125){
            pchar->crop_stage = 1;	/* planting */
        }else if(loct->doy>=125 && loct->doy<185){
            pchar->crop_stage = 2;	/* vegetative growth */
        }else if(loct->doy>=185 && loct->doy<215){
            pchar->crop_stage = 3;	/* reproductive growoth */
        }else if(loct->doy>=215 && loct->doy<250){
            if(loct->dvi>=1.0){
                pchar->crop_stage = 4;	/* maturation */
            }else{
                pchar->crop_stage = 3;	/* reproductive growoth */
            }
        }else if(loct->doy>=250 && loct->doy<255){
            pchar->crop_stage = 5;	/* harvesting */
        }else if(loct->doy>=255){
            pchar->crop_stage = 6;	/* fallow */
        }
    }else{
        pchar->crop_stage = 0;
    }
}
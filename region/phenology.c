/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
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

/* growth period and phenological process *************************/
void growthperiod(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar,
	struct Pmas *mass
){
	/** cumulative temperature, growing degree days (GDD) *******/	
	if(loct->tair_dayav > pchar->ctmp_lfdsp){
		mass->gdd += (loct->tair_dayav - pchar->ctmp_lfdsp);
	}
	/* initialization of warmth */
	if( (loct->tair_dayav < -7.0) ){
		/* initialization by very cold spell */
		mass->gdd = 0.0;
	}
    if(grid->lat >= 0.0 && loct->doy ==0){
        mass->gdd = 0.0;
    }else if(grid->lat < 0.0 && loct->doy ==182){
        mass->gdd = 0.0;
    }
    
	/* cumulative coldness (CDD) *******************/
	if(loct->tair_dayav <= pchar->ctmp_lfshd){
    
        if((grid->lat >= 0.0 && loct->doy >= 210) ||
            (grid->lat < 0.0 && loct->doy >= 30)){
            
            mass->cdd += (loct->tair_dayav - pchar->ctmp_lfshd);
        }
	}
    
	/* initialization of coldness */
    if((grid->lat >= 0.0 && loct->doy == 209) ||
                (grid->lat < 0.0 && loct->doy == 29)){
        /* summer */
        mass->cdd = 0.0;
    }
	if( (loct->tair_dayav > 18.0) ){
		/* very warm */
		mass->cdd = 0.0;
	}
	
	/** growing period length, days **/
	if(loct->doy==0) mass->grwpd = 0.0;
	if(loct->tair_dayav > 0.0){
		mass->grwpd += 1.0;
	}

	/** phenology scheme *********************************/
	switch(pchar->phenoltype){
		case 0:
			/* no vegetation, bareland */
			phenology_bareland(grid, pchar, mass); 
			break;
		case 1:
			/* evergreen */
			phenology_evergreen(grid, loct, pchar, mass); 
			break;
		case 2:
			/* cold deciduous */
			phenology_colddeciduous(grid, loct, pchar, mass); 
			break;
	}
	
	/* pchar->season = 1; */ /****************************/
}

/* bare land ************************************************/
void phenology_bareland(
	struct Grid *grid,
	struct Pchar *pchar, 
	struct Pmas *mass
){
	mass->season = -1; /* no growth */
	mass->day_frush = 0;
}

/* evergreen biomes **********************************/
void phenology_evergreen(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	if(loct->tmp_sfc >= pchar->tmin){ 
		/* warmer months */
		mass->season = 1; /* growing */
		mass->day_frush = 0;
	}else{
		/* cooler months */
		mass->season = 0; /* no growth */
		mass->day_frush = 0;
	}
}

/* deciduous biomes in higher latitudes *************/
void phenology_colddeciduous(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float crit_gdd_b, crit_cdd_b;
	float crit_gdd_e, crit_cdd_e;
	
	if(grid->veg_type == 4 || grid->veg_type == 5){
	/* deciduous broadleaf */
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 260.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 580.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0; /*	-170.0	*/
	}else if(grid->veg_type == 3){
	/* deciduous needle-leaf */
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 160.0;
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 480.0;
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0;
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -200.0; /*	-170.0	*/
	}else if(grid->veg_type == 10){
	/* grassland */
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 120.0;
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 770.0;
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -50.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -200.0; /*	-170.0	*/
	}else if(grid->veg_type == 12 || grid->veg_type == 14){
    /* crop */
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 150.0;
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 590.0;
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0;
	}else{
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 260.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 580.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0;
	}
		
	if(mass->gdd < crit_gdd_b){
		/* dormancy */
		mass->season = 0;
		mass->day_frush = 0;
		mass->day_shed += 1.0;
		//mass->gdd = 0.0;	/* 2009/04/30 by A.Ito */
	}else if( (mass->gdd >= crit_gdd_b)&&(mass->gdd < crit_gdd_e) ){
		/* leaf display */
		mass->season = 2;
		mass->cdd = 0.0;
		mass->day_frush += 1.0;
		mass->day_shed = 0;
	}else if(mass->gdd >= crit_gdd_e){
		/* growing period */
		mass->season = 1;
		mass->day_frush += 1.0;
		
		if( (mass->cdd < crit_cdd_b)&&(mass->cdd >= crit_cdd_e) ){
			/* leaf shedding */
			mass->season = 3;
			mass->day_shed += 1.0;
			mass->day_frush = 0;	/* 2009/04/30 by A.Ito */
			//mass->gdd = 0.0;
		}else if(mass->cdd < crit_cdd_e){
			/* dormancy */
			mass->season = 0;
			mass->day_shed += 1.0;
			mass->gdd = 0.0;
			mass->day_frush = 0;	/* 2009/04/30 by A.Ito */
		}else{
			//mass->day_shed = 0;		/* 2009/04/30 by A.Ito */
			//mass->cdd = 0.0;
		}
	}
}

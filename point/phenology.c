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

/* growth period and phenological process *************************/
void f_growth_period(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	/** cumulative temperature, growing degree days (GDD) **/	
	if(loct->tmp_2m > pchar->ctmp_lfdsp){
		pchar->gdd += loct->tmp_2m - pchar->ctmp_lfdsp;		
	}
	/* initialization of warmth */
	if( (loct->tmp_2m < -7.0) ){
		/* very cold */
		pchar->gdd = 0.0;
	}
	/* cumulative coldness (CDD) */
	if(loct->tmp_2m <= pchar->ctmp_lfshd && loct->doy>=210){
		pchar->cdd += loct->tmp_2m - pchar->ctmp_lfshd;
	}else{
		//pchar->cdd = 0.0;
	}
	/* initialization of coldness */
	if( (loct->tmp_2m > 18.0)){
		/* very warm */
		pchar->cdd = 0.0;
	}
	
	/** growing period length, days **/
	if(grid->lat>=0.0 && loct->doy==0){
		pchar->grw_pd = 0.0;
	}if(grid->lat<0.0 && loct->doy==182){
		pchar->grw_pd = 0.0;
	}
	if(loct->tmp_2m > 0.0){
		pchar->grw_pd++;
	}
	
	/* revised phenology scheme: 2012/01/06 (A.Ito) */
	if(grid->lat>=0.0){
		if(loct->doy==0){
			pchar->gddjan = 0.0;
		}
		pchar->gddjan += loct->tmp_2m - 0.0;
	}else if(grid->lat<0.0){
		if(loct->doy==182){
			pchar->gddjan = 0.0;
		}
		pchar->gddjan += loct->tmp_2m - 0.0;
	}

	/** phenology **/
	switch(pchar->phenoltype){
		case 0:
			/* no vegetation, bareland */
			phenology_bareland(grid, pchar); 
			break;
		case 1:
			/* evergreen */
			phenology_evergreen(grid, loct, pchar); 
			break;
		case 2:
			/* cold deciduous */
			phenology_colddeciduous(grid, loct, pchar); 
			break;
        case 3:
            /* crop lands */
			phenology_crop(grid, loct, pchar); 
			break;
	}
}

/* bare land ************************************************/
void phenology_bareland(
	struct Grid *grid,
	struct Pchar *pchar
){
	pchar->season = -1; /* no growth */
	pchar->day_flush = 0;
}

/* evergreen biomes ******************************************/
void phenology_evergreen(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	if(loct->tmp_sfc >= pchar->tmin){ 
		/* warmer months */
		pchar->season = 1; /* growing */
		pchar->day_flush = 0;
	}else{
		/* cooler months */
		pchar->season = 0; /* no growth */
		pchar->day_flush = 0;
	}
}

/* deciduous biomes in higher latitudes *********************/
void phenology_colddeciduous(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	double tindex, cindex;
	double crit_gdd_b, crit_cdd_b;
	double crit_gdd_e, crit_cdd_e;
	
	if((strcmp(grid->site_id, "TKY")==0)){
	/* Takayama */
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 260.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 580.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0; /*	-170.0	*/
	}else if((strcmp(grid->site_id, "TMK")==0)){
	/* Tomakomai */
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 120.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 440.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0; /*	-170.0	*/
	}else if((strcmp(grid->site_id, "TSE")==0)){
	/* Teshio */
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 60.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 320.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0; /*	-170.0	*/
	}else if((strcmp(grid->site_id, "QHB")==0)){
	/* Tibet */
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 90.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 750.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -50.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -200.0; /*	-170.0	*/
	}else if((strcmp(grid->site_id, "KSR")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 120.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 400.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -200.0;
	}else if((strcmp(grid->site_id, "KBU")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 120.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 300.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -150.0;
	}else if((strcmp(grid->site_id, "Tongyu")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 120.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 300.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -300.0;
	}else if((strcmp(grid->site_id, "CEAMIP_GDK")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 200.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 420.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -300.0;
	}else if((strcmp(grid->site_id, "CEAMIP_TUR")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 100.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 300.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -300.0;
	}else if((strcmp(grid->site_id, "CEAMIP_YLF")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 120.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 280.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -300.0;
	}else if((strcmp(grid->site_id, "CEAMIP_TMK")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 120.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 360.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -300.0;
	}else if((strcmp(grid->site_id, "CEAMIP_LSH")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 220.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 400.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -300.0;
	}else if((strcmp(grid->site_id, "CEAMIP_SKT")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 130.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 330.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -30.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -300.0;
	}else if((strcmp(grid->site_id, "CEAMIP_MMF")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 260.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 580.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0;
	}else if((strcmp(grid->site_id, "CEAMIP_MBF")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 260.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 580.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0;
	}else if((strcmp(grid->site_id, "CEAMIP_SMF")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 160.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 380.0; 
		/* critical GDD for the beginning of leaf shedding */
		crit_cdd_b = -20.0; 
		/* critical GDD for the end of leaf shedding */
		crit_cdd_e = -190.0;
	}else if((strcmp(grid->site_id, "LUCMIP0")==0)){
		/* critical GDD for the beginning of leaf flush */
		crit_gdd_b = 160.0; 
		/* critical GDD for the end of leaf flush */
		crit_gdd_e = 380.0; 
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
	
	/* crit_gdd_b = 200.0; 
	crit_gdd_e = 450.0; 
	crit_cdd_b = -20.0; 
	crit_cdd_e = -200.0;  */
	
	//tindex = pchar->gddjan;
	tindex = pchar->gdd;
	cindex = pchar->cdd;
	
	/* crit_gdd_b = pchar->gddth;
	crit_gdd_e = pchar->gddth + 300.0;
	crit_cdd_b = -20.0; 
	crit_cdd_e = -200.0; */
	
		
	if(tindex < crit_gdd_b){
		/* dormancy */
		pchar->season = 0;
		pchar->day_flush = 0;
		pchar->day_shed ++;
        //mass->gdd = 0.0;
	}else if( (tindex >= crit_gdd_b)&&(tindex < crit_gdd_e) ){
		/* leaf display */
		pchar->season = 2;
		pchar->cdd = 0.0;
		pchar->day_flush ++;
		pchar->day_shed = 0;
	}else if(tindex >= crit_gdd_e){
		/* growing period */
		pchar->season = 1;
		pchar->day_flush ++;
		
		if( (cindex < crit_cdd_b)&&(cindex >= crit_cdd_e) ){
			/* leaf shedding */
			pchar->season = 3;
			pchar->day_shed ++;
			pchar->day_flush = 0;	/* 2009/04/30 by A.Ito */
			//pchar->gdd = 0.0;
		}else if(cindex < crit_cdd_e){
			/* dormancy */
			pchar->season = 0;
			pchar->gdd = 0.0;
			pchar->day_shed ++;
			pchar->day_flush = 0;	/* 2009/04/30 by A.Ito */
		}else{
			//pchar->day_shed = 0;		/* 2009/04/30 by A.Ito */
			//pchar->cdd = 0.0;
        }
	}
}

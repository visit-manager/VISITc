/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/*   initialize environmental locations  */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

#define TERM_HYD  0.1

/* location conditions derived from the primary data (Secondary data1) *************/
void f_loct_init(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	short f, h, nn, ndy, nyr, nlyr;
	double  tem_grow;
	double atmp2m_a[366];
	double stmp10cm_a[366];
	double prate_sfc_a[366];
	
	/* species composition (if available) ************/
	if(WMODE == 1 && (strcmp(grid->site_id, "TKY")==0)){
		loct->comp_over1 = 0.407;		/* Quercus */
		loct->comp_over2 = 0.593;		/* Betula */
	}else{
		loct->comp_over1 = 1.0;
		loct->comp_over2 = 0.0;
	}
    
    loct->veg_state = 0; /* natural vegetation */
	
	loct->npp_max = 1.0;
	
	/* initial atmospheric CO2 level at Year *********/
	loct->CO2y = 1947; 
	loct->doy = 364;
	atmco2_trend(grid, loct);
	
	/*************************************************/
	/* long-term average climatology */
	for(h=0;h<366;h++){
		atmp2m_a[h] = stmp10cm_a[h] = prate_sfc_a[h] = 0.0;
	}

	nyr = nlyr = 0;
	for(f=BYR;f<=EYR;f++){
		ndy = (f%4==0)?366:365;
		nyr++;
		if(f%4==0) nlyr++;
		for(h=0;h<ndy;h++){
			atmp2m_a[h] += tmp_2m_d[f-BYR][h]-ZAT;
			stmp10cm_a[h] += tmp10_soil_d[f-BYR][h]-ZAT;
			prate_sfc_a[h] += prate_sfc_d[f-BYR][h];
		}
	}
	for(h=0;h<365;h++){
		atmp2m_a[h] *= 1.0/(double)nyr;
		stmp10cm_a[h] *= 1.0/(double)nyr;
		prate_sfc_a[h] *= 1.0/(double)nyr;
	}
	atmp2m_a[365] *= 1.0/(double)nlyr;
	stmp10cm_a[365] *= 1.0/(double)nlyr;
	prate_sfc_a[365] *= 1.0/(double)nlyr;
	
	/* 1980-1999 average */
	for(f=1980;f<2010;f++){
		for(h=0;h<365;h++){
			tmp_sfc_dav[h] += tmp_sfc_d[f-BYR][h]/30.0;
			tmp_2m_dav[h] += tmp_2m_d[f-BYR][h]/30.0;
			tmp10_soil_dav[h] += tmp10_soil_d[f-BYR][h]/30.0;
			tmp200_soil_dav[h] += tmp200_soil_d[f-BYR][h]/30.0;
			prate_sfc_dav[h] += prate_sfc_d[f-BYR][h]/30.0;
			dswrf_sfc_dav[h] += dswrf_sfc_d[f-BYR][h]/30.0;
			spfh_2m_dav[h] += spfh_2m_d[f-BYR][h]/30.0;
			tcdc_clm_dav[h] += tcdc_clm_d[f-BYR][h]/30.0;
			wind_10m_dav[h] += wind_10m_d[f-BYR][h]/30.0;
		}
	}
    tmp_sfc_dav[365] = tmp_sfc_dav[364];
    tmp_2m_dav[365] = tmp_2m_dav[364];
    tmp10_soil_dav[365] = tmp10_soil_dav[364];
    tmp200_soil_dav[365] = tmp200_soil_dav[364];
    prate_sfc_dav[365] = prate_sfc_dav[364];
    dswrf_sfc_dav[365] = dswrf_sfc_dav[364];
    spfh_2m_dav[365] = spfh_2m_dav[364];
    tcdc_clm_dav[365] = tcdc_clm_dav[364];
    wind_10m_dav[365] = wind_10m_dav[364];

	grid->atmp2m_av = grid->stmp10cm_av = 0.0;
	grid->atmp2m_max = -100.0;
	grid->atmp2m_min = 100.0;
	grid->gp_atem = grid->prate_sfc_ann = tem_grow = 0.0;
	loct->tmp_soil_mean = 0.0;
	nn = 0;
	for(h=0;h<365;h++){
		/*  annual mean temperature  */
		grid->atmp2m_av += atmp2m_a[h]/365.0;  
		grid->stmp10cm_av += stmp10cm_a[h]/365.0;  
		/*  annual  maximum  */
		grid->atmp2m_max = (atmp2m_a[h]>grid->atmp2m_max)?atmp2m_a[h]:grid->atmp2m_max;  
		/*  annual  minimum  */
		grid->atmp2m_min = (atmp2m_a[h]<grid->atmp2m_min)?atmp2m_a[h]:grid->atmp2m_min;  
		/*  annual  total  precipitation  */
		grid->prate_sfc_ann += prate_sfc_a[h];  
		
		/* growing period condition */
		if(atmp2m_a[h]>5.0){
			nn++;
			grid->gp_atem += (atmp2m_a[h]-5.0);
			tem_grow += atmp2m_a[h];
			grid->gp_pre += prate_sfc_a[h]; 
		}
		
		loct->tmp_soil_mean += (tmp10_soil_dav[h]-ZAT)/365.0;
	}
	grid->gp_tem = (nn>0)?tem_grow/(double)nn:0.0;

	/* long-term static radiation components *****************/
	for(f=0;f<366;f++){
		loct->doy = f;
		/* solar declination */
		loct->soldec[f] = f_soldec(grid, loct); 
		/* day length */
		loct->daylen[f] = f_daylen(grid, loct);	
		
		loct->day_srad[f] = 0.0;
		for(h=0;h<DSTEP;h++){
			loct->hour = h;
			/* hour angle */
			loct->hangle = -180.0+(double)h*7.5;
			/* solar height */
			loct->solhgt_h[f][h] = f_solhgt(grid, loct);
			/* atmosphere-top solar radiation */
			loct->toprad_h[f][h] = f_toprad(grid, loct); 
			
			loct->day_srad[f] += loct->toprad_h[f][h]*1800.0/1000000.0;
		}
	}
	
	/* revised phenology scheme: 2012/01/06 (A.Ito) *******************/
	(echar->tree).ncd = 0.0;
	(echar->c3).ncd = 0.0;
	(echar->c4).ncd = 0.0;
	if(grid->lat>0.0){
		/* Northern Hemisphere */
		for(f=300;f<365;f++){
			loct->tmp_2m = tmp_2m_dav[f]-ZAT;
			
			if(WMODE == 1 && (strcmp(grid->site_id, "TKY")==0)){ 
				f_clim_correct_TKY(grid, loct);		
			}
			
			if(loct->tmp_2m < 5.0){
				(echar->tree).ncd += 1.0;
				(echar->c3).ncd += 1.0;
				(echar->c4).ncd += 1.0;
			}
		}
		for(f=0;f<180;f++){
			loct->tmp_2m = tmp_2m_dav[f]-ZAT;
			
			if(WMODE == 1 && (strcmp(grid->site_id, "TKY")==0)){ 
				f_clim_correct_TKY(grid, loct);		
			}
			
			if(loct->tmp_2m < 5.0){
				(echar->tree).ncd += 1.0;
				(echar->c3).ncd += 1.0;
				(echar->c4).ncd += 1.0;
			}
		}
	}else{
		/* Southern Hemisphere */
		for(f=120;f<365;f++){
			loct->tmp_2m = tmp_2m_dav[f]-ZAT;
			
			if(loct->tmp_2m < 5.0){
				(echar->tree).ncd += 1.0;
				(echar->c3).ncd += 1.0;
				(echar->c4).ncd += 1.0;
			}
		}
	}
	
	if((echar->tree).phenoltype == 2){
		(echar->tree).gddth = -68.0 + 638.0 * exp(-0.010 * (echar->tree).ncd);
		printf("NCD  : %lf\n", (echar->tree).ncd);
		printf("GDDTH: %lf\n", (echar->tree).gddth);
	}
	if((echar->c3).phenoltype == 2){
		(echar->c3).gddth = -68.0 + 638.0 * exp(-0.010 * (echar->c3).ncd);
		printf("NCD  : %lf\n", (echar->c3).ncd);
		printf("GDDTH: %lf\n", (echar->c3).gddth);
	}
	if((echar->c4).phenoltype == 2){
		(echar->c4).gddth = -68.0 + 638.0 * exp(-0.010 * (echar->c4).ncd);
		printf("NCD  : %lf\n", (echar->c4).ncd);
		printf("GDDTH: %lf\n", (echar->c4).gddth);
	}
	
	/* initialization of water pools ****************/
	mass->sww = 0.5*grid->fieldcap;     
	mass->sww = (mass->sww>0.0)?mass->sww:0.0;

	mass->sw30 = 0.5*grid->fieldcap30;
	mass->sw30 = (mass->sw30>0.0)?mass->sw30:0.0;
	
	mass->snwa = 0.0;
	mass->snwa = (mass->snwa>0.0)?mass->snwa:0.0;
	
	/* initial local condition **********************/
	loct->climy = 1948;
	loct->CO2y = 1948;
	for(h=0;h<365;h++){
		loct->doy = h;
		f_doyTmody(1948, h, &(loct->month), &(loct->mday));
		
		/* portion of tree and sasa plants *******/
		/* local environmental conditions */
		f_loct_proc(grid, loct,  echar, mass, flux);
		
		/* phenology *****************************/
		f_growth_period(grid, loct, &(echar->tree));
		f_growth_period(grid, loct, &(echar->c3));
		f_growth_period(grid, loct, &(echar->c4));
	}
	
	/* initial soil CH4 concentration (assumed as atmpspheric level) **********/
	for(h=0;h<=SOIL_LAYER;h++){
		loct->prof_ch4[h] = loct->atm_ch4_a1[loct->CO2y - 1750]/1000.0 
			* loct->air_prsr / (UGC * (loct->tmp10_soil + ZAT));
	}
}

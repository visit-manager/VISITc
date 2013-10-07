/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/*********** functions of light environment ************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

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

/* solar declination, degree *******************************/
double f_soldec(
	struct Grid *grid, 
	struct Loct *loct
){
	double ge, aa, bb, cc, dlt;
	/* day of the year */

	ge = 2.0*PI/365.0*loct->doy;
	
	aa = -0.399912*cos(ge) + 0.070257*sin(ge);
	bb = -0.006758*cos(2.0*ge) + 0.000907*sin(2.0*ge);
	cc = -0.002697*cos(3.0*ge) + 0.00148*sin(3.0*ge);
	dlt = (0.006918 + aa + bb + cc);
	
	return (dlt);	
}

/* solar hight, degree ********************************/
double f_solhgt(
	struct Grid *grid, 
	struct Loct *loct
){
	double bbb;
	
	loct->sinb_h[loct->doy][loct->hour] = sin(grid->lat*dTr)*sin(loct->soldec[loct->doy]) 
		+cos(grid->lat*dTr)*cos(loct->soldec[loct->doy])*cos((double)(loct->hangle)*dTr);
	
	loct->sinb_h[loct->doy][loct->hour] = 
		(loct->sinb_h[loct->doy][loct->hour]<=1.0)?loct->sinb_h[loct->doy][loct->hour]:1.0; 
	loct->sinb_h[loct->doy][loct->hour] = 
		(loct->sinb_h[loct->doy][loct->hour]>=-1.0)?loct->sinb_h[loct->doy][loct->hour]:-1.0;
		
	bbb = asin(loct->sinb_h[loct->doy][loct->hour])*rTd;
	
	return (bbb);
}

/* day length, or photoperiod of vegetation production, hour *****************/
double f_daylen(
	struct Grid *grid, 
	struct Loct *loct
){
	double ho, sr, ss, dl, ha;
	
	/** hour-angle when the sun-angle equals to zero **/
	ho = -sin(grid->lat*dTr)*sin(loct->soldec[loct->doy])/
		(cos(grid->lat*dTr)*cos(loct->soldec[loct->doy])); 
	ho = (ho<=1.0)?ho:1.0; ho=(ho>=-1.0)?ho:-1.0;
	ha = acos(ho);
	
	sr = 12.0 - ha*rTd/15.0; /** time of sunrise **/
	ss = 12.0 + ha*rTd/15.0; /** time of sunset **/
	
	/** day length **/
	dl = ss - sr; 
	
	return (dl);
}

/* solar radiation at the atmosphere-top ***********************************/
double f_toprad(
	struct Grid *grid, 
	struct Loct *loct
){
	double ge, dtc, ho, solconst;
	double aa, bb;
	
	ge = 2.0*PI/365.0*loct->doy;
		
	/** relative distance between the sun and the earth **/
	/* Mizu-Knakyo-no-Kisyougaku (J. Kondo) pp 57, eq.4.5 */
	aa = 0.034221*cos(ge)+0.00128*sin(ge);
	bb = 0.000719*cos(2.0*ge)+0.000077*sin(2.0*ge);
	dtc = 1.00011+aa+bb;
		 
	/* solar constant = 4.921 MJ/m2 =1367 W/m2 = 1.96 cal/cm2/min
		 = 6151.5 micro mol photons/m2/s */	
	solconst = 1367.0; 
	
	/* holizontally incident radiation at the top of the atmosphere */
	/* Mizu-Knakyo-no-Kisyougaku (J. Kondo) pp 57, eq.4.1 */
	/* sin(solar altitude) = cos(zenith angle; tencho-kaku) */
	ho = solconst * dtc * loct->sinb_h[loct->doy][loct->hour]; 
	ho = (ho>=0.0)?ho:0.0;
	
	return (ho);
}

/* surface incident SW, PAR (W / m2), & PPFD (micro mol photon / m2 / s) ******************/
void f_ppfd(
	struct Grid *grid, 
	struct Loct *loct
){
	double kt, hd, dd;
	double e2p_b, e2p_d;
	double cloudiness, jj;
	
	/* cloudiness factor */
	cloudiness = loct->tcdc_clm;
	
	/* Black's equation */
	jj = 0.803 - 0.34*cloudiness - 0.458*cloudiness*cloudiness;
	/* jj = 0.8964-0.5392*cloudiness; */ /* new regression based on NCEP/NCAR data*/
	jj = (jj<=1.0)?jj:1.0; 
	jj = (jj>=0.0)?jj:0.0;
	
	/* if(SENSANS_SRAD==2){
		if(loct->day_srad[loct->doy] > 0.0){
			jj = loct->dswrf_sfc / (loct->day_srad[loct->doy]*1000000.0/24.0/3600.0);
			jj = (jj<=1.0)?jj:1.0; 
			jj = (jj>=0.0)?jj:0.0;
		}else{
			jj = 0.0;
		}
	} */
    
    /* case: no cloud data *********/
    /*
    if(loct->day_srad[loct->doy] > 0.0){
        jj = loct->dswrf_sfc / (loct->day_srad[loct->doy]*1000000.0/24.0/3600.0);
        jj = (jj<=1.0)?jj:1.0; 
        jj = (jj>=0.0)?jj:0.0;
    }else{
        jj = 0.0;
    }*/
	
	/* shortwave radiation ************************************/
	/* radiation at the ground, including cloud effect based on the empirical Equbal's equation */
	loct->sfcrad_h[loct->hour] = jj * loct->toprad_h[loct->doy][loct->hour];
	
	/* clear-sky cloud-free 0.803 */
	loct->cssfcrad_h[loct->hour] = 0.803 * loct->toprad_h[loct->doy][loct->hour];
	
	/* photosynthetically active radiation (par) in the global radiation, 
		based on the empirical Tooming's equation *****************************/	
	if(loct->toprad_h[loct->doy][loct->hour] > 0.0){
		/* beam / diffuse coefficients */
		kt = loct->sfcrad_h[loct->hour]/loct->toprad_h[loct->doy][loct->hour];
		hd = loct->sfcrad_h[loct->hour]*(0.958 - 0.982*kt); 
		
		/* beam */
		loct->parb_h[loct->hour] = 0.43*(loct->sfcrad_h[loct->hour] - hd);
		/* diffuse */
		loct->pard_h[loct->hour] = 0.57*hd;
		
		/* total PPFD with conventional conversion factors */
		/* constant, after McCree (1971) */
		/* e2p_b = 4.6; */ /* W/m2 to micro-mol photon /m2/s for beam radiation*/
		/* e2p_d = 4.2; */ /* W/m2 to micro-mol photon /m2/s for diffuse radiation*/
		loct->ppfd_c_h[loct->hour] = loct->parb_h[loct->hour]*4.6 + loct->pard_h[loct->hour]*4.2; 
		
		/* total PAR ***************************************/
		loct->par_h[loct->hour] = loct->parb_h[loct->hour] + loct->pard_h[loct->hour]; 

		/* updated PPFD ************************************/
		/* PAR->PPFD conversion factor after Dye (2004) */
		/* diffuse / total */
		dd = 0.958-0.982*kt;
		dd = (dd>0.01)?dd:0.01;		/* minimum 1 percent */
		
		/* diffuse*/
		e2p_d = 4.5886*dd/(0.010773+dd);
		e2p_d = (e2p_d>4.2)?e2p_d:4.2;
		loct->ppfdd_h[loct->hour] = loct->pard_h[loct->hour]*e2p_d;
		/* beam */
		e2p_b = 4.576 - 0.033144*dd;
		loct->ppfdb_h[loct->hour] = loct->parb_h[loct->hour]*e2p_b;

		/* total PPFD with Dye's conversion factors (default PPFD in VISIT) */
		loct->ppfd_h[loct->hour] = loct->ppfdb_h[loct->hour] + loct->ppfdd_h[loct->hour]; 
		
		/* cs: clear sky */
		loct->cspar_h[loct->hour] = loct->cssfcrad_h[loct->hour]*(0.43*(1.0 - 0.169454) + 0.57*0.169454); 
		loct->csparb_h[loct->hour] = loct->cssfcrad_h[loct->hour]*0.43*(1.0 - 0.169454);
		loct->cspard_h[loct->hour] = loct->cssfcrad_h[loct->hour]*0.57*0.169454;
		loct->csppfdb_h[loct->hour] = loct->cssfcrad_h[loct->hour]*0.43*(1.0 - 0.169454)*4.570383617;
		loct->csppfdd_h[loct->hour] = loct->cssfcrad_h[loct->hour]*0.57*0.169454*4.31431819;
	}else{
		loct->par_h[loct->hour] = 0.0;
		loct->ppfd_h[loct->hour] = 0.0;
		loct->ppfd_c_h[loct->hour] = 0.0;

		loct->parb_h[loct->hour] = 0.0;
		loct->pard_h[loct->hour] = 0.0;
		loct->ppfdb_h[loct->hour] = 0.0;
		loct->ppfdd_h[loct->hour] = 0.0;

		/* cs: clear sky */
		loct->cspar_h[loct->hour] = 0.0;
		loct->csparb_h[loct->hour] = 0.0;
		loct->cspard_h[loct->hour] = 0.0;
		loct->csppfdb_h[loct->hour] = 0.0;
		loct->csppfdd_h[loct->hour] = 0.0;
	}
}

/* net radiation of canopy and ground surface ******************************/
void f_net_rad(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass
){
	long f;
	double transmit;
	double in_short;
	double fsw_tree, fsw_c3, fsw_c4, fsw_ground;
	double aaa = 1.0;
	double emissiv_sfc, emissiv_atm;
	double lw_up, lw_down;
	double dav_dswrf;
	
	/* clear */
	loct->rn_long_eco = 0.0;
	loct->rn_long_tree = 0.0;
	loct->rn_long_c3 = 0.0;
	loct->rn_long_c4 = 0.0;
	loct->rn_long_ground = 0.0;
	loct->rn_short_eco = 0.0;
	loct->rn_short_tree = 0.0;
	loct->rn_short_c3 = 0.0;
	loct->rn_short_c4 = 0.0;
	loct->rn_short_ground = 0.0;
	loct->rn_eco = 0.0;
	loct->rn_tree = 0.0;
	loct->rn_c3 = 0.0;
	loct->rn_c4 = 0.0;
	loct->rn_ground = 0.0;

	/* canopy transmittance */
	transmit = 0.1;
	
	/* longwave budget *******************************************/
	/* revised 2008/06/09 */
	/* surface emissivity */
	emissiv_sfc = 0.95;
	/* upward, W/m2 */
	lw_up = aaa * emissiv_sfc * pow((loct->tmp_sfc+ZAT), 4.0)*SB_CONST;
	
	/* atmospheric emissivity: Monteith (1973) */
	emissiv_atm = 0.53 + 0.06*sqrt(loct->vp);
	/* downward, W/m2 */
	lw_down = (1.0 - loct->tcdc_clm) * emissiv_atm * pow((loct->tmp_2m+ZAT), 4.0)*SB_CONST 
			+ loct->tcdc_clm * (pow((loct->tmp_2m+ZAT), 4.0)*SB_CONST - 9.0);
	
	/* net surface (downward in positive), W/m2 */
	loct->rn_long_eco = lw_down - lw_up;
	//loct->rn_long_eco = lw_down;
	
	/* separation into surface cover */
	loct->rn_long_tree = loct->rn_long_eco * loct->fcover_tree;
	loct->rn_long_c3 = loct->rn_long_eco * loct->fcover_c3;
	loct->rn_long_c4 = loct->rn_long_eco * loct->fcover_c4;
	loct->rn_long_ground = loct->rn_long_eco * loct->fcover_ground;
	
	/* shortwave budget *****************************************/
	/* hourly solar radiation */
	dav_dswrf = 0.0;
	for(f=0;f<DSTEP;f++){
		dav_dswrf += loct->sfcrad_h[f]/(double)DSTEP;
	}
	
	in_short = dav_dswrf;
	/* in_short = loct->sfcrad_h[loct->hour]; */
	if(in_short < 0.0){
		in_short = 0.0;
	}
	
	/* fractional absorptance */
	fsw_tree = 1.0 - exp(-(echar->tree).eK*(mass->tree).lai * (1.0 - transmit));
	fsw_c3 = (1.0-fsw_tree)*loct->funder_c3*(1.0 - exp(-(echar->c3).eK*(mass->c3).lai) * (1.0 - transmit));
	fsw_c4 = (1.0-fsw_tree)*loct->funder_c4*(1.0 - exp(-(echar->c4).eK*(mass->c4).lai) * (1.0 - transmit));
	fsw_ground = 1.0 - (fsw_tree + fsw_c3 + fsw_c4);

	loct->albedo_sfc = fsw_ground*(echar->soil).albedo + fsw_tree*(echar->tree).albedo
		 + fsw_c3*(echar->c3).albedo + fsw_c4*(echar->c4).albedo;
	
	loct->rn_short_tree = in_short * fsw_tree * (1.0-(echar->tree).albedo);
	loct->rn_short_c3 = in_short * fsw_c3 * (1.0-(echar->c3).albedo);
	loct->rn_short_c4 = in_short * fsw_c4 * (1.0-(echar->c4).albedo);
	loct->rn_short_ground = in_short * fsw_ground * (1.0-(echar->soil).albedo);
	loct->rn_short_eco = loct->rn_short_tree + loct->rn_short_c3 + 
		loct->rn_short_c4 + loct->rn_short_ground;
	
	/* net radiation *******************************/
	loct->rn_eco = loct->rn_short_eco - loct->rn_long_eco; 
	loct->rn_tree = loct->rn_short_tree - loct->rn_long_tree; 
	loct->rn_c3 = loct->rn_short_c3 - loct->rn_long_c3; 
	loct->rn_c4 = loct->rn_short_c4 - loct->rn_long_c4; 
	loct->rn_ground = loct->rn_short_ground - loct->rn_long_ground; 

	/* canopy-top and -under PAR */
	(echar->tree).ppfd_t = loct->ppfd_c_h[loct->hour];
	(echar->tree).ppfd_b = (echar->tree).ppfd_t * (1.0 - fsw_tree);
	(echar->c3).ppfd_t = (echar->c4).ppfd_t = (echar->tree).ppfd_b; 
	(echar->c3).ppfd_b = (echar->c3).ppfd_t * (1.0 - (fsw_c3+fsw_c4));
	(echar->c4).ppfd_b = (echar->c4).ppfd_t * (1.0 - (fsw_c3+fsw_c4));
}

/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

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

/* DOY to Month and Day ****************************************/
void doyTmody(
	long year,
	long doy, 
	long *month, 
	long *day
){
	if(year%4==0){
		if(doy>=0 && doy<31){
			*month = 0;
			*day = doy-0;
		}else if(doy>=31 && doy<60){
			*month = 1;
			*day = doy-31;
		}else if(doy>=60 && doy<91){
			*month = 2;
			*day = doy-60;
		}else if(doy>=91&&doy<121){
			*month = 3;
			*day = doy-91;
		}else if(doy>=121&&doy<152){
			*month = 4;
			*day = doy-121;
		}else if(doy>=152&&doy<182){
			*month = 5;
			*day = doy-152;
		}else if(doy>=182&&doy<213){
			*month = 6;
			*day = doy-182;
		}else if(doy>=213&&doy<244){
			*month = 7;
			*day = doy-213;
		}else if(doy>=244&&doy<274){
			*month = 8;
			*day = doy-244;
		}else if(doy>=274&&doy<305){
			*month = 9;
			*day = doy-274;
		}else if(doy>=305&doy<335){
			*month = 10;
			*day = doy-305;
		}else if(doy>=335&&doy<366){
			*month = 11;
			*day = doy-335;
		}
	}else{
		if(doy>=0 && doy<31){
			*month = 0;
			*day = doy-0;
		}else if(doy>=31 && doy<59){
			*month = 1;
			*day = doy-31;
		}else if(doy>=59 && doy<90){
			*month = 2;
			*day = doy-59;
		}else if(doy>=90&&doy<120){
			*month = 3;
			*day = doy-90;
		}else if(doy>=120&&doy<151){
			*month = 4;
			*day = doy-120;
		}else if(doy>=151&&doy<181){
			*month = 5;
			*day = doy-151;
		}else if(doy>=181&&doy<212){
			*month = 6;
			*day = doy-181;
		}else if(doy>=212&&doy<243){
			*month = 7;
			*day = doy-212;
		}else if(doy>=243&&doy<273){
			*month = 8;
			*day = doy-243;
		}else if(doy>=273&&doy<304){
			*month = 9;
			*day = doy-273;
		}else if(doy>=304&doy<334){
			*month = 10;
			*day = doy-304;
		}else if(doy>=334&&doy<365){
			*month = 11;
			*day = doy-334;
		}
	}
}

/* Month and Day to DOY **************************************************/
void modyTdoy(
	long month, 
	long day, 
	long *doy
){
	long f;
	long mm[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	
	*doy = 0;
	for(f=0;f<month;f++){
		*doy += mm[month];
	}
	*doy += day;
}

/******************************************************/
void f_terminate(
	FILE *fp_r[N_OFILE]
){
	long h;
	
	if(NOTICE==1){
		printf("Close output files...");
	}
	
	for(h=0;h<N_OFILE;h++){
		fclose(fp_r[h]); 
	}

	if(NOTICE==1){
		printf("done\n");
	}
}

/* *********************************************************/
/* area (ha) of grid lat1-lat2 x lon1-lon2 (in degree)    */
float grid_area(
	float lat1,
	float lat2,
	float lon1,
	float lon2
){
	float area;
	float e_rad = 6378137.0; /* m, GRS-80(revised) */
	float e_exc, e_flat = 298.257;
	float l_lat, l_lon1, l_lon2;
	float m_lat, aa;
	
	if(lat1>90.0) lat1 = 90.0;
	if(lat2<-90.0) lat2 = -90.0;
	
	m_lat = (lat1 + lat2)/2.0 * PI/180.0;
	e_exc = sqrt(2.0/e_flat - 1.0/(e_flat*e_flat));
	
	aa = 1.0 - e_exc*e_exc*sin(m_lat)*sin(m_lat);
	l_lat = PI/180.0*e_rad*(1.0 - e_exc*e_exc)/pow(aa, 1.5) * fabs(lat1 - lat2);
	
	/* bug fixed 2007/12/25 E.Kato and A.Ito */
	aa = 1.0 - e_exc * e_exc*sin(lat1*PI/180.0) * sin(lat1*PI/180.0);
	l_lon1 = PI/180.0 * e_rad*cos(lat1*PI/180.0)/sqrt(aa) * fabs(lon1-lon2);
	aa = 1.0 - e_exc * e_exc*sin(lat2*PI/180.0) * sin(lat2*PI/180.0);
	l_lon2 = PI/180.0 * e_rad*cos(lat2*PI/180.0)/sqrt(aa) * fabs(lon1-lon2);
	
	area = (l_lon1 + l_lon2)*l_lat/2.0;
	
	return (area/1000000.0); /* km2 */
}

/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* Selected references
Ito, A., and Oikawa, T.: A simulation model of the carbon cycle in land ecosystems (Sim-CYCLE): 
  A description based on dry-matter production theory and plot-scale validation, 
  Ecological Modelling, 151, 147-179, 2002.
Ito, A., Saigusa, N., Murayama, S., and Yamamoto, S.: Modeling of gross and net carbon dioxide 
  exchange over a cool-temperate deciduous broad-leaved forest in Japan: Analysis of seasonal 
  and interannual change, Agricultural and Forest Meteorology, 134, 122-134, 2005.
Inatomi, M., Ito, A., Ishijima, K., and Murayama, S.: Greenhouse gas budget of a cool temperate 
  deciduous broadleaved forest in Japan estimated with a process-based model, 
  Ecosystems, 13, 472-483, 10.1007/s10021-010-9332-7, 2010.
*/

/* DOY to Month and Day ****************************************/
void f_doyTmody(
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
void f_modyTdoy(
	long month, 
	long day, 
	long *doy
){
	long f;
	long mm[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	long mm_l[12]={31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	
	*doy = 0;
	for(f=0;f<month;f++){
		*doy += mm[month];
	}
	*doy += day;
}

/***********************************************************/
void f_clim_correct_TKY(
	struct Grid *grid, 
	struct Loct *loct
){
	/*  temperature  */
	loct->tmp_2m = loct->tmp_2m * 1.1079 - 8.3206;
	loct->tmp_sfc = loct->tmp_sfc * 1.1079 - 8.3206;
	loct->tmp10_soil = loct->tmp10_soil * 1.1079 - 8.3206;
	loct->tmp200_soil = loct->tmp200_soil * 1.1079 - 8.3206;
	loct->tmp300_soil = loct->tmp300_soil * 1.1079 - 8.3206;	
	
	/* precipitation */
	switch(loct->month){
		case 0: loct->prate_sfc *= 2.3636; break;
		case 1: loct->prate_sfc *= 1.9079; break;
		case 2: loct->prate_sfc *= 1.3716; break;
		case 3: loct->prate_sfc *= 1.3368; break;
		case 4: loct->prate_sfc *= 1.1556; break;
		case 5: loct->prate_sfc *= 1.1687; break;
		case 6: loct->prate_sfc *= 0.9626; break;
		case 7: loct->prate_sfc *= 0.9846; break;
		case 8: loct->prate_sfc *= 2.0962; break;
		case 9: loct->prate_sfc *= 1.6507; break;
		case 10: loct->prate_sfc *= 1.7851; break;
		case 11: loct->prate_sfc *= 1.8242; break;
	}
}

/***********************************************************/
void f_clim_correct_QHB(
	struct Grid *grid, 
	struct Loct *loct
){
	/*  temperature  */
	loct->tmp_2m = (loct->tmp_2m + 2.9893) / 1.1212;
	loct->tmp_sfc = (loct->tmp_sfc + 2.9893) / 1.1212;
	loct->tmp10_soil = (loct->tmp10_soil + 2.9893) / 1.1212;
	loct->tmp200_soil = (loct->tmp200_soil + 2.9893) / 1.1212;
	loct->tmp300_soil = (loct->tmp300_soil + 2.9893) / 1.1212;
	
	/* precipitation */
	switch(loct->month){
		case 0: loct->prate_sfc *= 1.0000; break;
		case 1: loct->prate_sfc *= 1.0000; break;
		case 2: loct->prate_sfc *= 3.6717; break;
		case 3: loct->prate_sfc *= 2.1927; break;
		case 4: loct->prate_sfc *= 1.1906; break;
		case 5: loct->prate_sfc *= 2.1835; break;
		case 6: loct->prate_sfc *= 1.2303; break;
		case 7: loct->prate_sfc *= 1.4216; break;
		case 8: loct->prate_sfc *= 1.7679; break;
		case 9: loct->prate_sfc *= 3.8342; break;
		case 10: loct->prate_sfc *= 1.0000; break;
		case 11: loct->prate_sfc *= 1.0000; break;
	}
}

/***********************************************************/
void f_clim_correct_TMK(
	struct Grid *grid, 
	struct Loct *loct
){
	/*  temperature  */
	loct->tmp_2m = loct->tmp_2m * 0.915755053 -0.01030835;
	loct->tmp_sfc = loct->tmp_sfc * 0.915755053 -0.01030835;
	loct->tmp10_soil = loct->tmp10_soil * 0.915755053 -0.01030835;
	loct->tmp200_soil = loct->tmp200_soil * 0.915755053 -0.01030835;
	loct->tmp300_soil = loct->tmp300_soil * 0.915755053 -0.01030835;	
	
	/* precipitation */
	switch(loct->month){
		case 0: loct->prate_sfc *= 0.355032452; break;
		case 1: loct->prate_sfc *= 0.565984941; break;
		case 2: loct->prate_sfc *= 0.966911385; break;
		case 3: loct->prate_sfc *= 0.393909238; break;
		case 4: loct->prate_sfc *= 1.033119121; break;
		case 5: loct->prate_sfc *= 1.070011338; break;
		case 6: loct->prate_sfc *= 1.012129501; break;
		case 7: loct->prate_sfc *= 2.029436429; break;
		case 8: loct->prate_sfc *= 1.657190412; break;
		case 9: loct->prate_sfc *= 1.362057047; break;
		case 10: loct->prate_sfc *= 0.976141485; break;
		case 11: loct->prate_sfc *= 0.349070662; break;
	}
}

/***********************************************************/
void f_clim_correct_LSH(
	struct Grid *grid, 
	struct Loct *loct
){
	/*  temperature  */
	loct->tmp_2m = loct->tmp_2m * 0.954712713 + 1.480632709;
	loct->tmp_sfc = loct->tmp_sfc * 0.954712713 + 1.480632709;
	loct->tmp10_soil = loct->tmp10_soil * 0.954712713 + 1.480632709;
	loct->tmp200_soil = loct->tmp200_soil * 0.954712713 + 1.480632709;
	loct->tmp300_soil = loct->tmp300_soil * 0.954712713 + 1.480632709;	
	
	/* precipitation */
	switch(loct->month){
		case 0: loct->prate_sfc *= 0.395212285; break;
		case 1: loct->prate_sfc *= 0.453657175; break;
		case 2: loct->prate_sfc *= 0.681015961; break;
		case 3: loct->prate_sfc *= 0.053505452; break;
		case 4: loct->prate_sfc *= 1.023523524; break;
		case 5: loct->prate_sfc *= 0.203850553; break;
		case 6: loct->prate_sfc *= 1.213680203; break;
		case 7: loct->prate_sfc *= 1.176242353; break;
		case 8: loct->prate_sfc *= 0.911712151; break;
		case 9: loct->prate_sfc *= 2.326388889; break;
		case 10: loct->prate_sfc *= 1.574293488; break;
		case 11: loct->prate_sfc *= 1.338528414; break;
	}
}

/***********************************************************/
void f_clim_correct_TSE(
	struct Grid *grid, 
	struct Loct *loct
){
	/*  temperature  */
	loct->tmp_2m = loct->tmp_2m * 1.013362273 -1.612987859;
	loct->tmp_sfc = loct->tmp_sfc * 1.013362273 -1.612987859;
	loct->tmp10_soil = loct->tmp10_soil * 1.013362273 -1.612987859;
	loct->tmp200_soil = loct->tmp200_soil * 1.013362273 -1.612987859;
	loct->tmp300_soil = loct->tmp300_soil * 1.013362273 -1.612987859;	
	
	/* precipitation */
	switch(loct->month){
		case 0: loct->prate_sfc *= 0.1228533; break;
		case 1: loct->prate_sfc *= 0.450010628; break;
		case 2: loct->prate_sfc *= 1.178791378; break;
		case 3: loct->prate_sfc *= 1.212195924; break;
		case 4: loct->prate_sfc *= 1.242690058; break;
		case 5: loct->prate_sfc *= 1.134526207; break;
		case 6: loct->prate_sfc *= 1.442785638; break;
		case 7: loct->prate_sfc *= 1.042956257; break;
		case 8: loct->prate_sfc *= 2.656665104; break;
		case 9: loct->prate_sfc *= 1.06587352; break;
		case 10: loct->prate_sfc *= 0.721634718; break;
		case 11: loct->prate_sfc *= 0.206349206; break;
	}
}

/***********************************************************/
void f_clim_correct_SKT(
	struct Grid *grid, 
	struct Loct *loct
){
	/*  temperature  */
	loct->tmp_2m = loct->tmp_2m * 0.873677879 -0.143452211;
	loct->tmp_sfc = loct->tmp_sfc * 0.873677879 -0.143452211;
	loct->tmp10_soil = loct->tmp10_soil * 0.873677879 -0.143452211;
	loct->tmp200_soil = loct->tmp200_soil * 0.873677879 -0.143452211;
	loct->tmp300_soil = loct->tmp300_soil * 0.873677879 -0.143452211;	
	
	/* precipitation */
	switch(loct->month){
		case 0: loct->prate_sfc *= 0; break;
		case 1: loct->prate_sfc *= 0; break;
		case 2: loct->prate_sfc *= 0.009645062; break;
		case 3: loct->prate_sfc *= 0.049042687; break;
		case 4: loct->prate_sfc *= 0.547092014; break;
		case 5: loct->prate_sfc *= 0.698796229; break;
		case 6: loct->prate_sfc *= 3.104236563; break;
		case 7: loct->prate_sfc *= 1.580782089; break;
		case 8: loct->prate_sfc *= 1.10410575; break;
		case 9: loct->prate_sfc *= 0.568668379; break;
		case 10: loct->prate_sfc *= 0; break;
		case 11: loct->prate_sfc *= 0; break;
	}
}

/***********************************************************/
void f_clim_correct_MKL(
	struct Grid *grid, 
	struct Loct *loct
){
	/*  temperature  */
	loct->tmp_2m = loct->tmp_2m * 0.647725457 +8.990984673;
	loct->tmp_sfc = loct->tmp_sfc * 0.647725457 +8.990984673;
	loct->tmp10_soil = loct->tmp10_soil * 0.647725457 +8.990984673;
	loct->tmp200_soil = loct->tmp200_soil * 0.647725457 +8.990984673;
	loct->tmp300_soil = loct->tmp300_soil * 0.647725457 +8.990984673;	
	
	/* precipitation */
	switch(loct->month){
		case 0: loct->prate_sfc *= 0.019938112; break;
		case 1: loct->prate_sfc *= 3.262190762; break;
		case 2: loct->prate_sfc *= 0.834615657; break;
		case 3: loct->prate_sfc *= 0.617358108; break;
		case 4: loct->prate_sfc *= 0.701498769; break;
		case 5: loct->prate_sfc *= 0.443937088; break;
		case 6: loct->prate_sfc *= 0.487206867; break;
		case 7: loct->prate_sfc *= 0.28702829; break;
		case 8: loct->prate_sfc *= 0.370097667; break;
		case 9: loct->prate_sfc *= 0.1417892; break;
		case 10: loct->prate_sfc *= 0; break;
		case 11: loct->prate_sfc *= 0; break;
	}
}

/***********************************************************/
void f_clim_correct_KBU(
   struct Grid *grid, 
   struct Loct *loct
){
	double t_2m_ncep_k, t_2m_site_k;
	
	//t_2m_ncep_k = loct->tmp_2m + ZAT;
	
	//t_2m_site_k = t_2m_ncep_k * 0.9552 + 12.7404;
	
	/*  temperature  */
	//loct->tmp_2m += (t_2m_site_k - t_2m_ncep_k);
	//loct->tmp_2m = t_2m_site_k - ZAT;
	//loct->tmp_sfc += (t_2m_site_k - t_2m_ncep_k);
	//loct->tmp10_soil += (t_2m_site_k - t_2m_ncep_k);
	//loct->tmp200_soil += (t_2m_site_k - t_2m_ncep_k);
	//loct->tmp300_soil += (t_2m_site_k - t_2m_ncep_k);
	
	/* precipitation */
	/* if(loct->climy>=1951 && loct->climy<=2007){
		loct->prate_sfc = prec_altdata[loct->climy-1951][loct->doy];
	}else{
		loct->prate_sfc *= 1.8118;
	} */
}

/***********************************************************/
void f_clim_correct_Tongyu(
	struct Grid *grid, 
	struct Loct *loct
){
	double t_2m_ncep_k, t_2m_site_k;
	
	t_2m_ncep_k = loct->tmp_2m + ZAT;
	
	//t_2m_site_k = t_2m_ncep_k * 0.9666 + 10.4151;
	
	/*  temperature  */
	/* loct->tmp_2m += (t_2m_site_k - t_2m_ncep_k);
	loct->tmp_sfc += (t_2m_site_k - t_2m_ncep_k);
	loct->tmp10_soil += (t_2m_site_k - t_2m_ncep_k);
	loct->tmp200_soil += (t_2m_site_k - t_2m_ncep_k);
	loct->tmp300_soil += (t_2m_site_k - t_2m_ncep_k); */
	
	/* precipitation */
	/* if(loct->climy>=1951 && loct->climy<=2007){
		loct->prate_sfc = prec_altdata[loct->climy-1951][loct->doy];
	}else{
		loct->prate_sfc *= 0.644;
	} */
}

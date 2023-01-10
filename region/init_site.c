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

/**  initialize environmental characteristics of each site   **/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* climatology: *_d[] means the time-series during 1948 to 2004- */	
extern float tmp_sfc_d[PERIOD][366];		/* ground surface temperature, degree Celcius */
extern float tmp_2m_d[PERIOD][366];		/* 2m air temperature, degree Celcius */
extern float tmp10_soil_d[PERIOD][366];	/* soil temperature at 10 cm depth, degree Celcius */
extern float tmp200_soil_d[PERIOD][366];	/* soil temperature at 200 cm depth, degree Celcius */
extern float dswrf_sfc_d[PERIOD][366];		/* downward shortwave radiation at the surface, W m-2 */
extern float tcdc_clm_d[PERIOD][366];		/* total cloudiness, fraction */
extern float prate_sfc_d[PERIOD][366];		/* precipitation, mm mon-1 */
extern float spfh_2m_d[PERIOD][366];		/* specific humidity, kg kg-1 */
extern float wind_10m_d[PERIOD][366];		/* u-wind velocity, m s-1 */
extern float vpd_d[PERIOD][366];   		/* VPD, hPa */
extern long	month_day[12];
extern long	WMODE;

/*  initialize  grid  conditions   **************************************************/
/**  initialize climate and soil parameters by using the prepared data files  **/
void init_site(
	struct Grid *grid
){
	long e, f, year, day, ndy;
	float bbb, uwnd, vwnd;
	float tair, vsat, vpres;
	char filename[128];
	char verid[128];
	char dateid[128];
	FILE *fp_site, *fp_clim;
		
	/* input data is available: site */
	grid->flag_datavl = 1;

	/* site condition ***********************************************************/
	/* site description file */
	strcpy(filename, grid->file_site);
	if((fp_site = fopen(filename,"rt")) == NULL){
		printf("No site data !:%s\n", filename);
		exit (1);
	}

	/* version of parameter file */
	fscanf(fp_site,"%s", verid);
	
	/* updated date of parameter file */
	fscanf(fp_site,"%s", dateid);

	/* site ID */
	fscanf(fp_site,"%s", &(grid->site_id));
	printf("%s ", grid->site_id);
	/*
		TKY : Takayama temperate deciduous broadleaved forest
		TKC : Takayama conifer forest
		TMK : Tomakomai larch forest
		SGD : Sugadaira experimental grassland
		FJY : Fujiyoshida pine forest
		PSH : Pasoh tropical rain forest
		QHB : Qinghai - Tibetan alpine meadow
		MSE : Mase paddy field
		LSH : Laoshan larch forest
	*/
	
	/* vegetation type */
	fscanf(fp_site,"%d", &(grid->veg_type));
	printf("(SIGBP type:%d)", grid->veg_type);
	/*  grid->veg_type = 4;  */
	/*  veg type by MOD12
	0: Water
	1: Evergreen needleleaf forest
	2: Evergreen broadleaf forest
	3: Deciduous needleleaf forest
	4: Deciduous broadleaf forest
	5: Mixed forest
	6: Closed shrublands
	7: Open shrublands
	8: Woody savanna
	9: Savanna
	10: Grassland
	11: Permanent wetland
	12: Cropland
	13: Urban and build-up
	14: Cropland/natural vegetation mosaic
	15: Snow and ice
	16: Barren or sparsely vegetated
	17: Water? 
	*/
	
	/* area, ha */	
	grid->area = 1.0; 			
	/* latitude, degree */
	fscanf(fp_site,"%f", &(grid->lat));
	/* longitude. degree */
	fscanf(fp_site,"%f", &(grid->lon));
	/*  topography, m above sea level  */
	fscanf(fp_site,"%f", &(grid->topo));

	/* input  soil properties *********************/
	/* soil texture for N2O model  */
	fscanf(fp_site,"%ld", &(grid->stexture));
		/* 0: sandy */
		/* 1: medium */
		/* 2: fine */
	/*  water hoding capacity (WHC),  30cm soil  */
	fscanf(fp_site,"%f", &(grid->fieldcap30));
	/*  water hoding capacity (WHC),  whole  */
	fscanf(fp_site,"%f", &(grid->fieldcap));
	/*  soil rooting depth  */
	fscanf(fp_site,"%f", &(grid->soildpth));
	/*  hydraulic conductivity (k_s), in m/s  */
	fscanf(fp_site,"%f", &(grid->hyd_cond));
	/*  bulk density  */
	fscanf(fp_site,"%f", &(grid->bulkdens));
	/* soil pH */
	fscanf(fp_site,"%f", &(grid->soil_ph));
	/* sand fraction */
	fscanf(fp_site,"%f", &(grid->sand_frac));
	/* clay fraction */
	fscanf(fp_site,"%f", &(grid->clay_frac));
	
	/* check */
	fscanf(fp_site, "%f", &bbb);
	if((long)bbb!=9999999){
		printf("!! ERROR: broken site file !!\n");
		exit(1);
	}

	fclose(fp_site);
	
	/**  input daily climate condition  **********************************/	
	/* NCEP/NCAR reanalysis file */
	if((fp_clim = fopen(grid->file_ncep,"rt")) == NULL){
		printf("No input NCEP data !\n");
		exit (1);
	}
	
	for(f=BYR;f<=EYR;f++){
		ndy = (f%4==0)?366:365;
		for(e=0;e<ndy;e++){
			/* year and day */
			fscanf(fp_clim,"%ld %ld", &year, &day);  
			/* air temperature at 2m */
			fscanf(fp_clim,"%f",&bbb);  tmp_2m_d[f-BYR][e] = bbb;
			/* maximum temperature */
			fscanf(fp_clim,"%f",&bbb);  
			/* minimum temperature */
			fscanf(fp_clim,"%f",&bbb);  
			/* specific humidity at 2m */
			fscanf(fp_clim,"%f",&bbb);  spfh_2m_d[f-BYR][e] = bbb;
			/* precipitation */
			fscanf(fp_clim,"%f",&bbb);  prate_sfc_d[f-BYR][e] = bbb;
			/* downward shortwave radiation */
			fscanf(fp_clim,"%f",&bbb);  dswrf_sfc_d[f-BYR][e] = bbb;
			/* cloud cover */
			fscanf(fp_clim,"%f",&bbb);  tcdc_clm_d[f-BYR][e] = bbb;
			/* surface temperature */
			fscanf(fp_clim,"%f",&bbb);  tmp_sfc_d[f-BYR][e] = bbb;
			/* soil temperature 0-10 cm */
			fscanf(fp_clim,"%f",&bbb);  tmp10_soil_d[f-BYR][e] = bbb;
			/* soil temperature 10-200 cm */
			fscanf(fp_clim,"%f",&bbb);  tmp200_soil_d[f-BYR][e] = bbb;
			/* soil temperature 300 cm */
			fscanf(fp_clim,"%f",&bbb);
			/* soil water content, 0-10 cm */
			fscanf(fp_clim,"%f",&bbb);  
			/* soil water content, 10-200 cm */
			fscanf(fp_clim,"%f",&bbb);  
			/* u wind */
			fscanf(fp_clim,"%f",&uwnd); 
			/* v wind */
			fscanf(fp_clim,"%f",&vwnd);  
			wind_10m_d[f-BYR][e] = sqrt(uwnd*uwnd + vwnd*vwnd);
			/* snow depth */
			fscanf(fp_clim,"%f",&bbb);  
			/* air pressure (Pa) */
			fscanf(fp_clim,"%f",&bbb);  
			
			/* vapor pressure deficit */
			tair = tmp_2m_d[f-BYR][e] -ZAT;  /* deg C */
			bbb *= 0.01;	/* hPa */
			if(tmp_2m_d[f-BYR][e]>=0.0){
				vsat = 6.1078*pow(10.0, 7.5*tair/(237.3 + tair));
			}else{
				vsat = 6.1078*pow(10.0, 9.5*tair/(265.3 + tair));
			}
			vpres = bbb*spfh_2m_d[f-BYR][e]/(0.622 + 0.378*spfh_2m_d[f-BYR][e]);
			
			vpd_d[f-BYR][e] = vsat - vpres;
		}
	}
	fclose(fp_clim);

	/* soil parameterizatrions using Saxton (1986) *********/
	f_soil_saxton(grid);
	
	grid->fieldcap30 = grid->field_cap;
	grid->fieldcap = grid->fieldcap - grid->fieldcap30;
	
	if(NOTICE==1){
		printf("\nSAND: %f\n", grid->sand_frac);
		printf("CLAY: %f\n", grid->clay_frac);
		printf("A: %f\n", grid->a_sw);
		printf("B: %f\n", grid->b_sw);
		printf("FC: %f\n", grid->field_cap);
		printf("PS: %f\n", grid->pore_cap);
	}
}

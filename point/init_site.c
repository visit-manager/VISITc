/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/**  initialize environmental characteristics of each site   **/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* initialize  grid  conditions **************************************************/
/* initialize climate and soil parameters by using the prepared data files  **/
void f_init_site(
	struct Grid *grid
){
	long e, f, g, year, day, ndy;
	double bbb, uwnd, vwnd, apres;
	double tair, vsat, vpres, cloudiness;
	char filename[128];
	char verid[128];
	char dateid[128];
	FILE *fp_site, *fp_clim;
		
	/* input data is available */
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
	fscanf(fp_site,"%ld", &(grid->veg_type));
	printf("(SIGBP type:%ld)", grid->veg_type);
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
	fscanf(fp_site,"%lf", &(grid->lat));
	/* longitude. degree */
	fscanf(fp_site,"%lf", &(grid->lon));
	/*  topography, m above sea level  */
	fscanf(fp_site,"%lf", &(grid->topo));

	/* input  soil properties *********************/
	/* soil texture for N2O model  */
	fscanf(fp_site,"%ld", &(grid->stexture));
		/* 0: sandy */
		/* 1: medium */
		/* 2: fine */
	/*  water hoding capacity (WHC),  30cm soil  */
	fscanf(fp_site,"%lf", &(grid->fieldcap30));
	/*  water hoding capacity (WHC),  whole  */
	fscanf(fp_site,"%lf", &(grid->fieldcap));
	/*  soil rooting depth  */
	fscanf(fp_site,"%lf", &(grid->soildpth));
	/*  hydraulic conductivity (k_s), in m/s  */
	fscanf(fp_site,"%lf", &(grid->hyd_cond));
	/*  bulk density  */
	fscanf(fp_site,"%lf", &(grid->bulkdens));
	/* soil pH */
	fscanf(fp_site,"%lf", &(grid->soil_ph));
	/* sand fraction */
	fscanf(fp_site,"%lf", &(grid->sand_frac));
	/* clay fraction */
	fscanf(fp_site,"%lf", &(grid->clay_frac));
	
	/* check */
	fscanf(fp_site, "%lf", &bbb);
	if((long)bbb!=9999999){
		printf("!! ERROR: broken site file !!\n");
		exit(1);
	}

	fclose(fp_site);
	
	/**  input daily climate condition  **********************************/	
	if((fp_clim = fopen(grid->file_clim,"rt")) == NULL){
		printf("No input CLIMATE data !\n");
		exit (1);
	}
	
	if(strcmp(grid->site_id, "CEAMIP_TKY")==0 || strcmp(grid->site_id, "CEAMIP_TUR")==0 
        || strcmp(grid->site_id, "CEAMIP_TKC")==0|| strcmp(grid->site_id, "CEAMIP_SKR")==0
        || strcmp(grid->site_id, "CEAMIP_TMK")==0|| strcmp(grid->site_id, "CEAMIP_LSH")==0
        || strcmp(grid->site_id, "CEAMIP_MMF")==0|| strcmp(grid->site_id, "CEAMIP_MBF")==0
        || strcmp(grid->site_id, "CEAMIP_MKL")==0|| strcmp(grid->site_id, "CEAMIP_YLF")==0
        || strcmp(grid->site_id, "CEAMIP_YPF")==0|| strcmp(grid->site_id, "CEAMIP_SKT")==0
        || strcmp(grid->site_id, "CEAMIP_GDK")==0|| strcmp(grid->site_id, "CEAMIP_SMF")==0
        || strcmp(grid->site_id, "CEAMIP_TSE")==0|| strcmp(grid->site_id, "CEAMIP_CBS")==0
        || strcmp(grid->site_id, "CEAMIP_QYZ")==0|| strcmp(grid->site_id, "CEAMIP_PDF")==0
        || strcmp(grid->site_id, "CEAMIP_KBU")==0|| strcmp(grid->site_id, "CEAMIP_QHB")==0
        || strcmp(grid->site_id, "CEAMIP_HBG")==0|| strcmp(grid->site_id, "CEAMIP_YCS")==0
        || strcmp(grid->site_id, "CEAMIP_HFK")==0|| strcmp(grid->site_id, "CEAMIP_MSE")==0){
		/* CEAMIP file */
		for(f=BYR;f<=EYR;f++){
			for(e=0;e<365;e++){
				/* year and day */
				fscanf(fp_clim,"%ld %ld",&year,&day);  
				/* air temperature at 2m: degC => degK */
				fscanf(fp_clim,"%lf",&bbb);  tmp_2m_d[f-BYR][e] = bbb + ZAT;
				/* maximum temperature */
				fscanf(fp_clim,"%lf",&bbb);  
				/* minimum temperature */
				fscanf(fp_clim,"%lf",&bbb);  
				/* precipitation */
				fscanf(fp_clim,"%lf",&bbb);  prate_sfc_d[f-BYR][e] = bbb;
				
				/* relative humidity at 2m */
				fscanf(fp_clim,"%lf",&bbb);  spfh_2m_d[f-BYR][e] = bbb;
				/* VPD at 2m */
				fscanf(fp_clim,"%lf",&bbb);  vpd_d[f-BYR][e] = bbb/100.0;
				
				/* downward shortwave radiation */
				fscanf(fp_clim,"%lf",&bbb);  dswrf_sfc_d[f-BYR][e] = bbb;
				/* downward longwave radiation */
				fscanf(fp_clim,"%lf",&bbb);
				
				/* wind */
				fscanf(fp_clim,"%lf",&bbb);  wind_10m_d[f-BYR][e] = bbb;
								
				tcdc_clm_d[f-BYR][e] = 50.0;
				tmp_sfc_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
				tmp10_soil_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
				tmp200_soil_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
			}
			tmp_2m_d[f-BYR][365] = tmp_2m_d[f-BYR][364];
			prate_sfc_d[f-BYR][365] = prate_sfc_d[f-BYR][364];
			spfh_2m_d[f-BYR][365] = spfh_2m_d[f-BYR][364];
			vpd_d[f-BYR][365] = vpd_d[f-BYR][364];
			dswrf_sfc_d[f-BYR][365] = dswrf_sfc_d[f-BYR][364];
			wind_10m_d[f-BYR][365] = wind_10m_d[f-BYR][364];
			tcdc_clm_d[f-BYR][365] = tcdc_clm_d[f-BYR][364];
			tmp_sfc_d[f-BYR][365] = tmp_sfc_d[f-BYR][364];
			tmp10_soil_d[f-BYR][365] = tmp10_soil_d[f-BYR][364];
			tmp200_soil_d[f-BYR][365] = tmp200_soil_d[f-BYR][364];
		}
	}else if(strcmp(grid->site_id, "LUCMIP0")==0 || strcmp(grid->site_id, "LUCMIP1")==0 ||
            strcmp(grid->site_id, "LUCMIP2")==0 || strcmp(grid->site_id, "LUCMIP3")==0){
        
		for(f=BYR;f<=EYR;f++){
			for(e=0;e<365;e++){
            
                dswrf_sfc_d[f-BYR][e] = 0.0;
                tmp_2m_d[f-BYR][e] = 0.0;
                prate_sfc_d[f-BYR][e] = 0.0;
                spfh_2m_d[f-BYR][e] = 0.0;
                wind_10m_d[f-BYR][e] = 0.0;
                vpd_d[f-BYR][e] = 0.0;
                apres = 0.0;
                for(g=0;g<4;g++){
                    /* year and day */
                    fscanf(fp_clim,"%ld %ld",&year,&day);
                    
                    /* downward longwave */
                    fscanf(fp_clim,"%lf",&bbb);  
                    
                    /* downward shortwave radiation: W/m2 */
                    fscanf(fp_clim,"%lf",&bbb);  dswrf_sfc_d[f-BYR][e] += bbb / 4.0;
                    
                    /* air temperature at 2m: K */
                    fscanf(fp_clim,"%lf",&bbb);  tmp_2m_d[f-BYR][e] += bbb / 4.0;

                    /* precipitation */
                    fscanf(fp_clim,"%lf",&bbb);  prate_sfc_d[f-BYR][e] += bbb;

                    /* specific humidity at 2m */
                    fscanf(fp_clim,"%lf",&bbb);  spfh_2m_d[f-BYR][e] += bbb / 4.0;

                    /* air pressure at 2m */
                    fscanf(fp_clim,"%lf",&bbb);        apres += bbb / 4.0;

                    /* wind */
                    fscanf(fp_clim,"%lf",&uwnd);
                    fscanf(fp_clim,"%lf",&vwnd);
                    
                    wind_10m_d[f-BYR][e] += sqrt(uwnd*uwnd + vwnd*vwnd)/4.0;
                }
            
                /* temperatures */
                tmp_sfc_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
                tmp10_soil_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
                tmp200_soil_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
                
                /* cloudiness */
                tcdc_clm_d[f-BYR][e] = 50.0; /* re-calculated after */

                /* vapor pressure deficit */
                tair = tmp_2m_d[f-BYR][e] -ZAT;  /* deg C */
                apres *= 0.01;	/* hPa */
                if(tair >= 0.0){
                    vsat = 6.1078 * pow(10.0, 7.5*tair/(237.3 + tair));
                }else{
                    vsat = 6.1078 * pow(10.0, 9.5*tair/(265.3 + tair));
                }
                vpres = apres * spfh_2m_d[f-BYR][e]/(0.622 + 0.378*spfh_2m_d[f-BYR][e]);
                vpd_d[f-BYR][e] += (vsat - vpres)/4.0;
            }
        }
    
    }else if(EX_ASIAMIP == 1){
        
        for(f=BYR;f<=EYR;f++){
			for(e=0;e<365;e++){
				/* year and day */
				fscanf(fp_clim,"%ld %ld",&year,&day);  
				/* air temperature at 2m: degC => degK */
				fscanf(fp_clim,"%lf",&bbb);  tmp_2m_d[f-BYR][e] = bbb + ZAT;
				/* maximum temperature */
				fscanf(fp_clim,"%lf",&bbb);  
				/* minimum temperature */
				fscanf(fp_clim,"%lf",&bbb);  
				/* precipitation */
				fscanf(fp_clim,"%lf",&bbb);  prate_sfc_d[f-BYR][e] = bbb;
				
				/* relative humidity at 2m */
				fscanf(fp_clim,"%lf",&bbb);  spfh_2m_d[f-BYR][e] = bbb;
				/* VPD at 2m */
				fscanf(fp_clim,"%lf",&bbb);  vpd_d[f-BYR][e] = bbb/100.0;
				
				/* downward shortwave radiation, W/m2 */
				fscanf(fp_clim,"%lf",&bbb);  dswrf_sfc_d[f-BYR][e] = bbb;
				/* downward longwave radiation, W/m2 */
				fscanf(fp_clim,"%lf",&bbb);
				
				/* wind.¥,  */
				fscanf(fp_clim,"%lf",&bbb);  wind_10m_d[f-BYR][e] = bbb;
								
				tcdc_clm_d[f-BYR][e] = 50.0;
				tmp_sfc_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
				tmp10_soil_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
				tmp200_soil_d[f-BYR][e] = tmp_2m_d[f-BYR][e];
			}
			tmp_2m_d[f-BYR][365] = tmp_2m_d[f-BYR][364];
			prate_sfc_d[f-BYR][365] = prate_sfc_d[f-BYR][364];
			spfh_2m_d[f-BYR][365] = spfh_2m_d[f-BYR][364];
			vpd_d[f-BYR][365] = vpd_d[f-BYR][364];
			dswrf_sfc_d[f-BYR][365] = dswrf_sfc_d[f-BYR][364];
			wind_10m_d[f-BYR][365] = wind_10m_d[f-BYR][364];
			tcdc_clm_d[f-BYR][365] = tcdc_clm_d[f-BYR][364];
			tmp_sfc_d[f-BYR][365] = tmp_sfc_d[f-BYR][364];
			tmp10_soil_d[f-BYR][365] = tmp10_soil_d[f-BYR][364];
			tmp200_soil_d[f-BYR][365] = tmp200_soil_d[f-BYR][364];
        }
    
    }else{
		/* NCEP/NCAR reanalysis file */
		for(f=BYR;f<=EYR;f++){
			ndy = (f%4==0)?366:365;
			/* if(strcmp(grid->site_id, "KBU")==0){ ndy = 365; } */
            
			for(e=0;e<ndy;e++){
				/* year and day */
				fscanf(fp_clim,"%ld %ld",&year,&day);  
				/* air temperature at 2m */
				fscanf(fp_clim,"%lf",&bbb);  tmp_2m_d[f-BYR][e] = bbb;
				/* maximum temperature */
				fscanf(fp_clim,"%lf",&bbb);  
				/* minimum temperature */
				fscanf(fp_clim,"%lf",&bbb);  
				/* specific humidity at 2m */
				fscanf(fp_clim,"%lf",&bbb);  spfh_2m_d[f-BYR][e] = bbb;
				/* precipitation */
				fscanf(fp_clim,"%lf",&bbb);  prate_sfc_d[f-BYR][e] = bbb*24.0*3600.0;
                if(strcmp(grid->site_id, "PSO")==0){
                    prate_sfc_d[f-BYR][e] *= 0.45;
                }
                
				/* downward shortwave radiation */
				fscanf(fp_clim,"%lf",&bbb);  dswrf_sfc_d[f-BYR][e] = bbb;
				/* cloud cover */
				fscanf(fp_clim,"%lf",&bbb);  tcdc_clm_d[f-BYR][e] = bbb;
				/* surface temperature */
				fscanf(fp_clim,"%lf",&bbb);  tmp_sfc_d[f-BYR][e] = bbb;
				/* soil temperature 0-10 cm */
				fscanf(fp_clim,"%lf",&bbb);  tmp10_soil_d[f-BYR][e] = bbb;
				/* soil temperature 10-200 cm */
				fscanf(fp_clim,"%lf",&bbb);  tmp200_soil_d[f-BYR][e] = bbb;
				/* soil temperature 300 cm */
				fscanf(fp_clim,"%lf",&bbb);
				/* soil water content, 0-10 cm */
				fscanf(fp_clim,"%lf",&bbb);  
				/* soil water content, 10-200 cm */
				fscanf(fp_clim,"%lf",&bbb);  
				/* u wind */
				fscanf(fp_clim,"%lf",&uwnd); 
				/* v wind */
				fscanf(fp_clim,"%lf",&vwnd);  
				wind_10m_d[f-BYR][e] = sqrt(uwnd*uwnd + vwnd*vwnd);
				/* snow depth */
				fscanf(fp_clim,"%lf",&bbb);  
				/* air pressure (Pa) */
				fscanf(fp_clim,"%lf",&bbb);  apres = bbb;
				
                if(strcmp(grid->site_id, "GSM")==0){
                    tmp_2m_d[f-BYR][e] -= 300.0 *0.006;
                    tmp_sfc_d[f-BYR][e] -= 300.0 *0.006;
                    tmp10_soil_d[f-BYR][e] -= 300.0 *0.006;
                    tmp200_soil_d[f-BYR][e] -= 300.0 *0.006;
                }

				/* vapor pressure deficit */
				tair = tmp_2m_d[f-BYR][e] -ZAT;  /* deg C */
				apres *= 0.01;	/* hPa */
				if(tair>=0.0){
					vsat = 6.1078*pow(10.0, 7.5*tair/(237.3+tair));
				}else{
					vsat = 6.1078*pow(10.0, 9.5*tair/(265.3+tair));
				}
				vpres = apres*spfh_2m_d[f-BYR][e]/(0.622+0.378*spfh_2m_d[f-BYR][e]);
				
				vpd_d[f-BYR][e] = vsat - vpres;
                
			}
            /* if(strcmp(grid->site_id, "KBU")==0){ 
                tmp_2m_d[f-BYR][365] = tmp_2m_d[f-BYR][364];
                prate_sfc_d[f-BYR][365] = prate_sfc_d[f-BYR][364];
                spfh_2m_d[f-BYR][365] = spfh_2m_d[f-BYR][364];
                vpd_d[f-BYR][365] = vpd_d[f-BYR][364];
                dswrf_sfc_d[f-BYR][365] = dswrf_sfc_d[f-BYR][364];
                wind_10m_d[f-BYR][365] = wind_10m_d[f-BYR][364];
                tcdc_clm_d[f-BYR][365] = tcdc_clm_d[f-BYR][364];
                tmp_sfc_d[f-BYR][365] = tmp_sfc_d[f-BYR][364];
                tmp10_soil_d[f-BYR][365] = tmp10_soil_d[f-BYR][364];
                tmp200_soil_d[f-BYR][365] = tmp200_soil_d[f-BYR][364];
            } */
		}
	}
	fclose(fp_clim);

	/* soil parameterizatrions using Saxton (1986) *********/
	f_soil_saxton(grid);
	
	grid->fieldcap30 = grid->field_cap;
	grid->fieldcap = grid->fieldcap - grid->fieldcap30;
    //grid->fieldcap = grid->field_cap * 1.0/0.3;
	
	if(NOTICE==1){
		printf("\nSAND: %lf\n", grid->sand_frac);
		printf("CLAY: %lf\n", grid->clay_frac);
		printf("A: %lf\n", grid->a_sw);
		printf("B: %lf\n", grid->b_sw);
		printf("FC: %lf\n", grid->field_cap);
		printf("PS: %lf\n", grid->pore_cap);
	}
}

/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/** initialize environmental characteristics of each 1deg x 1deg grid **/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* initialize 1deg x 1deg region conditions ***********************************/
/** initialize climate and soil parameters by using the prepared data files **/
void f_init_region(
	struct Grid grid[DROW][DCOL]
){
	long f, g, i, j, year, month, day, doy, ndy, nceprow, ncepcol;
	double ndat[10];
	double tair, vsat, vpres, dist, dist_min, glat, glon, aprsr;
	char gridid[20];
	char ord[]="-";
	char clat[10], clon[10];
	double nlat[16] = {54.2846, 52.3799, 50.4752, 48.5705, 46.6658, 44.7611, 
		42.8564, 40.9517, 39.0470, 37.1422, 35.2375, 33.3328, 31.4281, 
		29.5234, 27.6186, 25.7139};
	double nlon[19] = {120.9375, 122.8125, 124.6875, 126.5625, 128.4375, 
		130.3125, 132.1875, 134.0625, 135.9375, 137.8125, 139.6875, 
		141.5625, 143.4375, 145.3125, 147.1875, 149.0625, 150.9375, 
		152.8125, 154.6875};
	double prec_wclim[12], tmin_wclim[12], tmax_wclim[12];
	FILE *fp_s1dat, *fp_ncepclim[10];
	
	/* S1 data ***************************************************/
	strcpy(gridid,"s1dat2b_");
	snprintf(clat, 10, "%ld", WROW);
	snprintf(clon, 10, "%ld", WCOL);
	strcat(gridid, clat); 
	strcat(gridid, ord); 
	strcat(gridid, clon);
	
	if((fp_s1dat = fopen(gridid,"rt")) == NULL){
		printf("!! NO s1 input data: %s!\n",gridid);
		exit (1);
	}
	
	strcpy(gridid,"RCB_");
	strcat(gridid, clat); 
	strcpy(ord,"-");
	strcat(gridid, ord); 
	strcat(gridid, clon);
	for(i=0; i<NROW; i++){ /* 120 */
		for(j=0; j<NCOL; j++){
			strcpy(grid[i][j].site_id, gridid);
			
			fscanf(fp_s1dat,"%lf %lf", &grid[i][j].lat, &grid[i][j].lon);
			
			if(grid[i][j].lat>(89.5-(double)WROW+0.5) || grid[i][j].lat<(89.5-(double)WROW-0.5)){
				printf("!! BAD latitude in s1 input data: %s - %lf!\n", gridid, grid[i][j].lat);
				exit (1);
			}
			if(grid[i][j].lon>((double)WCOL+1.0) || grid[i][j].lon<((double)WCOL)){
				printf("!! BAD longitude in s1 input data: %s: %lf!\n", gridid, grid[i][j].lon);
				exit (1);
			}
			
			/* for debug */
			/* printf("%lf %lf\n", grid[i][j].lat, grid[i][j].lon); */

			/* altitude by WorldClim */
			fscanf(fp_s1dat,"%lf", &grid[i][j].topo);
			grid[i][j].topo = (grid[i][j].topo>0.0)?grid[i][j].topo:0.0;
			
			/* new landcover data */
			fscanf(fp_s1dat,"%ld", &grid[i][j].veg_nlc);
			/* 0: Water (and Goode's interrupted space) */
			/* 1: Evergreen Needleleaf Forest */
			/* 2: Evergreen Broadleaf Forest */
			/* 3: Deciduous Needleleaf Forest */
			/* 4: Deciduous Broadleaf Forest */
			/* 5: Mixed Forest */
			/* 6: Woodland */
			/* 7: Wooded Grassland */
			/* 8: Closed Shrubland */
			/* 9: Open Shrubland */
			/* 10: Grassland */
			/* 11: Cropland */
			/* 12: Bare Ground */
			/* 13: Urban and Built-up */
			/* 14: Wetland */
			/* 15: Snow & Ice */

			/* MOD12 landcover data */
			fscanf(fp_s1dat,"%ld", &grid[i][j].veg_mod12);
			/* IGBP Class Scheme
			 0 Water
			 1 Evergreen Needleleaf Forest
			 2 Evergreen Broadleaf Forest
			 3 Deciduous Needleleaf Forest  
			 4 Deciduous Broadleaf Forest  
			 5 Mixed Forests               
			 6 Closed Shrublands          
			 7 Open Shrublands        
			 8 Woody Savannas          
			 9 Savannas               
			10 Grasslands                  
			11 Permanent Wetlands           
			12 Croplands                   
			13 Urban and Built-Up          
			14 Cropland/Natural Vegetation Mosaic
			15 Snow and Ice               
			16 Barren or Sparsely Vegetated 
			17 (IGBP Water Bodies, recoded to 0 for MODIS Land Product consistency.)
			254 Unclassified
			255 Fill Value
			*/
			fscanf(fp_s1dat,"%ld", &grid[i][j].veg_mod12a);
			/* UMD Class Scheme
			0 Water
			1 Evergreen Needleleaf Forest
			2 Evergreen Broadleaf Forest
			3 Deciduous Needleleaf Forest
			4 Deciduous Broadleaf Forest
			5 Mixed Forests
			6 Closed Shrublands
			7 Open Shrubland
			8 Woody Savannas
			9 Savannas
			10 Grasslands
			12 Croplands
			13 Urban and Built-Up
			16 Barren or Sparsely Vegetated
			254 Unclassified
			255 Fill Value
			*/
			fscanf(fp_s1dat,"%ld", &grid[i][j].veg_mod12b);
			/* LAI/fPAR Class Scheme
			0 Water
			1 Grasses/Cereal Crops
			2 Shrubs
			3 Broadleaf Crops
			4 Savannah
			5 Broadleaf Forest
			6 Needleleaf Forest
			7 Unvegetated
			8 Urban
			254 Unclassified
			255 Fill Value
			*/
			fscanf(fp_s1dat,"%ld", &grid[i][j].veg_mod12c);
			/* BGC Class Scheme
			0 Water
			1 Evergreen Needleleaf Vegetation
			2 Evergreen Broadleaf Vegetation
			3 Deciduous Needleleaf Vegetation
			4 Deciduous Broadleaf Vegetation
			5 Annual Broadleaf Vegetation
			6 Annual Grass Vegetation
			7 Non-vegetated Land
			8 Urban 
			254 Unclassified
			255 Fill Value
			*/
			fscanf(fp_s1dat,"%ld", &grid[i][j].veg_mod12d);
			/* PFT Class Scheme
			0 Water
			1 Needleleaf evergreen tree
			2 Broadleaf evergreen tree
			3 Needleleaf deciduous tree
			4 Broadleaf deciduoud tree
			5 Shrub
			6 Grass
			7 Cereal crop
			8 Broadleaf crop
			9 Urban and built-up
			10 Snow and ice
			11 Barren or sparsely vegetated
			254 Unclassified
			*/
			fscanf(fp_s1dat,"%ld", &grid[i][j].veg_mod12e);
			/* secondary IGBP type */
			
			/* GLC2000 landcover data */
			fscanf(fp_s1dat,"%lf %ld", &dist, &grid[i][j].veg_glc);
			/*
			1	Tree Cover, broadleaved, evergreen
			2	Tree Cover, broadleaved, deciduous, closed 
			3	Tree Cover, broadleaved, deciduous, open
			4	Tree Cover, needle-leaved, evergreen
			5	Tree Cover, needle-leaved, deciduous
			6	Tree Cover, mixed leaf type
			7	Tree Cover, regularly flooded, fresh water (& brackish)
			8	Tree Cover, regularly flooded, saline water
			9	Mosaic: Tree cover / Other natural vegetation 
			10	Tree Cover, burnt
			11	Shrub Cover, closed-open, evergreen
			12	Shrub Cover, closed-open, deciduous 
			13	Herbaceous Cover, closed-open 
			14	Sparse Herbaceous or sparse Shrub Cover
			15	Regularly flooded Shrub and/or Herbaceous Cover
			16	Cultivated and managed areas
			17	Mosaic: Cropland / Tree Cover / Other natural vegetation
			18	Mosaic: Cropland / Shrub or Grass Cover 
			19	Bare Areas
			20	Water Bodies (natural & artificial)
			21	Snow and Ice (natural & artificial)
			22	Artificial surfaces and associated areas
			*/
			
			grid[i][j].veg_type = grid[i][j].veg_nlc;

			/********************************************************************/			
			/* USDA (FAO/UNESCO) soil type data */
			fscanf(fp_s1dat,"%lf %ld", &dist, &grid[i][j].soil_usda);
			/* IGBP-DIS soil property data */
			fscanf(fp_s1dat,"%lf %lf %lf %lf %lf %lf %lf %lf", 
				&dist, &grid[i][j].wiltpoint, &grid[i][j].total_n, &grid[i][j].thermcap, 
				&grid[i][j].total_c, &grid[i][j].plf_water, &grid[i][j].fieldcap, 
				&grid[i][j].bulkdens);
			grid[i][j].fieldcap30 = 0.3 * grid[i][j].fieldcap;
			
			grid[i][j].hyd_cond = 0.004313;
			
			switch(grid[i][j].soil_usda){
				case 0:
					grid[i][j].soil_type = 0; /* ocean */
					break;
				case 1: case 2: case 3:
					grid[i][j].soil_type = 1; /* bare ground */
					break;
				case 5: case 6: case 7:
					grid[i][j].soil_type = 2; 
					/* Gelisols (Soil Taxonomy) */
					break;
				case 10: case 11: case 12: case 13:	
					grid[i][j].soil_type = 3; 
					/* Histosols (Soil Taxonomy) */
					break;
				case 15: case 16: case 17: case 18: case 19: 
					grid[i][j].soil_type = 4; 
					/* Spodosols */
					break;
				case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: 
					grid[i][j].soil_type = 5; 
					/* Andisols */
					break;
				case 30: case 31: case 32: case 33: case 34: 
					grid[i][j].soil_type = 6; 
					/* Oxisols */
					break;
				case 40: case 41: case 42: case 43: case 44: case 45: 
					grid[i][j].soil_type = 7; 
					/* Vertisols */
					break;
				case 50: case 51: case 52: case 53: case 54: case 55: case 56: 
					grid[i][j].soil_type = 8; 
					/* Aridisols */
					break;
				case 60: case 61: case 62: case 63: case 64: 
					grid[i][j].soil_type = 9; 
					/* Ultisols */
					break;
				case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: 
					grid[i][j].soil_type = 10; 
					/* Mollisols */
					break;
				case 80: case 81: case 82: case 83: case 84: 
					grid[i][j].soil_type = 11; 
					/* Alfisols */
					break;
				case 90: case 91: case 92: case 93: case 94: case 85: case 86: 
					grid[i][j].soil_type = 12; 
					/* Inceptisols */
					break;
				case 95: case 96: case 97: case 98: case 99:
					grid[i][j].soil_type = 13; 
					/* Entisols */
					break;
				default:
					grid[i][j].soil_type = 1; 
					break;
			}
			
			/********************************************************************/			
			/* WorldClim climatology data */
			for(f=0;f<12;f++){
				fscanf(fp_s1dat,"%lf", &prec_wclim[f]);
				grid[i][j].prec_wclim[f] = prec_wclim[f];
			}
			for(f=0;f<12;f++){
				fscanf(fp_s1dat,"%lf", &tmin_wclim[f]);
			}
			for(f=0;f<12;f++){
				fscanf(fp_s1dat,"%lf", &tmax_wclim[f]);
				
				grid[i][j].temp_wclim[f] = (tmin_wclim[f] + tmax_wclim[f])/2.0;
			}
		}
	}
	
	if(NOTICE==1){
		printf("Data availability\n");
	}
	for(i=0; i<DROW; i++){ /* 120 */
		for(j=0; j<DCOL; j++){
			if(grid[i][j].veg_type==0){
				grid[i][j].flag_datavl = -1; /* no land */
			}else if(grid[i][j].prec_wclim[0]<-900.0){
				grid[i][j].flag_datavl = -2; /* no climate data */
			}else if(grid[i][j].wiltpoint<-900.0){
				grid[i][j].flag_datavl = -3; /* no soil data */
			}else{
				grid[i][j].flag_datavl = 1;
			}
			
			if(NOTICE==1){
				printf("%2d", grid[i][j].flag_datavl);
				if((j+1)%12==0){
					printf("\n");
				}
			}
		}
		if(NOTICE==1){
			printf("\n");
		}
	}
	if(NOTICE==1){
		printf("done\n");
	}

	/**********************************/
	dist_min = 100.0;
	glat = 89.5 - (double)WROW;
	glon = (double)WCOL + 0.5;
	for(i=0;i<16;i++){
		for(j=0;j<19;j++){
			dist = (nlat[i] - glat)*(nlat[i] - glat) + (nlon[j] - glon)*(nlon[j] - glon);
			if(dist_min>dist){
				dist_min = dist;
				nceprow = i;
				ncepcol = j;
			}
		}
	}
	if(NOTICE==1){
		printf("NCEP data extraction: %ld %ld: %lf...", nceprow, ncepcol, dist_min);
	}
	
	if(( fp_ncepclim[0] = fopen("dswrf_sfc_s1ea.dat","rt")) == NULL){
		printf("No climate data: dswrf_sfc_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[1] = fopen("tcdc_eatm_s1ea.dat","rt")) == NULL){
		printf("No climate data: tcdc_eatm_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[2] = fopen("air_2m_s1ea.dat","rt")) == NULL){
		printf("No climate data: air_2m_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[3] = fopen("skt_sfc_s1ea.dat","rt")) == NULL){
		printf("No climate data: skt_sfc_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[4] = fopen("tmp_0-10cm_s1ea.dat","rt")) == NULL){
		printf("No climate data: tmp_0-10cm_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[5] = fopen("tmp_10-200cm_s1ea.dat","rt")) == NULL){
		printf("No climate data: tmp_10-200cm_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[6] = fopen("shum_2m_s1ea.dat","rt")) == NULL){
		printf("No climate data: shum_2m_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[7] = fopen("prate_sfc_s1ea.dat","rt")) == NULL){
		printf("No climate data: prate_sfc_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[8] = fopen("uwnd_10m_s1ea.dat","rt")) == NULL){
		printf("No climate data: uwnd_10m_s1ea.dat !\n");
		exit (1);
	}
	if(( fp_ncepclim[9] = fopen("vwnd_10m_s1ea.dat","rt")) == NULL){
		printf("No climate data: vwnd_10m_s1ea.dat !\n");
		exit (1);
	}
	
	/* clear */
	for(g=0;g<366;g++){
		tmp_2m_ncep_dav[g] = 0.0;
	}
	for(g=0;g<12;g++){
		prate_sfc_ncep_mav[g] = 0.0;
	}
	
	for(f=1948;f<=2005;f++){
		ndy = (f%4==0)?366:365;
		for(g=0;g<ndy;g++){
			f_doyTmody(f, g, &month, &day);
			fscanf(fp_ncepclim[0], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[1], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[2], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[3], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[4], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[5], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[6], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[7], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[8], "%ld %ld", &year, &doy);
			fscanf(fp_ncepclim[9], "%ld %ld", &year, &doy);
			
			for(i=0;i<16;i++){
				for(j=0;j<19;j++){
					fscanf(fp_ncepclim[0], "%lf", &ndat[0]);
					fscanf(fp_ncepclim[1], "%lf", &ndat[1]);
					fscanf(fp_ncepclim[2], "%lf", &ndat[2]);
					fscanf(fp_ncepclim[3], "%lf", &ndat[3]);
					fscanf(fp_ncepclim[4], "%lf", &ndat[4]);
					fscanf(fp_ncepclim[5], "%lf", &ndat[5]);
					fscanf(fp_ncepclim[6], "%lf", &ndat[6]);
					fscanf(fp_ncepclim[7], "%lf", &ndat[7]);
					fscanf(fp_ncepclim[8], "%lf", &ndat[8]);
					fscanf(fp_ncepclim[9], "%lf", &ndat[9]);
										
					if(i==nceprow && j==ncepcol){
						dswrf_sfc_d[f-BYR][g] = ndat[0];
						tcdc_clm_d[f-BYR][g] = ndat[1];						
						tmp_2m_d[f-BYR][g] = ndat[2];
						tmp_sfc_d[f-BYR][g] = ndat[3];
						tmp10_soil_d[f-BYR][g] = ndat[4];
						tmp200_soil_d[f-BYR][g] = ndat[5];
						spfh_2m_d[f-BYR][g] = ndat[6];
						prate_sfc_d[f-BYR][g] = ndat[7];
						wind_10m_d[f-BYR][g] = sqrt(ndat[8]*ndat[8] + ndat[9]*ndat[9]);
						
						/* when no soil temperature data */
						if(tmp10_soil_d[f-BYR][g]<170.0 || tmp10_soil_d[f-BYR][g]>330.0){
							tmp10_soil_d[f-BYR][g] = 0.796342*tmp_2m_d[f-BYR][g] + 55.345676;
						}
						if(tmp200_soil_d[f-BYR][g]<170.0 || tmp200_soil_d[f-BYR][g]>330.0){
							tmp200_soil_d[f-BYR][g] = 0.606448*tmp_2m_d[f-BYR][g] + 107.496955;
						}
						
						/* VPD */
						tair = tmp_2m_d[f-BYR][g] - ZAT;
						aprsr = 1013.25*exp(-1.0*(28.964*0.001)*9.8*500.0/(UGC*(tair+ZAT)));	/* hPa */
						if(tair>=0.0){
							vsat = 6.1078*pow(10.0, 7.5*tair/(237.3 + tair));
						}else{
							vsat = 6.1078*pow(10.0, 9.5*tair/(265.3 + tair));
						}
						vpres = aprsr * spfh_2m_d[f-BYR][g]/(0.622 + 0.378*spfh_2m_d[f-BYR][g]);
						
						vpd_d[f-BYR][g] = vsat - vpres;

						/* NCEP climatology */
						if(f>=1950 && f<=2000 && g<365){
							tmp_2m_ncep_dav[g] += (ndat[2]-ZAT)/51.0;
							prate_sfc_ncep_mav[month] += (ndat[7]*24.0*3600.0)/51.0;
						}
						
						/* for debug */
						/* printf("%ld %ld %lf %lf\n",f,g,ndat[0],ndat[4]); */
					}
				}
			}
		}
	}
	
	if(NOTICE==1){
		for(g=0;g<12;g++){
			printf("%ld %lf\n", g, prate_sfc_ncep_mav[g]);
		}
		for(g=0;g<365;g++){
			printf("%ld %lf\n", g, tmp_2m_ncep_dav[g]);
		}
		
		printf("done\n");
	}
	
	fclose(fp_s1dat);
	fclose(fp_ncepclim[0]);
	fclose(fp_ncepclim[1]);
	fclose(fp_ncepclim[2]);
	fclose(fp_ncepclim[3]);
	fclose(fp_ncepclim[4]);
	fclose(fp_ncepclim[5]);
	fclose(fp_ncepclim[6]);
	fclose(fp_ncepclim[7]);
	fclose(fp_ncepclim[8]);
	fclose(fp_ncepclim[9]);	
}

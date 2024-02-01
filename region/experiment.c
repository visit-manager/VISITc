/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: GHG schemes addition by Motoko Inatomi
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2007 / 05 / 30 by Akihiko Ito: Version 0.2.0.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

#ifdef _OPENMP
#include<omp.h>
#endif

extern long	month_day[12];
extern long	WMODE;
extern long	WGRIDS;
extern long NBIOME;
extern float fdat[NROW*NCOL];
extern int idat[NROW*NCOL];
extern struct Echar echar_type[MAX_BIOME];

double mean_veg[20][32];

/* SIMULATION ***********************************************************/
/* area */
void f_experiment(
  struct Grid grid[],
  struct Loct loct[],
  struct Echar echar[],
  struct Mass mass[],
  struct Flux flux[]
){
	long e, f, h, i, j;
	char num[8], filename[128];
	long ndy, end_year, calc_flag, pstart, pend;
	float gpp_a, npp_a, nep_a, lai_a, plant_a, soil_a, xx1_a, xx2_a, xx3_a, ch4_a, nn;
    float gpp_ga, npp_ga, nep_ga, plant_ga, soil_ga, prec_ga, rdata, vps;
    struct Grid grid0;
    struct Loct loct0;
    struct Echar echar0;
    struct Mass mass0;
    struct Flux flux0;
    float* out_ann;
    float* out_hr;
	/* FILE *fp_o; */
	FILE *fp_clim[N_CLIMD];
	FILE *fp_outann, *fp_outhr[N_OUTHR], *fp_log;
    FILE *fp_ss_grid, *fp_ss_loct, *fp_ss_mass, *fp_ss_flux, *fp_ss_echar;
    FILE *fp_monitor, *fp_veg;
    FILE *fp_lai;
    
    if(P_MONI == 0){
        pstart = 0;
        pend = WGRIDS -1;
    }else{
        pstart = P_MONI;
        pend = P_MONI;
    }
	
    /* experimental phase */
    for(i=pstart; i<=pend; i++){
        loct[i].phase = 2;
    }
    
	if(NOTICE == 1){
		printf("Start experimental phase\n");
	}
	
    for(i=0;i<20;i++){
        for(j=0;j<32;j++){
            mean_veg[i][j] = 0.0;
        }
    }
    
    out_ann = (float *)malloc(sizeof(float) * N_OUTANN*(pend-pstart+1));
    out_hr = (float *)malloc(sizeof(float) * N_OUTHR*(pend-pstart+1));

	/* **************************************************************************/
    if(USE_RESTART == 1){
        if(NOTICE == 1){
            printf("Read restart file...");
        }

        strcpy(filename, grid[0].area_id);
        strcat(filename, "_restart_grid.flt");
        if((fp_ss_grid = fopen(filename,"rb")) == NULL){
            printf("!! NO %s\n",filename);
            exit (1);
        }

        strcpy(filename, grid[0].area_id);
        strcat(filename, "_restart_loct.flt");
        if((fp_ss_loct = fopen(filename,"rb")) == NULL){
            printf("!! NO %s\n",filename);
            exit (1);
        }

        strcpy(filename, grid[0].area_id);
        strcat(filename, "_restart_mass.flt");
        if((fp_ss_mass = fopen(filename,"rb")) == NULL){
            printf("!! NO %s\n",filename);
            exit (1);
        }

        strcpy(filename, grid[0].area_id);
        strcat(filename, "_restart_flux.flt");
        if((fp_ss_flux = fopen(filename,"rb")) == NULL){
            printf("!! NO %s\n",filename);
            exit (1);
        }

        strcpy(filename, grid[0].area_id);
        strcat(filename, "_restart_echar.flt");
        if((fp_ss_echar = fopen(filename,"rb")) == NULL){
            printf("!! NO %s\n",filename);
            exit (1);
        }

        for(f=0;f<NCOL*NROW;f++){
            fread(&grid[f], sizeof(grid0), 1, fp_ss_grid);
            fwrite(&loct[f], sizeof(loct0), 1, fp_ss_loct);
            fread(&mass[f], sizeof(mass0), 1, fp_ss_mass);
            fread(&flux[f], sizeof(flux0), 1, fp_ss_flux);
            fread(&echar[f], sizeof(echar0), 1, fp_ss_echar);
        }
        
        fclose(fp_ss_grid);
        fclose(fp_ss_loct);
        fclose(fp_ss_mass);
        fclose(fp_ss_flux);
        fclose(fp_ss_echar);
        
        if(NOTICE == 1){
            printf("done\n");
        }
    }
    
	/* annual mean output */
	strcpy(filename, grid[0].area_id);
    strcat(filename, "_ann_");
	strcat(filename, grid[0].file_name);
	strcat(filename, ".flt");
    fp_outann = fopen(filename,"wb");
    
     /* log file */
    fp_log = fopen("log_ext.txt","wt");
    fp_monitor = fopen("monitor.txt","wt");
    fp_veg = fopen("veg_mean.txt","wt");
    
    end_year = EYR;
    
    /* number of pararelization threds with OpenMP */
	#ifdef _OPENMP
	omp_set_num_threads(48); /* cores */
	#endif
    
	/* roop for experimental stage ************************************************/
	for(e=BYR; e<=end_year; e++){ /**/
        for(i=pstart; i<=pend; i++){
            loct[i].climy = e;
            loct[i].adyear = e;

            loct[i].CO2y = e;
            if(FIX_CO2 == 1){
                loct[i].CO2y = BYR;
            }
        }
		
		/* number of seasonal roop */
        if(LEAP_YR ==1){
            ndy = (e%4==0)?366:365; /* */
        }else{
            ndy = 365; /* no leap year */
        }
		
		/* string of year ADXXXX */
		snprintf(num, 8, "%04d", (short)e);

        f_create_dfile_global(loct[0].phase, loct[0].adyear, &grid[0], filename, fp_outhr);
        
        /* fix canopy N exp: 2014/08/17 by A.Ito *************************/
        if((FIX_LAI == 2 || FIX_LAI == 3) && CONST_KN==0){
            strcpy(filename,"./data/CFSR_V04_");
            snprintf(num,8,"%4ld",e);
            strcat(filename, num);
            /* strcat(filename, "_140806.flt"); */
            strcat(filename, "_170124.flt");
            
            fp_lai = fopen(filename,"rb");
        }
        if((FIX_LAI == 2 || FIX_LAI == 3) && CONST_KN==1){
            strcpy(filename,"./data/CFSR_V04_");
            snprintf(num,8,"%4ld",e);
            strcat(filename, num);
            strcat(filename, "_170220.flt");
        
            fp_lai = fopen(filename,"rb");
        }
        if((FIX_LAI == 2 || FIX_LAI == 3) && CONST_KN==3){
            strcpy(filename,"./data/CFSR_V04_");
            snprintf(num,8,"%4ld",e);
            strcat(filename, num);
            //strcat(filename, "_161227.flt");
            strcat(filename, "_170321.flt");
            
            fp_lai = fopen(filename,"rb");
        }
        
        /* open meteorological data */
        if(strcmp(grid[0].area_id, "BAMIYAN")==0){
            f_open_bamiyan_clim(loct[0].phase, (short)e, fp_clim);
        }

		/* ************************************************/
        gpp_ga = npp_ga = nep_ga = plant_ga = soil_ga = prec_ga = 0.0;
        
        for(i=0; i<(N_OUTANN*(pend-pstart+1)); i++){
            out_ann[i] = 0.0;
        }
                 
        for(h=0; h<ndy; h++){
            doyTmody(e, h, &(loct[0].month), &(loct[0].mday));
            for(i=pstart; i<=pend; i++){
                loct[i].doy = h;
                loct[i].hour = 0;
                loct[i].month = loct[0].month;
                loct[i].mday = loct[0].mday;
            }
            /* doyTmody(1979, h, &(loct->month), &(loct->mday)); */  /* no leap year */
            
            /* open meteorological data */
            if(loct[0].mday == 0 && strcmp(grid[0].area_id, "GLOBAL")==0){
                f_open_global_clim(loct[0].phase, (short)e, (short)(loct[0].month)+1, fp_clim);
            }
            if(loct[0].mday == 0 && strcmp(grid[0].area_id, "EASIA")==0){
                f_open_easia_clim(loct[0].phase, (short)e, (short)(loct[0].month)+1, fp_clim);
            }
            if(loct[0].mday == 0 && strcmp(grid[0].area_id, "PAWCs") == 0){
                f_open_pawcs_clim(loct[0].phase, (short)e, (short)(loct[0].month)+1, (short)(loct[0].mday)+1, (short)(loct[0].hour)+1, fp_clim);
            }

            for(f=0; f<DSTEP; f++){
                for(i=pstart; i<=pend; i++){
                    loct[i].hour = f;
                    for(j=0;j<N_OUTHR;j++){
                        out_hr[j*(pend-pstart+1) + i] = 0.0;
                    }
                }

                if(strcmp(grid[0].area_id, "JAPAN") == 0){
                    //f_open_japan_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_japan_clim(loct[0].phase, (short)e, (short)(loct[0].month)+1, (short)(loct[0].mday)+1, (short)(loct[0].hour)+1, fp_clim);
                }
                if(strcmp(grid[0].area_id, "BB") == 0){
                    //f_open_bb_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_bb_clim(loct[0].phase, (short)e, (short)(loct[0].month)+1, (short)(loct[0].mday)+1, (short)(loct[0].hour)+1, fp_clim);
                }
                if(strcmp(grid[0].area_id, "JAPANc") == 0){
                    //f_open_japanc_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_japanc_clim(loct[0].phase, (short)e, (short)(loct[0].month)+1, (short)(loct[0].mday)+1, (short)(loct[0].hour)+1, fp_clim);
                }
                if(strcmp(grid[0].area_id, "JAPANh") == 0){
                    //f_open_japanh_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_japanh_clim(loct[0].phase, (short)e, (short)(loct[0].month)+1, (short)(loct[0].mday)+1, (short)(loct[0].hour)+1, fp_clim);
                }
                if(strcmp(grid[0].area_id, "JAPANk") == 0){
                    //f_open_japank_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_japank_clim(loct[0].phase, (short)e, (short)(loct[0].month)+1, (short)(loct[0].mday)+1, (short)(loct[0].hour)+1, fp_clim);
                }

                printf("%4ld %3ld %2ld: ", e, h, f);
                fprintf(fp_log,"%ld %ld %ld ", e, h, f);

                /* ********************************************************/
                if(strcmp(grid[0].area_id, "GLOBAL")==0){
                    /* read regional climate data ****/
                    
                    /* precipitation */
                    fread(fdat, 4, WGRIDS, fp_clim[0]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        grid[i].prec_region = fdat[i] * 3600.0; /* precipitation, mm/h */
                    }
                    
                    /* shortwave radiation */
                    fread(fdat, 4, WGRIDS, fp_clim[1]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        grid[i].srad_region = fdat[i] / 3600.0;      /* ssrd, W/m2 */
                    }

                    /* longwave radiation */
                    fread(fdat, 4, WGRIDS, fp_clim[2]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        grid[i].trad_region = fdat[i] / 3600.0;      /* strd, W/m2 */
                    }

                    /* temperature */
                    fread(fdat, 4, WGRIDS, fp_clim[3]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].tmax_region = fdat[i] - ZAT; /* temp max, K */
                        grid[i].tmin_region = fdat[i] - ZAT; /* temp min, K */
                    }

                    /* humidity */
                    fread(fdat, 4, WGRIDS, fp_clim[4]);
                    for(i=pstart; i<=pend; i++){
                        fdat[i] -= ZAT;
                        /* dew point, K => vapor pressure */
                        if(fdat[i] > 0.0){ /* at water surface */
                            vps = 6.1078 * pow(10.0,(7.5*fdat[i])/(237.3 + fdat[i]));
                        }else{ /* at ice surface */
                            vps = 6.1078 * pow(10.0,(9.5*fdat[i])/(265.3 + fdat[i]));
                        }
                        vps = (vps>=0.0)?vps:0.0;

                        grid[i].humd_region = vps;
                    }

                    /* wind */
                   fread(fdat, 4, WGRIDS, fp_clim[5]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].wind_region = fdat[i];  /* u-wind at 10m, m/s */
                    }
                    fread(fdat, 4, WGRIDS, fp_clim[6]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].wind_region = sqrt(fdat[i]*fdat[i] + grid[i].wind_region*grid[i].wind_region);  /* v-wind at 10m, m/s */
                    }

                    if(FIX_LAI == 2 || FIX_LAI == 3){
                        fread(fdat, 4, WGRIDS, fp_lai);
                        for(i=pstart; i<=pend; i++){
                            grid[i].p_lai = fdat[i]; /* LAI by fix-co2 experiment */
                        }
                    }
                }else if(strcmp(grid[0].area_id, "BAMIYAN")==0){
                    /* BAMIYAN: 2017/08/01 by A.Ito */
                    /* solar radiation */
                    fread(idat, 4, WGRIDS, fp_clim[1]);
                    for(i=pstart; i<=pend; i++){ 
                        if(idat[i] < 0){
                            idat[i] = 0;
                        }
                        grid[i].srad_region = (float)idat[i];      /* dsw rad, W/m2 */
                    }
                    
                    /* temperature */
                    if(f==0){
                        fread(fdat, 4, WGRIDS, fp_clim[2]);
                        for(i=pstart; i<=pend; i++){ 
                            grid[i].tmax_region = fdat[i] - ZAT; /* temp max, K */
                            grid[i].tmin_region = fdat[i] - ZAT; /* temp min, K */
                        }
                    }

                    /* precipitation */
                    if(f==0){
                        /* fread(fdat, 4, WGRIDS, fp_clim[0]); */
                        fscanf(fp_clim[0],"%f", &rdata);
                        
                        for(i=pstart; i<=pend; i++){
                            if(rdata < 0.0){
                                rdata = 0.0;
                            }
                            grid[i].prec_region = rdata * 3600.0/(float)DSTEP; /* precipitation, mm/h */
                        }
                    }
                    
                    /* humidity */
                    if(f==0){
                        fread(fdat, 4, WGRIDS, fp_clim[3]);
                        for(i=pstart; i<=pend; i++){ 
                            grid[i].humd_region = fdat[i];	/* specific humidity, g/g */
                        }
                    }
                    
                    /* wind */
                    if(f==0){
                        /* fread(fdat, 4, WGRIDS, fp_clim[4]); */
                        fscanf(fp_clim[4],"%f", &rdata);
                        
                        for(i=pstart; i<=pend; i++){
                            if(rdata < 0.0){
                                rdata = 0.0;
                            }
                            grid[i].wind_region = rdata; /* wind at 10m, m/s */
                        }
                    }
                }else if(strcmp(grid[0].area_id, "EASIA") == 0){
                    fread(fdat, 4, WGRIDS, fp_clim[0]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        grid[i].prec_region = fdat[i]; /* precipitation, mm/h */
                    }
                    fread(fdat, 4, WGRIDS, fp_clim[1]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        grid[i].srad_region = fdat[i];      /* dsw rad, W/m2 */
                    }
                    fread(fdat, 4, WGRIDS, fp_clim[2]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].tmax_region = fdat[i] - ZAT; /* temp max, K */
                        grid[i].tmin_region = fdat[i] - ZAT; /* temp min, K */
                    }
                    fread(fdat, 4, WGRIDS, fp_clim[3]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].humd_region = fdat[i];  /* dew-point temperature, K */
                    }
                }else if(strcmp(grid[0].area_id, "JAPAN")==0 || strcmp(grid[0].area_id, "BB")==0
                        || strcmp(grid[0].area_id, "JAPANc")==0|| strcmp(grid[0].area_id, "JAPANh")==0
                        || strcmp(grid[0].area_id, "JAPANk")==0|| strcmp(grid[0].area_id, "PAWCs")==0){
                    /* read regional climate data ****/
                    /* precipitation */
                    fread(fdat, 4, NROW*NCOL, fp_clim[0]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        /* grid[i].prec_region = 3600.0*fdat[i]; */ /* precipitation, mm/h */
                        grid[i].prec_region = fdat[i] * 3600.0; /* precipitation, mm/h */
                    }
                    
                    /* solar radiation */
                    fread(fdat, 4, NROW*NCOL, fp_clim[1]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        grid[i].srad_region = fdat[i] / 3600.0;      /* dsw rad, W/m2 */
                    }
                    
                    /* temperature */
                    fread(fdat, 4, NROW*NCOL, fp_clim[2]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].tmax_region = fdat[i] - ZAT; /* temp max, K */
                        grid[i].tmin_region = fdat[i] - ZAT; /* temp min, K */
                    }
                    
                    /* specific humidity */
                    fread(fdat, 4, NROW*NCOL, fp_clim[3]);
                    for(i=pstart; i<=pend; i++){
                        fdat[i] -= ZAT;
                        /* dew point, K => vapor pressure */
                        if(fdat[i] > 0.0){ /* at water surface */
                            vps = 6.1078 * pow(10.0,(7.5*fdat[i])/(237.3 + fdat[i]));
                        }else{ /* at ice surface */
                            vps = 6.1078 * pow(10.0,(9.5*fdat[i])/(265.3 + fdat[i]));
                        }
                        vps = (vps>=0.0)?vps:0.0;

                        grid[i].humd_region = vps;
                    }
                    
                    /* wind */
                    fread(fdat, 4, NROW*NCOL, fp_clim[4]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].wind_region = fdat[i]; /* wind at 10m, m/s */
                    }
                }
                
                /* initial assumption ******/
                /* grid[i].wind_region = 1.5;	*/
                
                /* ****************************************************/
                /* OpenMP */
                #pragma omp parallel for schedule(dynamic, 4) private(calc_flag)

                for(i=pstart; i<=pend; i++){ 
                    calc_flag = 1;
                    if(grid[i].flag_datavl!=1){
                        calc_flag = 0;
                    }
                    if(grid[i].prec_region<-0.1 || grid[i].prec_region>1500.0){
                        calc_flag = 0;
                    }
                    if(grid[i].srad_region<-1.0 || grid[i].srad_region>2000.0){
                        calc_flag = 0;
                    }
                    if(grid[i].tmax_region<-90.0 || grid[i].tmax_region>90.0){
                        calc_flag = 0;
                    }
                    if(grid[i].tmin_region<-90.0 || grid[i].tmin_region>90.0){
                        calc_flag = 0;
                    }
                    if(grid[i].humd_region<-1.0 || grid[i].humd_region>100.0){
                        calc_flag = 0;
                    }
                    if(grid[i].wind_region<-120.0 || grid[i].wind_region>120.0){
                        calc_flag = 0;
                    }
                    
                    grid[i].calc_flag = calc_flag;
                    
                    /* data for calculation available **********************/
                    if(calc_flag==1  && (i+0)%RUNFAST==0){
                    /* if(calc_flag==1  && (i+0)%RUNFAST==0 && grid[i].veg_type ==2){ */
                        grid[i].time++;
                        grid[i].age_stand += 1.0;
                        loct[i].hour = f;
                        loct[i].CO2y = e;
                        
                        if(FIX_CO2 == 1){
                            loct[i].CO2y = BYR;
                        }

                        /* land-use */
                        /* flux2[i].lue_conv = 0.0;
                        flux2[i].lue_detr_1 = 0.0;
                        flux2[i].lue_detr_10 = 0.0;
                        flux2[i].lue_detr_100 = 0.0;
                        f_conv_landuse(&(grid[i]),&(loct2[i]));
                        f_emit_landuse(&(grid[i]), &(loct2[i]), &(mass[i]), &(flux2[i])); */
                        
                        /* local condition ************************************/
                        f_init_cond(&(grid[i]), &(loct[i]), &(echar[i]), &(mass[i]), &(flux[i]));
                        
                        /* disturbance *****/
                        /* loct->climy = e; */
                        /*   disturbance_regime(e, &(grid[i]), &(loct2[i]), &(mass[i]), &(flux2[i])); */
                        
                        /* BAMIYAN grazing scenarios: 2017/09/03 by A.Ito */
                        if(strcmp(grid[0].area_id, "BAMIYAN")==0){
                            grazing_event(&(grid[i]), &(loct[i]), &(mass[i]));
                        }

                        /* basic scheme ****************************************/
                        f_daily_scheme(&(grid[i]), &(loct[i]), &(echar[i]), &(mass[i]), &(flux[i]));
                    }
                    
                    /* erosion */
                    /* f_erosion_rusle(&grid[i], loct, &mass[i], flux);
                    (mass[i].soil).msl_a -= (flux->soil).erosion_carbon; */
                    
                    /* *********************************************/
                    /* monitoring */
                    if(WMODE == 2){
                        /* printf("%4d ", grid[i].calc_flag); */
                    }
                }
                /* end of daily roop ***************************/
                #pragma omp barrier    
               
                /* monitoring on console ******************/
                gpp_a = npp_a = nep_a = lai_a = plant_a = soil_a = 0.0;
                xx1_a = xx2_a = xx3_a = ch4_a = nn = 0.0;
                for(i=pstart; i<=pend; i++){ 
                    if(grid[i].calc_flag == 1){
                        nn += 1.0;
                        gpp_a += (flux[i].tree).gpp + loct[i].funder_c3 * (flux[i].c3).gpp +
                                    loct[i].funder_c4 * (flux[i].c4).gpp;
                        npp_a += (flux[i].tree).npp + loct[i].funder_c3 * (flux[i].c3).npp +
                                   loct[i].funder_c4 * (flux[i].c4).npp;

                        /* //npp_a += (flux2[i].tree).rpg + loct2[i].funder_c3 * (flux2[i].c3).rpg +
                        //            loct2[i].funder_c4 * (flux2[i].c4).rpg; */

                        nep_a += flux[i].nep;
                        lai_a += ((mass[i].tree).lai + loct[i].funder_c3 * (mass[i].c3).lai +
                                  loct[i].funder_c4 * (mass[i].c4).lai);
                        plant_a += ((mass[i].tree).plant + loct[i].funder_c3 * (mass[i].c3).plant +
                                    loct[i].funder_c4 * (mass[i].c4).plant); /* */
                        soil_a += (mass[i].soil).soil;
                        
                        xx1_a += loct[i].tmp_2m;
                        xx2_a += loct[i].prate_sfc;
                        xx3_a += loct[i].ppfd_h;
                        
                        ch4_a += (flux[i].soil).ch4_wh;
                        
                        gpp_ga += grid[i].area * ((flux[i].tree).gpp + loct[i].funder_c3 * (flux[i].c3).gpp
                            + loct[i].funder_c4 * (flux[i].c4).gpp) / 1000000.0;
                        npp_ga += grid[i].area * ((flux[i].tree).npp + loct[i].funder_c3 * (flux[i].c3).npp
                            + loct[i].funder_c4 * (flux[i].c4).npp) / 1000000.0;
                        nep_ga += grid[i].area * flux[i].nep / 1000000.0;
                        plant_ga += grid[i].area * ((mass[i].tree).plant + loct[i].funder_c3 * (mass[i].c3).plant +
                                    loct[i].funder_c4 * (mass[i].c4).plant) / (float)(DSTEP*ndy) / 1000000.0;
                        soil_ga += grid[i].area * (mass[i].soil).soil / (float)(DSTEP*ndy) / 1000000.0;
                        prec_ga += grid[i].area * loct[i].prate_sfc / 100.0;

                        /* summation for output *********/
                        /* hourly outputs */
                        
                        /* CFSR */
                        /* //outdat01[i] = flux2[i].gpp;
                        //outdat02[i] = flux2[i].er;
                        //outdat03[i] = flux2[i].nep;
                        //outdat04[i] = flux2[i].sr;
                        //outdat04[i] = loct2[i].lai; */

                        /* BAMIYAN */
                        /* GLOBAL */
                        out_hr[0*(pend-pstart+1) + i] = flux[i].gpp;
                        out_hr[1*(pend-pstart+1) + i] = flux[i].npp;
                        out_hr[2*(pend-pstart+1) + i] = flux[i].nep;
                        out_hr[3*(pend-pstart+1) + i] = (flux[i].soil).ch4_wh;
                        
                        out_hr[4*(pend-pstart+1) + i] = grid[i].srad_region;
                        out_hr[5*(pend-pstart+1) + i] = loct[i].sfcrad_h;
                        out_hr[6*(pend-pstart+1) + i] = loct[i].ppfd_h;

                        
                        
                        
                        
                        //outdat03[i] = ((mass[i].tree).plant + loct2[i].funder_c3 * (mass[i].c3).plant + loct2[i].funder_c4 * (mass[i].c4).plant);
                        //outdat04[i] = loct2[i].lai;
                        
                        /* //outdat01[i] = (mass[i].tree).gdd;
                        //outdat02[i] = (mass[i].tree).lai;
                        //outdat03[i] = (float)(mass[i].tree).season;
                        //outdat04[i] = loct2[i].tair_dayav; */
                        
                        /* annual mean/sum outputs ***************/
                        out_ann[0*(pend-pstart+1) + i] += flux[i].gpp;
                        out_ann[1*(pend-pstart+1) + i] += flux[i].npp;
                        out_ann[2*(pend-pstart+1) + i] += flux[i].nep;
                        out_ann[3*(pend-pstart+1) + i] += (flux[i].soil).ch4_wh;

                        out_ann[4*(pend-pstart+1) + i] += loct[i].lai/(double)ndy/(double)DSTEP;
                        out_ann[5*(pend-pstart+1) + i] += ((mass[i].tree).plant + loct[i].funder_c3 * (mass[i].c3).plant
                                        + loct[i].funder_c4 * (mass[i].c4).plant)/(double)ndy/(double)DSTEP;
                        out_ann[6*(pend-pstart+1) + i] += (mass[i].soil).soil/(double)ndy;
                        
                        out_ann[7*(pend-pstart+1) + i] += loct[i].tmp_2m/(double)ndy/(double)DSTEP;
                        out_ann[8*(pend-pstart+1) + i] += loct[i].ppfd_h/(double)ndy/(double)DSTEP;
                        out_ann[9*(pend-pstart+1) + i] += loct[i].vpd/(double)ndy/(double)DSTEP;
                        out_ann[10*(pend-pstart+1) + i] += loct[i].prate_sfc;
                        out_ann[11*(pend-pstart+1) + i] += loct[i].aet;
                    }

                    mean_veg[grid[i].veg_type][0] += (double)(grid[i].area / 24.0 / (double)ndy);
                    mean_veg[grid[i].veg_type][1] += (double)(grid[i].area * flux[i].gpp);
                    mean_veg[grid[i].veg_type][2] += (double)(grid[i].area * flux[i].npp);
                    mean_veg[grid[i].veg_type][3] += (double)(grid[i].area * flux[i].nep);
                    mean_veg[grid[i].veg_type][4] += (double)(grid[i].area * loct[i].lai / 24.0 / (double)ndy);
                    mean_veg[grid[i].veg_type][5] += (double)(grid[i].area * ((mass[i].tree).plant + loct[i].funder_c3 * (mass[i].c3).plant
                        + loct[i].funder_c4 * (mass[i].c4).plant) / 24.0 / (double)ndy);
                    mean_veg[grid[i].veg_type][6] += (double)(grid[i].area * (mass[i].soil).soil / 24.0 / (double)ndy);
                }
                printf("%8.4f %8.4f %8.4f : %12.2f %12.2f %12.2f : %12.2f %12.2f %12.2f: %12.3f",
                       100.0*gpp_a/nn, 100.0*npp_a/nn, 100.0*nep_a/nn, 
                       lai_a/nn, plant_a/nn, soil_a/nn, xx1_a/nn, xx2_a/nn, xx3_a/nn, ch4_a/nn);
                printf("\n");

                fprintf(fp_log,"%.0f %f %f %f %f %f %f %f %f %f ",
                       nn, gpp_a, npp_a, nep_a, lai_a, plant_a, soil_a, xx1_a, xx2_a, xx3_a);
                fprintf(fp_log,"%12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f\n",
                       100.0*gpp_a/nn, 100.0*npp_a/nn, 100.0*nep_a/nn, 
                       lai_a/nn, plant_a/nn, soil_a/nn, ch4_a/nn);
                
                /* *****************/
                fprintf(fp_monitor,"%ld %ld %ld ",loct[P_MONI].adyear,loct[P_MONI].doy,loct[P_MONI].hour);
                
                /* fprintf(fp_monitor,"%f %f %f ",loct2[P_MONI].tmp_2m,loct2[P_MONI].dswrf_sfc,loct2[P_MONI].prate_sfc);
                fprintf(fp_monitor,"%ld %f ",(mass[P_MONI].tree).season,(mass[P_MONI].tree).gdd);
                fprintf(fp_monitor,"%f %f %f ",flux2[P_MONI].gpp,flux2[P_MONI].npp,flux2[P_MONI].nep);
                fprintf(fp_monitor,"%f %f %f ",(mass[P_MONI].tree).fol,(mass[P_MONI].tree).stm,(mass[P_MONI].tree).rot); */
                
                fprintf(fp_monitor,"%f %f ",mass[P_MONI].sw30, mass[P_MONI].sww);
                fprintf(fp_monitor,"%f %f ",loct[P_MONI].dswrf_sfc, loct[P_MONI].tcdc_clm);
                fprintf(fp_monitor,"%f %f ",loct[P_MONI].rn_short_eco, loct[P_MONI].rn_long_eco);
                fprintf(fp_monitor,"%f %f %f %f ",loct[P_MONI].prate_sfc, loct[P_MONI].pet, loct[P_MONI].aet, loct[P_MONI].ro2);
                fprintf(fp_monitor,"\n");

                if(e>=BYR && e<=EYR){
                    for(j=0;j<N_OUTHR;j++){
                        for(i=pstart; i<=pend; i++){
                            fdat[i] = out_hr[j*(pend-pstart+1) + i];
                        }
                        fwrite(fdat, sizeof(float), (pend-pstart+1), fp_outhr[j]);
                    }
                }
                
                if(strcmp(grid[0].area_id, "JAPAN")==0 || strcmp(grid[0].area_id, "BB")==0
                    || strcmp(grid[0].area_id, "JAPANc")==0 || strcmp(grid[0].area_id, "JAPANh")==0
                    || strcmp(grid[0].area_id, "JAPANk")==0){
                    fclose(fp_clim[0]);
                    fclose(fp_clim[1]);
                    fclose(fp_clim[2]);
                    fclose(fp_clim[3]);
                    fclose(fp_clim[4]);
                }
            }
            
            /* close climate files */
            j = 0;
            switch(loct[0].month){
                case 0: if(loct[0].mday==30){ j = 1; } break;
                case 1: if((e%4 == 0 && loct[0].mday==28) || (e%4 != 0 && loct[0].mday==27)){ j = 1; } break;
                //case 1: if(loct->mday==27){ j = 1; } break;
                case 2: if(loct[0].mday==30){ j = 1; } break;
                case 3: if(loct[0].mday==29){ j = 1; } break;
                case 4: if(loct[0].mday==30){ j = 1; } break;
                case 5: if(loct[0].mday==29){ j = 1; } break;
                case 6: if(loct[0].mday==30){ j = 1; } break;
                case 7: if(loct[0].mday==30){ j = 1; } break;
                case 8: if(loct[0].mday==29){ j = 1; } break;
                case 9: if(loct[0].mday==30){ j = 1; } break;
                case 10: if(loct[0].mday==29){ j = 1; } break;
                case 11: if(loct[0].mday==30){ j = 1; } break;
            }
            
            if(j == 1 && (strcmp(grid[0].area_id, "EASIA")==0 || strcmp(grid[0].area_id, "PAWCs")==0) ){
                fclose(fp_clim[0]);
                fclose(fp_clim[1]);
                fclose(fp_clim[2]);
                fclose(fp_clim[3]);
                fclose(fp_clim[4]);
            }

            if(j == 1 && strcmp(grid[0].area_id, "GLOBAL")==0 ){
                fclose(fp_clim[0]);
                fclose(fp_clim[1]);
                fclose(fp_clim[2]);
                fclose(fp_clim[3]);
                fclose(fp_clim[4]);
                fclose(fp_clim[5]);
                fclose(fp_clim[6]);
            }

            /* if(loct->climy == 2015 && loct->month==3 && loct->mday==29){
                break;
            } */
        }
		
        if(FIX_LAI == 2 || FIX_LAI == 3){
            fclose(fp_lai);
        }
        
		/* output annual maps */
		if(e>=BYR && e<=EYR){
			fwrite(out_ann, sizeof(float), N_OUTANN*(pend-pstart+1), fp_outann);
		}
        
        for(i=0;i<N_OUTHR;i++){
            fclose(fp_outhr[i]);
        }

        printf("ANNUAL %ld %f %f %f %f %f %f\n",loct[0].adyear, gpp_ga, npp_ga, nep_ga, plant_ga, soil_ga, prec_ga);

        /* if(loct->climy == 2015 && loct->month==3 && loct->mday==29){
            break;
        } */
        
        if(strcmp(grid[0].area_id, "BAMIYAN")==0){
            fclose(fp_clim[0]);
            fclose(fp_clim[1]);
            fclose(fp_clim[2]);
            fclose(fp_clim[3]);
            fclose(fp_clim[4]);
        }
	}

    for(i=0;i<20;i++){
        fprintf(fp_veg,"%ld ",i);
        for(j=0;j<32;j++){
            fprintf(fp_veg,"%lf ",mean_veg[i][j]);
        }
        fprintf(fp_veg,"\n");
    }
    
    fclose(fp_outann);
    fclose(fp_log);
    fclose(fp_monitor);
    fclose(fp_veg);
		
	if(NOTICE == 1){
		printf("done\n");
	}
}

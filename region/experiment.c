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
extern float tmp10_soil_d[PERIOD][366];
extern float fdat[NROW*NCOL];
extern int idat[NROW*NCOL];

extern struct Loct loct2[NROW*NCOL];
extern struct Echar echar2[NROW*NCOL];
extern struct Flux flux2[NROW*NCOL];

#if OUT_DAY==1
extern float outdat01[NROW*NCOL];
extern float outdat02[NROW*NCOL];
extern float outdat03[NROW*NCOL];
extern float outdat04[NROW*NCOL];
#endif

/* annual maps **********/
extern float out_a[N_ANNRES][NROW*NCOL];

double mean_veg[20][32];

/* SIMULATION ***********************************************************/
/* area */
void f_experiment(
	struct Grid grid[], 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Echar echar_type[MAX_BIOME], 
	struct Mass mass[], 
	struct Flux *flux, 
	FILE *fp_r[N_OFILE]
){
	long e, f, h, i, j;
	char num[8], filename[128];
	long ndy, end_year, calc_flag, pstart, pend;
	float gpp_a, npp_a, nep_a, lai_a, plant_a, soil_a, xx1_a, xx2_a, xx3_a, ch4_a, nn;
    float gpp_ga, npp_ga, nep_ga, plant_ga, soil_ga, prec_ga, rdata, vps;
	/* FILE *fp_o; */
	FILE *fp_restart;
	FILE *fp_clim[N_CLIMD];
	FILE *fp_out[1+24], *fp_log;
    FILE *fp_monitor, *fp_veg;
    FILE *fp_lai;
    
    /* experimental phase */
    loct->phase = 2;
    
    if(P_MONI == 0){
        pstart = 0;
        pend = WGRIDS -1;
    }else{
        pstart = P_MONI;
        pend = P_MONI;
    }
	
	if(NOTICE == 1){
		printf("Start experimental phase\n");
	}
	
	if(WMODE == 1){
		for(e=1;e<N_OFILE;e++){
			fprintf(fp_r[e],"%s %s\n", echar->para_ver_id, echar->para_date_id);
		}
	}else if(WMODE == 2){
		for(e=1;e<N_OFILE;e++){
			fprintf(fp_r[e],"%s %s\n", echar_type[0].para_ver_id, echar_type[0].para_date_id);
		}
	}
    
    for(i=0;i<20;i++){
        for(j=0;j<32;j++){
            mean_veg[i][j] = 0.0;
        }
    }
			
	/***************************************************************************/
	if(USE_RESTART == 1){
		strcpy(filename, grid[0].site_id);
		strcat(filename, "_restart.txt");
		if( (fp_restart = fopen(filename,"rt")) == NULL){
			printf("! No restart data !\n");
			exit (1);
		}
		for(i=pstart; i<=pend; i++){ 
			fscanf(fp_restart,"%f", &(mass[i].snwa));
			fscanf(fp_restart,"%f", &(mass[i].sw30));
			fscanf(fp_restart,"%f", &(mass[i].sww));

			fscanf(fp_restart,"%f", &((mass[i].tree).gdd));
			fscanf(fp_restart,"%f", &((mass[i].tree).cdd));
			fscanf(fp_restart,"%f", &((mass[i].tree).grwpd));
			fscanf(fp_restart,"%ld", &((mass[i].tree).day_frush));
			fscanf(fp_restart,"%ld", &((mass[i].tree).day_shed));
            
			fscanf(fp_restart,"%f", &((mass[i].c3).gdd));
			fscanf(fp_restart,"%f", &((mass[i].c3).cdd));
			fscanf(fp_restart,"%f", &((mass[i].c3).grwpd));
			fscanf(fp_restart,"%ld", &((mass[i].c3).day_frush));
			fscanf(fp_restart,"%ld", &((mass[i].c3).day_shed));
            
			fscanf(fp_restart,"%f", &((mass[i].c4).gdd));
			fscanf(fp_restart,"%f", &((mass[i].c4).cdd));
			fscanf(fp_restart,"%f", &((mass[i].c4).grwpd));
			fscanf(fp_restart,"%ld", &((mass[i].c4).day_frush));
			fscanf(fp_restart,"%ld", &((mass[i].c4).day_shed));
        
			fscanf(fp_restart,"%f", &((mass[i].tree).fol));
			fscanf(fp_restart,"%f", &((mass[i].tree).stm));
			fscanf(fp_restart,"%f", &((mass[i].tree).rot));
            
			fscanf(fp_restart,"%f", &((mass[i].c3).fol));
			fscanf(fp_restart,"%f", &((mass[i].c3).stm));
			fscanf(fp_restart,"%f", &((mass[i].c3).rot));
            
			fscanf(fp_restart,"%f", &((mass[i].c4).fol));
			fscanf(fp_restart,"%f", &((mass[i].c4).stm));
			fscanf(fp_restart,"%f", &((mass[i].c4).rot));
            
            fscanf(fp_restart,"%f", &((mass[i].tree).nsch_storage));
            fscanf(fp_restart,"%f", &((mass[i].c3).nsch_storage));
            fscanf(fp_restart,"%f", &((mass[i].c4).nsch_storage));

			fscanf(fp_restart,"%f", &((mass[i].soil).ltr_tf));
			fscanf(fp_restart,"%f", &((mass[i].soil).ltr_tc));
			fscanf(fp_restart,"%f", &((mass[i].soil).ltr_tr));
			fscanf(fp_restart,"%f", &((mass[i].soil).ltr_gf));
			fscanf(fp_restart,"%f", &((mass[i].soil).ltr_gc));
			fscanf(fp_restart,"%f", &((mass[i].soil).ltr_gr));
			fscanf(fp_restart,"%f", &((mass[i].soil).msl_a));
			fscanf(fp_restart,"%f", &((mass[i].soil).msl_i));
			fscanf(fp_restart,"%f", &((mass[i].soil).msl_p));
		}
	}
    
	/* annual mean output */
	strcpy(filename, grid[0].site_id);
	strcat(filename, grid[0].file_name);
	strcat(filename, "_ann");
	strcat(filename, ".flt");
	fp_out[0] = fopen(filename,"wb");
    
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
		loct->climy = e;
		loct->adyear = e;
		
		/* number of seasonal roop */
        if(LEAP_YR ==1){
            ndy = (e%4==0)?366:365; /* */
        }else{
            ndy = 365; /* no leap year */
        }
		
        loct->CO2y = e;
        if(FIX_CO2 == 1){
            loct->CO2y = BYR;
        }
		
		/* string of year ADXXXX */
		snprintf(num, 8, "%04d", (short)e);

        #if OUT_DAY==1
            f_create_dfile_global(loct->phase, loct->adyear, &grid[0], filename, fp_out);
        #endif
        
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
        if(strcmp(grid[0].site_id, "BAMIYAN")==0){
            f_open_bamiyan_clim(loct->phase, (short)e, fp_clim);
        }

		/* ************************************************/
        gpp_ga = npp_ga = nep_ga = plant_ga = soil_ga = prec_ga = 0.0;
        
        for(h=0; h<ndy; h++){
            loct->doy = h;
            doyTmody(e, h, &(loct->month), &(loct->mday));
            /* doyTmody(1979, h, &(loct->month), &(loct->mday)); */  /* no leap year */
            
            /* open meteorological data */
            if(loct->mday == 0 && strcmp(grid[0].site_id, "GLOBAL")==0){
                f_open_global_clim(loct->phase, (short)e, (short)(loct->month)+1, fp_clim);
            }
            if(loct->mday == 0 && strcmp(grid[0].site_id, "EASIA")==0){
                f_open_easia_clim(loct->phase, (short)e, (short)(loct->month)+1, fp_clim);
            }

            for(f=0; f<DSTEP; f++){
                loct->hour = f;

                if(strcmp(grid[0].site_id, "JAPAN") == 0){
                    //f_open_japan_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_japan_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid[0].site_id, "BB") == 0){
                    //f_open_bb_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_bb_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid[0].site_id, "JAPANc") == 0){
                    //f_open_japanc_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_japanc_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid[0].site_id, "JAPANh") == 0){
                    //f_open_japanh_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_japanh_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid[0].site_id, "JAPANk") == 0){
                    //f_open_japank_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday), (short)(loct->hour), fp_clim);
                    f_open_japank_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid[0].site_id, "PAWCs") == 0){
                    f_open_pawcs_clim(loct->phase, (short)e, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }

                printf("%4ld %3ld %2ld: ", e, h, f);
                fprintf(fp_log,"%ld %ld %ld ", e, h, f);

                /*********************************************************/
                if(strcmp(grid[0].site_id, "GLOBAL")==0){
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
                }else if(strcmp(grid[0].site_id, "BAMIYAN")==0){
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
                }else if(strcmp(grid[0].site_id, "EASIA")==0){
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
                }else if(strcmp(grid[0].site_id, "JAPAN")==0 || strcmp(grid[0].site_id, "BB")==0
                        || strcmp(grid[0].site_id, "JAPANc")==0|| strcmp(grid[0].site_id, "JAPANh")==0
                        || strcmp(grid[0].site_id, "JAPANk")==0|| strcmp(grid[0].site_id, "PAWCs")==0){
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
                
                for(i=pstart; i<=pend; i++){  
                    loct2[i].adyear = loct->adyear;
                    loct2[i].climy = loct->climy;
                    loct2[i].doy = loct->doy;
                    loct2[i].month = loct->month;
                    loct2[i].mday = loct->mday;
                    loct2[i].hour = loct->hour;
                    
                    if(h == 0){
                        for(j=0;j<N_ANNRES;j++){
                            out_a[j][i] = 0.0;
                        }
                    }
                }
                         
                /*****************************************************/
                /* OpenMP */
                  #pragma omp parallel for schedule(dynamic, 4) private(calc_flag)

                for(i=pstart; i<=pend; i++){ 
                    #if OUT_DAY==1
                    outdat01[i] = 0.0;
                    outdat02[i] = 0.0;
                    outdat03[i] = 0.0;
                    outdat04[i] = 0.0;
                    #endif

                    if(WMODE==2){
                        echar2[i] = echar_type[grid[i].veg_type];
                    }
                    
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
                        loct2[i].hour = f;
                        loct2[i].CO2y = e;
                        
                        if(FIX_CO2 == 1){
                            loct2[i].CO2y = BYR;
                        }

                        /* land-use */
                        /* flux2[i].lue_conv = 0.0;
                        flux2[i].lue_detr_1 = 0.0;
                        flux2[i].lue_detr_10 = 0.0;
                        flux2[i].lue_detr_100 = 0.0;
                        f_conv_landuse(&(grid[i]),&(loct2[i]));
                        f_emit_landuse(&(grid[i]), &(loct2[i]), &(mass[i]), &(flux2[i])); */
                        
                        /* local condition ************************************/
                        f_init_cond(&(grid[i]), &(loct2[i]), &(echar2[i]), &(mass[i]), &(flux2[i]));
                        
                        /* disturbance *****/
                        /* loct->climy = e; */
                        /*   disturbance_regime(e, &(grid[i]), &(loct2[i]), &(mass[i]), &(flux2[i])); */
                        
                        /* BAMIYAN grazing scenarios: 2017/09/03 by A.Ito */
                        if(strcmp(grid[0].site_id, "BAMIYAN")==0){
                            grazing_event(&(grid[i]), &(loct2[i]), &(mass[i]));
                        }

                        /* basic scheme ****************************************/
                        f_daily_scheme(&(grid[i]), &(loct2[i]), &(echar2[i]), &(mass[i]), &(flux2[i]));
                                            
                        /* summation for output *********/
                        #if OUT_DAY==1		
                        /* daily outputs */
                        
                        /* CFSR */
                        /* //outdat01[i] = flux2[i].gpp;
                        //outdat02[i] = flux2[i].er;
                        //outdat03[i] = flux2[i].nep;
                        //outdat04[i] = flux2[i].sr;
                        //outdat04[i] = loct2[i].lai; */

                        /* BAMIYAN */
                        outdat01[i] = flux2[i].gpp;
                        outdat02[i] = flux2[i].npp;
                        //outdat03[i] = ((mass[i].tree).plant + loct2[i].funder_c3 * (mass[i].c3).plant + loct2[i].funder_c4 * (mass[i].c4).plant);
                        //outdat04[i] = loct2[i].lai;
                        
                        outdat03[i] = flux2[i].nep;
                        outdat04[i] = (flux2[i].soil).ch4_wh;
                        
                        /* //outdat01[i] = (mass[i].tree).gdd;
                        //outdat02[i] = (mass[i].tree).lai;
                        //outdat03[i] = (float)(mass[i].tree).season;
                        //outdat04[i] = loct2[i].tair_dayav; */
                        #endif
                        
                        /* annual mean/sum outputs ***************/
                        out_a[0][i] += flux2[i].gpp;
                        out_a[1][i] += flux2[i].npp;
                        out_a[2][i] += flux2[i].nep;
                        //out_a[3][i] += flux2[i].nep;
                        out_a[3][i] += (flux2[i].soil).ch4_wh;

                        out_a[4][i] += loct2[i].lai/(double)ndy/(double)DSTEP;
                        out_a[5][i] += ((mass[i].tree).plant + loct2[i].funder_c3 * (mass[i].c3).plant
                                        + loct2[i].funder_c4 * (mass[i].c4).plant)/(double)ndy/(double)DSTEP;
                        out_a[6][i] += (mass[i].soil).soil/(double)ndy;
                        
                        out_a[7][i] += loct2[i].tmp_2m/(double)ndy/(double)DSTEP;
                        out_a[8][i] += loct2[i].ppfd_h/(double)ndy/(double)DSTEP;
                        out_a[9][i] += loct2[i].vpd/(double)ndy/(double)DSTEP;
                        out_a[10][i] += loct2[i].prate_sfc;
                        out_a[11][i] += loct2[i].aet;
                    }
                    
                    /* erosion */
                    /* f_erosion_rusle(&grid[i], loct, &mass[i], flux);
                    (mass[i].soil).msl_a -= (flux->soil).erosion_carbon; */
                    
                    /**********************************************/
                    /* monitoring */
                    if(WMODE==2){
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
                        gpp_a += (flux2[i].tree).gpp + loct2[i].funder_c3 * (flux2[i].c3).gpp +
                                    loct2[i].funder_c4 * (flux2[i].c4).gpp;
                        npp_a += (flux2[i].tree).npp + loct2[i].funder_c3 * (flux2[i].c3).npp +
                                   loct2[i].funder_c4 * (flux2[i].c4).npp;

                        /* //npp_a += (flux2[i].tree).rpg + loct2[i].funder_c3 * (flux2[i].c3).rpg +
                        //            loct2[i].funder_c4 * (flux2[i].c4).rpg; */

                        nep_a += flux2[i].nep;
                        lai_a += ((mass[i].tree).lai + loct2[i].funder_c3 * (mass[i].c3).lai +
                                  loct2[i].funder_c4 * (mass[i].c4).lai);
                        plant_a += ((mass[i].tree).plant + loct2[i].funder_c3 * (mass[i].c3).plant +
                                    loct2[i].funder_c4 * (mass[i].c4).plant); /* */
                        soil_a += (mass[i].soil).soil;
                        
                        xx1_a += loct2[i].tmp_2m;
                        xx2_a += loct2[i].prate_sfc;
                        xx3_a += loct2[i].ppfd_h;
                        
                        ch4_a += (flux2[i].soil).ch4_wh;
                        
                        gpp_ga += grid[i].area * ((flux2[i].tree).gpp + loct2[i].funder_c3 * (flux2[i].c3).gpp
                            + loct2[i].funder_c4 * (flux2[i].c4).gpp) / 1000000.0;
                        npp_ga += grid[i].area * ((flux2[i].tree).npp + loct2[i].funder_c3 * (flux2[i].c3).npp
                            + loct2[i].funder_c4 * (flux2[i].c4).npp) / 1000000.0;
                        nep_ga += grid[i].area * flux2[i].nep / 1000000.0;
                        plant_ga += grid[i].area * ((mass[i].tree).plant + loct2[i].funder_c3 * (mass[i].c3).plant +
                                    loct2[i].funder_c4 * (mass[i].c4).plant) / (float)(DSTEP*ndy) / 1000000.0;
                        soil_ga += grid[i].area * (mass[i].soil).soil / (float)(DSTEP*ndy) / 1000000.0;
                        prec_ga += grid[i].area * loct2[i].prate_sfc / 100.0;
                    }
                    /* //outdat01[i] = loct2[i].hangle;
                    //outdat02[i] = loct2[i].solhgt_h;
                    //outdat03[i] = loct2[i].toprad_h;
                    //outdat04[i] = loct2[i].ppfd_h; */

                    mean_veg[grid[i].veg_type][0] += (double)(grid[i].area / 24.0 / (double)ndy);
                    mean_veg[grid[i].veg_type][1] += (double)(grid[i].area * flux2[i].gpp);
                    mean_veg[grid[i].veg_type][2] += (double)(grid[i].area * flux2[i].npp);
                    mean_veg[grid[i].veg_type][3] += (double)(grid[i].area * flux2[i].nep);
                    mean_veg[grid[i].veg_type][4] += (double)(grid[i].area * loct2[i].lai / 24.0 / (double)ndy);
                    mean_veg[grid[i].veg_type][5] += (double)(grid[i].area * ((mass[i].tree).plant + loct2[i].funder_c3 * (mass[i].c3).plant
                        + loct2[i].funder_c4 * (mass[i].c4).plant) / 24.0 / (double)ndy);
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
                fprintf(fp_monitor,"%ld %ld %ld ",loct2[P_MONI].adyear,loct2[P_MONI].doy,loct2[P_MONI].hour);
                
                /* fprintf(fp_monitor,"%f %f %f ",loct2[P_MONI].tmp_2m,loct2[P_MONI].dswrf_sfc,loct2[P_MONI].prate_sfc);
                fprintf(fp_monitor,"%ld %f ",(mass[P_MONI].tree).season,(mass[P_MONI].tree).gdd);
                fprintf(fp_monitor,"%f %f %f ",flux2[P_MONI].gpp,flux2[P_MONI].npp,flux2[P_MONI].nep);
                fprintf(fp_monitor,"%f %f %f ",(mass[P_MONI].tree).fol,(mass[P_MONI].tree).stm,(mass[P_MONI].tree).rot); */
                
                fprintf(fp_monitor,"%f %f ",mass[P_MONI].sw30, mass[P_MONI].sww);
                fprintf(fp_monitor,"%f %f ",loct2[P_MONI].dswrf_sfc, loct2[P_MONI].tcdc_clm);
                fprintf(fp_monitor,"%f %f ",loct2[P_MONI].rn_short_eco, loct2[P_MONI].rn_long_eco);
                fprintf(fp_monitor,"%f %f %f %f ",loct2[P_MONI].prate_sfc, loct2[P_MONI].pet, loct2[P_MONI].aet, loct2[P_MONI].ro2);
                
                fprintf(fp_monitor,"\n");

                #if OUT_DAY==1
                if(e>=BYR && e<=EYR){
                    fwrite(outdat01, sizeof(float), NROW*NCOL, fp_out[1]);
                    fwrite(outdat02, sizeof(float), NROW*NCOL, fp_out[2]);
                    fwrite(outdat03, sizeof(float), NROW*NCOL, fp_out[3]);
                    fwrite(outdat04, sizeof(float), NROW*NCOL, fp_out[4]);
                }
                #endif
                
                if(strcmp(grid[0].site_id, "JAPAN")==0 || strcmp(grid[0].site_id, "BB")==0
                    || strcmp(grid[0].site_id, "JAPANc")==0 || strcmp(grid[0].site_id, "JAPANh")==0
                    || strcmp(grid[0].site_id, "JAPANk")==0){
                    fclose(fp_clim[0]);
                    fclose(fp_clim[1]);
                    fclose(fp_clim[2]);
                    fclose(fp_clim[3]);
                    fclose(fp_clim[4]);
                }
                if(strcmp(grid[0].site_id, "PAWCs")==0 && loct->mday==(month_day[loct->month]-1) && loct->hour==(DSTEP-1)){
                    fclose(fp_clim[0]);
                    fclose(fp_clim[1]);
                    fclose(fp_clim[2]);
                    fclose(fp_clim[3]);
                    fclose(fp_clim[4]);
                }
            }
            
            /* close climate files */
            j = 0;
            switch(loct->month){
                case 0: if(loct->mday==30){ j = 1; } break;
                case 1: if((e%4 == 0 && loct->mday==28) || (e%4 != 0 && loct->mday==27)){ j = 1; } break;
                //case 1: if(loct->mday==27){ j = 1; } break;
                case 2: if(loct->mday==30){ j = 1; } break;
                case 3: if(loct->mday==29){ j = 1; } break;
                case 4: if(loct->mday==30){ j = 1; } break;
                case 5: if(loct->mday==29){ j = 1; } break;
                case 6: if(loct->mday==30){ j = 1; } break;
                case 7: if(loct->mday==30){ j = 1; } break;
                case 8: if(loct->mday==29){ j = 1; } break;
                case 9: if(loct->mday==30){ j = 1; } break;
                case 10: if(loct->mday==29){ j = 1; } break;
                case 11: if(loct->mday==30){ j = 1; } break;
            }
            
            if(j == 1 && strcmp(grid[0].site_id, "EASIA")==0 ){
                fclose(fp_clim[0]);
                fclose(fp_clim[1]);
                fclose(fp_clim[2]);
                fclose(fp_clim[3]);
                fclose(fp_clim[4]);
            }

            if(j == 1 && strcmp(grid[0].site_id, "GLOBAL")==0 ){
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
		
		/* output annual maps */
		if(e>=BYR && e<=EYR){
			fwrite(out_a, sizeof(float), N_ANNRES*NROW*NCOL, fp_out[0]);
		}
        
        #if OUT_DAY==1
		for(j=1;j<=4;j++){
            fclose(fp_out[j]);
        }
        #endif	
        
        if(FIX_LAI == 2 || FIX_LAI == 3){
            fclose(fp_lai);
        }

        printf("ANNUAL %ld %f %f %f %f %f %f\n",loct->adyear, gpp_ga, npp_ga, nep_ga, plant_ga, soil_ga, prec_ga);

        /* if(loct->climy == 2015 && loct->month==3 && loct->mday==29){
            break;
        } */
        
        if(strcmp(grid[0].site_id, "BAMIYAN")==0){
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
    
	fclose(fp_out[0]);
    fclose(fp_log);
    fclose(fp_monitor);
    fclose(fp_veg);
		
	if(USE_RESTART == 1){
		fclose(fp_restart);
	}

	if(NOTICE==1){
		printf("done\n");
	}
}

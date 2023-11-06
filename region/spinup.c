/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
	modify: 2008 / 02 / 06 by Akihiko Ito
    modify: 2012 / 02 / 08 by Akihiko Ito
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

extern long	WMODE;
extern long WGRIDS;
extern long SPUPT;
extern float fdat[NROW * NCOL];
extern int idat[NROW * NCOL];
extern long    month_day[12];

extern struct Loct loct2[NROW*NCOL];
extern struct Echar echar2[NROW*NCOL];
extern struct Flux flux2[NROW*NCOL];

#if OUT_DAY==1
extern float outdat01[NROW*NCOL];
extern float outdat02[NROW*NCOL];
extern float outdat03[NROW*NCOL];
extern float outdat04[NROW*NCOL];
#endif

/* spin-up ************************************************************************/
void f_spinup(
	struct Grid grid[], 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Echar echar_type[MAX_BIOME], 
	struct Mass mass[], 
	struct Flux *flux, 
	FILE *fp_spinup
){
	long e, f, h, i, j, calc_flag, sasu_flag, climyr, end_climyr;
	char filename[128];
	char num[8];
	long ndy, pstart, pend;
	float gpp_a, npp_a, nep_a, lai_a, plant_a, soil_a, xx[20], ch4_a, nn, mm;
    float gpp_ga, npp_ga, nep_ga, plant_ga, soil_ga;
    float wi, ti, ll, rdata, vps;
	FILE *fp_clim[N_CLIMD], *fp_error;
	FILE *fp_out[1+4], *fp_log;
    FILE *fp_ss_grid, *fp_ss_loct, *fp_ss_mass, *fp_ss_flux, *fp_ss_echar;
	struct Grid grid0;
	struct Loct loct0;
	struct Echar echar0;
	struct Mass mass0;
	struct Flux flux0;
    
    loct->phase = 1;
	
	/* roop to stable stage ***********/
    if(P_MONI == 0){
        pstart = 0;
        pend = WGRIDS -1;
    }else{
        pstart = P_MONI;
        pend = P_MONI;
    }
    
	/* strcpy(filename, grid[0].site_id);
	strcat(filename, "_restart.txt");
	fp_restart = fopen(filename,"wt"); */
    
	strcpy(filename, grid[0].site_id);
	strcat(filename, "_restart_grid.flt");
	fp_ss_grid = fopen(filename,"wb");
    
	strcpy(filename, grid[0].site_id);
	strcat(filename, "_restart_loct.flt");
	fp_ss_loct = fopen(filename,"wb");
    
	strcpy(filename, grid[0].site_id);
	strcat(filename, "_restart_mass.flt");
	fp_ss_mass = fopen(filename,"wb");
    
	strcpy(filename, grid[0].site_id);
	strcat(filename, "_restart_flux.flt");
	fp_ss_flux = fopen(filename,"wb");
    
	strcpy(filename, grid[0].site_id);
	strcat(filename, "_restart_echar.flt");
	fp_ss_echar = fopen(filename,"wb");

	
    /* error file */
	fp_error = fopen("log_error.txt","wt");
    
    /* log file */
	fp_log = fopen("log_spinup.txt","wt");
    
	strcpy(filename, grid[0].site_id);
    strcat(filename, grid[0].file_name);
	strcat(filename, "_watch.flt");
	fp_out[0] = fopen(filename,"wb");

	if(NOTICE==1){
		printf("Start spin-up phase\n");
	}
	
	if(WMODE == 1){
		fprintf(fp_spinup,"%s %s\n", echar->para_ver_id, echar->para_date_id);
	}else if(WMODE == 2){
		fprintf(fp_spinup,"%s %s\n", echar_type[0].para_ver_id, echar_type[0].para_date_id);
	}
	
	for(i=pstart; i<=pend; i++){ 
		grid[i].time = 0;
		grid[i].age_stand = 0.0;
        
        //clear_b(&loct2[i], &echar2[i], &flux2[i]);
	}
	
    /* number of pararelization threds with OpenMP */
	#ifdef _OPENMP
	omp_set_num_threads(48); /* cores */
	#endif
    
    mm = 1.0;
	
	/* repetition *******************************************************/
    for(e=0 ; e<SPUPT ; e++){
        
        for(i=pstart; i<=pend; i++){ 
            if(grid[i].flag_datavl==1){
                /* stand age, year */
                grid[i].age_stand += 1.0;
            }
        }
        loct->adyear = BYR - SPUPT + e;
        
        /* randomized climate-data-year to remove trends and periodicity *******/
        //if(EYR>=2011){
        //    end_climyr = 2010; /* use data before 2011 */
        //}else{
            end_climyr = EYR;
        //}
        climyr = (long)(((double)rand()/(double)RAND_MAX)*(double)(end_climyr - BYR+1) + (double)BYR);
        
        /* for test */
        if(strcmp(grid[0].site_id, "JAPAN")==0 || strcmp(grid[0].site_id, "BB")==0
                || strcmp(grid[0].site_id, "JAPANc")==0|| strcmp(grid[0].site_id, "JAPANh")==0
                || strcmp(grid[0].site_id, "JAPANk")==0){
            climyr = 2020;
        }
        if(strcmp(grid[0].site_id, "PAWCs")==0){
            climyr = 2003;
        }

        loct->climy = climyr;
        
        /* open regional climate data ****/
        
        /* string of year ADXXXX */
        snprintf(num, 8, "%04d", (short)e);

        #if OUT_DAY==1
            //f_create_dfile_global(loct->phase, &grid[0], filename, fp_out);
        #endif
        
        if(e>=6 && (e%2==0)){
            for(i=pstart; i<=pend; i++){ 
                (mass[i].soil).sasu_li_tf = 0.0;
                (mass[i].soil).sasu_li_tc = 0.0;
                (mass[i].soil).sasu_li_tr = 0.0;
                (mass[i].soil).sasu_li_gf = 0.0;
                (mass[i].soil).sasu_li_gc = 0.0;
                (mass[i].soil).sasu_li_gr = 0.0;

                (mass[i].soil).sasu_hf_ha = 0.0;
                (mass[i].soil).sasu_hf_hi = 0.0;
                (mass[i].soil).sasu_hf_hp = 0.0;

                (mass[i].soil).sasu_mr_tf = 0.0;
                (mass[i].soil).sasu_mr_tc = 0.0;
                (mass[i].soil).sasu_mr_tr = 0.0;
                (mass[i].soil).sasu_mr_gf = 0.0;
                (mass[i].soil).sasu_mr_gc = 0.0;
                (mass[i].soil).sasu_mr_gr = 0.0;

                (mass[i].soil).sasu_mr_ha = 0.0;
                (mass[i].soil).sasu_mr_hi = 0.0;
                (mass[i].soil).sasu_mr_hp = 0.0;
            }
        }
        
        if(e>=6 && (e%2==0)){
            sasu_flag = 1;
        }else if(e>=7 && (e%2==1)){
            sasu_flag = 2;
        }else{
            sasu_flag = 0;
        }
        
        if((e>=(SPUPT - 5))){
            sasu_flag = 0;
        }
        
        /* update soil carbon pools by SASU estimates */
        if(USE_SASU == 1 && sasu_flag == 2){
            
            for(i=pstart; i<=pend; i++){ 
                if((mass[i].soil).sasu_mr_tf > 0.0){
                    (mass[i].soil).ltr_tf = (mass[i].soil).sasu_li_tf / (mass[i].soil).sasu_mr_tf;
                }else{
                    (mass[i].soil).ltr_tf = INIT_C;
                }
                
                if((mass[i].soil).sasu_mr_tc > 0.0){
                    (mass[i].soil).ltr_tc = (mass[i].soil).sasu_li_tc / (mass[i].soil).sasu_mr_tc;
                }else{
                    (mass[i].soil).ltr_tc = INIT_C;
                }
                
                if((mass[i].soil).sasu_mr_tr > 0.0){
                    (mass[i].soil).ltr_tr = (mass[i].soil).sasu_li_tr / (mass[i].soil).sasu_mr_tr;
                }else{
                    (mass[i].soil).ltr_tr = INIT_C;
                }
                
                if((mass[i].soil).sasu_mr_gf > 0.0){
                    (mass[i].soil).ltr_gf = (mass[i].soil).sasu_li_gf / (mass[i].soil).sasu_mr_gf;
                }else{
                    (mass[i].soil).ltr_gf = INIT_C;
                }
                
                if((mass[i].soil).sasu_mr_gc > 0.0){
                    (mass[i].soil).ltr_gc = (mass[i].soil).sasu_li_gc / (mass[i].soil).sasu_mr_gc;
                }else{
                    (mass[i].soil).ltr_gc = INIT_C;
                }
                
                if((mass[i].soil).sasu_mr_gr > 0.0){
                    (mass[i].soil).ltr_gr = (mass[i].soil).sasu_li_gr / (mass[i].soil).sasu_mr_gr;
                }else{
                    (mass[i].soil).ltr_gr = INIT_C;
                }
                
                if((mass[i].soil).sasu_mr_ha > 0.0){
                    (mass[i].soil).msl_a = (mass[i].soil).sasu_hf_ha / (mass[i].soil).sasu_mr_ha;
                }else{
                    (mass[i].soil).msl_a = INIT_C;
                }
                
                if((mass[i].soil).sasu_mr_hi > 0.0){
                    (mass[i].soil).msl_i = (mass[i].soil).sasu_hf_hi / (mass[i].soil).sasu_mr_hi;
                }else{
                    (mass[i].soil).msl_i = INIT_C;
                }
                
                if((mass[i].soil).sasu_mr_hp > 0.0){
                    (mass[i].soil).msl_p = (mass[i].soil).sasu_hf_hp / (mass[i].soil).sasu_mr_hp;
                }else{
                    (mass[i].soil).msl_p = INIT_C;
                }
            }
        }

        /* open meteorological data */
        if(strcmp(grid[0].site_id, "BAMIYAN")==0){
            f_open_bamiyan_clim(loct->phase, (short)climyr, fp_clim);
        }

        /*********************************************************/
        gpp_ga = npp_ga = nep_ga = plant_ga = soil_ga = 0.0;
        
        if(LEAP_YR ==1){
            ndy = (e%4==0)?366:365; /* */
        }else{
            ndy = 365; /* no leap year */
        }

        for(f=0;f<ndy;f++){
        
            loct->doy = f;
            doyTmody(climyr, f, &(loct->month), &(loct->mday));
            
            /* open meteorological data */
            if(loct->mday == 0 && strcmp(grid->site_id, "GLOBAL")==0){
                f_open_global_clim(loct->phase, (short)climyr, (short)(loct->month)+1, fp_clim);
            }
            if(loct->mday == 0 && strcmp(grid->site_id, "EASIA")==0){
                f_open_easia_clim(loct->phase, (short)climyr, (short)(loct->month)+1, fp_clim);
            }

            /* hourly roop ********************************************************/
            for(h=0;h<DSTEP;h++){

                loct->hour = h;
                printf("%4ld %4ld %3ld %3ld: ", loct->adyear, loct->climy, f, h);
                fprintf(fp_log,"%ld %ld %ld %ld ", loct->adyear, loct->climy, f, h);

                if(strcmp(grid->site_id, "JAPAN")==0){
                    f_open_japan_clim(loct->phase, (short)climyr, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid->site_id, "BB")==0){
                    f_open_bb_clim(loct->phase, (short)climyr, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid->site_id, "JAPANc")==0){
                    f_open_japanc_clim(loct->phase, (short)climyr, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid->site_id, "JAPANh")==0){
                    f_open_japanh_clim(loct->phase, (short)climyr, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid->site_id, "JAPANk")==0){
                    f_open_japank_clim(loct->phase, (short)climyr, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }
                if(strcmp(grid->site_id, "PAWCs")==0){
                    f_open_pawcs_clim(loct->phase, (short)climyr, (short)(loct->month)+1, (short)(loct->mday)+1, (short)(loct->hour)+1, fp_clim);
                }

                /* read regional climate data ***********/
                if(strcmp(grid[0].site_id, "GLOBAL")==0){
                    
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
                    if(h==0){
                        fread(fdat, 4, WGRIDS, fp_clim[2]);
                        for(i=pstart; i<=pend; i++){ 
                            grid[i].tmax_region = fdat[i] - ZAT; /* temp max, K */
                            grid[i].tmin_region = fdat[i] - ZAT; /* temp min, K */
                        }
                    }

                    /* precipitation */
                    if(h==0){
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
                    if(h==0){
                        fread(fdat, 4, WGRIDS, fp_clim[3]);
                        for(i=pstart; i<=pend; i++){ 
                            grid[i].humd_region = fdat[i];	/* specific humidity, g/g */
                        }
                    }
                    
                    /* wind */
                    if(h==0){
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
                        || strcmp(grid[0].site_id, "JAPANc")==0 || strcmp(grid[0].site_id, "JAPANh")==0
                        || strcmp(grid[0].site_id, "JAPANk")==0 || strcmp(grid[0].site_id, "PAWCs")==0){
                    fread(fdat, 4, NROW * NCOL, fp_clim[0]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        grid[i].prec_region = fdat[i] * 3600.0; /* precipitation, mm/h */
                    }
                    
                    fread(fdat, 4, NROW * NCOL, fp_clim[1]);
                    for(i=pstart; i<=pend; i++){
                        if(fdat[i] < 0.0){
                            fdat[i] = 0.0;
                        }
                        grid[i].srad_region = fdat[i] / 3600.0;      /* dsw rad, W/m2 */
                    }
                    
                    fread(fdat, 4, NROW * NCOL, fp_clim[2]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].tmax_region = fdat[i] - ZAT; /* temp max, K */
                        grid[i].tmin_region = fdat[i] - ZAT; /* temp min, K */
                    }

                    fread(fdat, 4, NROW * NCOL, fp_clim[3]);
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

                    fread(fdat, 4, NROW * NCOL, fp_clim[4]);
                    for(i=pstart; i<=pend; i++){
                        grid[i].wind_region = fdat[i];  /* wind, m -s */
                    }
                }
                    
                for(i=pstart; i<=pend; i++){ 
                    loct2[i].adyear = loct->adyear;
                    loct2[i].climy = loct->climy;
                    loct2[i].doy = loct->doy;
                    loct2[i].month = loct->month;
                    loct2[i].mday = loct->mday;
                    loct2[i].hour = loct->hour;
                }
                
                if(e == 0){
                    for(i=pstart; i<=pend; i++){
                        grid[i].tmp2m_ann += (grid[i].tmax_region + grid[i].tmin_region)/2.0/(float)ndy/(float)DSTEP;
                        grid[i].prec_ann += grid[i].prec_region;
                        grid[i].tsoil_annav = (grid[i].tsoil_annav*(mm-1.0) + (grid[i].tmax_region + grid[i].tmin_region)/2.0)/mm;
                    }
                }
                
                /* ****************************************************/
                /* OpenMP */
                #pragma omp parallel for schedule(dynamic, 4) private(calc_flag)
                
                for(i=pstart; i<=pend; i++){
                    /* printf("%ld\n", i); */
                    
                    /* ecophysiological parameter */
                    if(WMODE == 2){
                        echar2[i] = echar_type[grid[i].veg_type];
                    }
                    
                    calc_flag = 1; 
                    if(grid[i].flag_datavl!=1){
                        calc_flag = 0;
                    }else{ 	; }
                    if(grid[i].prec_region<-0.1 || grid[i].prec_region>1500.0){
                        calc_flag = 0;
                        printf("P %d\n", grid[i].flag_datavl);
                    }else{ ; }
                    if(grid[i].srad_region<-1.0 || grid[i].srad_region>2000.0){
                        calc_flag = 0;
                        printf("R %d\n", grid[i].flag_datavl);
                    }else{ ; }
                    if(grid[i].tmax_region<-90.0 || grid[i].tmax_region>90.0){
                        calc_flag = 0;
                        printf("TX %d\n", grid[i].flag_datavl);
                    }else{ ; }
                    if(grid[i].tmin_region<-90.0 || grid[i].tmin_region>90.0){
                        calc_flag = 0;
                        printf("TN %d\n", grid[i].flag_datavl);
                    }else{ ; }
                    if(grid[i].humd_region<-1.0 || grid[i].humd_region>100.0){
                        calc_flag = 0;
                        printf("H %d\n", grid[i].flag_datavl);
                    }else{ ; }
                    if(grid[i].wind_region<-120.0 || grid[i].wind_region>120.0){
                        calc_flag = 0;
                        printf("W %d\n", grid[i].flag_datavl);
                    }else{ ; } /* */
                    /*
                     if(calc_flag == 0) insufficient land data and no simulation 
                     */

                    grid[i].calc_flag = calc_flag;

                    if(calc_flag == 1 && (i+0)%RUNFAST==0){
                    //if(calc_flag == 1 && (i+0)%RUNFAST==0 && grid[i].veg_type ==2){
                        
                        loct2[i].CO2y = BYR; /* initial fixed */
                        grid[i].time++;
                                                        
                        /*** set environmental condition ***/
                        f_init_cond(&grid[i], &loct2[i], &echar2[i], &mass[i], &flux2[i]);
                        
                        /* if(nn==nn){
                            printf("Temp: %f %f %f\n", grid[i].tmax_region, grid[i].tmin_region, loct2[i].tmp_sfc);
                        } */
                        // printf("Q10: %d %f\n", grid[i].veg_type, (echar2[i].tree).qTf);
                        
                        if(loct->adyear && f==0){
                            flux->lue_conv = 0.0; 
                            flux->lue_detr_1 = 0.0; 
                            flux->lue_detr_10 = 0.0; 
                            flux->lue_detr_100 = 0.0; 
                            f_conv_landuse(&(grid[i]),loct);
                            f_emit_landuse(&(grid[i]), loct, &(mass[i]), flux);
                        }

                        /* **** disturbance *****/
                        /* loct->climy = e; */
                        /* disturbance_regime(loct->adyear, &grid[i], &loct2[i], &mass[i], &flux2[i]);  */
                        
                        /* BAMIYAN grazing scenarios: 2017/09/03 by A.Ito */
                        if(strcmp(grid[0].site_id, "BAMIYAN")==0){
                            grazing_event(&(grid[i]), &(loct2[i]), &(mass[i]));
                        }

                        /* loct2[i].climy = BYR; */
                        
                        /* basic scheme ******************************************/
                        f_daily_scheme(&grid[i], &loct2[i], &echar2[i], &mass[i], &flux2[i]);
                        
                        #if OUT_DAY==1		/* daily outputs */
                        /* CFSR run */
                        outdat01[i] = flux2[i].gpp;
                        outdat02[i] = flux2[i].er;
                        outdat03[i] = flux2[i].nep;
                        outdat04[i] = loct2[i].rn_eco;
                        #endif
                        
                        /* cumulative parameters for SASU */
                        if(USE_SASU == 1 && sasu_flag == 1){
                            (mass[i].soil).sasu_li_tf += (flux2[i].soil).li_tf;
                            (mass[i].soil).sasu_li_tc += (flux2[i].soil).li_tc;
                            (mass[i].soil).sasu_li_tr += (flux2[i].soil).li_tr;
                            
                            (mass[i].soil).sasu_li_gf += (flux2[i].soil).li_gf;
                            (mass[i].soil).sasu_li_gc += (flux2[i].soil).li_gc;
                            (mass[i].soil).sasu_li_gr += (flux2[i].soil).li_gr;
                            
                            (mass[i].soil).sasu_mr_tf += (echar2[i].soil).sr_lf/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_l;
                            (mass[i].soil).sasu_mr_tc += (echar2[i].soil).sr_lc/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_l;
                            (mass[i].soil).sasu_mr_tr += (echar2[i].soil).sr_lr/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_l;
                            
                            (mass[i].soil).sasu_mr_gf += (echar2[i].soil).sr_lf/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_l;
                            (mass[i].soil).sasu_mr_gc += (echar2[i].soil).sr_lc/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_l;
                            (mass[i].soil).sasu_mr_gr += (echar2[i].soil).sr_lr/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_l;

                            (mass[i].soil).sasu_hf_ha += (flux2[i].soil).hf_tfa + (flux2[i].soil).hf_tca + (flux2[i].soil).hf_tra +
                                                        (flux2[i].soil).hf_gfa + (flux2[i].soil).hf_gca + (flux2[i].soil).hf_gra;
                            (mass[i].soil).sasu_hf_hi += (flux2[i].soil).hf_tfi + (flux2[i].soil).hf_tci + (flux2[i].soil).hf_tri +
                                                        (flux2[i].soil).hf_gfi + (flux2[i].soil).hf_gci + (flux2[i].soil).hf_gri;
                            (mass[i].soil).sasu_hf_hp += (flux2[i].soil).hf_tfp + (flux2[i].soil).hf_tcp + (flux2[i].soil).hf_trp +
                                                        (flux2[i].soil).hf_gfp + (flux2[i].soil).hf_gcp + (flux2[i].soil).hf_grp;

                            (mass[i].soil).sasu_mr_ha += (echar2[i].soil).sr_ha/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_h;
                            (mass[i].soil).sasu_mr_hi += (echar2[i].soil).sr_hi/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_h;
                            (mass[i].soil).sasu_mr_hp += (echar2[i].soil).sr_hp/1000.0 / (double)DSTEP * (echar2[i].soil).f_tm_h;
                        }

                        if(flux->npp > grid->npp_max){
                            grid->npp_max = flux->npp;
                        }
                        if(grid->npp_max < 1.0){
                            grid->npp_max = 1.0;
                        }
                        
                        loct->m_casa_pre = loct->m_casa;
                        loct->vmc_pre = loct->vmc; /* */
                    }else{
                        grid[i].calc_flag = 0;
                        loct2[i].ppfd_h = 0.0;
                    }
                    /* erosion */
                    /* f_erosion_rusle(&grid[i], &loct2[i], &mass[i], &flux2[i]); */
                    /* if(grid[i].time > 100){
                        (mass[i].soil).msl_a -= (flux2[i].soil).erosion_carbon;
                    } */
                }
                #pragma omp barrier
                
                //printf("MONITOR: %ld %6.1f: ", (mass[P_MONI].tree).season, (mass[P_MONI].tree).gdd); /**********************/
                //printf("%6.1f %6.1f %6.1f: ", (mass[P_MONI].tree).lai, (mass[P_MONI].c3).lai, (mass[P_MONI].c4).lai);
                //printf("%6.1f %6.1f: ", loct2[P_MONI].tair_dayav, loct2[P_MONI].tair_dayav_c);
                
                gpp_a = npp_a = nep_a = lai_a = plant_a = soil_a = ch4_a = 0.0;
                nn = 0.0;
                for(i=0; i<20; i++){
                    xx[i] = 0.0;
                }
                
                for(i=pstart; i<=pend; i++){
                    
                    grid[i].n_clim += 1.0;
                    grid[i].tmp_clim += (grid[i].tmax_region + grid[i].tmin_region)/2.0/grid[i].n_clim;
                    grid[i].prec_clim += 365.0 * grid[i].prec_region/grid[i].n_clim;
                    grid[i].aet_clim += 365.0 * loct2[i].aet/grid[i].n_clim;
                    grid[i].pet_clim += 365.0 * loct2[i].pet/grid[i].n_clim;
               
                    if(grid[i].calc_flag == 1){
                        nn += 1.0;
                        gpp_a += (flux2[i].tree).gpp + loct2[i].funder_c3 * (flux2[i].c3).gpp + loct2[i].funder_c4 * (flux2[i].c4).gpp;
                        npp_a += (flux2[i].tree).npp + loct2[i].funder_c3 * (flux2[i].c3).npp + loct2[i].funder_c4 * (flux2[i].c4).npp;
                        nep_a += flux2[i].nep;
                        lai_a += ((mass[i].tree).lai + loct2[i].funder_c3 * (mass[i].c3).lai +
                                  loct2[i].funder_c4 * (mass[i].c4).lai);
                        plant_a += ((mass[i].tree).plant + loct2[i].funder_c3 * (mass[i].c3).plant +
                                    loct2[i].funder_c4 * (mass[i].c4).plant); /* */
                        //plant_a += ((mass[i].tree).plant);
                        soil_a += (mass[i].soil).soil;
                        
                        ch4_a += (flux2[i].soil).ch4_wh;
                        
                        /* for debugging: 2011/05/18 by A.Ito */
                        if((mass[i].soil).soil>=0.0 && (mass[i].soil).soil<=1000.0){
                            ;
                        }else{
                            fprintf(fp_error,"%ld %ld %ld %lf\n",e,f,i,(mass[i].soil).soil);
                            //fprintf(fp_error,"%.1lf %.1lf %.1lf ",(mass[i].soil).msl_a,(mass[i].soil).msl_i,(mass[i].soil).msl_p);
                            fprintf(fp_error,"\n");
                        }
                        
                        xx[0] += loct2[i].tmp_2m;
                        xx[1] += loct2[i].tmp_sfc;
                        //xx[2] += loct2[i].vpd;
                        xx[2] += loct2[i].dswrf_sfc;
                        xx[3] += loct2[i].ppfd_h;
                        xx[4] += (echar2[i].c3).psat;
                        xx[5] += (mass[i].c3).gdd;
                        xx[6] += loct2[i].soilwtr_l;
                        xx[7] += loct2[i].soilwtr_h;
                        xx[8] += loct2[i].aet;
                        xx[9] += loct2[i].ro2;

                        xx[10] += loct2[i].rn_eco;
                        xx[11] += (mass[i].tree).lai;
                        xx[12] += (mass[i].c3).lai;
                        xx[13] += (mass[i].c4).lai;
                        xx[14] += loct2[i].evpr;
                        xx[15] += loct2[i].trnsp;
                        
                        if(e>=1){
                           /* LAI max constrained by Iio et al. (2014): 2016/08/15 by A.Ito */
                            if(grid[i].pet_clim > 0.0){
                                wi = grid[i].prec_clim / grid[i].pet_clim;
                            }else{
                                wi = 0.01;
                            }
                            wi = log10(wi);
                            
                            ti = grid[i].tmp_clim;
                            
                            ll = 0.532191 + 0.593486*wi - 0.441269*wi*wi*wi - 0.245266*wi*wi + 0.0164922*ti - 0.000620463*ti*ti;
                            
                            /* correct: 2016/11/15 by A.Ito */
                            /* ll = log10(ll); */
                            ll = powf(10.0, ll);
                            
                            if(ll>1.0 && ll<12.0){
                                ;
                            }else{
                                if(ll<=1.0){
                                    ll = 1.0;
                                }else if(ll>=12.0){
                                    ll = 12.0;
                                }else{
                                    ;
                                }
                            }
                            
                            loct2[i].laimax_meta = ll;
                        }else{
                            loct2[i].laimax_meta = 3.0;
                        }
                        
                        xx[7] += loct2[i].laimax_meta;
                        
                        //xx[0] += loct2[i].tmp_2m;
                        //xx[0] += (echar2[i].tree).opt_lai;
                        //xx[1] += loct2[i].prate_sfc;
                        //xx[1] += (echar2[i].tree).malloc_f;
                        //xx[1] += (echar2[i].tree).psat;
                        //xx[2] += loct2[i].ppfd_h;
                        //xx[3] += loct2[i].soilwtr_l;
                        //xx[4] += loct2[i].soilwtr_h;
                        //xx[5] += loct2[i].vpd;
                        //xx[6] += (mass[i].tree).lai;
                        //xx[7] += (mass[i].tree).plant;
                        
                        gpp_ga += grid[i].area * ((flux2[i].tree).gpp + loct2[i].funder_c3 * (flux2[i].c3).gpp + loct2[i].funder_c4 * (flux2[i].c4).gpp) / 1000000.0;
                        npp_ga += grid[i].area * ((flux2[i].tree).npp + loct2[i].funder_c3 * (flux2[i].c3).npp + loct2[i].funder_c4 * (flux2[i].c4).npp) / 1000000.0;
                        nep_ga += grid[i].area * flux2[i].nep / 1000000.0;
                        plant_ga += grid[i].area * ((mass[i].tree).plant + loct2[i].funder_c3 * (mass[i].c3).plant +
                                    loct2[i].funder_c4 * (mass[i].c4).plant) / (float)(DSTEP*ndy) / 1000000.0;
                        soil_ga += grid[i].area * (mass[i].soil).soil / (float)(DSTEP*ndy) / 1000000.0;
                    }
                    #if OUT_DAY==1
                    outdat04[i] = loct2[i].ppfd_h;
                    #endif
                }
                printf("%8.4f %8.4f %8.4f : %10.2f %10.2f %10.2f: %10.3f: ",
                       100.0*gpp_a/nn, 100.0*npp_a/nn, 100.0*nep_a/nn, 
                       lai_a/nn, plant_a/nn, soil_a/nn, ch4_a/nn);
                /* printf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.4f: ",
                       xx[0]/nn, xx[1]/nn, xx[2]/nn, xx[3]/nn, xx[4]/nn, xx[5]/nn, xx[6]/nn, xx[7]/nn, xx[8]/nn, xx[9]/nn);
                printf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f:  ",
                       xx[10]/nn, xx[11]/nn, xx[12]/nn, xx[13]/nn, xx[14]/nn, xx[15]/nn); */
                printf("%.0f\n", nn);

                fprintf(fp_log,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
                       nn, gpp_a, npp_a, nep_a, lai_a, plant_a, soil_a, ch4_a,
                       xx[0], xx[1], xx[2], xx[3], xx[4], xx[5], xx[6], xx[7], xx[8], xx[9]);
                /* fprintf(fp_log,"%12.2f %12.2f %12.2f %12.2f %12.2f %12.2f\n",
                       100.0*gpp_a/nn, 100.0*npp_a/nn, 100.0*nep_a/nn, 
                       lai_a/nn, plant_a/nn, soil_a/nn); */
                fprintf(fp_log,"\n");
                
                /*****************************/
                #if OUT_DAY==1
                if(e==BYR && f==0){
                    fwrite(outdat04,sizeof(float),WGRIDS, fp_out[0]);
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
            }
            if(strcmp(grid[0].site_id, "PAWCs")==0 && loct->mday==(month_day[loct->month]-1) && loct->hour==(DSTEP-1)){
                fclose(fp_clim[0]);
                fclose(fp_clim[1]);
                fclose(fp_clim[2]);
                fclose(fp_clim[3]);
                fclose(fp_clim[4]);
            }

            j = 0;
            switch(loct->month){
                case 0: if(loct->mday==30){ j = 1; } break;
                case 1: if( (climyr%4 != 0 && loct->mday==27) || (climyr%4 == 0 && loct->mday==28) ){ j = 1; } break;
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
            
            if(j == 1 && (strcmp(grid[0].site_id, "EASIA")==0)){
                fclose(fp_clim[0]);
                fclose(fp_clim[1]);
                fclose(fp_clim[2]);
                fclose(fp_clim[3]);
                fclose(fp_clim[4]);
            }
            
            if(j == 1 && (strcmp(grid[0].site_id, "GLOBAL")==0)){
                fclose(fp_clim[0]);
                fclose(fp_clim[1]);
                fclose(fp_clim[2]);
                fclose(fp_clim[3]);
                fclose(fp_clim[4]);
                fclose(fp_clim[5]);
                fclose(fp_clim[6]);
            }

            mm += 1.0;
        }
        
        printf("ANNUAL %ld %f %f %f %f %f\n",loct->adyear, gpp_ga, npp_ga, nep_ga, plant_ga, soil_ga);
        
        if(strcmp(grid[0].site_id, "BAMIYAN")==0){
            fclose(fp_clim[0]);
            fclose(fp_clim[1]);
            fclose(fp_clim[2]);
            fclose(fp_clim[3]);
            fclose(fp_clim[4]);
        }
    }

    fclose(fp_out[0]);
    
    /* RESTART FILES ***********/
    /* for(i=pstart; i<=pend; i++){
        fprintf(fp_restart,"%f ", mass[i].snwa);
        fprintf(fp_restart,"%f ", mass[i].sw30);
        fprintf(fp_restart,"%f ", mass[i].sww);

        fprintf(fp_restart,"%f ", (mass[i].tree).gdd);
        fprintf(fp_restart,"%f ", (mass[i].tree).cdd);
        fprintf(fp_restart,"%f ", (mass[i].tree).grwpd);
        fprintf(fp_restart,"%ld ", (mass[i].tree).day_frush);
        fprintf(fp_restart,"%ld ", (mass[i].tree).day_shed);
        fprintf(fp_restart,"%f ", (mass[i].c3).gdd);
        fprintf(fp_restart,"%f ", (mass[i].c3).cdd);
        fprintf(fp_restart,"%f ", (mass[i].c3).grwpd);
        fprintf(fp_restart,"%ld ", (mass[i].c3).day_frush);
        fprintf(fp_restart,"%ld ", (mass[i].c3).day_shed);
        fprintf(fp_restart,"%f ", (mass[i].c4).gdd);
        fprintf(fp_restart,"%f ", (mass[i].c4).cdd);
        fprintf(fp_restart,"%f ", (mass[i].c4).grwpd);
        fprintf(fp_restart,"%ld ", (mass[i].c4).day_frush);
        fprintf(fp_restart,"%ld ", (mass[i].c4).day_shed);

        fprintf(fp_restart,"%f ", (mass[i].tree).fol);
        fprintf(fp_restart,"%f ", (mass[i].tree).stm);
        fprintf(fp_restart,"%f ", (mass[i].tree).rot);
        fprintf(fp_restart,"%f ", (mass[i].c3).fol);
        fprintf(fp_restart,"%f ", (mass[i].c3).stm);
        fprintf(fp_restart,"%f ", (mass[i].c3).rot);
        fprintf(fp_restart,"%f ", (mass[i].c4).fol);
        fprintf(fp_restart,"%f ", (mass[i].c4).stm);
        fprintf(fp_restart,"%f ", (mass[i].c4).rot);

        fprintf(fp_restart,"%f ", (mass[i].tree).nsch_storage);
        fprintf(fp_restart,"%f ", (mass[i].c3).nsch_storage);
        fprintf(fp_restart,"%f ", (mass[i].c4).nsch_storage);

        fprintf(fp_restart,"%f ", (mass[i].soil).ltr_tf);
        fprintf(fp_restart,"%f ", (mass[i].soil).ltr_tc);
        fprintf(fp_restart,"%f ", (mass[i].soil).ltr_tr);
        fprintf(fp_restart,"%f ", (mass[i].soil).ltr_gf);
        fprintf(fp_restart,"%f ", (mass[i].soil).ltr_gc);
        fprintf(fp_restart,"%f ", (mass[i].soil).ltr_gr);
        fprintf(fp_restart,"%f ", (mass[i].soil).msl_a);
        fprintf(fp_restart,"%f ", (mass[i].soil).msl_i);
        fprintf(fp_restart,"%f\n", (mass[i].soil).msl_p);
    } */
    
    for(f=0;f<NCOL*NROW;f++){
        fwrite(&grid[f], sizeof(grid0), 1, fp_ss_grid);
        fwrite(&loct2[f], sizeof(loct0), 1, fp_ss_loct);
        fwrite(&mass[f], sizeof(mass0), 1, fp_ss_mass);
        fwrite(&flux2[f], sizeof(flux0), 1, fp_ss_flux);
        fwrite(&echar2[f], sizeof(echar0), 1, fp_ss_echar);
        //fwrite(&grid[f], sizeof(&grid[f]), 1, fp_ss_grid);
        //fwrite(&loct2[f], sizeof(&loct2[f]), 1, fp_ss_loct);
        //fwrite(&mass[f], sizeof(&mass[f]), 1, fp_ss_mass);
        //fwrite(&flux2[f], sizeof(&flux2[f]), 1, fp_ss_flux);
        //fwrite(&echar2[f], sizeof(&echar2[f]), 1, fp_ss_echar);
    }

	if(NOTICE==1){
		printf("done\n");
	}
	
	fclose(fp_error);
    fclose(fp_log);

	fclose(fp_ss_grid);
	fclose(fp_ss_loct);
	fclose(fp_ss_mass);
	fclose(fp_ss_flux);
	fclose(fp_ss_echar);
}

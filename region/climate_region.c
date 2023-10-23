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

/* **********************************************/
void f_open_global_clim(
    short stage,
    short yr, 
    short mon,
    FILE *fp_clim[N_CLIMD]
){
	char syear[8];
	char smon[8];
	char filename[128];
    short yrc;
    
    yrc = yr;
    
    /* ***********/
	snprintf(syear, 8, "%04d", yr);
	snprintf(smon, 8, "%02d", mon);
    
    /* prec */
    strcpy(filename, "./data/TP/");
    strcat(filename, syear);
    strcat(filename, "/ERA5_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_TP.flt");
    if((fp_clim[0] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }

    /* srad */
    strcpy(filename, "./data/SSRD/");
    strcat(filename, syear);
    strcat(filename, "/ERA5_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_SSRD.flt");
    if((fp_clim[1] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }

    /* trad */
    strcpy(filename, "./data/STRD/");
    strcat(filename, syear);
    strcat(filename, "/ERA5_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_STRD.flt");
    if((fp_clim[2] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }

    /* temp */
    strcpy(filename, "./data/2T/");
    strcat(filename, syear);
    strcat(filename, "/ERA5_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_2T.flt");
    if((fp_clim[3] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }

    /* humidity */
    strcpy(filename, "./data/2D/");
    strcat(filename, syear);
    strcat(filename, "/ERA5_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_2D.flt");
    if((fp_clim[4] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }

    /* wind */
    strcpy(filename, "./data/10U/");
    strcat(filename, syear);
    strcat(filename, "/ERA5_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_10U.flt");
    if((fp_clim[5] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }

    strcpy(filename, "./data/10V/");
    strcat(filename, syear);
    strcat(filename, "/ERA5_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_10V.flt");
    if((fp_clim[6] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
}

/***********************************************/
void f_open_bamiyan_clim(
    short stage,
    short yr, 
    FILE *fp_clim[6]
){
	char syear[8], syear_r[8];
	char filename[128];
    short yrc;
    
    yrc = yr;
    
    /************/
	snprintf(syear, 8, "%04d", yr);
    snprintf(syear_r, 8, "%04d", 2042);
    
    /****************/
    yr = yrc;
    snprintf(syear, 8, "%04d", yr);

    /* humidity */
    if(EX_BAMIYAN==0){
        if(yr>=1979 && yr <=2016){
            strcpy(filename, "./data/ecmwf/hum");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==1){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/hum");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp26/hum");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==2){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/hum");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp85/hum");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }
    
    strcat(filename, syear);
    strcat(filename, ".bin");
    if((fp_clim[3] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* prec */
    if(EX_BAMIYAN==0){
        if(yr>=1979 && yr <=2016){
            strcpy(filename, "./data/ecmwf/rain");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==1){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/rain");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp26/rain");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==2){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/rain");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp85/rain");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }

    strcat(filename, syear);
    strcat(filename, ".txt");
    if((fp_clim[0] = fopen(filename,"rt")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* srad */
    if(EX_BAMIYAN==0){
        if(yr>=1979 && yr <=2016){
            strcpy(filename, "./data/ecmwf/rad");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==1){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/rad");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp26/rad");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==2){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/rad");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp85/rad");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }
    
    if(yr>=2043 && (EX_BAMIYAN==1 || EX_BAMIYAN==2)){
        strcat(filename, syear_r);
    }else{
        strcat(filename, syear);
    }
    strcat(filename, ".bin");
    if((fp_clim[1] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* temp */
    if(EX_BAMIYAN==0){
        if(yr>=1979 && yr <=2016){
            strcpy(filename, "./data/ecmwf/temp");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==1){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/temp");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp26/temp");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==2){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/temp");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp85/temp");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }

    strcat(filename, syear);
    strcat(filename, ".bin");
    if((fp_clim[2] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* wind */
    if(EX_BAMIYAN==0){
        if(yr>=1979 && yr <=2016){
            strcpy(filename, "./data/ecmwf/win");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==1){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/win");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp26/wind");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }else if(EX_BAMIYAN==2){
        if(yr>=1979 && yr <=2005){
            strcpy(filename, "./data/ecmwf/win");
        }else if(yr>=2006 && yr <=2055){
            strcpy(filename, "./data/rcp85/wind");
        }else{
            printf("BAD setting period\n");
            exit(1);
        }
    }

    strcat(filename, syear);
    strcat(filename, ".txt");
    if((fp_clim[4] = fopen(filename,"rt")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* replace CFSR2 ***********/
    /* if(yr>=2011){
        yr = 2006 + yr%5;
    }
	
	snprintf(syear, 8, "%04d", yr);
	snprintf(smon, 8, "%02d", mon); */
    
	/* printf("!! %s data\n", filename); */
	
}

void f_open_easia_clim(
    short stage,
    short yr,
    short mon,
    FILE *fp_clim[6]
){
    char syear[8];
    char smon[8];
    char filename[128];
    short yrc;
    
    yrc = yr;
    
    /************/
    snprintf(syear, 8, "%04d", yr);
    snprintf(smon, 8, "%02d", mon);
    
    /* humidity */
    strcpy(filename, "./data/d2m_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_easia_qd.flt2");
    if((fp_clim[3] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* prec */
    strcpy(filename, "./data/tp_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_easia_qd.flt2");
    if((fp_clim[0] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* srad */
    strcpy(filename, "./data/ssrd_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_easia_qd.flt2");
    if((fp_clim[1] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* temp */
    strcpy(filename, "./data/t2m_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_easia_qd.flt2");
    if((fp_clim[2] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* wind */
    /* strcpy(filename, "./data/wind10m_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, "_easia_qd.flt2"); */
    if((fp_clim[4] = fopen("dummy.flt","rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
}

void f_open_japan_clim(
    short stage,
    short yr,
    short mon,
    short day,
    short hour,
    FILE *fp_clim[6]
){
    char syear[8],smon[8],sday[8],shour[8];
    char filename[128];
    short yrc;
    
    yrc = yr;
    
    /************/
    snprintf(syear, 8, "%04d", yr);
    snprintf(smon, 8, "%02d", mon);
    snprintf(sday, 8, "%02d", day);
    snprintf(shour, 8, "%02d", hour-1);
    
    /* prec */
    strcpy(filename, "./data/TP/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JAPAN_ERA5int_TP_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[0] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* srad */
    strcpy(filename, "./data/SSRD/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JAPAN_ERA5int_SSRD_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[1] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* temp */
    strcpy(filename, "./data/2T/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JAPAN_ERA5int_2T_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[2] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* humidity */
   strcpy(filename, "./data/2D/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JAPAN_ERA5int_2D_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[3] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* wind */
   strcpy(filename, "./data/10U/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JAPAN_ERA5int_10U_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[4] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
}

void f_open_bb_clim(
    short stage,
    short yr,
    short mon,
    short day,
    short hour,
    FILE *fp_clim[6]
){
    char syear[8],smon[8],sday[8],shour[8];
    char filename[128];
    short yrc;
    
    yrc = yr;
    
    /* ***********/
    snprintf(syear, 8, "%04d", yr);
    snprintf(smon, 8, "%02d", mon);
    snprintf(sday, 8, "%02d", day);
    snprintf(shour, 8, "%02d", hour-1);
    
    /* prec */
    strcpy(filename, "./data/TP/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "BB1km_ERA5int_TP_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[0] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* srad */
    strcpy(filename, "./data/SSRD/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "BB1km_ERA5int_SSRD_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[1] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* temp */
    strcpy(filename, "./data/2T/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "BB1km_ERA5int_2T_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[2] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* humidity */
   strcpy(filename, "./data/2D/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "BB1km_ERA5int_2D_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[3] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* wind */
   strcpy(filename, "./data/10U/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "BB1km_ERA5int_10U_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[4] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
}

void f_open_japanc_clim(
    short stage,
    short yr,
    short mon,
    short day,
    short hour,
    FILE *fp_clim[6]
){
    char syear[8],smon[8],sday[8],shour[8];
    char filename[128];
    short yrc;
    
    yrc = yr;
    
    /* ***********/
    snprintf(syear, 8, "%04d", yr);
    snprintf(smon, 8, "%02d", mon);
    snprintf(sday, 8, "%02d", day);
    snprintf(shour, 8, "%02d", hour-1);
    
    /* prec */
    strcpy(filename, "./data/TP/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNc1km_ERA5_TP_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[0] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* srad */
    strcpy(filename, "./data/SSRD/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNc1km_ERA5_SSRD_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[1] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* temp */
    strcpy(filename, "./data/2T/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNc1km_ERA5_2T_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[2] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* humidity */
   strcpy(filename, "./data/2D/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNc1km_ERA5_2D_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[3] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* wind */
   strcpy(filename, "./data/10U/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNc1km_ERA5_10U_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[4] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
}

void f_open_japanh_clim(
    short stage,
    short yr,
    short mon,
    short day,
    short hour,
    FILE *fp_clim[6]
){
    char syear[8],smon[8],sday[8],shour[8];
    char filename[128];
    short yrc;
    
    yrc = yr;
    
    /* ***********/
    snprintf(syear, 8, "%04d", yr);
    snprintf(smon, 8, "%02d", mon);
    snprintf(sday, 8, "%02d", day);
    snprintf(shour, 8, "%02d", hour-1);
    
    /* prec */
    strcpy(filename, "./data/TP/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNh1km_ERA5_TP_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[0] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* srad */
    strcpy(filename, "./data/SSRD/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNh1km_ERA5_SSRD_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[1] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* temp */
    strcpy(filename, "./data/2T/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNh1km_ERA5_2T_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[2] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* humidity */
   strcpy(filename, "./data/2D/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNh1km_ERA5_2D_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[3] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* wind */
   strcpy(filename, "./data/10U/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNh1km_ERA5_10U_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[4] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
}

void f_open_japank_clim(
    short stage,
    short yr,
    short mon,
    short day,
    short hour,
    FILE *fp_clim[6]
){
    char syear[8],smon[8],sday[8],shour[8];
    char filename[128];
    short yrc;
    
    yrc = yr;
    
    /* ***********/
    snprintf(syear, 8, "%04d", yr);
    snprintf(smon, 8, "%02d", mon);
    snprintf(sday, 8, "%02d", day);
    snprintf(shour, 8, "%02d", hour-1);
    
    /* prec */
    strcpy(filename, "./data/TP/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNk1km_ERA5_TP_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[0] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* srad */
    strcpy(filename, "./data/SSRD/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNk1km_ERA5_SSRD_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[1] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* temp */
    strcpy(filename, "./data/2T/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNk1km_ERA5_2T_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[2] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* humidity */
   strcpy(filename, "./data/2D/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNk1km_ERA5_2D_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[3] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
    
    /* wind */
   strcpy(filename, "./data/10U/");
    strcat(filename, syear);
    strcat(filename, "/");
    strcat(filename, "JPNk1km_ERA5_10U_");
    strcat(filename, syear);
    strcat(filename, smon);
    strcat(filename, sday);
    strcat(filename, shour);
    strcat(filename, ".flt");
    if((fp_clim[4] = fopen(filename,"rb")) == NULL){
        printf("!! %s data\n", filename);
        exit (1);
    }
}

void f_open_pawcs_clim(
    short stage,
    short yr,
    short mon,
    short day,
    short hour,
    FILE *fp_clim[6]
){
    char syear[8],smon[8],sday[8],shour[8];
    char filename[128];
    short yrc;
    
    yrc = yr;
    
    if(day==1 && hour==1){
        
        /* ***********/
        snprintf(syear, 8, "%04d", yr);
        snprintf(smon, 8, "%02d", mon);
        snprintf(sday, 8, "%02d", day);
        snprintf(shour, 8, "%02d", hour-1);
        
        /* prec */
        strcpy(filename, "./data/TP/");
        strcat(filename, syear);
        strcat(filename, "/");
        strcat(filename, "PAWCs_TP_");
        strcat(filename, syear);
        strcat(filename, smon);
        strcat(filename, ".flt");
        if((fp_clim[0] = fopen(filename,"rb")) == NULL){
            printf("!! %s data\n", filename);
            exit (1);
        }
        
        /* srad */
        strcpy(filename, "./data/SSRD/");
        strcat(filename, syear);
        strcat(filename, "/");
        strcat(filename, "PAWCs_SSRD_");
        strcat(filename, syear);
        strcat(filename, smon);
        strcat(filename, ".flt");
        if((fp_clim[1] = fopen(filename,"rb")) == NULL){
            printf("!! %s data\n", filename);
            exit (1);
        }
        
        /* temp */
        strcpy(filename, "./data/2T/");
        strcat(filename, syear);
        strcat(filename, "/");
        strcat(filename, "PAWCs_2T_");
        strcat(filename, syear);
        strcat(filename, smon);
        strcat(filename, ".flt");
        if((fp_clim[2] = fopen(filename,"rb")) == NULL){
            printf("!! %s data\n", filename);
            exit (1);
        }
        
        /* humidity */
        strcpy(filename, "./data/2D/");
        strcat(filename, syear);
        strcat(filename, "/");
        strcat(filename, "PAWCs_2D_");
        strcat(filename, syear);
        strcat(filename, smon);
        strcat(filename, ".flt");
        if((fp_clim[3] = fopen(filename,"rb")) == NULL){
            printf("!! %s data\n", filename);
            exit (1);
        }
        
        /* wind */
        strcpy(filename, "./data/10U/");
        strcat(filename, syear);
        strcat(filename, "/");
        strcat(filename, "PAWCs_10U_");
        strcat(filename, syear);
        strcat(filename, smon);
        strcat(filename, ".flt");
        if((fp_clim[4] = fopen(filename,"rb")) == NULL){
            printf("!! %s data\n", filename);
            exit (1);
        }
    }else{
        ;
    }
}

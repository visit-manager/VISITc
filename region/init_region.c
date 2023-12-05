/*    VISIT: Vegetation Integrative SImulator for Tarce gases                */
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems        */
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC                    */
/* Carbon cycle, erosion, biomass burning, land-use change,                */
/* CH4 emission and oxidation, N2O emission,,,,,                        */
/*    version b 1.0.0    cerated in November 02, 2007                        */

/* [history]
    create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
    modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
    modify: 2007 / 12 / 01 by Akihiko Ito: Version 0.2.0.
*/

/** initialize environmental characteristics of each 1deg x 1deg grid **/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* climatology: *_d[] means the time-series during 1948 to 2004- */
extern float tmp_2m_ncep_dav[366];            /* 2m air temperature, degree Celcius */
extern float prate_sfc_ncep_mav[12];        /* precipitation, mm mon-1 */
extern long    month_day[12];
extern long    WGRIDS;
extern float fdat[NROW*NCOL];
extern unsigned char cdat[NROW*NCOL];

/* initialize ERA5 global conditions ****************************/
/** initialize climate and soil parameters by using the prepared data files **/
void f_init_global_run(
    struct Grid grid[]
){
    long i, j;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    if(NOTICE==1){
        printf(" Setting latitude...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = 90.0 - ((float)(i/NCOL) + 0.5) * (180.0 / (float)NROW);
    }
    if(NOTICE==1){
        printf("done\n");
    }
 
    /* longitude */
    if(NOTICE==1){
        printf(" Setting longitude...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = -180.0 - ((float)(i%NCOL) + 0.5) * (360.0 / (float)NCOL);
    }
    if(NOTICE==1){
        printf("done\n");
    }
 
    /* grid area */
    for(i=0; i<NROW; i++){
        garea = grid_area(90.0 - ((float)(i/NCOL)) * (180.0 / (float)NROW), 90.0 - ((float)(i/NCOL) + 1.0) * (180.0 / (float)NROW),  0.0, 360.0/(double)NCOL);
        
        for(j=0;j<NCOL;j++){
            grid[i*NCOL+j].area = garea;
        }
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/global_srtm_qd.flt","rb")) == NULL){
        printf("!! NO global_srtm_qd.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE == 1){
        printf(" Reading ./data/global_srtm_qd.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
        //grid[i].topo = 100.0;
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
 
    /* crop fraction *****************************/
    if((fp_dat = fopen("./data/global_crop2015_qd.flt","rb")) == NULL){
        printf("!! NO global_crop2015_qd.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    for(i=0; i<WGRIDS; i++){
        grid[i].agrPt = fdat[i];
    }
    if(NOTICE == 1){
        printf(" Reading ./data/global_crop2015_qd.flt...");
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* land cover ********/
    if((fp_dat = fopen("./data/global_synmap_qd.flt","rb")) == NULL){
        printf("!! NO global_synmap_qd.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE == 1){
        printf(" Reading ./data/global_synmap_qd.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].veg_type = (short)fdat[i];
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         */
        
        switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } /* */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /* Urban greening: 2022/01/14 */
        /* if(EX_URBAN_FOREST == 1 && grid[i].veg_type == 13){
            grid[i].veg_type = 2;
        } */
        
        /* ****************************************** for control: 2009/04/05 *****/
        if(grid[i].veg_type != 0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    if((fp_dat = fopen("./data/global_soiltaxonomy_order_qd.flt","rb")) == NULL){
        printf("!! NO global_soiltaxonomy_order_qd.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    for(i=0; i<WGRIDS; i++){
        grid[i].soil_type = (short)fdat[i];
    }
    if(NOTICE == 1){
        printf(" Reading ./data/global_soiltaxonomy_order_qd.flt...");
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* clay fraction *********************/
    if((fp_dat = fopen("./data/global_hwsd_soilmap_tclay_qd.flt","rb")) == NULL){
        printf("!! NO global_hwsd_soilmap_tclay_qd.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading global_hwsd_soilmap_tclay_qd.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
        if(grid[i].clay_frac < 0.0){
            grid[i].field_cap = 8; /* no soil data */
        }
       if(grid[i].clay_frac < 0.01){
            grid[i].clay_frac = 0.011;
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* sand fraction *********************/
    if((fp_dat = fopen("./data/global_hwsd_soilmap_tsand_qd.flt","rb")) == NULL){
        printf("!! NO global_hwsd_soilmap_tsand_qd.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading global_hwsd_soilmap_tsand_qd.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
        if(grid[i].sand_frac < 0.0){
            grid[i].field_cap = 8; /* no soil data */
        }
        if(grid[i].sand_frac < 0.01){
            grid[i].sand_frac = 0.01;
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n  ");
    }
    
    for(i=0; i<WGRIDS; i++){
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }

        if(grid[i].veg_type < 0 || grid[i].veg_type > 15){
            printf("%ld %d %d\n", i,grid[i].veg_type,grid[i].flag_datavl);
        }
    }
    
    /* *******************************************************/
    for(i=0; i<WGRIDS; i++){
        /* data available */
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no vegetated land */
            //printf("%ld %ld %ld\n", i,grid[i].veg_type,grid[i].flag_datavl);
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        
        if(NOTICE == 1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%10==0 && ((i/NCOL)%1==0)) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0 && ((i/NCOL)%1==0)){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
       if(grid[i].field_cap <= 1.0){
            grid[i].field_cap = 1.0;
        }
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;

    }
}

/* initialize BAMIYAN  conditions ****************************/
/* * initialize climate and soil parameters by using the prepared data files **/
void f_init_bamiyan_run(
    struct Grid grid[]
){
    long i;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    /*
    if((fp_dat = fopen("./data/bamiyan_lat_3sec.flt","rb")) == NULL){
        printf("!! NO bamiyan_lat_3sec.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bamiyan_lat_3sec.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    } */
    
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = 35.00351986 - 1.0/1200.0 * (double)i/(double)NCOL;
    }
    
    /* longitude */
    /* if((fp_dat = fopen("./data/bamiyan_lon_3sec.flt","rb")) == NULL){
        printf("!! NO bamiyan_lon_3sec.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bamiyan_lon_3sec.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    } */
    
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = 67.33698109 + 1.0/1200.0 * (double)(i%NCOL);
    }
    
    /* grid area */
    for(i=0; i<WGRIDS; i++){
        garea = grid_area(grid[i].lat+0.5/1200.0, grid[i].lat-0.5/1200.0,  0.0, 1.0/1200.0);
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/bamiyan_alt_3sec.flt","rb")) == NULL){
        printf("!! NO bamiyan_alt_3sec.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bamiyan_alt_3sec.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
        //grid[i].topo = 100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* land cover */
    if((fp_dat = fopen("./data/bamiyan_landcover_3sec.flt","rb")) == NULL){
        printf("!! NO bamiyan_landcover_3sec.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE == 1){
        printf(" Reading ./data/bamiyan_landcover_3sec.flt...");
    }
    for(i=0; i<WGRIDS; i++){

        grid[i].veg_type = (short)fdat[i];
        
        /* TEST: 2017/08/09 by A.Ito */
        grid[i].veg_type = 10;
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         
         */
                
        /* switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /******************************************* for control: 2009/04/05 *****/
        if(grid[i].veg_type!=0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    for(i=0; i<WGRIDS; i++){
        /* grid[i].soil_usda = (short)cdat[i]; */
        
        grid[i].soil_usda = 80;
        
        switch(grid[i].soil_usda){
            case 0:
                grid[i].soil_type = 0; /* ocean */
                break;
            case 1: case 2: case 3:
                grid[i].soil_type = 1; /* bare ground */
                break;
            case 5: case 6: case 7:
                grid[i].soil_type = 2;
                /* Gelisols (Soil Taxonomy) */
                break;
            case 10: case 11: case 12: case 13:
                grid[i].soil_type = 3;
                /* Histosols (Soil Taxonomy) */
                break;
            case 15: case 16: case 17: case 18: case 19:
                grid[i].soil_type = 4;
                /* Spodosols */
                break;
            case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27:
                grid[i].soil_type = 5;
                /* Andisols */
                break;
            case 30: case 31: case 32: case 33: case 34:
                grid[i].soil_type = 6;
                /* Oxisols */
                break;
            case 40: case 41: case 42: case 43: case 44: case 45:
                grid[i].soil_type = 7;
                /* Vertisols */
                break;
            case 50: case 51: case 52: case 53: case 54: case 55: case 56:
                grid[i].soil_type = 8;
                /* Aridisols */
                break;
            case 60: case 61: case 62: case 63: case 64:
                grid[i].soil_type = 9;
                /* Ultisols */
                break;
            case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77:
                grid[i].soil_type = 10;
                /* Mollisols */
                break;
            case 80: case 81: case 82: case 83: case 84:
                grid[i].soil_type = 11;
                /* Alfisols */
                break;
            case 90: case 91: case 92: case 93: case 94: case 85: case 86:
                grid[i].soil_type = 12;
                /* Inceptisols */
                break;
            case 95: case 96: case 97: case 98: case 99:
                grid[i].soil_type = 13;
                /* Entisols */
                break;
            default:
                grid[i].soil_type = 1;
                break;
        }
    }
    
    /* fclose(fp_dat);
     if(NOTICE==1){
     printf("done\n");
     } */
    
    /* clay fraction */
    /* if((fp_dat = fopen("./data/bamiyan_clay_3sec.flt","rb")) == NULL){
        printf("!! NO ./data/bamiyan_clay_3sec.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bamiyan_clay_3sec.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    } */

    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = 0.23;
    }
    
    /* sand fraction */
    /* if((fp_dat = fopen("./data/bamiyan_sand_3sec.flt","rb")) == NULL){
        printf("!! NO ./data/bamiyan_sand_3sec.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bamiyan_sand_3sec.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n  ");
    } */

    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = 0.43;
    }
    
    /********************************************************/
    for(i=0; i<WGRIDS; i++){
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no land */
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        if(grid[i].fieldcap < 0.0){
            grid[i].flag_datavl = 8; /* no soil data */
        }
        
        if(NOTICE==1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%15==0) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;
    }
}

/* initialize East Asia conditions ****************************/
/** initialize climate and soil parameters by using the prepared data files **/
void f_init_easia_run(
    struct Grid grid[]
){
    long i, j;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    if((fp_dat = fopen("./data/lat_easia_qd.flt","rb")) == NULL){
        printf("!! NO cfsr_lat_t382.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/lat_easia_qd.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* longitude */
    if((fp_dat = fopen("./data/lon_easia_qd.flt","rb")) == NULL){
        printf("!! NO cfsr_lon_t382.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/lon_easia_qd.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* grid area */
    for(i=0; i<NROW; i++){
    
        garea = grid_area(55.0 - ((float)i)/4.0,   55.0 - ((float)(i+1))/4.0,  0.0, 0.25);
        
        for(j=0;j<NCOL;j++){
            grid[i*NCOL+j].area = garea;
        }
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/etopo1_easia_qd.flt","rb")) == NULL){
        printf("!! NO etopo1_easia_qd.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/etopo1_easia_qd.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* land cover */
    if((fp_dat = fopen("./data/synmaps_easia_qd.flt","rb")) == NULL){
        printf("!! NO synmaps_easia_qd.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/synmaps_easia_qd.flt...");
    }
    for(i=0; i<WGRIDS; i++){

        grid[i].veg_type = (short)fdat[i];
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         
         */
        
        /* switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /******************************************* for control: 2009/04/05 *****/
        if(grid[i].veg_type!=0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    for(i=0; i<WGRIDS; i++){
        /* grid[i].soil_usda = (short)cdat[i]; */
        
        grid[i].soil_usda = 80;
        
        switch(grid[i].soil_usda){
            case 0:
                grid[i].soil_type = 0; /* ocean */
                break;
            case 1: case 2: case 3:
                grid[i].soil_type = 1; /* bare ground */
                break;
            case 5: case 6: case 7:
                grid[i].soil_type = 2;
                /* Gelisols (Soil Taxonomy) */
                break;
            case 10: case 11: case 12: case 13:
                grid[i].soil_type = 3;
                /* Histosols (Soil Taxonomy) */
                break;
            case 15: case 16: case 17: case 18: case 19:
                grid[i].soil_type = 4;
                /* Spodosols */
                break;
            case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27:
                grid[i].soil_type = 5;
                /* Andisols */
                break;
            case 30: case 31: case 32: case 33: case 34:
                grid[i].soil_type = 6;
                /* Oxisols */
                break;
            case 40: case 41: case 42: case 43: case 44: case 45:
                grid[i].soil_type = 7;
                /* Vertisols */
                break;
            case 50: case 51: case 52: case 53: case 54: case 55: case 56:
                grid[i].soil_type = 8;
                /* Aridisols */
                break;
            case 60: case 61: case 62: case 63: case 64:
                grid[i].soil_type = 9;
                /* Ultisols */
                break;
            case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77:
                grid[i].soil_type = 10;
                /* Mollisols */
                break;
            case 80: case 81: case 82: case 83: case 84:
                grid[i].soil_type = 11;
                /* Alfisols */
                break;
            case 90: case 91: case 92: case 93: case 94: case 85: case 86:
                grid[i].soil_type = 12;
                /* Inceptisols */
                break;
            case 95: case 96: case 97: case 98: case 99:
                grid[i].soil_type = 13;
                /* Entisols */
                break;
            default:
                grid[i].soil_type = 1;
                break;
        }
    }
    
    /* fclose(fp_dat);
     if(NOTICE==1){
     printf("done\n");
     } */
    
    /* clay fraction */
    if((fp_dat = fopen("./data/soilmap_easia_025d_tclay.flt","rb")) == NULL){
        printf("!! NO ./data/soilmap_easia_025d_tclay.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/soilmap_easia_025d_tclay.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* sand fraction */
    if((fp_dat = fopen("./data/soilmap_easia_025d_tsand.flt","rb")) == NULL){
        printf("!! NO ./data/soilmap_easia_025d_tsand.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/soilmap_easia_025d_tsand.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n  ");
    }
    
    /********************************************************/
    for(i=0; i<WGRIDS; i++){
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no land */
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        
        if(NOTICE==1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%18==0) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;

        if(grid[i].fieldcap < 0.0){
            grid[i].flag_datavl = 8; /* no soil data */
        }
    }
}

/* initialize Japan conditions ****************************/
/** initialize climate and soil parameters by using the prepared data files **/
void f_init_japan_run(
    struct Grid grid[]
){
    long i, j;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    if((fp_dat = fopen("./data/japan_lat_hm.flt","rb")) == NULL){
        printf("!! NO japan_lat_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japan_lat_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* longitude */
    if((fp_dat = fopen("./data/japan_lon_hm.flt","rb")) == NULL){
        printf("!! NO japan_lon_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japan_lon_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* grid area */
    for(i=0; i<NROW; i++){
    
        garea = grid_area(46.0 - ((float)i)/120.0,   46.0 - ((float)(i+1))/120.0,  0.0, 1.0/120.0);
        
        for(j=0;j<NCOL;j++){
            grid[i*NCOL+j].area = garea;
        }
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/japan_srtm_hm.flt","rb")) == NULL){
        printf("!! NO japan_srtm_hm.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japan_srtm_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* land cover */
    if((fp_dat = fopen("./data/japan-ex2_synmap_hm.flt","rb")) == NULL){
        printf("!! NO japan_synmap_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japan-ex2_synmap_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){

        grid[i].veg_type = (short)fdat[i];
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         
         */
        
        switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } /* */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /* Urban greening: 2022/01/14 */
        if(EX_URBAN_FOREST == 1 && grid[i].veg_type == 13){
            grid[i].veg_type = 2;
        }
        
        /* ****************************************** for control: 2009/04/05 *****/
        if(grid[i].veg_type!=0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    for(i=0; i<WGRIDS; i++){
        /* grid[i].soil_usda = (short)cdat[i]; */
        
        grid[i].soil_usda = 80;
        
        switch(grid[i].soil_usda){
            case 0:
                grid[i].soil_type = 0; /* ocean */
                break;
            case 1: case 2: case 3:
                grid[i].soil_type = 1; /* bare ground */
                break;
            case 5: case 6: case 7:
                grid[i].soil_type = 2;
                /* Gelisols (Soil Taxonomy) */
                break;
            case 10: case 11: case 12: case 13:
                grid[i].soil_type = 3;
                /* Histosols (Soil Taxonomy) */
                break;
            case 15: case 16: case 17: case 18: case 19:
                grid[i].soil_type = 4;
                /* Spodosols */
                break;
            case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27:
                grid[i].soil_type = 5;
                /* Andisols */
                break;
            case 30: case 31: case 32: case 33: case 34:
                grid[i].soil_type = 6;
                /* Oxisols */
                break;
            case 40: case 41: case 42: case 43: case 44: case 45:
                grid[i].soil_type = 7;
                /* Vertisols */
                break;
            case 50: case 51: case 52: case 53: case 54: case 55: case 56:
                grid[i].soil_type = 8;
                /* Aridisols */
                break;
            case 60: case 61: case 62: case 63: case 64:
                grid[i].soil_type = 9;
                /* Ultisols */
                break;
            case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77:
                grid[i].soil_type = 10;
                /* Mollisols */
                break;
            case 80: case 81: case 82: case 83: case 84:
                grid[i].soil_type = 11;
                /* Alfisols */
                break;
            case 90: case 91: case 92: case 93: case 94: case 85: case 86:
                grid[i].soil_type = 12;
                /* Inceptisols */
                break;
            case 95: case 96: case 97: case 98: case 99:
                grid[i].soil_type = 13;
                /* Entisols */
                break;
            default:
                grid[i].soil_type = 1;
                break;
        }
    }
    
    /* fclose(fp_dat);
     if(NOTICE==1){
     printf("done\n");
     } */
    
    /* clay fraction */
    if((fp_dat = fopen("./data/japan_hwsd_soilmap_tclay_hm.flt","rb")) == NULL){
        printf("!! NO ./data/japan_hwsd_soilmap_tclay_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japan_hwsd_soilmap_tclay_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* sand fraction */
    if((fp_dat = fopen("./data/japan_hwsd_soilmap_tsand_hm.flt","rb")) == NULL){
        printf("!! NO ./data/japan_hwsd_soilmap_tsand_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japan_hwsd_soilmap_tsand_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n  ");
    }
    
    /********************************************************/
    for(i=0; i<WGRIDS; i++){
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no land */
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        if(grid[i].clay_frac < 0.0 || grid[i].sand_frac < 0.0){
            grid[i].flag_datavl = 8; /* no soil */
        }
        
        if(NOTICE == 1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%48==0) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;

        if(grid[i].fieldcap <= 0.0){
            grid[i].flag_datavl = 8; /* no soil data */
        }
    }
}

/* initialize Bhutan-Bangradesh conditions ****************************/
/* * initialize climate and soil parameters by using the prepared data files **/
void f_init_bb_run(
    struct Grid grid[]
){
    long i, j;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    if((fp_dat = fopen("./data/bb_lat_hm.flt","rb")) == NULL){
        printf("!! NO bb_lat_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bb_lat_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* longitude */
    if((fp_dat = fopen("./data/bb_lon_hm.flt","rb")) == NULL){
        printf("!! NO bb_lon_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bb_lon_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* grid area */
    for(i=0; i<NROW; i++){
    
        garea = grid_area(46.0 - ((float)i)/120.0,   46.0 - ((float)(i+1))/120.0,  0.0, 1.0/120.0);
        
        for(j=0;j<NCOL;j++){
            grid[i*NCOL+j].area = garea;
        }
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/bb_srtm_hm.flt","rb")) == NULL){
        printf("!! NO bb_srtm_hm.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bb_srtm_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* land cover */
    if((fp_dat = fopen("./data/bb_synmap_hm.flt","rb")) == NULL){
        printf("!! NO bb_synmap_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bb_synmap_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){

        grid[i].veg_type = (short)fdat[i];
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         
         */
        
        switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } /* */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /* Urban greening: 2022/01/14 */
        /* if(EX_URBAN_FOREST == 1 && grid[i].veg_type == 13){
            grid[i].veg_type = 2;
        } */
        
        /* ****************************************** for control: 2009/04/05 *****/
        if(grid[i].veg_type!=0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    for(i=0; i<WGRIDS; i++){
        /* grid[i].soil_usda = (short)cdat[i]; */
        
        grid[i].soil_usda = 80;
        
        switch(grid[i].soil_usda){
            case 0:
                grid[i].soil_type = 0; /* ocean */
                break;
            case 1: case 2: case 3:
                grid[i].soil_type = 1; /* bare ground */
                break;
            case 5: case 6: case 7:
                grid[i].soil_type = 2;
                /* Gelisols (Soil Taxonomy) */
                break;
            case 10: case 11: case 12: case 13:
                grid[i].soil_type = 3;
                /* Histosols (Soil Taxonomy) */
                break;
            case 15: case 16: case 17: case 18: case 19:
                grid[i].soil_type = 4;
                /* Spodosols */
                break;
            case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27:
                grid[i].soil_type = 5;
                /* Andisols */
                break;
            case 30: case 31: case 32: case 33: case 34:
                grid[i].soil_type = 6;
                /* Oxisols */
                break;
            case 40: case 41: case 42: case 43: case 44: case 45:
                grid[i].soil_type = 7;
                /* Vertisols */
                break;
            case 50: case 51: case 52: case 53: case 54: case 55: case 56:
                grid[i].soil_type = 8;
                /* Aridisols */
                break;
            case 60: case 61: case 62: case 63: case 64:
                grid[i].soil_type = 9;
                /* Ultisols */
                break;
            case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77:
                grid[i].soil_type = 10;
                /* Mollisols */
                break;
            case 80: case 81: case 82: case 83: case 84:
                grid[i].soil_type = 11;
                /* Alfisols */
                break;
            case 90: case 91: case 92: case 93: case 94: case 85: case 86:
                grid[i].soil_type = 12;
                /* Inceptisols */
                break;
            case 95: case 96: case 97: case 98: case 99:
                grid[i].soil_type = 13;
                /* Entisols */
                break;
            default:
                grid[i].soil_type = 1;
                break;
        }
    }
    
    /* fclose(fp_dat);
     if(NOTICE==1){
     printf("done\n");
     } */
    
    /* clay fraction */
    if((fp_dat = fopen("./data/bb_hwsd_soilmap_tclay_hm.flt","rb")) == NULL){
        printf("!! NO ./data/bb_hwsd_soilmap_tclay_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bb_hwsd_soilmap_tclay_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* sand fraction */
    if((fp_dat = fopen("./data/bb_hwsd_soilmap_tsand_hm.flt","rb")) == NULL){
        printf("!! NO ./data/bb_hwsd_soilmap_tsand_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/bb_hwsd_soilmap_tsand_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n  ");
    }
    
    /* *******************************************************/
    for(i=0; i<WGRIDS; i++){
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no land */
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        if(grid[i].clay_frac < 0.0 || grid[i].sand_frac < 0.0){
            grid[i].flag_datavl = 8; /* no soil */
        }
        
        if(NOTICE == 1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%48==0) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;

        if(grid[i].fieldcap <= 0.0){
            grid[i].flag_datavl = 8; /* no soil data */
        }
    }
}

/* initialize Japan center conditions ****************************/
/* * initialize climate and soil parameters by using the prepared data files **/
void f_init_japanc_run(
    struct Grid grid[]
){
    long i, j;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    if((fp_dat = fopen("./data/japanc_lat_hm.flt","rb")) == NULL){
        printf("!! NO japanc_lat_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanc_lat_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* longitude */
    if((fp_dat = fopen("./data/japanc_lon_hm.flt","rb")) == NULL){
        printf("!! NO japanc_lon_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanc_lon_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* grid area */
    for(i=0; i<NROW; i++){
    
        garea = grid_area(46.0 - ((float)i)/120.0,   46.0 - ((float)(i+1))/120.0,  0.0, 1.0/120.0);
        
        for(j=0;j<NCOL;j++){
            grid[i*NCOL+j].area = garea;
        }
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/japanc_srtm_hm.flt","rb")) == NULL){
        printf("!! NO japanc_srtm_hm.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanc_srtm_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* land cover */
    if((fp_dat = fopen("./data/japanc_synmap_hm.flt","rb")) == NULL){
        printf("!! NO japanc_synmap_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanc_synmap_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){

        grid[i].veg_type = (short)fdat[i];
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         
         */
        
        switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } /* */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /* Urban greening: 2022/01/14 */
        /* if(EX_URBAN_FOREST == 1 && grid[i].veg_type == 13){
            grid[i].veg_type = 2;
        } */
        
        /* ****************************************** for control: 2009/04/05 *****/
        if(grid[i].veg_type!=0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    for(i=0; i<WGRIDS; i++){
        /* grid[i].soil_usda = (short)cdat[i]; */
        
        grid[i].soil_usda = 80;
        
        switch(grid[i].soil_usda){
            case 0:
                grid[i].soil_type = 0; /* ocean */
                break;
            case 1: case 2: case 3:
                grid[i].soil_type = 1; /* bare ground */
                break;
            case 5: case 6: case 7:
                grid[i].soil_type = 2;
                /* Gelisols (Soil Taxonomy) */
                break;
            case 10: case 11: case 12: case 13:
                grid[i].soil_type = 3;
                /* Histosols (Soil Taxonomy) */
                break;
            case 15: case 16: case 17: case 18: case 19:
                grid[i].soil_type = 4;
                /* Spodosols */
                break;
            case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27:
                grid[i].soil_type = 5;
                /* Andisols */
                break;
            case 30: case 31: case 32: case 33: case 34:
                grid[i].soil_type = 6;
                /* Oxisols */
                break;
            case 40: case 41: case 42: case 43: case 44: case 45:
                grid[i].soil_type = 7;
                /* Vertisols */
                break;
            case 50: case 51: case 52: case 53: case 54: case 55: case 56:
                grid[i].soil_type = 8;
                /* Aridisols */
                break;
            case 60: case 61: case 62: case 63: case 64:
                grid[i].soil_type = 9;
                /* Ultisols */
                break;
            case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77:
                grid[i].soil_type = 10;
                /* Mollisols */
                break;
            case 80: case 81: case 82: case 83: case 84:
                grid[i].soil_type = 11;
                /* Alfisols */
                break;
            case 90: case 91: case 92: case 93: case 94: case 85: case 86:
                grid[i].soil_type = 12;
                /* Inceptisols */
                break;
            case 95: case 96: case 97: case 98: case 99:
                grid[i].soil_type = 13;
                /* Entisols */
                break;
            default:
                grid[i].soil_type = 1;
                break;
        }
    }
    
    /* fclose(fp_dat);
     if(NOTICE==1){
     printf("done\n");
     } */
    
    /* clay fraction */
    if((fp_dat = fopen("./data/japanc_hwsd_soilmap_tclay_hm.flt","rb")) == NULL){
        printf("!! NO ./data/japanc_hwsd_soilmap_tclay_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanc_hwsd_soilmap_tclay_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* sand fraction */
    if((fp_dat = fopen("./data/japanc_hwsd_soilmap_tsand_hm.flt","rb")) == NULL){
        printf("!! NO ./data/japanc_hwsd_soilmap_tsand_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanc_hwsd_soilmap_tsand_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n  ");
    }
    
    /* *******************************************************/
    for(i=0; i<WGRIDS; i++){
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no land */
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        if(grid[i].clay_frac < 0.0 || grid[i].sand_frac < 0.0){
            grid[i].flag_datavl = 8; /* no soil */
        }
        
        if(NOTICE == 1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%48==0) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;

        if(grid[i].fieldcap <= 0.0){
            grid[i].flag_datavl = 8; /* no soil data */
        }
    }
}

/* initialize Japan Hokkaido conditions ****************************/
/* * initialize climate and soil parameters by using the prepared data files **/
void f_init_japanh_run(
    struct Grid grid[]
){
    long i, j;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    if((fp_dat = fopen("./data/japanh_lat_hm.flt","rb")) == NULL){
        printf("!! NO japanh_lat_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanh_lat_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* longitude */
    if((fp_dat = fopen("./data/japanh_lon_hm.flt","rb")) == NULL){
        printf("!! NO japanh_lon_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanh_lon_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* grid area */
    for(i=0; i<NROW; i++){
    
        garea = grid_area(46.0 - ((float)i)/120.0,   46.0 - ((float)(i+1))/120.0,  0.0, 1.0/120.0);
        
        for(j=0;j<NCOL;j++){
            grid[i*NCOL+j].area = garea;
        }
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/japanh_srtm_hm.flt","rb")) == NULL){
        printf("!! NO japanh_srtm_hm.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanh_srtm_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* land cover */
    if((fp_dat = fopen("./data/japanh_synmap_hm.flt","rb")) == NULL){
        printf("!! NO japanh_synmap_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japanh_synmap_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){

        grid[i].veg_type = (short)fdat[i];
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         */
        
        switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } /* */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /* Urban greening: 2022/01/14 */
        /* if(EX_URBAN_FOREST == 1 && grid[i].veg_type == 13){
            grid[i].veg_type = 2;
        } */
        
        /* ****************************************** for control: 2009/04/05 *****/
        if(grid[i].veg_type!=0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    for(i=0; i<WGRIDS; i++){
        /* grid[i].soil_usda = (short)cdat[i]; */
        
        grid[i].soil_usda = 80;
        
        switch(grid[i].soil_usda){
            case 0:
                grid[i].soil_type = 0; /* ocean */
                break;
            case 1: case 2: case 3:
                grid[i].soil_type = 1; /* bare ground */
                break;
            case 5: case 6: case 7:
                grid[i].soil_type = 2;
                /* Gelisols (Soil Taxonomy) */
                break;
            case 10: case 11: case 12: case 13:
                grid[i].soil_type = 3;
                /* Histosols (Soil Taxonomy) */
                break;
            case 15: case 16: case 17: case 18: case 19:
                grid[i].soil_type = 4;
                /* Spodosols */
                break;
            case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27:
                grid[i].soil_type = 5;
                /* Andisols */
                break;
            case 30: case 31: case 32: case 33: case 34:
                grid[i].soil_type = 6;
                /* Oxisols */
                break;
            case 40: case 41: case 42: case 43: case 44: case 45:
                grid[i].soil_type = 7;
                /* Vertisols */
                break;
            case 50: case 51: case 52: case 53: case 54: case 55: case 56:
                grid[i].soil_type = 8;
                /* Aridisols */
                break;
            case 60: case 61: case 62: case 63: case 64:
                grid[i].soil_type = 9;
                /* Ultisols */
                break;
            case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77:
                grid[i].soil_type = 10;
                /* Mollisols */
                break;
            case 80: case 81: case 82: case 83: case 84:
                grid[i].soil_type = 11;
                /* Alfisols */
                break;
            case 90: case 91: case 92: case 93: case 94: case 85: case 86:
                grid[i].soil_type = 12;
                /* Inceptisols */
                break;
            case 95: case 96: case 97: case 98: case 99:
                grid[i].soil_type = 13;
                /* Entisols */
                break;
            default:
                grid[i].soil_type = 1;
                break;
        }
    }
    
    /* fclose(fp_dat);
     if(NOTICE==1){
     printf("done\n");
     } */
    
    /* clay fraction */
    if((fp_dat = fopen("./data/japanh_hwsd_soilmap_tclay_hm.flt","rb")) == NULL){
        printf("!! NO ./data/japanh_hwsd_soilmap_tclay_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE == 1){
        printf(" Reading ./data/japanh_hwsd_soilmap_tclay_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* sand fraction */
    if((fp_dat = fopen("./data/japanh_hwsd_soilmap_tsand_hm.flt","rb")) == NULL){
        printf("!! NO ./data/japanh_hwsd_soilmap_tsand_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE == 1){
        printf(" Reading ./data/japanh_hwsd_soilmap_tsand_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n  ");
    }
    
    /* *******************************************************/
    for(i=0; i<WGRIDS; i++){
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no land */
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        if(grid[i].clay_frac < 0.0 || grid[i].sand_frac < 0.0){
            grid[i].flag_datavl = 8; /* no soil */
        }
        
        if(NOTICE == 1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%48==0) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;

        if(grid[i].fieldcap <= 0.0){
            grid[i].flag_datavl = 8; /* no soil data */
        }
    }
}

/* initialize Japan center conditions ****************************/
/* * initialize climate and soil parameters by using the prepared data files **/
void f_init_japank_run(
    struct Grid grid[]
){
    long i, j;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    if((fp_dat = fopen("./data/japank_lat_hm.flt","rb")) == NULL){
        printf("!! NO japank_lat_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japank_lat_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* longitude */
    if((fp_dat = fopen("./data/japank_lon_hm.flt","rb")) == NULL){
        printf("!! NO japank_lon_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japank_lon_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* grid area */
    for(i=0; i<NROW; i++){
    
        garea = grid_area(46.0 - ((float)i)/120.0,   46.0 - ((float)(i+1))/120.0,  0.0, 1.0/120.0);
        
        for(j=0;j<NCOL;j++){
            grid[i*NCOL+j].area = garea;
        }
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/japank_srtm_hm.flt","rb")) == NULL){
        printf("!! NO japank_srtm_hm.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE == 1){
        printf(" Reading ./data/japank_srtm_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* land cover */
    if((fp_dat = fopen("./data/japank_synmap_hm.flt","rb")) == NULL){
        printf("!! NO japank_synmap_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE == 1){
        printf(" Reading ./data/japank_synmap_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){

        grid[i].veg_type = (short)fdat[i];
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         */
        
        switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } /* */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /* Urban greening: 2022/01/14 */
        /* if(EX_URBAN_FOREST == 1 && grid[i].veg_type == 13){
            grid[i].veg_type = 2;
        } */
        
        /* ****************************************** for control: 2009/04/05 *****/
        if(grid[i].veg_type != 0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    for(i=0; i<WGRIDS; i++){
        /* grid[i].soil_usda = (short)cdat[i]; */
        
        grid[i].soil_usda = 80;
        
        switch(grid[i].soil_usda){
            case 0:
                grid[i].soil_type = 0; /* ocean */
                break;
            case 1: case 2: case 3:
                grid[i].soil_type = 1; /* bare ground */
                break;
            case 5: case 6: case 7:
                grid[i].soil_type = 2;
                /* Gelisols (Soil Taxonomy) */
                break;
            case 10: case 11: case 12: case 13:
                grid[i].soil_type = 3;
                /* Histosols (Soil Taxonomy) */
                break;
            case 15: case 16: case 17: case 18: case 19:
                grid[i].soil_type = 4;
                /* Spodosols */
                break;
            case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27:
                grid[i].soil_type = 5;
                /* Andisols */
                break;
            case 30: case 31: case 32: case 33: case 34:
                grid[i].soil_type = 6;
                /* Oxisols */
                break;
            case 40: case 41: case 42: case 43: case 44: case 45:
                grid[i].soil_type = 7;
                /* Vertisols */
                break;
            case 50: case 51: case 52: case 53: case 54: case 55: case 56:
                grid[i].soil_type = 8;
                /* Aridisols */
                break;
            case 60: case 61: case 62: case 63: case 64:
                grid[i].soil_type = 9;
                /* Ultisols */
                break;
            case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77:
                grid[i].soil_type = 10;
                /* Mollisols */
                break;
            case 80: case 81: case 82: case 83: case 84:
                grid[i].soil_type = 11;
                /* Alfisols */
                break;
            case 90: case 91: case 92: case 93: case 94: case 85: case 86:
                grid[i].soil_type = 12;
                /* Inceptisols */
                break;
            case 95: case 96: case 97: case 98: case 99:
                grid[i].soil_type = 13;
                /* Entisols */
                break;
            default:
                grid[i].soil_type = 1;
                break;
        }
    }
    
    /* fclose(fp_dat);
     if(NOTICE==1){
     printf("done\n");
     } */
    
    /* clay fraction */
    if((fp_dat = fopen("./data/japank_hwsd_soilmap_tclay_hm.flt","rb")) == NULL){
        printf("!! NO ./data/japank_hwsd_soilmap_tclay_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japank_hwsd_soilmap_tclay_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* sand fraction */
    if((fp_dat = fopen("./data/japank_hwsd_soilmap_tsand_hm.flt","rb")) == NULL){
        printf("!! NO ./data/japank_hwsd_soilmap_tsand_hm.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/japank_hwsd_soilmap_tsand_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n  ");
    }
    
    /* *******************************************************/
    for(i=0; i<WGRIDS; i++){
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no land */
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        if(grid[i].clay_frac < 0.0 || grid[i].sand_frac < 0.0){
            grid[i].flag_datavl = 8; /* no soil */
        }
        
        if(NOTICE == 1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%48==0) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;

        if(grid[i].fieldcap <= 0.0){
            grid[i].flag_datavl = 8; /* no soil data */
        }
    }
}

/* initialize PAWCs conditions ****************************/
/* * initialize climate and soil parameters by using the prepared data files **/
void f_init_pawcs_run(
    struct Grid grid[]
){
    long i, j;
    double garea;
    double rdepth[18]={0.0, 1.8, 3.0, 2.0, 2.0, 2.4, 2.5, 3.1, 1.7, 2.4, 1.5, 1.0, 1.5, 1.5, 1.5, 1.0, 4.0, 1.0};
    FILE *fp_dat;
    
    /* region data ***************************************************/
    /* latitude */
    if((fp_dat = fopen("./data/lat_pawcs_01d.flt","rb")) == NULL){
        printf("!! NO lat_pawcs_01d.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/lat_pawcs_hm.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lat = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* longitude */
    if((fp_dat = fopen("./data/lon_pawcs_01d.flt","rb")) == NULL){
        printf("!! NO lon_pawcs_01d.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/lon_pawcs_01d.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].lon = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* grid area */
    for(i=0; i<NROW; i++){
    
        garea = grid_area(46.0 - ((float)i)/120.0,   46.0 - ((float)(i+1))/120.0,  0.0, 1.0/120.0);
        
        for(j=0;j<NCOL;j++){
            grid[i*NCOL+j].area = garea;
        }
    }
    
    /* elevation */
    if((fp_dat = fopen("./data/srtm_pawcs_01d.flt","rb")) == NULL){
        printf("!! NO srtm_pawcs_01d.flt file\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/srtm_pawcs_01d.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].topo = fdat[i];
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* land cover */
    if((fp_dat = fopen("./data/synmap_pawcs_01d.flt","rb")) == NULL){
        printf("!! NO synmap_pawcs_01d.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/synmap_pawcs_01d.flt...");
    }
    for(i=0; i<WGRIDS; i++){

        grid[i].veg_type = (short)fdat[i];
        
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
    
        /*
         SYNMAP
         Value    Life forms          Tree leaf type    Tree leaf longevity    % land
         0 = water
         1        Trees               Needle            Evergreen              9.8
         2        Trees               Needle            Deciduous              1.7
         3        Trees               Needle            Mixed                  0.6
         4        Trees               Broad             Evergreen              8.2
         5        Trees               Broad             Deciduous              3
         6        Trees               Broad             Mixed                  0.5
         7        Trees               Mixed             Evergreen              0.25
         8        Trees               Mixed             Deciduous              0.03
         9        Trees               Mixed             Mixed                  3.5
         10       Trees & Shrubs      Needle            Evergreen              1.5
         11       Trees & Shrubs      Needle            Deciduous              0.2
         12       Trees & Shrubs      Needle            Mixed                  0.05
         13       Trees & Shrubs      Broad             Evergreen              0.3
         14       Trees & Shrubs      Broad             Deciduous              2.5
         15       Trees & Shrubs      Broad             Mixed                  0.2
         16       Trees & Shrubs      Mixed             Evergreen              0.03
         17       Trees & Shrubs      Mixed             Deciduous              0.04
         18       Trees & Shrubs      Mixed             Mixed                  0.3
         19       Trees & Grasses     Needle            Evergreen              0.2
         20       Trees & Grasses     Needle            Deciduous              0.03
         21       Trees & Grasses     Needle            Mixed                  0.01
         22       Trees & Grasses     Broad             Evergreen              0.3
         23       Trees & Grasses     Broad             Deciduous              2.15
         24       Trees & Grasses     Broad             Mixed                  0.15
         25       Trees & Grasses     Mixed             Evergreen              0.005
         26       Trees & Grasses     Mixed             Deciduous              0.02
         27       Trees & Grasses     Mixed             Mixed                  0.2
         28       Trees & Crops       Needle            Evergreen              0.3
         29       Trees & Crops       Needle            Deciduous              0.006
         30       Trees & Crops       Needle            Mixed                  0.003
         31       Trees & Crops       Broad             Evergreen              0.7
         32       Trees & Crops       Broad             Deciduous              1.1
         33       Trees & Crops       Broad             Mixed                  0.2
         34       Trees & Crops       Mixed             Evergreen              0.01
         35       Trees & Crops       Mixed             Deciduous              0.01
         36       Trees & Crops       Mixed             Mixed                  0.4
         37       Shrubs              -                 -                      4.5
         38       Shrubs & Grasses    -                 -                      8.3
         39       Shrubs & Crops      -                 -                      0.4
         40       Shrubs & Barren     -                 -                      10.5
         41       Grasses             -                 -                      8.3
         42       Grasses & Crops     -                 -                      1.5
         43       Grasses & Barren    -                 -                      0.3
         44       Crops               -                 -                      10.7
         45       Barren              -                 -                      11.7
         46       Urban               -                 -                      0.2
         47       Snow & Ice          -                 -                      5.2
         */
        
        switch( (short)fdat[i] ){
            case 0:     grid[i].veg_type = 0;       break;
            case 1:     grid[i].veg_type = 1;       break;
            case 2:     grid[i].veg_type = 3;       break;
            case 3:     grid[i].veg_type = 5;       break;
            case 4:     grid[i].veg_type = 2;       break;
            case 5:     grid[i].veg_type = 4;       break;
            case 6:     grid[i].veg_type = 5;       break;
            case 7:     grid[i].veg_type = 5;       break;
            case 8:     grid[i].veg_type = 5;       break;
            case 9:     grid[i].veg_type = 5;       break;
            case 10:     grid[i].veg_type = 6;       break;
            case 11:     grid[i].veg_type = 6;       break;
            case 12:     grid[i].veg_type = 6;       break;
            case 13:     grid[i].veg_type = 6;       break;
            case 14:     grid[i].veg_type = 6;       break;
            case 15:     grid[i].veg_type = 6;       break;
            case 16:     grid[i].veg_type = 6;       break;
            case 17:     grid[i].veg_type = 6;       break;
            case 18:     grid[i].veg_type = 6;       break;
            case 19:     grid[i].veg_type = 10;       break;
            case 20:     grid[i].veg_type = 10;       break;
            case 21:     grid[i].veg_type = 10;       break;
            case 22:     grid[i].veg_type = 10;       break;
            case 23:     grid[i].veg_type = 10;       break;
            case 24:     grid[i].veg_type = 10;       break;
            case 25:     grid[i].veg_type = 10;       break;
            case 26:     grid[i].veg_type = 10;       break;
            case 27:     grid[i].veg_type = 10;       break;
            case 28:     grid[i].veg_type = 12;       break;
            case 29:     grid[i].veg_type = 12;       break;
            case 30:     grid[i].veg_type = 12;       break;
            case 31:     grid[i].veg_type = 12;       break;
            case 32:     grid[i].veg_type = 12;       break;
            case 33:     grid[i].veg_type = 12;       break;
            case 34:     grid[i].veg_type = 12;       break;
            case 35:     grid[i].veg_type = 12;       break;
            case 36:     grid[i].veg_type = 12;       break;
            case 37:     grid[i].veg_type = 6;       break;
            case 38:     grid[i].veg_type = 6;       break;
            case 39:     grid[i].veg_type = 6;       break;
            case 40:     grid[i].veg_type = 6;       break;
            case 41:     grid[i].veg_type = 12;       break;
            case 42:     grid[i].veg_type = 12;       break;
            case 43:     grid[i].veg_type = 12;       break;
            case 44:     grid[i].veg_type = 12;       break;
            case 45:     grid[i].veg_type = 15;       break;
            case 46:     grid[i].veg_type = 13;       break;
            case 47:     grid[i].veg_type = 15;       break;
            default:     grid[i].veg_type = 15;       break;
        } /* */
        
        if(grid[i].veg_type < 0){
            grid[i].veg_type = 0;
        }if(grid[i].veg_type >= 16){
            grid[i].veg_type = 0;
        }
        
        /* Urban greening: 2022/01/14 */
        /* if(EX_URBAN_FOREST == 1 && grid[i].veg_type == 13){
            grid[i].veg_type = 2;
        } */
        
        /* ****************************************** for control: 2009/04/05 *****/
        if(grid[i].veg_type!=0){
            /* grid[i].veg_type = 2;*/
        }
    }
    fclose(fp_dat);
    if(NOTICE == 1){
        printf("done\n");
    }
    
    /* soil type *****************************/
    for(i=0; i<WGRIDS; i++){
        /* grid[i].soil_usda = (short)cdat[i]; */
        
        grid[i].soil_usda = 80;
        
        switch(grid[i].soil_usda){
            case 0:
                grid[i].soil_type = 0; /* ocean */
                break;
            case 1: case 2: case 3:
                grid[i].soil_type = 1; /* bare ground */
                break;
            case 5: case 6: case 7:
                grid[i].soil_type = 2;
                /* Gelisols (Soil Taxonomy) */
                break;
            case 10: case 11: case 12: case 13:
                grid[i].soil_type = 3;
                /* Histosols (Soil Taxonomy) */
                break;
            case 15: case 16: case 17: case 18: case 19:
                grid[i].soil_type = 4;
                /* Spodosols */
                break;
            case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27:
                grid[i].soil_type = 5;
                /* Andisols */
                break;
            case 30: case 31: case 32: case 33: case 34:
                grid[i].soil_type = 6;
                /* Oxisols */
                break;
            case 40: case 41: case 42: case 43: case 44: case 45:
                grid[i].soil_type = 7;
                /* Vertisols */
                break;
            case 50: case 51: case 52: case 53: case 54: case 55: case 56:
                grid[i].soil_type = 8;
                /* Aridisols */
                break;
            case 60: case 61: case 62: case 63: case 64:
                grid[i].soil_type = 9;
                /* Ultisols */
                break;
            case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77:
                grid[i].soil_type = 10;
                /* Mollisols */
                break;
            case 80: case 81: case 82: case 83: case 84:
                grid[i].soil_type = 11;
                /* Alfisols */
                break;
            case 90: case 91: case 92: case 93: case 94: case 85: case 86:
                grid[i].soil_type = 12;
                /* Inceptisols */
                break;
            case 95: case 96: case 97: case 98: case 99:
                grid[i].soil_type = 13;
                /* Entisols */
                break;
            default:
                grid[i].soil_type = 1;
                break;
        }
    }
    
    /* fclose(fp_dat);
     if(NOTICE==1){
     printf("done\n");
     } */
    
    /* clay fraction */
    if((fp_dat = fopen("./data/hwsd_soilmap_pawcs_01d_tclay.flt","rb")) == NULL){
        printf("!! NO ./data/hwsd_soilmap_pawcs_01d_tclay.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/hwsd_soilmap_pawcs_01d_tclay.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].clay_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n");
    }
    
    /* sand fraction */
    if((fp_dat = fopen("./data/hwsd_soilmap_pawcs_01d_tsand.flt","rb")) == NULL){
        printf("!! NO ./data/hwsd_soilmap_pawcs_01d_tsand.flt data\n");
        exit (1);
    }
    fread(fdat, sizeof(float), WGRIDS, fp_dat);
    if(NOTICE==1){
        printf(" Reading ./data/hwsd_soilmap_pawcs_01d_tsand.flt...");
    }
    for(i=0; i<WGRIDS; i++){
        grid[i].sand_frac = fdat[i]/100.0;
    }
    fclose(fp_dat);
    if(NOTICE==1){
        printf("done\n  ");
    }
    
    /* *******************************************************/
    for(i=0; i<WGRIDS; i++){
        grid[i].flag_datavl = 1;
        
        if(grid[i].veg_type <= 0 || grid[i].veg_type >= 15){
            grid[i].flag_datavl = 6; /* no land */
        }
        if(grid[i].topo < 0.0){
            grid[i].flag_datavl = 7; /* no land */
        }
        if(grid[i].clay_frac < 0.0 || grid[i].sand_frac < 0.0){
            grid[i].flag_datavl = 8; /* no soil */
        }
        
        if(NOTICE == 1){
            /* printf("%2d", grid[i].flag_datavl); */
            /* if((i+1)%30==0) printf("%.2f ", grid[i].lat); */
            if((i+1)%48==0) printf("%2d", grid[i].veg_type);
            if((i+1)%NCOL==0){
                printf("\n  ");
            }
        }
        
        if(grid[i].clay_frac <= 0.15){
            grid[i].stexture = 0;
        }else if(grid[i].clay_frac > 0.15 && grid[i].clay_frac <= 0.55){
            grid[i].stexture = 1;
        }else if(grid[i].clay_frac > 0.55){
            grid[i].stexture = 2;
        }
        
        /* soil parameterizatrions using Saxton (1986) *********/
        f_soil_saxton(&grid[i]);
        
        grid[i].fieldcap30 = grid[i].field_cap; /* top 30cm */
        //grid[i].fieldcap = grid[i].field_cap * 3.333; /* 100cm */
        grid[i].fieldcap = grid[i].field_cap * (rdepth[grid[i].veg_type] - 0.3)/0.3; /* 30-rooting_depth */
        grid[i].hyd_cond = 0.004313;

        if(grid[i].fieldcap <= 0.0){
            grid[i].flag_datavl = 8; /* no soil data */
        }
    }
}


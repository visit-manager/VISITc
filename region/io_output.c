/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2012 / 12 / 12 by Akihiko Ito
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/*******************************************************/
void f_create_dfile_global(
    short phase,
    short year, 
    struct Grid *grid,
    char filename[128],
    FILE *fp_out[N_OFILE]
){
    char num[8];
    
    snprintf(num, 8, "%04d", year);

    if(phase==1){
        strcpy(filename, grid->area_id);
        strcat(filename, "_su_V01_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[1] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        strcat(filename, "_su_V02_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[2] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        strcat(filename, "_su_V03_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[3] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        strcat(filename, "_su_V04_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[4] = fopen(filename,"wb");

    }else if(phase==2){
        strcpy(filename, grid->area_id);
        //strcat(filename, "_gddt_");
        strcat(filename, "_V01_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[1] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        //strcat(filename, "_lait_");
        strcat(filename, "_V02_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[2] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        //strcat(filename, "_ssnt_");
        strcat(filename, "_V03_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[3] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        //strcat(filename, "_tairav_");
        strcat(filename, "_V04_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[4] = fopen(filename,"wb");

    }
}

void f_create_dfile_bamiyan(
    short phase,
    short year, 
    struct Grid *grid,
    char filename[128],
    FILE *fp_out[N_OFILE]
){
    char num[8];
    
    snprintf(num, 8, "%04d", year);

    if(phase==1){
        strcpy(filename, grid->area_id);
        strcat(filename, "_su_V01_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[1] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        strcat(filename, "_su_V02_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[2] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        strcat(filename, "_su_V03_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[3] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        strcat(filename, "_su_V04_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[4] = fopen(filename,"wb");

    }else if(phase==2){
        strcpy(filename, grid->area_id);
        //strcat(filename, "_gddt_");
        strcat(filename, "_V01_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[1] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        //strcat(filename, "_lait_");
        strcat(filename, "_V02_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[2] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        //strcat(filename, "_ssnt_");
        strcat(filename, "_V03_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[3] = fopen(filename,"wb");
        
        strcpy(filename, grid->area_id);
        //strcat(filename, "_tairav_");
        strcat(filename, "_V04_");
        strcat(filename, num);
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_out[4] = fopen(filename,"wb");

    }
}

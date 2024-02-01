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
    FILE *fp_outhr[N_OUTHR]
){
    short i;
    char num[8], syr[8];
    
    snprintf(syr, 8, "%04d", year);
    
    for(i=0;i<N_OUTHR;i++){
        snprintf(num, 8, "%02d", i);
        
        strcpy(filename, grid->area_id);
        //strcat(filename, "_gddt_");
        strcat(filename, "_V");
        strcat(filename, num);
        strcat(filename, "_");
        strcat(filename, syr);
        strcat(filename, "_");
        strcat(filename, grid->file_name);
        strcat(filename, ".flt");
        fp_outhr[i] = fopen(filename,"wb");
    }
}

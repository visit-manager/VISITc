/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in November 02, 2007						*/

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

extern float fdat[NROW*NCOL];

/* read land use data *************************************/
void f_open_landuse_mstmipgl(
    FILE *fp_lu[N_LUDATA],
    long yr_offset,
    struct Grid grid[]
){    
    long f;
    long yr_start;
    
    printf("\nReading land-use data...\n");
    
    /* Hurtt's harminized data */
    yr_start = 1700;
    
    fp_lu[0] = fopen("./data/mstmip_hurtts_gcrop_hd.flt","rb");
    fp_lu[1] = fopen("./data/mstmip_hurtts_gflcp_hd.flt","rb");
    fp_lu[2] = fopen("./data/mstmip_hurtts_gflcs_hd.flt","rb");
    fp_lu[3] = fopen("./data/mstmip_hurtts_gflcu_hd.flt","rb");
    fp_lu[4] = fopen("./data/mstmip_hurtts_gflpc_hd.flt","rb");
    fp_lu[5] = fopen("./data/mstmip_hurtts_gflps_hd.flt","rb");
    fp_lu[6] = fopen("./data/mstmip_hurtts_gflpu_hd.flt","rb");
    fp_lu[7] = fopen("./data/mstmip_hurtts_gflsc_hd.flt","rb");
    fp_lu[8] = fopen("./data/mstmip_hurtts_gflsp_hd.flt","rb");
    fp_lu[9] = fopen("./data/mstmip_hurtts_gflsu_hd.flt","rb");
    fp_lu[10] = fopen("./data/mstmip_hurtts_gfluc_hd.flt","rb");
    fp_lu[11] = fopen("./data/mstmip_hurtts_gflup_hd.flt","rb");
    fp_lu[12] = fopen("./data/mstmip_hurtts_gflus_hd.flt","rb");
    fp_lu[13] = fopen("./data/mstmip_hurtts_gflvc_hd.flt","rb");
    fp_lu[14] = fopen("./data/mstmip_hurtts_gflvp_hd.flt","rb");
    fp_lu[15] = fopen("./data/mstmip_hurtts_gflvu_hd.flt","rb");
    fp_lu[16] = fopen("./data/mstmip_hurtts_gfsh1_hd.flt","rb");
    fp_lu[17] = fopen("./data/mstmip_hurtts_gfsh2_hd.flt","rb");
    fp_lu[18] = fopen("./data/mstmip_hurtts_gfsh3_hd.flt","rb");
    fp_lu[19] = fopen("./data/mstmip_hurtts_gfvh1_hd.flt","rb");
    fp_lu[20] = fopen("./data/mstmip_hurtts_gfvh2_hd.flt","rb");
    fp_lu[21] = fopen("./data/mstmip_hurtts_gothr_hd.flt","rb");
    fp_lu[22] = fopen("./data/mstmip_hurtts_gpast_hd.flt","rb");
    fp_lu[23] = fopen("./data/mstmip_hurtts_gsbh1_hd.flt","rb");
    fp_lu[24] = fopen("./data/mstmip_hurtts_gsbh2_hd.flt","rb");
    fp_lu[25] = fopen("./data/mstmip_hurtts_gsbh3_hd.flt","rb");
    fp_lu[26] = fopen("./data/mstmip_hurtts_gsecd_hd.flt","rb");
    fp_lu[27] = fopen("./data/mstmip_hurtts_gssma_hd.flt","rb");
    fp_lu[28] = fopen("./data/mstmip_hurtts_gssmb_hd.flt","rb");
    fp_lu[29] = fopen("./data/mstmip_hurtts_gsumm_hd.flt","rb");
    fp_lu[30] = fopen("./data/mstmip_hurtts_gurbn_hd.flt","rb");
    fp_lu[31] = fopen("./data/mstmip_hurtts_gvbh1_hd.flt","rb");
    fp_lu[32] = fopen("./data/mstmip_hurtts_gvbh2_hd.flt","rb");
    fp_lu[33] = fopen("./data/mstmip_hurtts_gzdis_hd.flt","rb");
    
    /**/
    for(f=0;f<N_LUDATA;f++){
        fseek(fp_lu[f], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    }

    /**/
    fread(fdat,sizeof(float),NROW*NCOL,fp_lu[13]);
    for(f=0;f<NROW*NCOL;f++){
        grid[f].gflvc_pi = fdat[f];
    }
    fseek(fp_lu[13], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    
    fread(fdat,sizeof(float),NROW*NCOL,fp_lu[14]);
    for(f=0;f<NROW*NCOL;f++){
        grid[f].gflvp_pi = fdat[f];
    }
    fseek(fp_lu[14], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    
    fread(fdat,sizeof(float),NROW*NCOL,fp_lu[7]);
    for(f=0;f<NROW*NCOL;f++){
        grid[f].gflsc_pi = fdat[f];
    }
    fseek(fp_lu[7], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    
    fread(fdat,sizeof(float),NROW*NCOL,fp_lu[8]);
    for(f=0;f<NROW*NCOL;f++){
        grid[f].gflsp_pi = fdat[f];
    }
    fseek(fp_lu[8], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    
    for(f=0;f<NROW*NCOL;f++){
        grid[f].deforest_pi = grid[f].gflvc_pi + grid[f].gflvp_pi 
                            + grid[f].gflsc_pi + grid[f].gflsp_pi;
    }
    
    printf("done\n");
}

void f_open_landuse_mstmipna(
    FILE *fp_lu[N_LUDATA],
    long yr_offset,
    struct Grid grid[]
){    
    long f;
    long yr_start;
    
    printf("\nReading land-use data...\n");
    
    /* Hurtt's harminized data */
    yr_start = 1700;
    
    fp_lu[0] = fopen("./data/mstmip_na_hurtts_gcrop_qd.flt","rb");
    fp_lu[1] = fopen("./data/mstmip_na_hurtts_gflcp_qd.flt","rb");
    fp_lu[2] = fopen("./data/mstmip_na_hurtts_gflcs_qd.flt","rb");
    fp_lu[3] = fopen("./data/mstmip_na_hurtts_gflcu_qd.flt","rb");
    fp_lu[4] = fopen("./data/mstmip_na_hurtts_gflpc_qd.flt","rb");
    fp_lu[5] = fopen("./data/mstmip_na_hurtts_gflps_qd.flt","rb");
    fp_lu[6] = fopen("./data/mstmip_na_hurtts_gflpu_qd.flt","rb");
    fp_lu[7] = fopen("./data/mstmip_na_hurtts_gflsc_qd.flt","rb");
    fp_lu[8] = fopen("./data/mstmip_na_hurtts_gflsp_qd.flt","rb");
    fp_lu[9] = fopen("./data/mstmip_na_hurtts_gflsu_qd.flt","rb");
    fp_lu[10] = fopen("./data/mstmip_na_hurtts_gfluc_qd.flt","rb");
    fp_lu[11] = fopen("./data/mstmip_na_hurtts_gflup_qd.flt","rb");
    fp_lu[12] = fopen("./data/mstmip_na_hurtts_gflus_qd.flt","rb");
    fp_lu[13] = fopen("./data/mstmip_na_hurtts_gflvc_qd.flt","rb");
    fp_lu[14] = fopen("./data/mstmip_na_hurtts_gflvp_qd.flt","rb");
    fp_lu[15] = fopen("./data/mstmip_na_hurtts_gflvu_qd.flt","rb");
    fp_lu[16] = fopen("./data/mstmip_na_hurtts_gfsh1_qd.flt","rb");
    fp_lu[17] = fopen("./data/mstmip_na_hurtts_gfsh2_qd.flt","rb");
    fp_lu[18] = fopen("./data/mstmip_na_hurtts_gfsh3_qd.flt","rb");
    fp_lu[19] = fopen("./data/mstmip_na_hurtts_gfvh1_qd.flt","rb");
    fp_lu[20] = fopen("./data/mstmip_na_hurtts_gfvh2_qd.flt","rb");
    fp_lu[21] = fopen("./data/mstmip_na_hurtts_gothr_qd.flt","rb");
    fp_lu[22] = fopen("./data/mstmip_na_hurtts_gpast_qd.flt","rb");
    fp_lu[23] = fopen("./data/mstmip_na_hurtts_gsbh1_qd.flt","rb");
    fp_lu[24] = fopen("./data/mstmip_na_hurtts_gsbh2_qd.flt","rb");
    fp_lu[25] = fopen("./data/mstmip_na_hurtts_gsbh3_qd.flt","rb");
    fp_lu[26] = fopen("./data/mstmip_na_hurtts_gsecd_qd.flt","rb");
    fp_lu[27] = fopen("./data/mstmip_na_hurtts_gssma_qd.flt","rb");
    fp_lu[28] = fopen("./data/mstmip_na_hurtts_gssmb_qd.flt","rb");
    fp_lu[29] = fopen("./data/mstmip_na_hurtts_gsumm_qd.flt","rb");
    fp_lu[30] = fopen("./data/mstmip_na_hurtts_gurbn_qd.flt","rb");
    fp_lu[31] = fopen("./data/mstmip_na_hurtts_gvbh1_qd.flt","rb");
    fp_lu[32] = fopen("./data/mstmip_na_hurtts_gvbh2_qd.flt","rb");
    fp_lu[33] = fopen("./data/mstmip_na_hurtts_gzdis_qd.flt","rb");
    
    /**/
    for(f=0;f<N_LUDATA;f++){
        fseek(fp_lu[f], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    }

    /**/
    fread(fdat,sizeof(float),NROW*NCOL,fp_lu[13]);
    for(f=0;f<NROW*NCOL;f++){
        grid[f].gflvc_pi = fdat[f];
    }
    fseek(fp_lu[13], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    
    fread(fdat,sizeof(float),NROW*NCOL,fp_lu[14]);
    for(f=0;f<NROW*NCOL;f++){
        grid[f].gflvp_pi = fdat[f];
    }
    fseek(fp_lu[14], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    
    fread(fdat,sizeof(float),NROW*NCOL,fp_lu[7]);
    for(f=0;f<NROW*NCOL;f++){
        grid[f].gflsc_pi = fdat[f];
    }
    fseek(fp_lu[7], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    
    fread(fdat,sizeof(float),NROW*NCOL,fp_lu[8]);
    for(f=0;f<NROW*NCOL;f++){
        grid[f].gflsp_pi = fdat[f];
    }
    fseek(fp_lu[8], sizeof(float)*NROW*NCOL*(yr_offset - yr_start),SEEK_SET);
    
    for(f=0;f<NROW*NCOL;f++){
        grid[f].deforest_pi = grid[f].gflvc_pi + grid[f].gflvp_pi 
                            + grid[f].gflsc_pi + grid[f].gflsp_pi;
    }
    
    printf("done\n");
}

/* close *************************************************/

/* read land use data *************************************/
void f_read_landuse_mstmipgl(
    short yr,
    struct Grid grid[],
    FILE *fp_lu[N_LUDATA]
){
    long f,g;
    
    for(f=0;f<N_LUDATA;f++){
        fseek(fp_lu[f], sizeof(float)*NROW*NCOL*(yr - 1700),SEEK_SET);
    
        fread(fdat, sizeof(float), NROW*NCOL, fp_lu[f]);
        
        for(g=0;g<NROW*NCOL;g++){
            switch(f){
                case 0: grid[g].gcrop = fdat[g]; break;
                case 1: grid[g].gflcp = fdat[g]; break;
                case 2: grid[g].gflcs = fdat[g]; break;
                case 3: grid[g].gflcu = fdat[g]; break;
                case 4: grid[g].gflpc = fdat[g]; break;
                case 5: grid[g].gflps = fdat[g]; break;
                case 6: grid[g].gflpu = fdat[g]; break;
                case 7: grid[g].gflsc = fdat[g]; break;
                case 8: grid[g].gflsp = fdat[g]; break;
                case 9: grid[g].gflsu = fdat[g]; break;
                case 10: grid[g].gfluc = fdat[g]; break;
                case 11: grid[g].gflup = fdat[g]; break;
                case 12: grid[g].gflus = fdat[g]; break;
                case 13: grid[g].gflvc = fdat[g]; break;
                case 14: grid[g].gflvp = fdat[g]; break;
                case 15: grid[g].gflvu = fdat[g]; break;
                case 16: grid[g].gfsh1 = fdat[g]; break;
                case 17: grid[g].gfsh2 = fdat[g]; break;
                case 18: grid[g].gfsh3 = fdat[g]; break;
                case 19: grid[g].gfvh1 = fdat[g]; break;
                case 20: grid[g].gfvh2 = fdat[g]; break;
                case 21: grid[g].gothr = fdat[g]; break;
                case 22: grid[g].gpast = fdat[g]; break;
                case 23: grid[g].gsbh1 = fdat[g]; break;
                case 24: grid[g].gsbh2 = fdat[g]; break;
                case 25: grid[g].gsbh3 = fdat[g]; break;
                case 26: grid[g].gsecd = fdat[g]; break;
                case 27: grid[g].gssma = fdat[g]; break;
                case 28: grid[g].gssmb = fdat[g]; break;
                case 29: grid[g].gsumm = fdat[g]; break;
                case 30: grid[g].gurbn = fdat[g]; break;
                case 31: grid[g].gvbh1 = fdat[g]; break;
                case 32: grid[g].gvbh2 = fdat[g]; break;
                case 33: grid[g].gzdis = fdat[g]; break;
                default:
                    printf("!! Bad access to land use data!!\n");
                    exit(1);
            }
        }
    }
}

/* land use change ************************************************************/
void f_conv_landuse(
    struct Grid *grid,
    struct Loct *loct
){
    
	if(LANDUSE == 0){
		/* Natural condition *************/
		grid->gm_crop_con = 0.0;
		grid->gm_past_con = 0.0;
	}else if(LANDUSE>=1 && LANDUSE<=5){
		/* SAGE land-use data:
		 Ramankutty, N., and J. A. Foley (1999), Estimating historical changes in global 
		 land cover: croplands from 1700 to 1992, Global Biogeochemical Cycles, 13(4), 997-1027.
		 */
		if(loct->climy <= 1990){
			/* SAGE, net land use, only cropland */
			grid->gm_crop_con = grid->gcrop;
			grid->gm_past_con = 0.0;
		}else{
			switch(LANDUSE){
				case 1:
					/* fixed 1900 */
					grid->gm_crop_con = grid->gm_crop_base; /*  */
					break;
				case 2:
					/* fixed 1990 */
					grid->gm_crop_con = grid->gm_crop_base; /*  */
					break;
				case 3:
					/* conventional scenario */
					grid->gm_crop_con = grid->gm_crop_base + grid->gm_crop_trend 
                                        *(float)(loct->climy - 1990);
					break;
				case 4:
					/* high scenario : ex e */
					if(grid->gm_crop_trend>0.0){
						grid->gm_crop_con = grid->gm_crop_base + grid->gm_crop_trend*(1.0+0.01*
                                            ((float)(loct->climy - 1990))) *(float)(loct->climy - 1990);
					}else{
						grid->gm_crop_con = grid->gm_crop_base;
					}
					break;
				case 5:
					/* IMAGE-based scenario */
					/* grid->gcrop_con = grid->gcrop_base + ((grid->fcrop3_image[grid->climy - 1990] 
                                + grid->fcrop4_image[grid->climy - 1990]) 
                                - (grid->fcrop3_image[0]+grid->fcrop4_image[0]))/100.0; */
                    grid->gm_crop_con = grid->gm_crop_base;
					break;
                default:
                    grid->gm_crop_con = 0.0;
                    break;
			}
			
			grid->gm_past_con = 0.0;
		}
	}else if(LANDUSE == 6){
		/* EOS-WEBSTER */
		/* Hurtt, G. C., S. Frolking, M. G. Fearon, B. Moore, E. Shevliakova, S. Malyshev, 
         S. W. Pacala, and R. A. Houghton. 2006. The underpinnings of land-use history: 
         three centuries of global gridded land-use transitions, wood-harvest activity, 
         and resulting secondary lands. Global Change Biology 12:1-22. */
        grid->gm_crop_con = grid->gcrop;
        grid->gm_past_con = grid->gpast;
	}else if(LANDUSE == 7){
		/* SAGE land-use data 1700-2007 (Revised: 2010/01/07):
		 Ramankutty, N., and J. A. Foley (1999), Estimating historical changes in global 
		 land cover: croplands from 1700 to 1992, Global Biogeochemical Cycles, 13(4), 997-1027.
         */
        grid->gm_crop_con = grid->gcrop;
        grid->gm_past_con = grid->gpast;
	}else if(LANDUSE == 8){
		/* Hurtt harmonized land-use change, 1700-2005 (added 2010/01/31) */
        grid->gm_crop_con = grid->gcrop;
        grid->gm_past_con = grid->gpast;
	}else{
		printf("Bad land-use setting ID\n");
		exit(1);
	}
	
	/*********************************/
	if(grid->gm_crop_con < 0.0){
		grid->gm_crop_con = 0.0;
	}
	if(grid->gm_past_con < 0.0){
		grid->gm_past_con = 0.0;
	}
	
	if((grid->gm_crop_con + grid->gm_past_con) > 0.999){
		grid->gm_crop_con = 0.999 / (grid->gm_crop_con + grid->gm_past_con) * grid->gm_crop_con;
		grid->gm_past_con = 0.999 / (grid->gm_crop_con + grid->gm_past_con) * grid->gm_past_con;
	}
    
	/*********************************/
	/* annual deforestation */
	if(loct->phase==1){
		/* spin-up */
		loct->deforest = grid->deforest_pi;
		/* 2008/08/20 corrected by A.Ito (thanks to E.Kato) */
        
	}else if(loct->phase==2){
        /* experiment */
		if(LANDUSE == 0){
			loct->deforest = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			loct->deforest = grid->gm_crop_con - grid->gm_crop_pre;
		}else if(LANDUSE==6 || LANDUSE==8){
			if(loct->climy<=1999){
				loct->deforest = grid->gflvc + grid->gflvp
                    + grid->gflsc + grid->gflsp;
				loct->deforest_v = grid->gflvc + grid->gflvp;
				loct->deforest_s = grid->gflsc + grid->gflsp;
			}else{
				loct->deforest = (grid->gm_crop_con - grid->gm_crop_pre) 
                    + (grid->gm_past_con - grid->gm_past_pre);
				loct->deforest_v = loct->deforest;
				loct->deforest_s = 0.0;
			}
		}else if(LANDUSE==7){
			loct->deforest = (grid->gm_crop_con - grid->gm_crop_pre) 
                + (grid->gm_past_con - grid->gm_past_pre);
		}
	}
	
	/* abandonment */
	if(loct->deforest < 0.0){
		loct->deforest = 0.0;
	}
	
	/* historical change in paddy field area: added by A.Ito (2011/2/28) ********/
	if(LANDUSE==7){
		if(loct->climy<=2005){
			if(grid->gm_paddy_base > 0.0 && grid->gm_crop_base > 0.0){
				grid->gm_paddy = grid->gm_paddy_base * 
                    (grid->gm_crop_con / grid->gm_crop_base);
			}
		}else if(loct->climy>=2006){
			if(grid->gm_paddy_base > 0.0 && grid->gm_crop_base > 0.0){
				grid->gm_paddy = grid->gm_paddy_base;
			}
		}
		
		if(grid->gm_paddy > 1.0){
			grid->gm_paddy = 1.0;
		}
		if(grid->gm_paddy < 0.0){
			grid->gm_paddy = 0.0;
		}
	}if(LANDUSE==8){
		if(loct->climy<=2005){
			if(grid->gm_paddy_base > 0.0 && grid->gm_crop_base > 0.0){
				grid->gm_paddy = grid->gm_paddy_base * 
                    (grid->gm_crop_con / grid->gm_crop_base);
			}
		}else if(loct->climy>=2006){
			if(grid->gm_paddy_base > 0.0 && grid->gm_crop_base > 0.0){
				grid->gm_paddy = grid->gm_paddy_base;
			}
		}
		
		if(grid->gm_paddy > 1.0){
			grid->gm_paddy = 1.0;
		}
		if(grid->gm_paddy < 0.0){
			grid->gm_paddy = 0.0;
		}
	}else if(LANDUSE==0){
		grid->gm_paddy = 0.0;
	}
}

/* emission from land use change ****************************************************/
void f_emit_landuse(
    struct Grid *grid, 
    struct Loct *loct,
    struct Mass *mass, 
    struct Flux *flux
){
	long f;
	float f_luc;
	float fe_conv;		/* fraction of conversion flux (1year) */
	float fe_10;		/* fraction of 10-year pool flux */
	float fe_100;		/* fraction of 100-year pool flux */
	float fe_detr;		/* fraction of detritus flux */
	float new_detr, new_conv, new_10, new_100;	/* added by A.Ito based on E.Kato (2009/03/30) */
	float eff_mass, f_mass_secfor;
    float pmas_fol, pmas_stm, pmas_rot;
	
	/* mass fraction of secondary forest */
	/* Nelson et al. (2000) Bioscience 50:419-431 */
	f_mass_secfor = 0.1;
	
	switch(grid->veg_type){
            /* detritus production by land-use change:
             McGuire, A. D., S. Sitch, J. S. Clein, R. Dargaville, G. Esser, J. Foley, 
             M. Heimann, F. Joos, J. Kaplan, D. W. Kicklighter, R. A. Meier, J. M. Melillo, 
             B. I. Moore, L. J. Williams, and U. Wittenberg, 2001: Carbon balance of the 
             terrestrial biosphere in the twentieth century: analysis of CO2, 
             climate and land use effects with four process-based ecosystem models. 
             Global Biogeochemical Cycles, 15, 183-206.
             */
            /* tropical forests */
		case 1: case 2:  
			fe_detr = 0.33;
			fe_conv = 0.40;
			fe_10 = 0.27;
			fe_100 = 0.0;
			break;
            /* temperate & boreal forests */
		case 3: case 4: case 5: case 6: case 7: case 8:
			fe_detr = 0.33;
			fe_conv = 0.40;
			fe_10 = 0.20;
			fe_100 = 0.07;
			break;
            /* grasslands & tundra */
		case 10: case 13: case 14: case 15:
			fe_detr = 0.50;
			fe_conv = 0.50;
			fe_10 = 0.0;
			fe_100 = 0.0;
			break;
            /* deserts & shrubs */
		case 9: case 11: case 12:
			fe_detr = 0.50;
			fe_conv = 0.40;
			fe_10 = 0.10;
			fe_100 = 0.0;
			break;
		default:
			fe_detr = 0.50;
			fe_conv = 0.40;
			fe_10 = 0.10;
			fe_100 = 0.0;
			break;
	}
    
    /* total vegetation carbon pools */
    pmas_fol = (mass->tree).fol + loct->funder_c3*(mass->c3).fol + loct->funder_c4*(mass->c4).fol;
    pmas_stm = (mass->tree).stm + loct->funder_c3*(mass->c3).stm + loct->funder_c4*(mass->c4).stm;
    pmas_rot = (mass->tree).rot + loct->funder_c3*(mass->c3).rot + loct->funder_c4*(mass->c4).rot;
	
	/* effective biomass */
	eff_mass = pmas_fol + pmas_stm + 0.8*pmas_rot;
	
    /* initialization added: A. Ito (with Hamada-san's comment) 2012/01/30 */
    f_luc = 0.0;
    
    /* LANDUSE: senstivity analysis fot estimation schemes */
    
	if(loct->phase==0){
		/* spin-up: fluxes for 1801-1900 *******************************/
		/* modified by A.Ito (2009/06/05: 2010/01/07) */
		if(LANDUSE == 0){
			f_luc = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			f_luc = grid->deforest_pi;
		}else if(LANDUSE==6){
			f_luc = (grid->gflvc_pi + grid->gflvp_pi)
                    + (grid->gflsc_pi + grid->gflsp_pi)*f_mass_secfor;
		}else if(LANDUSE==7){
			/* added 2010/01/07 (A.Ito) */
			f_luc = grid->deforest_pi;
		}
		
		/* modified by A.Ito (2009/08/19) */
        /* emission at land-use conversion */
		if(f_luc > 0.0){
            /* deforestation */
			mass->c_ludetr_1 = f_luc * 0.2*pmas_rot;
			flux->lue_conv = f_luc * eff_mass * fe_conv/(fe_conv + fe_10 + fe_100);
		}else{
			mass->c_ludetr_1 = 0.0;
			flux->lue_conv = 0.0;
		}
		
        /* emission from 10-yr pool */
		for(f=1891;f<=1900;f++){
			
			if(LANDUSE == 0){
				f_luc = 0.0;
			}else if(LANDUSE>=1 && LANDUSE<=5){
				f_luc = grid->deforest_pi;
			}else if(LANDUSE==6 || LANDUSE==8){
				f_luc = (grid->gflvc_pi + grid->gflvp_pi)
                        + (grid->gflsc_pi + grid->gflsp_pi)*f_mass_secfor;
			}else if(LANDUSE==7){
				/* added 2010/01/07 (A.Ito) */
				f_luc = grid->deforest_pi;
			}
			
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			if(f_luc > 0.0){
				new_10 = f_luc * eff_mass * fe_10/(fe_conv + fe_10 + fe_100);
				mass->c_ludetr_10 = new_10;
			}else{
                /* corrected: A. Ito (with Hamada-san's comment) 2012/01/30 */
				new_10 = 0.0;
				mass->c_ludetr_10 = 0.0;
			}
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			flux->lue_detr_10 += new_10/10.0;
		}
		
        /* emission from 100-yr pool */
		for(f=1801;f<=1900;f++){
			/* senstivity analysis */
			if(LANDUSE == 0){
				f_luc = 0.0;
			}else if(LANDUSE>=1 && LANDUSE<=5){
				f_luc = grid->deforest_pi;
			}else if(LANDUSE==6 || LANDUSE==8){
				f_luc = (grid->gflvc_pi + grid->gflvp_pi)
                    + (grid->gflsc_pi + grid->gflsp_pi)*f_mass_secfor;
			}else if(LANDUSE==7){
				/* added 2010/01/07 (A.Ito) */
				f_luc = grid->deforest_pi;
			}
			
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			if(f_luc > 0.0){
				new_100 = f_luc * eff_mass * fe_100/(fe_conv + fe_10 + fe_100);				
				mass->c_ludetr_100 = new_100;
			}else{
				new_100 = 0.0;
				mass->c_ludetr_100 = 0.0;
			}
			flux->lue_detr_100 += new_100/100.0;
		}
	}else{
		/* experiment: 1901 - 2000 - 2100 *****************************/
		
		/* emission from 10-year pool */
		/* corrected by A.Ito (2009/06/03) based on E.Kato (2008/11/21) */
		/* corrected: A.Ito and E.Kato (2009/08/16) */
        
		/* annual land use change */
		if(LANDUSE == 0){
			f_luc = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			f_luc = loct->deforest;
			/*  grid->f_crop_con - grid->f_crop_p;  */
		}else if(LANDUSE==6 || LANDUSE==8){
			/* assumption: biomass in secondary forest is half (0.5) of primary forest */
			f_luc = loct->deforest_v + loct->deforest_s * f_mass_secfor;
		}else if(LANDUSE==7){
			/* added 2010/01/07 (A.Ito) */
			f_luc = loct->deforest;
		}
		
        /* emission from 10-yr pool */
		flux->lue_detr_10 = mass->c_ludetr_10/10.0;
        mass->c_ludetr_10 -= flux->lue_detr_10;
        /* emission from 100-yr pool */
		flux->lue_detr_100 = mass->c_ludetr_100/100.0;
        mass->c_ludetr_100 -= flux->lue_detr_100;
        
		if(f_luc > 0.0){ /* deforestion occured */
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			new_detr = f_luc * 0.2*pmas_rot;
			new_conv = f_luc * eff_mass * fe_conv/(fe_conv + fe_10 + fe_100);
			new_10 = f_luc * eff_mass * fe_10/(fe_conv + fe_10 + fe_100);
			new_100 = f_luc * eff_mass * fe_100/(fe_conv + fe_10 + fe_100);
			
			/* emission from 1-yr or instantaneous pool */
			flux->lue_detr_1 = new_detr;
			flux->lue_conv = new_conv;
			
			mass->c_ludetr_10 += new_10;
			mass->c_ludetr_100 += new_100;
 		}else{
			flux->lue_conv = 0.0;
			flux->lue_detr_1 = 0.0;
		}
	}
}

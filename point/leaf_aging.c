/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

#include<stdio.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/*  leaf biochemical properties  *****************************/
/*
 Ito, A., H. Muraoka, H. Koizumi, N. Saigusa, S. Murayama, and S. Yamamoto. 2006. 
 Seasonal variation in leaf properties and ecosystem carbon budget in a cool-temperate 
 deciduous broad-leaved forest: simulation analysis at Takayama site, Japan. 
 Ecological Research 21:137-149.
*/
void f_leaf_aging(
	struct  Grid *grid, 
	struct  Loct  *loct,  
	struct  Echar *echar, 
	struct  Mass  *mass
){
	double jvrt1, rdt1, jvrt2, rdt2;
	
	/***********************************************************/		
	/* maximum carboxylation rate (Vcmax) at 25 deg C */
	if(AG_VMX == 1){
		if((echar->tree).day_flush<=0){
			(echar->overs1).vcmx_b = 0.0;
		}else if((echar->tree).day_flush>0&&(echar->tree).day_flush<15){
			(echar->overs1).vcmx_b = 0.0+3.0*((double)(echar->tree).day_flush);
		}else if((echar->tree).day_flush>=15&&(echar->tree).day_flush<140){
			(echar->overs1).vcmx_b = 45.0+0.29*((double)(echar->tree).day_flush-15.0+1.0);
		}else if((echar->tree).day_flush>=140&&(echar->tree).day_flush<180){
			(echar->overs1).vcmx_b = 81.25-2.0*((double)(echar->tree).day_flush-140.0+1.0);
			if((echar->overs1).vcmx_b<25.0) (echar->overs1).vcmx_b=25.0;
		}else if((echar->tree).day_flush>=180){
			(echar->overs1).vcmx_b = 25.0;
		}
		(echar->overs1).vcmx_b = ((echar->overs1).vcmx_b>0.0)?(echar->overs1).vcmx_b:0.0;

		/* (echar->overs1).vcmx_b = (echar->overs1).vcmx_b/(echar->overs1).vcmx_con; */
		(echar->overs1).vcmx_b = (echar->overs1).vcmx_b/81.25;
		(echar->overs1).vcmx = (echar->overs1).lai*(echar->overs1).vcmx_con*(echar->overs1).vcmx_b;
		
		if((echar->tree).day_flush<=0){
			(echar->overs2).vcmx_b = 0.0;
		}else if((echar->tree).day_flush>0&&(echar->tree).day_flush<15){
			(echar->overs2).vcmx_b = 0.0+4.5*((double)(echar->tree).day_flush);
		}else if((echar->tree).day_flush>=15&&(echar->tree).day_flush<140){
			(echar->overs2).vcmx_b = 63.0+0.12*((double)(echar->tree).day_flush-15.0+1.0);
		}else if((echar->tree).day_flush>=140&&(echar->tree).day_flush<180){
			(echar->overs2).vcmx_b = 78.0-4.5*((double)(echar->tree).day_flush-140.0+1.0);
		}else if((echar->tree).day_flush>=180){
			(echar->overs2).vcmx_b = 0.0;
		}
		(echar->overs2).vcmx_b = ((echar->overs2).vcmx_b>0.0)?(echar->overs2).vcmx_b:0.0;
		
		/* (echar->overs2).vcmx_b = (echar->overs2).vcmx_b/(echar->overs2).vcmx_con; */
		(echar->overs2).vcmx_b = (echar->overs2).vcmx_b/78.0;
		(echar->overs2).vcmx = (echar->overs2).lai*(echar->overs2).vcmx_con*(echar->overs2).vcmx_b;
	}else{
		(echar->overs1).vcmx_b = 1.0;
		(echar->overs1).vcmx = (echar->overs1).lai*(echar->overs1).vcmx_con*(echar->overs1).vcmx_b;  /* Quercus */

		(echar->overs2).vcmx_b = 1.0;
		(echar->overs2).vcmx = (echar->overs2).lai*(echar->overs2).vcmx_con*(echar->overs2).vcmx_b;  /* Betula */
	}
	(echar->unders).vcmx_b = 1.0;
	(echar->unders).vcmx = (echar->unders).lai*(echar->unders).vcmx_con*(echar->unders).vcmx_b;
	
	/***********************************************************/		
	/* Jmax/Vcmax at 25 deg C */	
	if(AG_JVR == 1){
		if((echar->tree).day_flush<=0){
			jvrt1 = 2.5;
		}else if((echar->tree).day_flush>0&&(echar->tree).day_flush<15){
			jvrt1 = 2.5-0.07*((double)(echar->tree).day_flush);
		}else if((echar->tree).day_flush>=15&&(echar->tree).day_flush<140){
			jvrt1 = 1.52+0.004*((double)(echar->tree).day_flush-15.0+1.0);
		}else if((echar->tree).day_flush>=140&&(echar->tree).day_flush<180){
			jvrt1 = 2.02+0.045*((double)(echar->tree).day_flush-140.0+1.0);
		}else if((echar->tree).day_flush>=180){
			jvrt1 = 3.865;
		}
		if((echar->tree).day_flush<=0){
			jvrt2 = 2.5;
		}else if((echar->tree).day_flush>0&&(echar->tree).day_flush<15){
			jvrt2 = 2.5-0.042*((double)(echar->tree).day_flush);
		}else if((echar->tree).day_flush>=15&&(echar->tree).day_flush<140){
			jvrt2 = 1.912+0.0016*((double)(echar->tree).day_flush-15.0+1.0);
		}else if((echar->tree).day_flush>=140&&(echar->tree).day_flush<180){
			jvrt2 = 2.112+0.045*((double)(echar->tree).day_flush-140.0+1.0);
		}else if((echar->tree).day_flush>=180){
			jvrt2 = 3.957;
		}
		(echar->overs1).jvr = jvrt1;
		(echar->overs2).jvr = jvrt2;
	}else{
		(echar->overs1).jvr = 1.772;
		(echar->overs2).jvr = 2.0128;
	}
	(echar->unders).jvr = 2.1;
	
	/***********************************************************/		
	/* dark respiration */	
	if(AG_RD == 1){
		if((echar->tree).day_flush<=0){
			rdt1 = 9.0;
		}else if((echar->tree).day_flush>0&&(echar->tree).day_flush<15){
			rdt1 = 9.0-0.49*((double)(echar->tree).day_flush);
		}else if((echar->tree).day_flush>=15&&(echar->tree).day_flush<140){
			rdt1 = 2.14-0.009*((double)(echar->tree).day_flush-15.0+1.0);
		}else if((echar->tree).day_flush>=140&&(echar->tree).day_flush<180){
			rdt1 = 1.015+0.1*((double)(echar->tree).day_flush-140.0+1.0);
		}else if((echar->tree).day_flush>=180){
			rdt1 = 5.0115;
		}
		if((echar->tree).day_flush<=0){
			rdt2 = 10.5;
		}else if((echar->tree).day_flush>0&&(echar->tree).day_flush<15){
			rdt2 = 10.5-0.58*((double)(echar->tree).day_flush);
		}else if((echar->tree).day_flush>=15&&(echar->tree).day_flush<140){
			rdt2 = 2.38-0.012*((double)(echar->tree).day_flush-15.0+1.0);
		}else if((echar->tree).day_flush>=140&&(echar->tree).day_flush<180){
			rdt2 = 0.88+0.1*((double)(echar->tree).day_flush-140.0+1.0);
		}else if((echar->tree).day_flush>=180){
			rdt2 = 4.98;
		}
		(echar->overs1).rd0 = rdt1;
		(echar->overs2).rd0 = rdt2;
	}else{
		(echar->overs1).rd0 = 3.115;
		(echar->overs2).rd0 = 2.98;
	}
	/*  (echar->unders).rd0 = 2.0;  */
	(echar->unders).rd0 = 0.3*1.198;
	/********************************************************************************************************/		
}

void f_sla_change(
    struct  Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass
){
	double lma1, lma2;
    
    if(strcmp(grid->site_id, "TKY")==0){
        if(AG_LMA==1){
            if((echar->tree).day_flush<=0 && loct->doy<210){
                lma1 = 50.0;
            }else if((echar->tree).day_flush>0&&(echar->tree).day_flush<15){
                lma1 = 50.0+2.0*((double)(echar->tree).day_flush);
            }else if((echar->tree).day_flush>=15&&(echar->tree).day_flush<180){
                lma1 = 78.0+0.24*((double)(echar->tree).day_flush-15.0+1.0);
                if((echar->tree).day_flush>=140&&((echar->tree).day_shed>=1&&(echar->tree).day_shed<=40)){
                    lma1 = 108.0-0.5*((double)(echar->tree).day_flush-140.0+1.0);
                }else if((echar->tree).day_shed>=40){
                    lma1 = 88.0;
                }
            }else{
                lma1 = 88.0;
            }
            if((echar->tree).day_flush<=0 && loct->doy<210){
                lma2 = 60.0;
            }else if((echar->tree).day_flush>0&&(echar->tree).day_flush<15){
                lma2 = 60.0+0.12*((double)(echar->tree).day_flush);
            }else if((echar->tree).day_flush>=15&&(echar->tree).day_flush<180){
                lma2 = 61.8+0.12*((double)(echar->tree).day_flush-15.0+1.0);
                if((echar->tree).day_flush>=140&&((echar->tree).day_shed>=1&&(echar->tree).day_shed<=40)){
                    lma2 = 76.68-0.5*((double)(echar->tree).day_flush-140.0+1.0);
                }else if((echar->tree).day_shed>=40){
                    lma2 = 56.68;
                }
            }else{
                lma2 = 56.68;
            }
        }else{
            lma1 = 93.12;
            lma2 = 69.24;
        }
        
        /* Quercus and Betula */
        (echar->overs1).lma = lma1;
        (echar->overs2).lma = lma2;
        (echar->overs1).lai = 10000.0/lma1*(mass->tree).fol*dmTc/100.0/2.0;
        (echar->overs2).lai = 10000.0/lma2*(mass->tree).fol*dmTc/100.0/2.0;
        (echar->tree).sla = 10000.0/(loct->comp_over1*lma1 + loct->comp_over2*lma2); 
        (echar->unders).lai = (mass->c3).lai;
    }else if(strcmp(grid->site_id, "CEAMIP_TSE")==0){
        lma1 = 93.12;
        lma2 = 69.24;
        /* Quercus and Betula */
        (echar->overs1).lma = lma1;
        (echar->overs2).lma = lma2;
        (echar->overs1).lai = 10000.0/lma1*(mass->tree).fol*dmTc/100.0/2.0;
        (echar->overs2).lai = 10000.0/lma2*(mass->tree).fol*dmTc/100.0/2.0;
        (echar->tree).sla = 10000.0/(loct->comp_over1*lma1 + loct->comp_over2*lma2); 
        (echar->unders).lai = (mass->c3).lai;
    }
}

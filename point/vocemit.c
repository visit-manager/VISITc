/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/* header files */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"structure.h"
#include"prototype.h"

/* VOC emission model by Guenther **********************************/
/*
Guenther, A., 1997. Seasonal and spatial variations in natural volatile 
organic compound emissions. Ecological Applications 7, 34-45.
*/
void f_voc_emit_guenther97(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	/* biogenic emission factors, micro g C/gdm/hour   */
	/* Latherie et al. (2006) */
	double emit_potent_isopr[18] = {0.0,
		24.0, 24.0, 16.0, 8.0, 45.0, 8.0, 8.0, 8.0, 
		24.0, 16.0, 16.0, 24.0, 16.0, 24.0, 16.0, 
		5.0, 5.0};
	double emit_potent_monotrp[18] = {0.0,
		0.8, 0.8, 1.2, 2.4, 0.8, 2.4, 2.4, 2.4, 
		1.2, 0.8, 0.8, 1.2, 0.8, 1.2, 0.8, 
		0.2, 0.2};
	double emit_potent_methanl[18] = {0.0,
		0.6, 0.6, 0.9, 1.8, 0.6, 1.8, 1.8, 1.8, 
		0.9, 0.6, 0.6, 0.9, 0.6, 0.9, 0.6, 
		2.0, 2.0};
	double emit_potent_acetone[18] = {0.0,
		0.29, 0.29, 0.43, 0.87, 0.29, 0.87, 0.87, 0.87, 
		0.43, 0.29, 0.29, 0.43, 0.29, 0.43, 0.29, 
		0.07, 0.07};
	double emit_potent_actaldhd[18] = {0.0,
		0.1, 0.1, 0.15, 0.3, 0.1, 0.3, 0.3, 0.3, 
		0.15, 0.1, 0.1, 0.15, 0.1, 0.15, 0.1, 
		0.025, 0.025};
	double emit_potent_frmardhd[18] = {0.0,
		0.07, 0.07, 0.1, 0.2, 0.07, 0.2, 0.2, 0.2, 
		0.1, 0.07, 0.07, 0.1, 0.07, 0.1, 0.07, 
		0.017, 0.017};
	double emit_potent_formacd[18] = {0.0,
		0.01, 0.01, 0.015, 0.03, 0.01, 0.03, 0.03, 0.03, 
		0.015, 0.01, 0.01, 0.015, 0.01, 0.015, 0.01, 
		0.0025, 0.0025};
	double emit_potent_acetacd[18] = {0.0,
		0.002, 0.002, 0.003, 0.006, 0.002, 0.006, 0.006, 0.006, 
		0.003, 0.002, 0.002, 0.003, 0.002, 0.003, 0.002, 
		0.0005, 0.0005};
	/* Tao and Jain (2005) **/
	double emit_potent_co[18] = {0.0,
		0.3, 0.3, 0.3, 0.3, 0.3, 0.36, 0.3, 0.3, 
		0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 
		0.3, 0.3};
	short f;
	double foliar_dens, leaf_temp, parday;
	double f_ppfd, f_temp_isopr, f_temp_monotrp, f_phenology;
	double aa, bb, cc, laiage[49], t_lai;
	
	/* veg_sage:
	0	Ocean
	1	Tropical Evergreen Forest/Woodland
	2	Tropical Deciduous Forest/Woodland
	3	Temperate Broadleaf Evergreen Forest/Woodland
	4	Temperate Needleleaf Evergreen Forest/Woodland
	5	Temperate Deciduous Forest/Woodland
	6	Boreal Evergreen Forest/Woodland
	7	Boreal Deciduous Forest/Woodland
	8	Evergreen/Deciduous Mixed Forest/Woodland
	9	Savanna
	10	Grassland/Steppe
	11	Dense Shrubland
	12	Open Shrubland
	13	Tundra
	14	Desert
	15	Polar Desert/Rock/Ice
	*/

	/* foliar density, g d.m. C / m2   */
	foliar_dens = ((mass->tree).fol + (mass->c3).fol*loct->funder_c3 + 
			(mass->c4).fol*loct->funder_c4) * 100.0 * dmTc;
	
	/* light factor */
	parday = loct->ppfd_h[loct->hour];
	f_ppfd = 0.0027 * 1.066 * parday / 
		sqrt(1.0 + 0.0027 * 0.0027 * parday * parday) * 0.5;
	
	/* temperature factor */
	leaf_temp =loct->tmp_sfc + ZAT;
	aa = exp(95000.0*(leaf_temp - 303.15)/(8.314 * leaf_temp * 303.15));
	bb = 0.961 + exp(230000.0*(leaf_temp - 314.0)/(8.314 * leaf_temp * 303.15));
	
	f_temp_isopr = aa/bb;
	
	f_temp_monotrp = exp(0.09*(leaf_temp - 303.15)); 
	
	/* other factor */
	f_phenology = 0.75;
	
	if((echar->tree).day_flush <=0){
		f_phenology = 0.1;
	}else if((echar->tree).day_flush >0 && (echar->tree).day_flush <=30){
		f_phenology = 0.3;
	}else if((echar->tree).day_flush >30 && (echar->tree).day_flush <=100){
		f_phenology = 1.1;
	}else if((echar->tree).day_flush >100 && (echar->tree).day_flush <=200){
		f_phenology = 0.4;
	}else if((echar->tree).day_flush >200){
		f_phenology = 0.1;
	}
	
/*	t_lai = 0.0;
	for(f=0;f<=48;f++){
		laiage[f] = (echar->c3).fleaf_age[f]*loct->C3ptn[grid->m] + (echar->c4).fleaf_age[f]*loct->C4ptn[grid->m];
		t_lai += laiage[f];
	}
	for(f=0;f<=48;f++){
		if(t_lai > 0.0){
			laiage[f] *= 1.0/t_lai;
		}else{
			laiage[f] = 0.0;
		}
	}
	
	switch(grid->veg_type){
		case 1: case 3: case 4: case 6: case 8: case 9: case 10: 
		case 11: case 12: case 13: case 14: case 15: 
			
			f_phenology = 0.05 * laiage[0];
			f_phenology += 0.5 * (laiage[1] + laiage[2]);
			for(f=3;f<=36;f++){
				f_phenology += 1.1 * laiage[f];
			}
			for(f=37;f<=48;f++){
				f_phenology += 0.4 * laiage[f];
			}
			
			break;
		case 2: case 5: case 7: 

			f_phenology = 0.05 * laiage[0];
			f_phenology += 0.5 * laiage[1];
			for(f=2;f<=8;f++){
				f_phenology += 1.1 * laiage[f];
			}
			for(f=9;f<=48;f++){
				f_phenology += 0.4 * laiage[f];
			}

			break;
		default:
			f_phenology = 0.5;
			break;
	}  */
	
	cc = foliar_dens * loct->daylen[loct->doy] * 10000.0 / 1000000.0;
	
	/* VOC emission, g C ha-1 day-1  */
	flux->voc_isopr_g97 = emit_potent_isopr[grid->veg_type] * cc * f_ppfd * f_temp_isopr * f_phenology;
	flux->voc_monotrp_g97 = emit_potent_monotrp[grid->veg_type] * cc * f_temp_monotrp * f_phenology;
	flux->voc_methanl_g97 = emit_potent_methanl[grid->veg_type] * cc * f_temp_monotrp * f_phenology;
	flux->voc_acetone_g97 = emit_potent_acetone[grid->veg_type] * cc * f_temp_monotrp * f_phenology;
	flux->voc_actaldhd_g97 = emit_potent_actaldhd[grid->veg_type] * cc * f_temp_monotrp * f_phenology;
	flux->voc_frmardhd_g97 = emit_potent_frmardhd[grid->veg_type] * cc *f_temp_monotrp * f_phenology;
	flux->voc_formacd_g97 = emit_potent_formacd[grid->veg_type] * cc * f_temp_monotrp * f_phenology;
	flux->voc_acetacd_g97 = emit_potent_acetacd[grid->veg_type] * cc * f_temp_monotrp * f_phenology;
	flux->voc_co_g97 = emit_potent_co[grid->veg_type] * cc * f_temp_monotrp * f_phenology;
}

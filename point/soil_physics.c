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
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* FC (field capacity) and PC (total pore space capacity) from SSSAJ86_Saxton */
/*
Saxton, K.E., Rawls, W.J., Romberger, J.S. and Papendick, R.I., 1986. 
Estimating generalized soil-water characteristics from texture. 
Soil Science Society of America Journal, 50:1031-1036.
*/
void f_soil_saxton(
	struct Grid *grid
){
	double crit_tension;
	double p_sand = grid->sand_frac*100.0;	/* percent sand fraction */
	double p_clay = grid->clay_frac*100.0;	/* percent clay fraction */
	
	switch(grid->stexture){
		/* 0: sandy */
		case 0: crit_tension = 10.0;	break;
		/* 1: medium */
		case 1: crit_tension = 33.0;	break;
		/* 2: fine */
		case 2: crit_tension = 33.0;	break;
	}
	
	/* parameters for soil water - potential equation */
	grid->a_sw = exp(-4.396 -0.0715*p_clay -4.488*0.0001*p_sand*p_sand 
					- 4.285*0.00001*p_sand*p_sand*p_clay)*100;
	grid->b_sw = -3.14 -0.00222*p_clay*p_clay -3.484*0.00001*p_sand*p_sand*p_clay;
	
	/* field capacity, mm for 30cm-soil */
	grid->field_cap = pow(crit_tension/grid->a_sw, 1.0/grid->b_sw);
	grid->field_cap *= 300.0;

	/* pore space capacity, mm for 30cm-soil */
	grid->pore_cap = 0.332 -7.251*0.0001*p_sand + 0.1276*log10(p_clay);	
	grid->pore_cap *= 300.0;
}

/* albedo of soil surface ************************/
double albedo_soil(
	struct Loct *loct,
	struct Schar *schar
){
	double albedo;
	
	/* a function of soil and snow coverage */
	albedo = schar->albedo0 + (0.7 - schar->albedo0)/(1.0 + exp(-0.05*(loct->snow_acc-70.0)));
	albedo = (albedo>0.05)?albedo:0.05; 
	albedo = (albedo<0.7)?albedo:0.7;
	
	return (albedo);
}

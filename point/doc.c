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
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/**** dissolved organic carbon (DOC) ****/
/*
Boyer, E.W., Hornberger, G.M., Bencala, K.E. and McKnight, D., 1996. 
Overview of a simple model describing variation of dissolved organic 
carbn in an upland catchment. Ecological Modelling, 86:183-188.
*/
void f_doc_boyer(
	struct Grid *grid,  
	struct Loct *loct,  
	struct Mass *mass, 
	struct Flux *flux
){
	double aa, kk;
	double stemp, doc_est;
	
	aa = 0.11;
	kk = 0.002;
	
	stemp = (loct->tmp10_soil + loct->tmp200_soil)/2.0;
	
	if(mass->sww > 0.0){
		doc_est = (mass->soil).doc + aa*pow(10.0, 0.04*stemp)
			 - (1.0 - exp(-kk*stemp))*(mass->soil).doc
			 - loct->ro2/mass->sww*(mass->soil).doc;
		if(doc_est <= 0.0){
			doc_est = 0.0;
		}
	}else{
		doc_est = 0.0;
		(flux->soil).doc_boyer = 0.0;
	}
	
	/* DOC mass, mg/L */
	(mass->soil).doc = doc_est;
	
	/* runoff DOC, mg/L -> gC/ha/day */
	(flux->soil).doc_boyer = (mass->soil).doc * loct->ro2 * 10000.0 / 1000.0;
}

/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/* [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: 2006 / 08 / 03 by Akihiko Ito: Version 0.2.0.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* from foliage, maintenance respiration **********************************/
float f_rfm(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float rfmt0, t0, rfm, ft;
	
	/* specific rate, at 15 deg C */
	rfmt0 = pchar->rmf/1000.0; 
	t0 = 15.0; /* base temperature */
	/* temperature dependence, exponential */
	//ft = exp(log(pchar->qTf)/10.0*(loct->tmp_sfc-t0));
	ft = exp(log(2.0)/10.0*(loct->tmp_sfc - t0));

	if(mass->fol >= 0.0){	
		rfm = mass->fol * rfmt0 * ft;
	}else{
		rfm = 0.0;
	}
    
    rfm *= 1.0/ (float)DSTEP;
	
	return (rfm);
}

/* from foliage, growth / construction respiration ********************/
float f_rfg(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pflx *flux
){
	float rfg;
	
	/* construction cost is proportional to biomass growth */
	rfg = pchar->rgf * flux->tpf; 
	/* flux->tpf: carbon translocation */
    
    rfg *= 1.0/ (float)DSTEP;
	
	return (rfg);
}

/* from stem and branch, maintenance respiration ***********************/
float f_rcm(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float rcmt0, t0, rcm, ft;
	
	/* specific rate, at 15 deg C */
	rcmt0 = pchar->rmc/1000.0;
	t0 = 15.0; 
	/* temperature dependence, exponential */
	ft = exp(log(pchar->qTc)/10.0*(loct->tmp_sfc - t0));
	ft = 1.0;
	
	if(mass->stm >= 0.0){	
		rcm = mass->stm * rcmt0 * ft;
	}else{
		rcm = 0.0;
	}
	
	/* rfm = rfmt0; */
    
    rcm *= 1.0/ (float)DSTEP;
	
	return (rcm);
}

/* from stem and branch, growth and construction respiration ********/
float f_rcg(
	struct Grid *grid, 
	struct Pchar *pchar,
	struct Pflx *flux
){
	float rcg;
	
	/* construction cost is proportional to biomass growth */
	rcg = pchar->rgc * flux->tpc;
    
    rcg *= 1.0/ (double)DSTEP;
	
	return (rcg);
}

/* from root, maintenance respiration ******************************/
float f_rrm(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float rrmt0, t0, rrm, ft;
	
	/* specific rate, at 15 deg C */
	rrmt0 = pchar->rmr/1000.0;
	t0 = 15.0; 

	/* temperature dependence, exponential */
	/* ft=exp(log(pchar->qTr)/10.0*(grid->tmp_sfc-t0)); */
	ft = exp(log(pchar->qTr)/10.0*(loct->tmp10_soil - t0));
	
	if(mass->rot >= 0.0){	
		rrm = mass->rot * rrmt0 * ft;
	}else{
		rrm = 0.0;
	}
	
	/* rfm = rfmt0; */
    
    rrm *= 1.0/ (double)DSTEP;
	
	return (rrm);
}

/* from root, growth and construction respiration **************/
float f_rrg(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pflx *flux
){
	float rrg;
	
	/* construction cost is proportional to biomass growth */
	rrg = pchar->rgr * flux->tpr;
    
    rrg *= 1.0/ (double)DSTEP;
	
	return (rrg);
}

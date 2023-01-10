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

/** This file gives description of functions releted to photosynthesis, or GPP **/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* formula of DAILY gross primary production **********************************/
/*
Monsi, M. and Saeki, T., 2005. On the factor light in plant communities and 
its importance for matter production. Annals of Botany, 95:549-567.
Kuroiwa, S., 1966. Dry matter production of plants. Ecology and Evolution. 
Iwanami Shoten, Tokyo, pp. 71-100.
*/
float f_gpp(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	float gpp;
	float cc1, cc2, cc3, bb;
	
	/* give daily GPP based on the Monsi-Saeki theory, or Kuroiwa's equation */
	if(pchar->psat > 0.0){
		cc1 = 2.0*pchar->psat*loct->daylen*lTs/pchar->eK; 
		bb = pchar->eK*pchar->lue*pchar->ppfd_top/pchar->psat;
		cc2 = 1.0 + sqrt(1.0 + bb);
		cc3 = 1.0 + sqrt(1.0 + bb*exp(-1.0*pchar->eK*mass->lai));
		gpp = cc1*log(cc2 / cc3);
	}else{
		gpp = 0.0;
	}
	
	/* gpp = pchar->lue; */

	return (gpp);
}

/* light-saturated photosynthetic rate ***************************************/
void f_pc_sat(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar
){
	float ftem, fstl, fnstl;
	float aa1, aa2, aa3;
	
	/* optimum temperature for photosynthesis ********************/
	if(pchar->phototype == 3){ /* C3 plants with change */
		pchar->topt = pchar->topt0 + 0.01*pchar->ci;	
	}else if(pchar->phototype == 4){ /* C4 plants without change */
		pchar->topt = pchar->topt0;
	}
	
	/** CO2 compensation point **********************************************/
	if(pchar->phototype == 3){ 
		/* C3 plants with change: Brooks & Farquhar(1985) */
		aa3 = 0.000347 * (loct->tmp_sfc - 20.0)*(loct->tmp_sfc - 20.0);
		aa1 = 1.0 + 0.0451*(loct->tmp_sfc - 20.0) + aa3;
		aa1 = (aa1 > 0.0)?aa1:0.0;
		pchar->cmpcd = pchar->cmpcd0 * aa1;	
	}else if(pchar->phototype == 4){ 
		/* C4 plants without change */
		pchar->cmpcd = pchar->cmpcd0;
	}
		
	/* temperature schalar on photosynthesis ***********************************/
	/* equation from Raich et al. (1991) Ecol. Appl (used in TEM) */
	aa1 = (loct->tmp_sfc - pchar->tmax)*(loct->tmp_sfc - pchar->tmin);
	aa2 = (loct->tmp_sfc - pchar->topt)*(loct->tmp_sfc - pchar->topt);
	ftem = aa1/(aa1 - aa2);
	ftem = (ftem<=1.0)?ftem:1.0; 
	ftem = (ftem>=0.0)?ftem:0.0;
	
	/* stomatal limitation on photosynthesis via intercellular CO2 concentration *****/
    fstl = 0.5;
	if(pchar->phototype == 3){ /* C3 plants */
		/*fstl=0.05+0.95*(pchar->ci-pchar->cmpcd)/(pchar->kmci+pchar->ci); */
		fstl = 0.30 + 0.70*(pchar->ci - pchar->cmpcd)/(pchar->kmci + pchar->ci); 
	}else if(pchar->phototype == 4){ /* C4 plants */
		fstl = 0.50 + 0.50*(pchar->ci - pchar->cmpcd)/(pchar->kmci + pchar->ci); 
	}
	fstl = (fstl<=1.0)?fstl:1.0; 
	fstl = (fstl>=0.0)?fstl:0.0;
	
	/* non-stomatal limitation (soil water effect) on photosynthesis ***************/
    fnstl = 0.5;
	if(pchar->phototype == 3){ /* C3 plants */
		fnstl = 0.95*loct->soilwtr_h/(loct->soilwtr_h + grid->fieldcap*pchar->km_nstl) + 0.05;
	}else if(pchar->phototype == 4){ /* C4 plants */
		fnstl = 0.89*loct->soilwtr_h/(loct->soilwtr_h + grid->fieldcap*pchar->km_nstl) + 0.11;
	}
	fnstl = (fnstl<=1.0)?fnstl:1.0; 
	fnstl = (fnstl>=0.0)?fnstl:0.0;
		
	/** light-saturated photosynthesis rate **/
	pchar->psat = pchar->pmax * ftem * fstl * fnstl;
	
	/* printf("%ld %.2lf %.2lf %.2lf\n", pchar->phototype, ftem, fstl, fnstl);*/
}

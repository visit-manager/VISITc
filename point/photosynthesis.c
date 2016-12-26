/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC                       */
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version in January 24, 2013                                         */

/** This file gives description of functions releted to photosynthesis, or GPP **/
#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"
#include"definition.h"

/* formula of DAILY gross primary production **********************************/
/*
Monsi, M. and Saeki, T., 2005. On the factor light in plant communities and 
  its importance for matter production. Annals of Botany, 95:549-567.
Kuroiwa, S., 1966. Dry matter production of plants. Ecology and Evolution. 
  Iwanami Shoten, Tokyo, pp. 71-100.
*/
double f_gpp(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double gpp;
	double cc1, cc2, cc3, bb;
		
	/* give daily GPP based on the Monsi-Saeki theory, or Kuroiwa's equation */
	if(pchar->psat>0.0){
		cc1 = 2.0*pchar->psat*loct->daylen[loct->doy]*lTs/pchar->eK; 
		bb = pchar->eK*pchar->lue*pchar->ppfd_t/pchar->psat;
		cc2 = 1.0 + sqrt(1.0 + bb);
		cc3 = 1.0 + sqrt(1.0 + bb*exp(-1.0*pchar->eK*mass->lai));
		gpp = cc1*log(cc2/cc3);
	}else{
		gpp = 0.0;
	}

	return (gpp);
}

/* light-saturated photosynthetic rate ***************************************/
void f_pc_sat(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar
){
	double ftem, fstl, fnstl;
	double aa1, aa2, aa3;
	
	/* optimum temperature for photosynthesis ********************/
	if(pchar->phototype==3){ /* C3 plants with change */
		pchar->topt = pchar->topt0 + 0.01*pchar->ci;	
	}else if(pchar->phototype==4){ /* C4 plants without change */
		pchar->topt = pchar->topt0;
	}
	
	/** CO2 compensation point ******************************************/
	if(pchar->phototype==3){ 
		/* C3 plants with change: Brooks & Farquhar(1985) */
		/*
		Brooks, A., and G. D. Farquhar. 1985. Effect of temperature on the 
		CO2/O2 specificity of ribulose-1,5-bisphosphate carboxylase/oxygenase 
		and the rate of respiration in the light. Planta 165:397-406.
		*/
		aa3 = 0.000347*(loct->tmp_sfc - 20.0)*(loct->tmp_sfc - 20.0);
		aa1 = 1.0 + 0.0451*(loct->tmp_sfc - 20.0)+aa3;
		aa1 = (aa1>0.0)?aa1:0.0;
		pchar->cmpcd = pchar->cmpcd0 * aa1;	
	}else if(pchar->phototype == 4){ 
		/* C4 plants without change */
		pchar->cmpcd = pchar->cmpcd0;
	}
		
	/* temperature schalar on photosynthesis ********************************/
	/* equation from Raich et al. (1991) Ecol. Appl (used in TEM) */
	aa1 = (loct->tmp_sfc - pchar->tmax)*(loct->tmp_sfc - pchar->tmin);
	aa2 = (loct->tmp_sfc - pchar->topt)*(loct->tmp_sfc - pchar->topt);
	ftem = aa1/(aa1 - aa2);
    
    /* frozen-soil constraint: 2016/01/27 by A.Ito */
    if(loct->tmp10_soil < 2.0){
        //ftem *= 0.1;
        ftem = 0.0;
    } /**/
    
	ftem = (ftem <= 1.0)?ftem:1.0; 
	ftem = (ftem >= 0.0)?ftem:0.0;
	
	/* stomatal limitation on photosynthesis via intercellular CO2 concentration *****/
	if(pchar->phototype==3){ /* C3 plants */
		/*fstl=0.05+0.95*(pchar->ci-pchar->cmpcd)/(pchar->kmci+pchar->ci); */
		fstl = 0.30 + 0.70*(pchar->ci - pchar->cmpcd)/(pchar->kmci + pchar->ci); 
	}else if(pchar->phototype==4){ /* C4 plants */
		fstl = 0.50 + 0.50*(pchar->ci - pchar->cmpcd)/(pchar->kmci + pchar->ci); 
	}
	fstl = (fstl<=1.0)?fstl:1.0; 
	fstl = (fstl>=0.0)?fstl:0.0;
	
	/* non-stomatal limitation (soil water effect) on photosynthesis ***************/
	if(pchar->phototype==3){ /* C3 plants */
		fnstl = 0.95*loct->soilwtr_h/(loct->soilwtr_h + grid->fieldcap*pchar->km_nstl) + 0.05;
	}else if(pchar->phototype==4){ /* C4 plants */
		fnstl = 0.89*loct->soilwtr_h/(loct->soilwtr_h + grid->fieldcap*pchar->km_nstl) + 0.11;
	}
	fnstl = (fnstl<=1.0)?fnstl:1.0; 
	fnstl = (fnstl>=0.0)?fnstl:0.0;
		
	/** light-saturated photosynthesis rate **/
	pchar->psat = pchar->pmax * ftem * fstl * fnstl; 
	
	/* printf("%ld %.2lf %.2lf %.2lf\n", pchar->phototype, ftem, fstl, fnstl);*/
}

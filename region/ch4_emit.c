/*	VISIT: Vegetation Integrative SImulator for Tarce gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version b 1.0.0	cerated in August 18, 2007							*/

/*  [history]
	create: 2006 / 07 / 09 by Akihiko Ito: Version 0.2.0.
	modify: GHG schemes addition by Motoko Inatomi
	modify: 2007 / 12 / 25 by Akihiko Ito: Version 0.2.0.
	ad W&H: 2008 / 05 / 31 by Akihiko Ito: Version 0.2.0.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"definition.h"
#include"setting.h"
#include"structure.h"
#include"prototype.h"

/* CH4 emission from plants */
/*
Keppler, F., Hamilton, J.T.G., Braß, M. and Röckmann, T., 2006. 
Methane emissions from terrestrial plants under aerobic conditions. 
Nature, 439:187-191.
Kirschbaum, M.U.F., Bruhn, D., Etheridge, D.M., Evans, J.R., Farquhar, G.D., 
Gifford, R.M., Paul, K.I. and Winters, A.J., 2006. 
A comment on the quantitative significance of aerobic methane releasse by plants. 
Functional Plant Biology, 33:521-530.
*/
void f_ch4emit_plant(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	float femit_sun, femit_shade;
	float sunshine;
	
	/* base emission rate */
	femit_sun = 374.0;		/* ng gdw-1 h-1 */
	femit_shade = 119.0;	/* ng gdw-1 h-1 */
	
	sunshine = loct->daylen * (1.0 - loct->tcdc_clm);
	if(sunshine<0.0){
		sunshine = 0.0;
	}
	
	/* MASS-based scaling up ***********************/
	/* Tree, g m-2 day-1 */
	if((mass->tree).season != 0){
		(flux->tree).emit_ch4_kirschbaum_mass = ((mass->tree).fol*dmTc*1000.0) * 
			(sunshine*femit_sun + (24.0 - sunshine)*femit_shade) * pow(10.0, -10.0);
	}else{
		(flux->tree).emit_ch4_kirschbaum_mass = 0.0;
	}

	/* C3, g m-2 day-1 */
	if((mass->c3).season != 0){
		(flux->c3).emit_ch4_kirschbaum_mass = ((mass->c3).fol*dmTc*1000.0) * 
			(sunshine*femit_sun + (24.0 - sunshine)*femit_shade) * pow(10.0, -10.0);
	}else{
		(flux->c3).emit_ch4_kirschbaum_mass = 0.0;
	}

	/* C4, g m-2 day-1 */
	if((mass->c4).season != 0){
		(flux->c4).emit_ch4_kirschbaum_mass = ((mass->c4).fol*dmTc*1000.0) * 
			(sunshine*femit_sun + (24.0 - sunshine)*femit_shade) * pow(10.0, -10.0);
	}else{
		(flux->c4).emit_ch4_kirschbaum_mass = 0.0;
	}
	
	/* PHOTO-based scaling up ***********************/
	if(sunshine > 0.0){
		if((flux->tree).npp>0.0){
			(flux->tree).emit_ch4_kirschbaum_photo = 2.0 * (16.0/12.0) * 
				((flux->tree).npp*1000.0) / 30000.0 * (1.0 + (24.0 - sunshine)/sunshine * 
				femit_shade / femit_sun);
		}else{
			(flux->tree).emit_ch4_kirschbaum_photo = 0.0;
		}
		if((flux->c3).npp>0.0){
			(flux->c3).emit_ch4_kirschbaum_photo = 2.0 * (16.0/12.0) * 
				((flux->c3).npp*1000.0) / 30000.0 * (1.0 + (24.0 - sunshine)/sunshine * 
				femit_shade / femit_sun);
		}else{
			(flux->c3).emit_ch4_kirschbaum_photo = 0.0;
		}
		if((flux->c4).npp>0.0){
			(flux->c4).emit_ch4_kirschbaum_photo = 2.0 * (16.0/12.0) * 
				((flux->c4).npp*1000.0) / 30000.0 * (1.0 + (24.0 - sunshine)/sunshine * 
				femit_shade / femit_sun);
		}else{
			(flux->c4).emit_ch4_kirschbaum_photo = 0.0;
		}
	}else{
		(flux->tree).emit_ch4_kirschbaum_photo = 0.0;
		(flux->c3).emit_ch4_kirschbaum_photo = 0.0;
		(flux->c4).emit_ch4_kirschbaum_photo = 0.0;
	}
}

/* Walter & Heimann CH4 emission scheme ********************************/
void f_ch4emit_walter(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Flux *flux
){
	long f, g, cumtime;
	float ww[SOILWET_LAYER+2];
	float df[SOILWET_LAYER+2], dpth[SOILWET_LAYER+2], t_mean;
	float ff[SOILWET_LAYER+2], q_ebull[SOILWET_LAYER+2], q_plant[SOILWET_LAYER+2], q_prod[SOILWET_LAYER+2], q_oxid[SOILWET_LAYER+2];
	float sdepth, wtdepth, b_thresh, tmp[SOILWET_LAYER+2], rdepth, poro;
	float f_in, f_org[SOILWET_LAYER+2], f_t[SOILWET_LAYER+2], f_grow, t_gr, t_mat;
	float t_veg, flux_ebull, flux_plant;
	float hh, rr, kk;
    float coeff_prod, q10_ch4prod;
	
	/* ***********************************************************************/
	sdepth = 2.0;		/* soil depth, m */
	kk = 0.01;			/* time step, hour */
	hh = sdepth/SOILWET_LAYER;	/* layer thickness, m */
    
    coeff_prod = 0.3;
    q10_ch4prod = 2.0;
	
	t_veg = 4.0;
	b_thresh = 500.0;	/* bubbling threshold, micro mol / liter */	
	rdepth = 0.3;		/* rooting depth, m */
	poro = 0.45*grid->sand_frac + 0.20*(1.0 - grid->sand_frac - grid->clay_frac)
			+0.14*grid->clay_frac + 0.45*0.03;		/* porosity */
    
	wtdepth = loct->water_table_depth;      /* water-table depth, m from surface */
	
	/* depth from the soil surface, m */
	dpth[0] = -0.04;
	for(f=1;f<=SOILWET_LAYER;f++){
		dpth[f] = 0.0 + ((float)f-0.5) * hh;
	}
	
	/* diffusion coefficient, m2 s-1 */
    df[0] = 0.2;
	for(f=1;f<=SOILWET_LAYER;f++){
		if(dpth[f]<wtdepth){
			df[f] = 0.2 * 0.66 * poro;  /* pow(poro, 2.3) */		
		}else{
			df[f] = 0.00002 * 0.66 * poro;  
		}
	}
	
	/* plant growth state parameters **************/
	/* W&H2000 pp.763 */
	if(grid->tsoil_annav < 5.0){
		t_gr = 2.0;
	}else{
		t_gr = 7.0;
	}
	t_mat = t_gr + 10.0;
    
    tmp[0] = loct->tmp_sfc;
	for(f=1;f<=SOILWET_LAYER;f++){
		tmp[f] = loct->tmp10_soil*(float)(SOILWET_LAYER - f)/(float)SOILWET_LAYER 
					+ (float)f/(float)SOILWET_LAYER*loct->tmp200_soil;
	}
	/* fgow: Eq. 20 */
    f_grow = 2.0;
	if(tmp[5] < t_gr){
		f_grow = 0.0;
	}else if(tmp[5] >= t_gr && tmp[5] <= t_mat){
		f_grow = 0.0 + 4.0 * (1.0 - pow((t_mat - tmp[5])/(t_mat - t_gr), 2.0));
	}else if(tmp[5] > t_mat){
		f_grow = 4.0;
	}
	
	t_mean = grid->tsoil_annav;
	
	/* organic matter factor **************************/
	f_org[0] = 1.0;
	for(f=1;f<=SOILWET_LAYER;f++){
		/* vegetated soil */
		if(dpth[f]<rdepth){
			f_org[f] = 1.0;
		}else{
			f_org[f] = exp( -fabs(dpth[f]*100.0 - rdepth*100.0) / 10.0 );
		}
	}
	
	/* soil temperature profile *******************/
    f_t[0] = 0.0;
	for(f=1;f<=SOILWET_LAYER;f++){		
		if(tmp[f] > 0.0){
			f_t[f] = 1.0;
		}else{
			f_t[f] = 0.0;
		}
	}

	/* carbon input factor: Eq.2 */
	if(grid->npp_max > 0.0){
		f_in = 1.0 + flux->npp/grid->npp_max;
	}else{
		f_in = 0.0;
	}
	
	/* *****************************************************************************/
	/* stability index, should be < 0.5 */
	rr = kk/hh/hh;
	if(rr*df[1] > 0.5){
		printf("BAD calculation condition: r = %f\n", rr);
		exit(1);
	}
	
	/* TIME *******************************************/
	if(grid->time == 1){
		cumtime = 10000;
	}else{
		cumtime = 100;
	}
    flux_ebull = flux_plant = 0.0;   
	for(g=1;g<=cumtime;g++){ 
		
		/* DEPTH **************************************/
		flux_ebull = flux_plant = 0.0;
		ff[0] = 0.0;
		for(f=1;f<=SOILWET_LAYER;f++){
			/* ebullition at high CH4 concentration */
			if( loct->prof_ch4[f] > b_thresh){
				q_ebull[f] = -1.0 * (loct->prof_ch4[f] - b_thresh);
				flux_ebull += -1000.0*hh * q_ebull[f];
			}else{
				q_ebull[f] = 0.0;
			}
			/* plant flux above rooting depth */
			if(dpth[f]<=rdepth){
				q_plant[f] = -0.01 * t_veg * 2.0*(rdepth - dpth[f])/rdepth * f_grow * loct->prof_ch4[f];
				flux_plant += -1000.0*hh * q_plant[f] * (1.0 - 0.5);
			}else{
				q_plant[f] = 0.0;
			}
			/* CH4 production and oxidation */
			if(dpth[f]>=wtdepth){
				/* below water table */
				q_prod[f] = coeff_prod * f_org[f] * f_in * f_t[f] * pow(q10_ch4prod, (tmp[f]-t_mean)/10.0);
				q_oxid[f] = 0.0;
			}else{
				/* above water table */
				q_prod[f] = 0.0;
				q_oxid[f] = -1.0*(20.0 * loct->prof_ch4[f])/(3.0 + loct->prof_ch4[f]) * pow(2.0, (tmp[f]-t_mean)/10.0);
			}
			
			/* total flux, micro M / hour */
			ff[f] = (q_oxid[f] + q_prod[f] + q_ebull[f] + q_plant[f]);
		}
		
		/* diffusion equation solved by explicit method */
		for(f=1;f<=(SOILWET_LAYER-1);f++){
			/* ww[f] = df[f]*rr*(loct->prof_ch4[f+1] + loct->prof_ch4[f-1]) + (1.0 - 2.0*df[f]*rr)*loct->prof_ch4[f] + kk*ff[f]; */
			
			ww[f] = loct->prof_ch4[f] + rr*(df[f]*(loct->prof_ch4[f+1] - 2.0*loct->prof_ch4[f] + loct->prof_ch4[f-1]) + ff[f]);
			
			if(ww[f]<=0.0){
				ww[f] = 0.0;
			}
		}
		/* boundary condition at the bottom: no downward flux */
		ww[SOILWET_LAYER] = loct->prof_ch4[SOILWET_LAYER] + rr*(df[SOILWET_LAYER]*(2.0*loct->prof_ch4[SOILWET_LAYER-1] - 2.0*loct->prof_ch4[SOILWET_LAYER]));
		if(ww[SOILWET_LAYER]<=0.0){
			ww[SOILWET_LAYER] = 0.0;
		}
					
		/* ww[LAYER] = 500.0*16.0/(8.3144*(tmp[0]+273.15))/1000000.0*apress; */
		
		/* concentration update */
		for(f=1;f<=SOILWET_LAYER;f++){
			loct->prof_ch4[f] = ww[f];
		}
	}
	
	/* flux: mg CH4 m-2 day-1 */
	(flux->soil).ch4_wh_plant =  flux_plant *(24.0/(float)DSTEP)*16.0/1000.0;
	(flux->soil).ch4_wh_ebull =  flux_ebull *(24.0/(float)DSTEP)*16.0/1000.0;
	(flux->soil).ch4_wh_diff =  df[1]/(dpth[1] - dpth[0])*(loct->prof_ch4[1] - loct->prof_ch4[0]) *(24.0/(float)DSTEP)*16.0/1000.0;
    
    (flux->soil).ch4_wh = (flux->soil).ch4_wh_plant + (flux->soil).ch4_wh_ebull + (flux->soil).ch4_wh_diff;
    //(flux->soil).ch4_wh = t_mean;
}

/* icc -o main_nc_output_hour main_nc_output_hour.c -L/usr/local/lib -lnetcdf */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>

#define NTIME 24
#define NLAT 1200
#define NLON 1200

float mread[NLAT*NLON];
float mans[NTIME*NLAT*NLON];

void handle_error(int status);
void f_doyTmody(long year, long doy, long *month, long *day);

int main(int argc, char *argv[])
{
	int status, ncid, dim_time, dim_lat, dim_lon, time_id, lat_id, lon_id, VAR_id;
	int time_dimids[1], lat_dimids[1], lon_dimids[1], VAR_dimids[3];
    float inlat[NLAT],inlon[NLON], intdata;
    double intime[NTIME];
    long year, doy, veg1, veg2;
	long f,g,h,i,j, month, day;
    short n_time[1];

	static char title_text[] = "VISITc BB map";
	static char stname_time[] = "time";			static char units_time[] = "time";
	static char stname_lat[] = "lat";			static char units_lat[] = "latitude, degree north (N->S) (float/real)";
	static char stname_lon[] = "lon";			static char units_lon[] = "longitude, degree east (W->E) (float/real)";
	static char stname_dat[] = "data";          static char units_dat[] = "data";
    
	static size_t start1[] = {0};
	static size_t count1[] = {1};
	static size_t start3[] = {0,0,0};
	static size_t count3[] = {NTIME,NLAT,NLON};
	float missing_value[] = {1.e+20f};
    char filename[128],num[8];
	FILE *fp_in;
    
    for(f=0;f<NTIME*NLAT*NLON;f++){
        mans[f] = 0.0;
    }
    
    fp_in = fopen("BB_V01_2020_20220824.flt","rb");
    for(f=0;f<24*119;f++){
        fread(mread,sizeof(float),NLAT*NLON, fp_in);
    }

    for(f=0;f<24;f++){
        fread(mread,sizeof(float),NLAT*NLON, fp_in);
        
        for(g=0;g<NLAT*NLON;g++){
            mans[f*NLAT*NLON + g] = mread[g];
        }
    }
    fclose(fp_in);
    
    /****************/
    for(f=0;f<NTIME;f++){
        intime[f] = (double)f;
    }
    for(f=0;f<NLAT;f++){
        //inlat[f] = 90.0 - ((float)f+0.5)*180.0/(double)NLAT;
        inlat[f] = 30.0 - ((float)f+0.5)/120.0;
    }
    for(f=0;f<NLON;f++){
        //inlon[f] = -180.0 + ((float)f+0.5)*360.0/(double)NLON;
        inlon[f] = 85.0 + ((float)f+0.5)/120.0;
    }
        
    printf("%s\n", nc_inq_libvers());
        
    strcpy(filename, "VISITc_BB_V01_2020_120_20220824_hour.nc");
    status = nc_create(filename, NC_CLOBBER, &ncid);
    if(status != NC_NOERR) handle_error(status); 

    /****************/
    status = nc_def_dim(ncid, "time", NTIME, &dim_time);
    if(status != NC_NOERR) handle_error(status); 
    status = nc_def_dim(ncid, "lat", NLAT, &dim_lat);
    if(status != NC_NOERR) handle_error(status); 
    status = nc_def_dim(ncid, "lon", NLON, &dim_lon);
    if(status != NC_NOERR) handle_error(status); 

    /****************/
    time_dimids[0] = dim_time;
    status = nc_def_var(ncid, "time", NC_DOUBLE, 1, time_dimids, &time_id);
    if(status != NC_NOERR) handle_error(status); 
    lat_dimids[0] = dim_lat;
    status = nc_def_var(ncid, "lat", NC_FLOAT, 1, lat_dimids, &lat_id);
    if(status != NC_NOERR) handle_error(status); 
    lon_dimids[0] = dim_lon;
    status = nc_def_var(ncid, "lon", NC_FLOAT, 1, lon_dimids, &lon_id);
    if(status != NC_NOERR) handle_error(status); 

    VAR_dimids[0] = dim_time;
    VAR_dimids[1] = dim_lat;
    VAR_dimids[2] = dim_lon;
    
    status = nc_def_var(ncid, "data", NC_FLOAT, 3, VAR_dimids, &VAR_id);
    if(status != NC_NOERR) handle_error(status); 

    /****************/
    status = nc_put_att_text(ncid, time_id, "standard_name", strlen(stname_time), stname_time);
    if(status != NC_NOERR) handle_error(status); 
    status = nc_put_att_text(ncid, time_id, "units", strlen(units_time), units_time);
    if(status != NC_NOERR) handle_error(status); 
    status = nc_put_att_text(ncid, lat_id, "standard_name", strlen(stname_lat), stname_lat);
    if(status != NC_NOERR) handle_error(status); 
    status = nc_put_att_text(ncid, lat_id, "units", strlen(units_lat), units_lat);
    if(status != NC_NOERR) handle_error(status); 
    status = nc_put_att_text(ncid, lon_id, "standard_name", strlen(stname_lon), stname_lon);
    if(status != NC_NOERR) handle_error(status); 
    status = nc_put_att_text(ncid, lon_id, "units", strlen(units_lon), units_lon);
    if(status != NC_NOERR) handle_error(status); 

    status = nc_put_att_text(ncid, VAR_id, "standard_name", strlen(stname_dat), stname_dat);
    if(status != NC_NOERR) handle_error(status); 
    status = nc_put_att_text(ncid, VAR_id, "units", strlen(units_dat), units_dat);
    if(status != NC_NOERR) handle_error(status); 

    /*****/
    status = nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(title_text), title_text);
    if(status != NC_NOERR) handle_error(status); 


    /****************/
    status = nc_enddef(ncid);
    if(status != NC_NOERR) handle_error(status); 

    /****************/
    count1[0] = NTIME;
    nc_put_vara_double(ncid, time_id, start1, count1, intime);
    if(status != NC_NOERR) handle_error(status);
    count1[0] = NLAT;
    status = nc_put_vara_float(ncid, lat_id, start1, count1, inlat);
    if(status != NC_NOERR) handle_error(status); 
    count1[0] = NLON;
    status = nc_put_vara_float(ncid, lon_id, start1, count1, inlon);
    if(status != NC_NOERR) handle_error(status); 
            
    count3[0] = NTIME;
    status = nc_put_vara_float(ncid, VAR_id, start3, count3, mans);
    if(status != NC_NOERR) handle_error(status); 
    
    /* close */
    status = nc_close(ncid);
    if(status != NC_NOERR) handle_error(status); 
	
	return 0;
}

void handle_error(int status){
	if (status != NC_NOERR) {
		fprintf(stderr, "%s\n", nc_strerror(status));
		exit(-1);
	}
}

/* DOY to Month and Day ****************************************/
void f_doyTmody(long year, long doy, long *month, long *day){
	if(year%4==0){
		if(doy>=0 && doy<31){
			*month = 0;
			*day = doy-0;
		}else if(doy>=31 && doy<60){
			*month = 1;
			*day = doy-31;
		}else if(doy>=60 && doy<91){
			*month = 2;
			*day = doy-60;
		}else if(doy>=91&&doy<121){
			*month = 3;
			*day = doy-91;
		}else if(doy>=121&&doy<152){
			*month = 4;
			*day = doy-121;
		}else if(doy>=152&&doy<182){
			*month = 5;
			*day = doy-152;
		}else if(doy>=182&&doy<213){
			*month = 6;
			*day = doy-182;
		}else if(doy>=213&&doy<244){
			*month = 7;
			*day = doy-213;
		}else if(doy>=244&&doy<274){
			*month = 8;
			*day = doy-244;
		}else if(doy>=274&&doy<305){
			*month = 9;
			*day = doy-274;
		}else if(doy>=305&doy<335){
			*month = 10;
			*day = doy-305;
		}else if(doy>=335&&doy<366){
			*month = 11;
			*day = doy-335;
		}
	}else{
		if(doy>=0 && doy<31){
			*month = 0;
			*day = doy-0;
		}else if(doy>=31 && doy<59){
			*month = 1;
			*day = doy-31;
		}else if(doy>=59 && doy<90){
			*month = 2;
			*day = doy-59;
		}else if(doy>=90&&doy<120){
			*month = 3;
			*day = doy-90;
		}else if(doy>=120&&doy<151){
			*month = 4;
			*day = doy-120;
		}else if(doy>=151&&doy<181){
			*month = 5;
			*day = doy-151;
		}else if(doy>=181&&doy<212){
			*month = 6;
			*day = doy-181;
		}else if(doy>=212&&doy<243){
			*month = 7;
			*day = doy-212;
		}else if(doy>=243&&doy<273){
			*month = 8;
			*day = doy-243;
		}else if(doy>=273&&doy<304){
			*month = 9;
			*day = doy-273;
		}else if(doy>=304&doy<334){
			*month = 10;
			*day = doy-304;
		}else if(doy>=334&&doy<365){
			*month = 11;
			*day = doy-334;
		}
	}
}



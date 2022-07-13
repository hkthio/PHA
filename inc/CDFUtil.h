/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This is an example which reads some surface pressure and
   temperatures. The data file read by this program is produced
   companion program sfc_pres_temp_wr.c. It is intended to illustrate
   the use of the netCDF C API.

   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id: sfc_pres_temp_rd.c,v 1.3 2006/06/25 11:40:00 ed Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include <math.h>

#define NDIMS 2
#define LAT_NAME "y"
#define LON_NAME "x"
#define DATA_NAME "z"
#define UNITS "units"
#define DEGREES_EAST "degrees_east"
#define DEGREES_NORTH "degrees_north"
#define RANGE_NAME "actual_range"


#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return -1;}

int read_grd_nf(char *, double **, double **, int *, int *, float ***);
int write_grd_nf(char *, double *, double *, int, int, float **);
int read_grd_hdr_cf(char *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *);
float *read_grd_data_cf( char *, float * );

#ifndef PI
#define PI      3.14159265358979323846
#define TPI     6.28318531717958647692
#define HPI     1.57079632579489661923 
#define CONV    0.017453292519943295
#endif

#ifndef Sphere
#define Sphere
#define REARTH  6371.0

void matmul();
void rotmat();

struct spherical
   {
   double phi;
   double the;
   double r;
   };
struct cartesian
   {
   double x;
   double y;
   double z;
   };
struct cubesphere 
   {
   double eta;
   double xi;
   double r;
   int reg;
   };
struct geographic 
   {
   double lon;
   double lat;
   double dep;
   };

struct Line
   {
   int np;
   struct geographic *pt;
   };

struct cubesphere sphtocub();
struct spherical cubtosph();
struct spherical cartosph();
struct cartesian sphtocart();
struct spherical geotosph();
struct geographic sphtogeo();
struct geographic newloc();
struct geographic rotgeo();
void distaz();
double distance();
void rottoequator();
void rottopole();

int csregion();
struct cartesian mult();
struct cartesian add();
struct cartesian vecprod();
double scaprod();
double norm();
int inpoly();
int inrect();

double distoline(struct geographic,struct geographic, struct geographic,struct geographic *,int *);
double distoplane(struct geographic,struct geographic,struct geographic,struct geographic,struct geographic *,int *,double *,double *);
#endif

#include "Sphere.h"

struct Station
   { 
   struct geographic lo;
   char name[8];
   char net[3];
   } ;

int StaRead();
int StaGetNetwork();
int StaGetDist();
int StaGetBox();
int StaGetPoly();
int StaMerge();

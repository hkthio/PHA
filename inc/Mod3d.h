#define ityp 139

struct model
  {
  float vp;
  float vs;
  float rh;
  float th;
  };

struct modelparams
  {
  int nlat;
  float lat0;
  float dlat;
  int nlon;
  float lon0;
  float dlon;
  int nlay;
  float lay0;
  };


struct modelparams mpar;
struct modelparams cpar;
/*struct modelparams 1Dpar = { 1,90.,181.,1,0.,361.,0,0.};*/

struct model mantle[90][180][25];
struct model crust[90][180][25];
struct model mod1D[1][1][1000];

struct model getcrust(float, float, float,struct modelparams);
struct model getvel(float, float, float,struct modelparams,struct model (*)[180][25]);
struct model getiasp(float);
float getdep(float float, int struct modelparams,struct model (*)[180][25]);

/* using Sphere.h */
struct 3dgrid
  {
  int nlat;
  int nlon;
  int npts; /* != nlat*nlon if only selected points ae available, in thas case use mask to find points */
  float dlat;
  float dlon;
  struct geographic UL;
  int nlay;
  float lay0; /*top of the first layer */
  };

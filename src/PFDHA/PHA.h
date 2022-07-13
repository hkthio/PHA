double pXceed(double, double, double, double);
float gammp(float, float);
float betai(float, float, float );


enum EquType{NORM, BETA, GAMMA};
enum EquSubType{BIL, QUAD, ELLI, BILN, QUADN, ELLIN, TAKAO, JNES, MOSS, YOUNGS};

/* structure for displacement relations */
struct FDEq
   {
   enum EquType type;
   enum EquSubType subtype;
   };


/* PHA-Disp.c*/
float SlipScaling(double, int , int , double *);
float AvSlip(float, int);
float MaxSlip(float, int);
struct FDEq *InitEq();
double probFD(float, double, double, double, double, int, double *);
double probFDDistr(float, double, double, double, double, int, double *, int);
double probSurfRup(float, int);
double probSurfRupArea(float,float,double,int,int);
double probSurfRupR(double,double,double,int, int,int);
double probSurfRupDistr(float ,float );

/* PHA-Math.c*/
double interp(double *x,double *y,double x0,int nx);
double interpl(double *x,double *y,double x0,int nx);
double interplilo(double *x,double *y,double x0,int nx);
double pXceed(double, double, double, double);
float gammp(float, float);
float betai(float, float, float);
double normCDF(double);
double normPDF(double, double, double);

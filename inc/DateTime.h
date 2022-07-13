#ifndef DTime
#define DTime

struct DateTime
   {
   int yr;
   int mo;
   int dy;
   int hr;
   int mi;
   double sc;
   };

static struct DateTime DATE_NULL= { 0, 0, 0, 0, 0, 0.0};
static struct DateTime t70 = {1970,1,1,0,0,0.0};


int julian();
int ndays();
float dtime();
struct DateTime addtime();
void jul2cal();
int  julmd();

#endif

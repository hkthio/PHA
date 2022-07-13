#ifndef CMTdef
#define CMTdef

struct MT
   {
   float m0;
   float mrr;
   float mtt;
   float mff;
   float mrt;
   float mrf;
   float mtf;
   int exp;
   };
struct FP
   {
   float str;
   float dip;
   float rak;
   };
struct AX
   {
   float mo;
   float pl;
   float az;
   };

static struct MT CMT_NULL = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float getm0();

#endif

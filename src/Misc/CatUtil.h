#include "Sphere.h"
#include "DateTime.h"
#include "CMT.h"

struct SubFault
   {
   struct geographic lo;
   struct FP fp;
   float D;
   float L;
   float W;
   int ix;
   int iy;
   int is;
   };

struct Event
   {
   struct geographic lo;
   struct geographic cm;
   struct DateTime ti;
   struct MT mt;
   struct FP fp1;
   struct FP fp2;
   float mag;
   char  magtype[2];
   float MS;
   float mb;
   float ML;
   int flag;
   float t70;
   };


void CMTWriteTbl();
void CMTWriteTblMag();
int CatInside();
int CMTReadTbl();
int CMTReadTblMag();
int SocalMechReadTbl();
int ISCMechReadTbl();
void ISCWriteTblMech();
int EventRead();
int CatRead();
int CatReadMag();
int CatLines();
int SegLines();
void CatWrite();
void CatWriteMag();
struct Event * CatSortDate();
void CatSort();
int CatTrim();
void LocWriteLine();
void MecWriteLine();
void EndLine();
void DateWriteLine();
void MTWriteLine();
int InvRead();


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PHA.h"

double GR(double, double, double);
double GR_tr(double, double, double, double, double);
double GR_tr_int(double, double, double, double, double, double);
double Mag2Mom(double);
double Mom2Mag(double);

int main()
   {
   /* conventions: amin = log Nmin     */
   /* Nmin = is rate per year M > Mmin */
   double L, W, A, D;
   double amin, Nmin, b, beta, M;
   double Mom, Mmin, Mmax, N0;
   /* integration parameters */
   double dMag, dMag2, dN;
   int i,nMag;
   double SlipRate,Mom_cumul, N_cumul;

   dMag=0.01;
   dMag2=dMag/2.;
   Mmin=4.5;

   /*SSC3-GAWA*/
   Mmax=7.6;
   L=300.e3;
   W=16.e3;
   Nmin=.3129;
   beta=2.439;
   /*SSC1-SAGO*/
   Mmax=7.5;
   L=120.e3;
   W=16.e3;
   Nmin=.1042;
   beta=2.382;
 
   amin=log10(Nmin);
   b=beta/log(10.);
   nMag=(Mmax-Mmin)/dMag-1;

   /* GR-relation */
   Mom_cumul=0;
   N_cumul=0.0;
   
   printf("  M         Mom     dN(M)        Mom_c\n");
   for(i=0;i<nMag;i++)
      {
      M=Mmin+dMag2+i*dMag;
      Mom=Mag2Mom(M);
      dN=GR_tr_int(amin,b,M-dMag2,M+dMag2,Mmin,Mmax);
      N_cumul+=dN;
      Mom_cumul+=Mom*dN;
      printf("%5.3f %12.6e %12.6e %12.6e\n",M,Mom,dN,Mom_cumul);
      }
   SlipRate=Mom_cumul/(L*W*Mu);
   printf("Moment rate   N check     Slip rate\n");
   printf("%12.6e % 9.6f %12.6e\n",Mom_cumul,N_cumul/Nmin,SlipRate);
   return 1;
   }

double GR_tr(double amin, double b, double M, double Mmin, double Mmax)
   {
   /* number of events between Mmin and M */
   double rate,N;
   N=pow(10.,amin);
   if(M > Mmax)
      {
      rate=N;
      }
   else if(M < Mmin)
      {
      rate = 0.0;
      }
   else
      {
      rate=N*(1-pow(10.,-b*(M-Mmin)))/(1-pow(10.,-b*(Mmax-Mmin)));
      }
   return rate;
   }

double GR_tr_int(double amin, double b, double M1, double M2, double Mmin, double Mmax)
   {
   double rate;
   if(M1 > M2)
      {
      printf("Error GR_tr_int: M1 should be smaller than M2\n");
      return -1;
      }
   else
      {
      rate=(GR_tr(amin,b,M2,Mmin,Mmax)-GR_tr(amin,b,M1,Mmin,Mmax));
      }
   return rate;
   }

double GR(double a, double b, double M)
   {
   double rate;
   rate=pow(10.,a-b*M);
   return rate;
   }

double Mag2Mom(double Mag)
   {
   double Mom;
   Mom=pow(10.,1.5*Mag+9.1);
   return Mom;
   }

double Mom2Mag(double Mom)
   {
   double Mag;
   Mag=(log10(Mom)-9.1)/1.5;
   return Mag;
   }

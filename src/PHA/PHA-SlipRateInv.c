/* Code for the inversion of sliprate distributions for PHA */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../inc/CatUtil.h"
int main()
   {
   /* Common variables */
   FILE *fp_in;
   int i,j,k;

   /* define a rupture: nx, ny(nx), n, sf(n)*/
   struct Rupture 
      {
      int n;
      int nx;
      int *ny;
      struct SubFault *sf;
      float Area;   /* actual area */
      float Length; /* nominal length */
      float Width;  /* nominal width */
      } Rup;

   float WtScale[3];

   /* constants */
   float Mu=32.e9;

   /*-- parameters for epistemic branches*/
   /* Scaling relations */
   int nScaling,*iMagArea,iScale;
   float *wtMagArea,Area;

   /*top of the rupture */
   int nTopRup, nBotRup;
   float *wtTopRup,*wtBotRup,*zTopRup,*zBotRup;
   int nW, nL;

   /* Rupture variables */
   int nSubFaults,*nSegments;
   float *wtSeg,*SlipRate;
   int **iStart, **iEnd;

   /* Segmentation models */
   int nSegMod;

   /* Logic tree setup */
   int nBranch[6],nLogicTree,i1,i2,i3,i4,i5,i6;
   int j0,j1;
   int k0;

   /* --- Aleatory braches ---*/

   /* Magnitude scaling */
   int nSigMag=3,iEps;
   float Mag,Mom;

   /* Slip */
   int nAsper;
   float Dav;

   /* distinguish between epistemic (i*), aleatory (j*), and segment (k*) loops */
   int n_epi_, n_ale_slp;
   /*--- finished declarations ---*/


   /* read i_invall with rates per subfault */
   nSubFaults=InvRead("i_invall",Rup.sf);

   /* read scaling relations */
   fscanf(fp_in,"%d",&nScaling);
   iMagArea  = (int   *) malloc(nScaling*sizeof(int));
   wtMagArea = (float *) malloc(nScaling*sizeof(float));
   for(j=0;j<nScaling;j++) fscanf(fp_in,"%d",(iMagArea+j));
   for(j=0;j<nScaling;j++) fscanf(fp_in,"%f",(wtMagArea+j));

   /* read in scenarios catalog for logic tree branch */

   /* read segmentation models */
   fscanf(fp_in,"%d",&nSegMod);
   wtSeg     = (float *) malloc(nSegMod*sizeof( float));
   nSegments = (int *)   malloc(nSegMod*sizeof( int));
   iStart =    (int **)  malloc(nSegMod*sizeof( int *));
   iEnd   =    (int **)  malloc(nSegMod*sizeof( int *));
   for(i=0;i<nSegMod;i++)
      {
      /* read segments */
      fscanf(fp_in,"%d %f",(nSegments+i),(wtSeg+i));
      iStart[i] =    (int *)   malloc(nSegments[i]*sizeof(int));
      iEnd[i]   =    (int *)   malloc(nSegments[i]*sizeof(int));
      SlipRate =  (float *) malloc(nSegments[i]*sizeof(float));
      for(j=0;j<nSegments[i];j++) fscanf(fp_in,"%d",(iStart[i]+j));
      for(j=0;j<nSegments[i];j++) fscanf(fp_in,"%d",(iEnd[i]+j));
      }


   /* loop over logic tree branches */
   for(i1=0;i1<nScaling;i1++)
      {
      /* branch 0 */
      /* area magnitude scaling relations for now */
      for(i2=0;i2<nTopRup;i2++)
         {
         /* branch 1 */
         /* top of rupture */
         for(i3=0;i3<nBotRup;i3++)
            {
            /* branch 2 */
            /* bottom of the rupture */
            for(i4=0;i4<nSegMod;i4++)
               {
               /* branch 4 */
               for(i5=0;i5<1;i5++)
                  {
                  /* branch 5 */
                  for(i6=0;i6<1;i6++)
                     {
                     /* branch 6 */
                     for(k0=0;k0<nSegments[i4];k0++)
                        {
                        /* all sources within this Tree */
                        /* determine the rupture extent for this segment */
                        for(i=iStart[i4][k0];i<=iEnd[i4][k0];i++)
                           {
                           for(j=0;j<nW[i];j++)
                              {
                              if(dep > zBot)
                                 {
                                 }
                              }
                           }
                        for(j0=1-nSigMag;j0<nSigMag;j0++)
                           /* integrate over aleatory magnitude variability */
                           {
                           Mag=MagArea(Area,iScale,iEps);
                           Mom=pow(10,1.5*Mag+16.1);
                           Dav=Mom/(Mu*Area);
                           for(j0=0;j0<nAsper;j0++)
                              {
                              /* Slip distribution aleatory branch */
                              /* slip location   */
                              for(j1=0;j1<nAleatory[1];j1++)
                                 {
                                 /* aleatory branch 1 */
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }





struct Rupture XtractRup(struct Rupture *Rup, int iS, int iE,float Top, float Bot)
   /* Set negative slip for araes of faults that non-active */
   {
   struct Rupture new;
   int i;

   newRup=Rup;
   for(i=iS;i<iE;i++)
      {
      newRup.sf[i]=0.0;
      }
   }

struct SubFault SlipAsper(SubFault *sf)
   {
   }


float MagArea(float Area, float epsilon, int iScale)
   {
   /* MagArea Scaling  Strasser, Papazachos, Murotani */
   float aM_A[3] =  [4.441, 3.279, 6.033];
   float bM_A[3] =  [0.846, 1.163, 0.667];
   float sM_A[3] =  [0.286, 0.250, 0.250];
   Mag=a+b*log10(Area)+sig*epsilon;
   return Mag;
   }

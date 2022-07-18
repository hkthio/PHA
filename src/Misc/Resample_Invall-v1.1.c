#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CatUtil.h"
/*    -- program to resample a set of subfaults to a coarser set ---------------
C     --
C     -- ired, jred - reduction factors in x (along strike) and y (along width) 
C     -- directions
C     -- ntaper - number of taper elements OUTSIDE the nominal area of the new
C     -- subfault
C     ------------------------------------------------------------------------*/

#define MAX_LEN   2600
#define MAX_WIDTH  600
#define min(x, y) (((x) < (y)) ? (x) : (y))

float ***alloc_data(size_t, size_t, size_t);
void free_data(float ***a, size_t, size_t);

int main(int argc, char *argv[])
   {
   /* input variables */
   int ired, jred, ntaper1;

   /* mask variables */
   int ntaper2, nsumx, nsumy;
   int i,j,k,l,m,n;
   float a[100][100],b[100][100];
   float scale;
   float ilo,ila;

   /* subfault variables */
   float alo,ala,adp,str,dip,rak,slp,aL,aW;
   /* output variables */
   float avlo,avla,avdp,avst,avdi,avrk,avsl,avLn,avWd,cumA;
   int kl,km,nl,ll,nlred,nmred,lprev,nsubs;
   FILE *fp;

   int minx=0, maxx=MAX_LEN;
   float maxdep=670.;
   
   char carg[80],line[80];
   int nm[MAX_LEN];
   float ***flo, ***fla;
   struct SubFault **sf,**sfr;

   /*printf("Initializing arrays\n");*/
   sf  = (struct SubFault **) calloc(MAX_LEN,sizeof( struct SubFault *));
   sfr = (struct SubFault **) calloc(MAX_LEN,sizeof( struct SubFault *));
   flo = alloc_data(MAX_LEN,MAX_WIDTH,4);
   fla = alloc_data(MAX_LEN,MAX_WIDTH,4);
   for(k=0;k<MAX_LEN;k++)
      {
      sf[k]  = (struct SubFault *) calloc (MAX_WIDTH,sizeof(struct SubFault));
      sfr[k] = (struct SubFault *) calloc (MAX_WIDTH,sizeof(struct SubFault));
      }

   sscanf(argv[1],"%d",&ired);
   sscanf(argv[2],"%d",&jred);
   sscanf(argv[3],"%d",&ntaper1);
   if(argc > 4)
      {
      sscanf(argv[4],"%f",&maxdep);
      }
   if(argc > 5)
      {
      sscanf(argv[5],"%d",&minx);
      sscanf(argv[6],"%d",&maxx);
      }
      
   /* set slip taper -- */
   ntaper2=2*ntaper1;
   nsumx=ired+ntaper2;
   nsumy=jred+ntaper2;
   scale=1./(1+ntaper2);
   /* determine the taper matrix */
   for(i=0;i<nsumx;i++)
      {
      if(i < ntaper2 )
         {
         l=ntaper2-i;
         }
      else if(i >= ired)
         {
         l=i-ired+1;
         }
      else
         {
         l=0;
         }
      for(j=0;j<nsumy;j++)
         {
         if(j < ntaper2 )
            {
            m=ntaper2-j;
            }
         else if(j >= jred) 
            {
            m=j-jred+1;
            }
         else
            {
            m=0;
            }
         a[i][j]=1.-sqrt(l*l+m*m)*scale;
         if(a[i][j] < 0.0) a[i][j]=0.0;
         }    
      /*for(j=0;j<nsumy;j++) printf("%6.4f ", a[i][j]);
      printf("\n");*/
      }
   /*printf("\n");*/
   for(i=0;i<nsumx;i++)
      {
      for(j=0;j<nsumy;j++)
         {
         if(i >= ntaper2 && i < ired)
            {
            l=i;
            if(j >= ntaper2 && j < jred) 
               {
               b[i][j]=a[i][j];
               }
            else
               {
               if(j < ntaper2)
                  {
                  m=j+jred;
                  }
               else
                  {
                  m=j-jred;
                  }
               b[i][j]=a[i][j]+a[l][m];
               }
            }
         else
            {
            if(i < ntaper2)
               {
               l=i+ired;
               }
            else
               {
               l=i-ired;
               }
            if(j >= ntaper2 && j < jred)
               {
               m=j;
               b[i][j]=a[i][j]+a[l][m];
               }
            else
               {
               if(j < ntaper2)
                  {
                  m=j+jred;
                  }
               else
                  {
                  m=j-jred;
                  }
               b[i][j]=a[i][j]+a[i][m]+a[l][j]+a[l][m];
               }
            }
         }
      /*for(j=0;j<nsumy;j++) printf("%6.4f ", b[i][j]);
      printf("\n");*/
      }
   /*printf("\n");*/
   for(i=0;i<nsumx;i++)
      {
      for(j=0;j<nsumy;j++)
         {
         a[i][j] /= b[i][j];
         }
      /*for(j=0;j<nsumy;j++) printf("%6.4f ", a[i][j]);
      printf("\n");*/
      }

   for(l=0;l<MAX_LEN;l++) nm[l]=0;

   fp=fopen("o_subfaults","r");
   /*printf("Reading subfaults\n");*/

   lprev=-1;
   nl=0;
   i=0;
   while(scanf("%f %f %f %f %f %f %f %f %f %d %d %d\n",&alo,&ala,&adp,&str,&dip,&rak,&slp,&aL,&aW,&l,&m,&n) == 12)
      {
      /*printf("%d %d %d\n",l,m,n);*/
      sf[l][m].lo.lon = alo;
      sf[l][m].lo.lat = ala;
      sf[l][m].lo.dep = adp;
      sf[l][m].fp.str = str;
      sf[l][m].fp.dip = dip;
      sf[l][m].fp.rak = rak;
      sf[l][m].D      = slp;
      sf[l][m].L      = aL;
      sf[l][m].W      = aW;
      sf[l][m].ix     = l;
      sf[l][m].iy     = m;
      sf[l][m].is     = n;
      nm[l]=nm[l]+1;
      if(l != lprev) 
         {
         lprev=l;
         nl++;
         }
      /* read subfaults */
      fscanf(fp, "%*[^\n]\n", NULL);
      for(k=0;k<4;k++) 
         {
         fscanf(fp,"%f %f\n",&ilo,&ila);
         flo[l][m][k]=ilo;
         fla[l][m][k]=ila;
         /*printf("%d %d %9.4f %9.4f\n",l,m,flo[l][m][k],fla[l][m][k]);*/
         }
      }
   fclose(fp);

   fp=fopen("o_subfaults_red","w");
   /* resampling */
   nlred=(nl-ntaper2)/ired;
   for(kl=0;kl<nlred;kl++)
      {
      nmred=9999;
      /* determine smallest width for this new row */
      for(i=0;i<ired;i++) nmred=min(nmred,(nm[kl*ired+i+ntaper1]-ntaper2)/jred);
      for(km=0;km<nmred;km++)
         {
         nsubs=0;
         avlo=avla=avdp=avst=avdi=avrk=avsl=avLn=avWd=cumA=0.0;
         for(i=0;i<ired;i++)
            {
            l=ntaper1+i+kl*ired;
            for(j=0;j<jred;j++)
               {
               m=ntaper1+j+km*jred;
               avlo +=  sf[l][m].lo.lon;
               avla +=  sf[l][m].lo.lat;
               avdp +=  sf[l][m].lo.dep;
               avst +=  sf[l][m].fp.str;
               avdi +=  sf[l][m].fp.dip;
               avrk +=  sf[l][m].fp.rak;
               avsl +=  sf[l][m].D;
               avLn +=  sf[l][m].L;
               avWd +=  sf[l][m].W;
               cumA +=  sf[l][m].L*sf[l][m].W;
               nsubs++;
               }
            }
         avlo /= nsubs;
         avla /= nsubs;
         avdp /= nsubs;
         avst /= nsubs;
         avdi /= nsubs;
         avrk /= nsubs;
         avsl /= nsubs;
         /* get average length, make sure LxW is consistent with area */
         avLn= avLn/jred;
         if(avdp < maxdep && kl >= minx && kl <= maxx)
            {
            printf("%9.4f %8.4f %6.2f %7.2f %7.2f %7.2f %5.1f %5.1f %5.1f %3d %3d %3d\n",avlo,avla,avdp,avst,avdi,avrk,avsl,avLn,cumA/avLn,kl+1,km+1,-1);
            fprintf(fp,"> %9.3f %8.3f %6.1f\n",avlo,avla,avdp);
            ll=ntaper1;
            fprintf(fp,"%10.4f %9.4f\n",flo[ll+    kl*ired][ll+    km*jred][0],fla[ll+    kl*ired][ll+    km*jred][0]);
            fprintf(fp,"%10.4f %9.4f\n",flo[ll+(kl+1)*ired-1][ll+    km*jred][1],fla[ll+(kl+1)*ired-1][ll+    km*jred][1]);
            fprintf(fp,"%10.4f %9.4f\n",flo[ll+(kl+1)*ired-1][ll+(km+1)*jred-1][2],fla[ll+(kl+1)*ired-1][ll+(km+1)*jred-1][2]);
            fprintf(fp,"%10.4f %9.4f\n",flo[ll+    kl*ired][ll+(km+1)*jred-1][3],fla[ll+    kl*ired][ll+(km+1)*jred-1][3]);
            for(i=0;i<nsumx;i++)
               {
               l=i+kl*ired;
               for(j=0;j<nsumy;j++)
                  {
                  m=j+km*jred;
                  if(a[i][j]*sf[l][m].D > 0.0) 
                     {
                     printf("%9.4f %8.4f %6.2f %7.2f %7.2f %7.2f %5.1f %5.1f %5.1f %3d %3d %3d\n",
                             sf[l][m].lo.lon,sf[l][m].lo.lat,sf[l][m].lo.dep,
                             sf[l][m].fp.str,sf[l][m].fp.dip,sf[l][m].fp.rak,
                             a[i][j]*sf[l][m].D,sf[l][m].L,sf[l][m].W,l+1,m+1,1);
                     }
                  }
               }
            }
         }
      }
   fclose(fp);
   return 0;
   }

float ***alloc_data(size_t xlen, size_t ylen, size_t zlen)
{
    float ***p;
    size_t i, j;

    if ((p = malloc(xlen * sizeof *p)) == NULL) {
        perror("malloc 1");
        return NULL;
    }

    for (i=0; i < xlen; ++i)
        p[i] = NULL;

    for (i=0; i < xlen; ++i)
        if ((p[i] = malloc(ylen * sizeof *p[i])) == NULL) {
            perror("malloc 2");
            free_data(p, xlen, ylen);
            return NULL;
        }

    for (i=0; i < xlen; ++i)
        for (j=0; j < ylen; ++j)
            p[i][j] = NULL;

    for (i=0; i < xlen; ++i)
        for (j=0; j < ylen; ++j)
            if ((p[i][j] = malloc(zlen * sizeof *p[i][j])) == NULL) {
                perror("malloc 3");
                free_data(p, xlen, ylen);
                return NULL;
            }

    return p;
}

void free_data(float ***data, size_t xlen, size_t ylen)
{
    size_t i, j;

    for (i=0; i < xlen; ++i) {
        if (data[i] != NULL) {
            for (j=0; j < ylen; ++j)
                free(data[i][j]);
            free(data[i]);
        }
    }
    free(data);
}


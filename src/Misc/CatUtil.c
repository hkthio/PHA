#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "CatUtil.h"
#include "DateTime.h"
/* Important: we need to introdcude a flag which specifies the order of the lon lat pairs
   suggest 1 - lon lat and -1 - lat lon */

int CatRead(char *filename, struct Event **tmp)
   {
   /* reading the base catalog */
   /* event array is intialized here, which is why we are using pointers to  */
   /* pointers */
   FILE *fp;
   int i,nevents,nl,nread;
   struct Event *ev;

   i=0;
   if((fp=fopen(filename,"r")) == NULL)
      {
      printf("Cannot open file %s\n",filename);
      return -1;
      }
   nl=CatLines(fp);
   ev = (struct Event *) malloc(nl*sizeof(struct Event));

   printf("Reading catalog %s - %9d\n",filename,nl);
   while((nread=fscanf(fp,"%d %d %d %d %d %lf %lf %lf %lf %f",&ev[i].ti.yr,&ev[i].ti.mo,
                    &ev[i].ti.dy,&ev[i].ti.hr,&ev[i].ti.mi,&ev[i].ti.sc,
                    &ev[i].lo.lon,&ev[i].lo.lat,&ev[i].lo.dep,&ev[i].mag)) >= 9 && i < nl)
      {
      if(nread==9) 
         {
         ev[i].mag=-9.;
         printf("No magnitude found: %4d\n",i);
         }
      ev[i].flag=0;
      ev[i].t70=dtime(ev[i].ti,t70);
      i++;
      }
   nevents=i;
   fclose(fp);
   printf("Read catalog %s - %4i events\n",filename,nevents);
   *tmp=ev;
   return nevents;
   }
int CatReadMag(char *filename, struct Event **tmp)
   {
   /* reading the base catalog */
   /* event array is intialized here, which is why we are using pointers to  */
   /* pointers */
   FILE *fp;
   char mtype[8],line[512];
   int i,nevents,nl,nread,pos,buf;
   int imb,iMS,iML,imag;
   struct Event *ev;
   float mag;

   if((fp=fopen(filename,"r")) == NULL)
      {
      printf("Cannot open file %s\n",filename);
      return -1;
      }
   nl=CatLines(fp);
   ev = (struct Event *) malloc(nl*sizeof(struct Event));

   printf("Reading catalog %s - %9d\n",filename,nl);
   i=0;
   while(fgets(line,512,fp))
      {
      buf=0;
      imb=iMS=iML=imag=0;
      nread=sscanf(line,"%d %d %d %d %d %lf %lf %lf %lf%n",&ev[i].ti.yr,&ev[i].ti.mo,
                       &ev[i].ti.dy,&ev[i].ti.hr,&ev[i].ti.mi,&ev[i].ti.sc,
                       &ev[i].lo.lon,&ev[i].lo.lat,&ev[i].lo.dep,&pos);
      buf+=pos;
      while((nread=sscanf(line+buf,"%f %s%n",&mag,mtype,&pos))==2)
         {
         buf+=pos;
         strncpy(ev[i].magtype,mtype,2);
         if(mtype[1]=='w'|| mtype[1] == 'W')
            {
            if(ev[i].mag <= 0.0) ev[i].mag= mag;
            } 
         else if(mtype[1]=='s'|| mtype[1] == 'S')
            {
            if(ev[i].MS <= 0.0) ev[i].MS = mag;
            if(ev[i].mag <= 0.0) 
               {
               ev[i].mag = mag;
               strncpy(ev[i].magtype,mtype,2);
               }
            } 
         else if(mtype[1]=='b'|| mtype[1] == 'B')
            {
            if(ev[i].mb <= 0.0) ev[i].mb = mag;
            if(ev[i].mag <= 0.0) 
               {
               ev[i].mag = mag;
               strncpy(ev[i].magtype,mtype,2);
               }
            }
         else if(mtype[1]=='l'|| mtype[1] == 'L')
            {
            if(ev[i].ML <= 0.0) ev[i].ML = mag;
            if(ev[i].mag <= 0.0) 
               {
               ev[i].mag = mag;
               strncpy(ev[i].magtype,mtype,2);
               }
            }
         }
      ev[i].flag=0;
      ev[i].t70=dtime(ev[i].ti,t70);
      i++;
      }
   nevents=i;
   fclose(fp);
   *tmp=ev;
   printf("Read catalog %s - %4i events\n",filename,nevents);
   return nevents;
   }

int CMTReadTbl(char *filename, struct Event **tmp)
   {
   /* reading CMT catalog in table format */
   /* event array is intialized here, which is why we are using pointers to  */
   /* pointers */
   /* format: 2017  6 30 13 33 44.7   1.79  58.92  10.0    1.90  59.00  12.0 23  1.370  0.157 -1.530  0.927 -0.164  0.221   1.731  17 45  130 146 57   57*/

   FILE *fp;
   int i,nevents,nl;
   char cdum[80],datestring[12];
   struct Event *ev;

   i=0;
   if((fp=fopen(filename,"r")) == NULL)
      {
      printf("Cannot open file %s\n",filename);
      return -1;
      }
   nl=CatLines(fp);
   ev = (struct Event *) malloc(nl*sizeof(struct Event));

   printf("Reading catalog %s - %9d\n",filename,nl);
   while(fscanf(fp,"%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %d %f %f %f %f %f %f %f %f %f %f %f %f %f", 
                &ev[i].ti.yr,&ev[i].ti.mo,&ev[i].ti.dy,&ev[i].ti.hr,&ev[i].ti.mi,&ev[i].ti.sc,
                &ev[i].lo.lon,&ev[i].lo.lat,&ev[i].lo.dep,
                &ev[i].cm.lon,&ev[i].cm.lat,&ev[i].cm.dep,
                &ev[i].mt.exp,
                &ev[i].mt.mrr,&ev[i].mt.mtt,&ev[i].mt.mff,
                &ev[i].mt.mrt,&ev[i].mt.mrf,&ev[i].mt.mtf,
                &ev[i].mt.m0,
                &ev[i].fp1.str,&ev[i].fp1.dip,&ev[i].fp1.rak,
                &ev[i].fp2.str,&ev[i].fp2.dip,&ev[i].fp2.rak)==26
                    && i < nl)
      {
      /*printf("%4d\n",ev[i].ti.yr);*/
      ev[i].flag=0;
      ev[i].ti.sc=0.0;
      ev[i].t70=dtime(ev[i].ti,t70);
      i++;
      }
   nevents=i;
   fclose(fp);
   printf("Read catalog %s - %4i events\n",filename,nevents);
   *tmp=ev;
   return nevents;
   }

int CMTReadTblMag(char *filename, struct Event **tmp)
   {
   FILE *fp;
   int i,nevents,nl;
   char cdum[80],datestring[12];
   struct Event *ev;

   i=0;
   if((fp=fopen(filename,"r")) == NULL)
      {
      printf("Cannot open file %s\n",filename);
      return -1;
      }
   nl=CatLines(fp);
   ev = (struct Event *) malloc(nl*sizeof(struct Event));

   printf("Reading catalog %s - %9d\n",filename,nl);
   while(fscanf(fp,"%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
                &ev[i].ti.yr,&ev[i].ti.mo,&ev[i].ti.dy,&ev[i].ti.hr,&ev[i].ti.mi,&ev[i].ti.sc,
                &ev[i].lo.lon,&ev[i].lo.lat,&ev[i].lo.dep,
                &ev[i].cm.lon,&ev[i].cm.lat,&ev[i].cm.dep,
                &ev[i].mt.exp,
                &ev[i].mt.mrr,&ev[i].mt.mtt,&ev[i].mt.mff,
                &ev[i].mt.mrt,&ev[i].mt.mrf,&ev[i].mt.mtf,
                &ev[i].mt.m0,
                &ev[i].fp1.str,&ev[i].fp1.dip,&ev[i].fp1.rak,
                &ev[i].fp2.str,&ev[i].fp2.dip,&ev[i].fp2.rak,
                &ev[i].mag,&ev[i].mb,&ev[i].MS,&ev[i].ML)==30
                    && i < nl)
      {
      /*printf("%f\n",ev[i].MS);*/
      ev[i].flag=0;
      ev[i].ti.sc=0.0;
      ev[i].t70=dtime(ev[i].ti,t70);
      i++;
      }
   nevents=i;
   fclose(fp);
   printf("Read catalog %s - %4i events\n",filename,nevents);
   *tmp=ev;
   return nevents;
   }

int SocalMechReadTbl(char *filename, struct Event **tmp)
   {
   /* reading SCSN mechanism catalog (Yang, W., E. Hauksson and P. M. Shearer, 
      Computing a large refined catalog of focal mechanisms for southern 
      California (1981 - 2010): Temporal Stability of the Style of Faulting, 
      Bull. Seismol. Soc. Am., June 2012, v. 102, p. 1179-1194, 
      doi:10.1785/0120110311, 2012) */
   /* event array is intialized here, which is why we are using pointers to  */
   /* pointers */
   /* format: 1981  1  1  4 13 55.710  3301565  33.25517 -115.96750   5.680  2.260  318  57 -168  37  39   18  0.17    0  0.00 C */

   FILE *fp;
   int i,nevents,nl,evid;
   char cdum[80],datestring[12];
   struct Event *ev;

   i=0;
   if((fp=fopen(filename,"r")) == NULL)
      {
      printf("Cannot open file %s\n",filename);
      return -1;
      }
   nl=CatLines(fp);
   ev = (struct Event *) malloc(nl*sizeof(struct Event));

   printf("Reading catalog %s - %9d\n",filename,nl);
   while(fscanf(fp,"%d %d %d %d %d %lf \
                    %d %lf %lf %lf \
                    %f %f %f %f \
                    %*d %*d %*d %*f %*d %*f %*s", 
                &ev[i].ti.yr,&ev[i].ti.mo,&ev[i].ti.dy,&ev[i].ti.hr,&ev[i].ti.mi,&ev[i].ti.sc,
                &evid,&ev[i].lo.lon,&ev[i].lo.lat,&ev[i].lo.dep,
                &ev[i].mag, &ev[i].fp1.str,&ev[i].fp1.dip,&ev[i].fp1.rak)==14 && i < nl)
      {
      printf("%4d\n",ev[i].ti.yr);
      ev[i].flag=0;
      ev[i].ti.sc=0.0;
      ev[i].t70=dtime(ev[i].ti,t70);
      i++;
      }
   nevents=i;
   fclose(fp);
   printf("Read catalog %s - %4i events\n",filename,nevents);
   *tmp=ev;
   return nevents;
   }

int ISCMechReadTbl(char *filename, struct Event **tmp)
   {
   /* reading the ISC mechanism catalog  */
   /* the CSV file is first reformatted to give yr mo dy hr mi sc lon lat dep str dip rak */
   /* event array is intialized here, which is why we are using pointers to  */
   /* pointers */
   /* format: 1981  1  1  4 13 55.710  33.25517 -115.96750  5.680 318  57 -168  */

   FILE *fp;
   int i,nevents,nl,evid;
   char cdum[80],datestring[12];
   struct Event *ev;

   i=0;
   if((fp=fopen(filename,"r")) == NULL)
      {
      printf("Cannot open file %s\n",filename);
      return -1;
      }
   nl=CatLines(fp);
   ev = (struct Event *) malloc(nl*sizeof(struct Event));

   printf("Reading catalog %s - %9d\n",filename,nl);
   while(fscanf(fp,"%d %d %d %d %d %lf \
                    %lf %lf %lf \
                    %f %f %f",\
                &ev[i].ti.yr,&ev[i].ti.mo,&ev[i].ti.dy,&ev[i].ti.hr,&ev[i].ti.mi,&ev[i].ti.sc,
                &ev[i].lo.lon,&ev[i].lo.lat,&ev[i].lo.dep,
                &ev[i].fp1.str,&ev[i].fp1.dip,&ev[i].fp1.rak)==12 && i < nl)
      {
      ev[i].flag=0;
      ev[i].ti.sc=0.0;
      ev[i].t70=dtime(ev[i].ti,t70);
      i++;
      }
   nevents=i;
   fclose(fp);
   printf("Read catalog %s - %4i events\n",filename,nevents);
   *tmp=ev;
   return nevents;
   }

int CMTRead(char *filename, struct Event **tmp)
   {
   /* reading CMT catalog in meca format */
   /* event array is intialized here, which is why we are using pointers to  */
   /* pointers */
   /* format: -71.10 -25.75 75 -0.34 -0.61 0.95 0.59 -1.76 0.46 24 X Y 201004101506A*/

   FILE *fp;
   int i,nevents,nl;
   char cdum[80],datestring[12];
   struct Event *ev;

   i=0;
   if((fp=fopen(filename,"r")) == NULL)
      {
      printf("Cannot open file %s\n",filename);
      return -1;
      }
   nl=CatLines(fp);
   ev = (struct Event *) malloc(nl*sizeof(struct Event));

   printf("Reading catalog %s - %9d\n",filename,nl);
   while(fscanf(fp,"%lf %lf %lf %f %f %f %f %f %f %d %s %s %s", 
                &ev[i].lo.lon,&ev[i].lo.lat,&ev[i].lo.dep,
                &ev[i].mt.mrr,&ev[i].mt.mtt,&ev[i].mt.mff,
                &ev[i].mt.mrt,&ev[i].mt.mrf,&ev[i].mt.mtf,
                &ev[i].mt.exp,cdum,cdum,datestring)== 13 
                    && i < nl)
      {
      ev[i].flag=0;
      sscanf(datestring,"%4d%2d%2d%2d%2d",&ev[i].ti.yr,&ev[i].ti.mo,&ev[i].ti.dy,
                                         &ev[i].ti.hr,&ev[i].ti.mi);
      ev[i].ti.sc=0.0;
      ev[i].t70=dtime(ev[i].ti,t70);
      i++;
      }
   nevents=i;
   fclose(fp);
   printf("Read catalog %s - %4i events\n",filename,nevents);
   *tmp=ev;
   return nevents;
   }

struct Event * CatSortDate(struct Event *ev, int n, int sign)
   {
   /* This version returns initializes and returns a new array */
   int i,j;
   struct Event temp, *evnew;
   evnew = (struct Event *) malloc(n*sizeof(struct Event));
   
   if(sign != -1 && sign !=1)
      {
      printf("Invalid value for sign, should be 1 or -1, not %d\n",sign);
      exit(-1);
      }
   for (i = 0; i < n; i++) 
      {
      /*printf("%d\n",ev[i].ti.yr);*/
      evnew[i]=ev[i];
      }
   for (i = 1; i < n; i++)
      {
      for (j = 0; j < n - i; j++)
         {
         if (sign*evnew[j].t70 < sign*evnew[j+1].t70)
            {
            temp = evnew[j];
            evnew[j] = evnew[j + 1];
            evnew[j + 1] = temp;
            }
         }
      }
   return evnew;
   }

void CatSort(struct Event *ev, int n)
   {
   int i,j;
   struct Event temp;
   for (i = 1; i < n; i++)
      {
      for (j = 0; j < n - i; j++)
         {
          if (ev[j].t70 < ev[j+1].t70)
            {
            temp = ev[j];
            ev[j] = ev[j + 1];
            ev[j + 1] = temp;
            }
         }
      }
   }

void CatWriteMag(struct Event ev, FILE *fp, int sign)
   {
   DateWriteLine(ev.ti,fp);
   LocWriteLine(ev.lo,fp);
   /*fprintf(fp,"%4.2f %4.2f %4.2f %4.2f",ev.mag,ev.mb,ev.MS,ev.ML);*/
   fprintf(fp,"%4.2f %s",ev.mag,ev.magtype);
   EndLine(fp);
   }

void CatWrite(struct Event ev, FILE *fp, int sign)
   {
   DateWriteLine(ev.ti,fp);
   LocWriteLine(ev.lo,fp);
   fprintf(fp,"%4.2f ",ev.mag);
   EndLine(fp);
   }

void CMTWriteTblMag(struct Event ev, FILE *fp, int sign)
   {
   if(sign < 0)
      {
      if(ev.lo.lon > 180. && ev.lo.lon <= 360.) ev.lo.lon -=360.;
      }
   else if(sign > 0)
      {
      if(ev.lo.lon >= -180. && ev.lo.lon < 0.0) ev.lo.lon +=360.;
      }
   DateWriteLine(ev.ti,fp);
   LocWriteLine(ev.lo,fp);
   LocWriteLine(ev.cm,fp);
   MTWriteLine(ev.mt,fp);
   MecWriteLine(ev.fp1,fp);
   MecWriteLine(ev.fp2,fp);
   fprintf(fp,"%4.2f %4.2f %4.2f %4.2f",ev.mag,ev.mb,ev.MS,ev.ML);
   EndLine(fp);
   }

void ISCWriteTblMech(struct Event ev, FILE *fp, int sign)
   {
   if(sign < 0)
      {
      if(ev.lo.lon > 180. && ev.lo.lon <= 360.) ev.lo.lon -=360.;
      }
   else if(sign > 0)
      {
      if(ev.lo.lon >= -180. && ev.lo.lon < 0.0) ev.lo.lon +=360.;
      }
   DateWriteLine(ev.ti,fp);
   LocWriteLine(ev.lo,fp);
   MecWriteLine(ev.fp1,fp);
   fprintf(fp,"%4.2f",ev.mag);
   EndLine(fp);
   }

void CMTWriteTbl(struct Event ev, FILE *fp, int sign)
   {
   if(sign < 0)
      {
      if(ev.lo.lon > 180. && ev.lo.lon <= 360.) ev.lo.lon -=360.;
      }
   else if(sign > 0)
      {
      if(ev.lo.lon >= -180. && ev.lo.lon < 0.0) ev.lo.lon +=360.;
      }
   DateWriteLine(ev.ti,fp);
   LocWriteLine(ev.lo,fp);
   LocWriteLine(ev.cm,fp);
   MTWriteLine(ev.mt,fp);
   MecWriteLine(ev.fp1,fp);
   MecWriteLine(ev.fp2,fp);
   EndLine(fp);
   }

void CMTWrite(struct Event ev, FILE *fp, int sign)
   {
   if(sign < 0)
      {
      if(ev.lo.lon > 180. && ev.lo.lon <= 360.) ev.lo.lon -=360.;
      }
   else if(sign > 0)
      {
      if(ev.lo.lon >= -180. && ev.lo.lon < 0.0) ev.lo.lon +=360.;
      }
   fprintf(fp,"%4i %02i %02i %02i %02i %04.1f %9.4f %8.4f %5.1f %4.2f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %2d\n",
               ev.ti.yr,ev.ti.mo,ev.ti.dy,ev.ti.hr,ev.ti.mi,ev.ti.sc,
               ev.lo.lon,ev.lo.lat,ev.lo.dep,ev.mag,
               ev.mt.mrr,ev.mt.mtt,ev.mt.mff,
               ev.mt.mrt,ev.mt.mrf,ev.mt.mtf,ev.mt.exp);
   }

void LocSign(struct geographic *Loc, int sign)
   /* if sign < 0, -180 < lon < 180 */
   {
   if(sign < 0)
      {
      if(Loc->lon > 180. && Loc->lon <= 360.) Loc->lon -=360.;
      }
   else if(sign >= 0)
      {
      if(Loc->lon >= -180. && Loc->lon < 0.0) Loc->lon +=360.;
      }
   }

/* write out indivdual structures */
int EventRead(FILE *fp, struct Event *ev)
   {
   if(fscanf(fp,"%d %d %d %d %d %lf %lf %lf %lf %f",&ev->ti.yr,&ev->ti.mo,
                    &ev->ti.dy,&ev->ti.hr,&ev->ti.mi,&ev->ti.sc,
                    &ev->lo.lat,&ev->lo.lon,&ev->lo.dep,&ev->mag)== 10)
      {
      return 1;
      }
   else
      {
      return -1;
      }
   }

void LocWriteLine(struct geographic Loc, FILE *fp)
   {
   fprintf(fp,"%9.4f %8.4f %5.1f ", Loc.lon,Loc.lat,Loc.dep);
   }
void DateWriteLine(struct DateTime Date, FILE *fp)
   {
   fprintf(fp,"%4i %02i %02i %02i %02i %04.1f ",
               Date.yr,Date.mo,Date.dy,Date.hr,Date.mi,Date.sc);
   }
void MTWriteLine(struct MT mt, FILE *fp)
   {
   fprintf(fp,"%2d %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f ",
               mt.exp,mt.mrr,mt.mtt,mt.mff,mt.mrt,mt.mrf,mt.mtf,mt.m0);
   }
void MecWriteLine(struct FP ft, FILE *fp)
   {
   fprintf(fp,"%6.1f %4.1f %6.1f ",
               ft.str,ft.dip,ft.rak);
   }
void EndLine(FILE *fp)
   {
   fprintf(fp,"\n");
   }
/* end individual output lines */

int CatTrim(struct Event *ev, int nev, struct Event **tmp, struct DateTime start, struct DateTime stop, struct geographic UL, struct geographic LR)
   {
   struct Event *evnew;
   struct geographic dum;
   evnew = (struct Event *) malloc(nev*sizeof(struct Event));

   int i,j,isign;
   j=0;
   if(UL.lon >= 0.0)
      {
      isign=1;
      if(LR.lon < 0.0) LR.lon += 360.;
      }
   else
      {
      isign=-1;
      if(LR.lon > 180.0) LR.lon -= 360.;
      }
   for(i=0;i<nev;i++)
      {
      if(dtime(ev[i].ti,start) > 0.0 && dtime(ev[i].ti,stop) < 0.0)
         {
         if(ev[i].lo.lat < UL.lat && ev[i].lo.lat > LR.lat)
            {
            dum=ev[i].lo;
            LocSign(&dum,isign);
            if(dum.lon > UL.lon && dum.lon < LR.lon)
               {
               evnew[j]=ev[i];
               j++;
               }
            }
         }
      }
   *tmp=evnew;
   return j;
   }

int CatLines(FILE * fp)
   {
   /* reads number of lines from current position to the end of file,
      and rewinds to current position */
   int lines=0;
   long int pos;
   int ch=0;
   
   pos = ftell(fp);
   while(!feof(fp))
      {
      ch = fgetc(fp);
      if(ch == '\n') lines++;
      }
   fseek(fp, pos, SEEK_SET);
   return lines;
   }

int SegLines(FILE * fp,char *seg, int *nseg, int *maxlines)
   {
   /* reads number of segments and lines  per segment from current position to the end 
      of file, and rewinds to current position, return value is total number of lines*/
   int lines=0,nlines;
   const char *cchr;
   long int pos;
   char *line_buf = NULL;
   size_t line_buf_size = 0;
   ssize_t line_size;
   
   pos = ftell(fp);
   nlines=0;
   *maxlines=0;
   *nseg=0;

   while((line_size = getline(&line_buf, &line_buf_size, fp)) >0)
      {
      lines++;
      cchr=line_buf;
      /*printf("%s %s\n",cchr,seg);*/
      if(strncmp(cchr,seg,1)==0)
         {
         if(nlines > *maxlines) *maxlines=nlines; 
         nlines=0;
         *nseg+=1;
         }
      else
         {
         nlines++;
         }
      }
   fseek(fp, pos, SEEK_SET);
   return lines;
   }

int CatInside(int nvert, struct geographic *vert, struct geographic point)
{
  int i, j, c = 0;
  /*printf("%12.6f %11.6f %12.6f %11.6f\n",point.lon,point.lat,vert[0].lon,vert[0].lat);*/
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((vert[i].lat>point.lat) != (vert[j].lat>point.lat)) &&
	 (point.lon < (vert[j].lon-vert[i].lon) * (point.lat-vert[i].lat) / (vert[j].lat-vert[i].lat) + vert[i].lon) )
       c = !c;
  }
  return c;
}

int InvRead(char *filename, struct SubFault **tmp)
   {
   /* reading i_invall format */
   /* event array is intialized here, which is why we are using pointers to  */
   /* pointers */

   /* subfault variables */
   float alo,ala,adp,str,dip,rak,slp,aL,aW;
   int i,l,m,n,nsub;
   FILE *fp;
   struct SubFault *sf;

   printf("Initializing arrays\n");
   fp=fopen(filename,"r");
   nsub=CatLines(fp);
   sf  = (struct SubFault *) calloc(nsub,sizeof( struct SubFault ));

   i=0;
   while(fscanf(fp,"%f %f %f %f %f %f %f %f %f %d %d %d\n",&alo,&ala,&adp,&str,&dip,&rak,&slp,&aL,&aW,&l,&m,&n) == 12)
      {
      sf[i].lo.lon = alo;
      sf[i].lo.lat = ala;
      sf[i].lo.dep = adp;
      sf[i].fp.str = str;
      sf[i].fp.dip = dip;
      sf[i].fp.rak = rak;
      sf[i].D      = slp;
      sf[i].L      = aL;
      sf[i].W      = aW;
      sf[i].ix     = l;
      sf[i].iy     = m;
      sf[i].is     = n;
      i++;
      }
   if(i != nsub)
      {
      printf("Incorrect entry: %d of %d lines\n",i,nsub);
      return -1;
      }
   fclose(fp);
   *tmp=sf;
   return nsub;
   }

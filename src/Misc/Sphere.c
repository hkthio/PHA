/* v1.0 20220226
Spherical coordinates routines:

Functions to translate from spherical to geographical coordinates and
cartesian coordinates.

Cubed-Sphere transform functions:
regions are defined as follows (folowing Ronchi et al.):
1 - centered on lat0/lon0
2 - east
3 - antipode
4 - west
5 - north
6 - south
*/

#include <math.h>
#include <stdio.h>
#include "Sphere.h"

struct geographic rotgeo(struct geographic old, double *rot) 
   {
   double a[3],b[3];
   struct geographic new;
   struct spherical s;
   struct cartesian c;
   c=sphtocart(geotosph(old));
   a[0]=c.x;
   a[1]=c.y;
   a[2]=c.z;
   matmul(rot,a,b,3,3,1);
   c.x=b[0];
   c.y=b[1];
   c.z=b[2];
   new=sphtogeo(cartosph(c));
   return new;
   }
   
void rotmat(struct cartesian e1, struct cartesian e2, struct cartesian e3, double *rot, double *irot)
   {
   double r1,r2,r3;
   r1=sqrt(e1.x*e1.x+e1.y*e1.y+e1.z*e1.z);
   r2=sqrt(e2.x*e2.x+e2.y*e2.y+e2.z*e2.z);
   r3=sqrt(e3.x*e3.x+e3.y*e3.y+e3.z*e3.z);
   *(irot+0)=*(rot+0)=e1.x/r1;
   *(irot+4)=*(rot+4)=e2.y/r2;
   *(irot+8)=*(rot+8)=e3.z/r3;
   *(irot+1)=*(rot+3)=e2.x/r2;
   *(irot+2)=*(rot+6)=e3.x/r3;
   *(irot+3)=*(rot+1)=e1.y/r1;
   *(irot+5)=*(rot+7)=e3.y/r3;
   *(irot+6)=*(rot+2)=e1.z/r1;
   *(irot+7)=*(rot+5)=e2.z/r2;
   }

void matmul(double *a, double *b, double *c, int l, int m, int n)
   {
   int i,j,k;
   for(i=0;i<m;i++)
      {
      for(k=0;k<n;k++)
         {
         *(c+i*n+k)=0.0;
         for(j=0;j<l;j++)
            {
            *(c+i*n+k)+= *(a+i*l+j) * *(b+j*n+k);
            }
         }
      }
   }

void distaz(struct geographic a,struct geographic b, double *dist, double *az, double *baz)
   {
   struct geographic p,q;
   struct cartesian c1, c2, c3;
   struct cartesian e1, e2, e3;
/* double rot_a[3][3],irot_a[3][3];
   double rot_b[3][3],irot_b[3][3];*/
   double rot_a[9],irot_a[9];
   double rot_b[9],irot_b[9];

   rottopole(a,rot_a,irot_a);
   rottopole(b,rot_b,irot_b);

   p=rotgeo(b,rot_a);
   *dist=90.0-p.lat;
   *az=p.lon;
   q=rotgeo(a,rot_b);
   *baz=q.lon;
   }

double distance(struct geographic a, struct geographic b)
   {
   double dist,az,baz;
   distaz(a,b,&dist,&az,&baz);
   return dist;
   }

void rottoequator(struct geographic a, double az, double *rot, double *irot)
   {
   struct spherical s1, s2, s3;
   struct cartesian c1,c2,c3;
   struct geographic b;

   b=newloc(a,az,90.);
   s1=geotosph(a);
   s3=geotosph(b);
   if(s3.the > PI)
      {
      s3.the=TPI-s3.the;
      s3.phi=s3.phi+PI;
      }
   c3=sphtocart(s3);
   c1=sphtocart(s1);
   c2=vecprod(c3,c1);
   /*c2.x=c3.y*c1.z-c3.z*c1.y;
   c2.y=c3.z*c1.x-c3.x*c1.z;
   c2.z=c3.x*c1.y-c3.y*c1.x;*/
   rotmat(c1,c2,c3,rot,irot);
   }

void rottopole(struct geographic a, double *rot, double *irot)
   {
   struct spherical s1, s2, s3;
   struct cartesian c1,c2,c3;
   s3=geotosph(a);
   s1.the=s3.the-HPI;
   s1.phi=s3.phi;
   s1.r=s3.r;
   if(s1.the < 0)
      {
      s1.the=-s1.the;
      s1.phi=s1.phi+PI;
      }
   c3=sphtocart(s3);
   c1=sphtocart(s1);
   c2.x=c1.y*c3.z-c1.z*c3.y;
   c2.y=c1.z*c3.x-c1.x*c3.z;
   c2.z=c1.x*c3.y-c1.y*c3.x;
   rotmat(c1,c2,c3,rot,irot);
   }

struct geographic newloc(struct geographic a, double az, double dist)
   {
   double sdis,saz,ct,st,cd,sd,cb,sb,sa,sc,stp,ctp,ath,ca,cc,cang,aang;
   struct spherical b;
   struct geographic c;
   b=geotosph(a);
   sdis=dist*CONV;
   saz=az*CONV;
   ct=cos(b.the);
   cd=cos(sdis);
   st=sin(b.the);
   sd=sin(sdis);
   cb=cos(saz);
   sb=sin(saz);
   
   ctp=ct*cd+st*sd*cb;
   b.the=acos(ctp);
   if(b.the<0.0) b.the+=PI;
   stp=sin(b.the);

   sc=st*sb/stp;
   sa=sd*sb/stp;

   if(b.the==HPI)
      {
      cc=-cd*sc*cb/sb;
      ca=-cb*cc+sb*sc*cd;
      }
   else
      {
      cc=(ct*sd-st*cd*cb)/stp;
      ca=(ctp*st-sd*cb)/(stp*ct);
      }
   cang=atan2(sc,cc);
   aang=atan2(sa,ca);
   b.phi=b.phi+aang;
   c=sphtogeo(b);
   c.dep=a.dep;
   return c;
   }
   
struct geographic sphtogeo(struct spherical a)
   {
   struct geographic b;
   b.dep=REARTH-a.r;
   b.lat=90.-a.the/CONV;
   b.lon=a.phi/CONV;
   return b;
   }
struct spherical geotosph(struct geographic a)
   {
   struct spherical b;
   b.r=REARTH-a.dep;
   b.the=(90.-a.lat)*CONV;
   b.phi=a.lon*CONV;
   return b;
   }
 struct spherical cartosph(struct cartesian a)
   {
   struct spherical b;
   b.r=sqrt(a.x*a.x+a.y*a.y+a.z*a.z)  ;
   b.phi=atan2(a.y,a.x);
   b.the=acos(a.z/b.r);
   return b;
   }
struct cartesian sphtocart(struct spherical a)
   {
   struct cartesian b;
   b.x=a.r*sin(a.the)*cos(a.phi);
   b.y=a.r*sin(a.the)*sin(a.phi);
   b.z=a.r*cos(a.the);
   return b;
   }

struct spherical cubtosph(struct cubesphere c)
   {
   struct spherical s;
   switch ( c.reg )
      {
      case 1:
         s.phi=c.xi;
         s.the=atan2(1.,cos(s.phi)*tan(c.eta));
         break;
      case 2:
/*         s.phi=atan2(-1.,tan(c.xi));*/
         s.phi=atan2(1.,-tan(c.xi));
         s.the=atan2(1.,sin(s.phi)*tan(c.eta));
         break;
      case 3:
         s.phi=atan2(-tan(c.xi),-1.);
         s.the=atan2(1.,-cos(s.phi)*tan(c.eta));
         break;
      case 4:
         s.phi=atan2(-1.,tan(c.xi));
         s.the=atan2(1.,-sin(s.phi)*tan(c.eta));
         break;
      case 5:
         s.phi=atan2(tan(c.xi),tan(-c.eta));
         s.the=atan(tan(c.xi)/sin(s.phi));
         break;
      case 6:
         s.phi=atan2(tan(c.xi),tan(c.eta));
         s.the=PI+atan(-tan(c.xi)/sin(s.phi));
         break;
      default:
         {
         printf("cubetosph: region not defined - c.reg = %d\n",c.reg);
         break;
         }
      }
/*   while(s.the < 0.)
      {
      s.the=PI-s.the;
      s.phi=s.phi+PI;
      }
   while(s.the > PI)
      {
      s.the=s.the-PI;
      s.phi=s.phi+PI;
      }
   while(s.phi > TPI) s.phi=s.phi-TPI;
   while(s.phi < 0.)  s.phi=s.phi+TPI; */
   s.r=c.r;
   return s;
   }
      
struct cubesphere cartocub(struct geographic a)
   {
   struct cubesphere c;
   struct cartesian b;
   b=sphtocart(geotosph(a));
   c.reg=csregion(b);
   switch ( c.reg )
      {
      case 1:
         c.xi=atan2(b.y,b.x);
         c.eta=atan2(b.z,b.x);
         break;
      case 2:
         c.xi=atan2(-b.x,b.y);
         c.eta=atan2(b.z,b.y);
         break;
      case 3:
         c.xi=atan2(-b.y,-b.x);
         c.eta=atan2(b.z,-b.x);
         break;
      case 4:
         c.xi=atan2(b.x,-b.y);
         c.eta=atan2(b.z,-b.y);
         break;
      case 5:
         c.xi=atan2(b.y,b.z);
         c.eta=atan2(-b.x,b.z);
         break;
      case 6:
         c.xi=atan2(b.y,-b.z);
         c.eta=atan2(b.x,-b.z);
         break;
      default:
         {
         printf("cubetosph: region not defined - c.reg = %d\n",c.reg);
         break;
         }
      }
   c.r=sqrt(b.x*b.x+b.y*b.y+b.z*b.z);
   return c;
   }

int csregion(struct cartesian b)
   {
   int reg;
   double ax,ay,az;
   ax=fabs(b.x);
   ay=fabs(b.y);
   az=fabs(b.z);
   if(ax>ay)
      {
      if(ax>az)
         {
         if(b.x > 0.0)
            {
            reg=1;
            }
         else
            {
            reg=3;
            }
         }
       else
         {
         if(b.z > 0.0)
            {
            reg=5;
            }
         else
            {
            reg=6;
            }
         }
      }
   else
      {
      if(b.y > 0.0)
         {
         reg=2;
         }
      else
         {
         reg=4;
         }
      }
   return reg;
   }
   
struct cartesian mult(struct cartesian a, double b)
   {
   struct cartesian c;
   c.x=b*a.x;
   c.y=b*a.y;
   c.z=b*a.z;
   return c;
   }
   
struct cartesian add(struct cartesian a, struct cartesian b)
   {
   struct cartesian c;
   c.x=a.x+b.x;
   c.y=a.y+b.y;
   c.z=a.z+b.z;
   return c;
   }
   
struct cartesian vecprod(struct cartesian a, struct cartesian b)
   {
   struct cartesian c;
   c.x=a.y*b.z-a.z*b.y;
   c.y=a.z*b.x-a.x*b.z;
   c.z=a.x*b.y-a.y*b.x;
   return c;
   }
   
double scaprod(struct cartesian a, struct cartesian b)
   {
   double c;
   c=a.x*b.x+a.y*b.y+a.z*b.z;
   return c;
   }
   
double norm(struct cartesian a)
   {
   double c;
   c=sqrt(scaprod(a,a));
   return c;
   }

int inrect(struct geographic a, struct geographic UL, struct geographic LR)
   {
   int c=0;
   if(UL.lon < 0.0) 
      {
      if(a.lon > 180.) a.lon-=360.;
      }
   else
      {
      if(a.lon < 0.) a.lon+=360.;
      }
   if(a.lon < LR.lon && a.lon > UL.lon && a.lat < UL.lat && a.lat > LR.lat) c=1;
   return c;
   }

int inpoly(struct geographic pt, struct geographic *vrtx, int nvert)
   {
   int i, j, c = 0,sign=1;
   for (i = 0;  i < nvert; i++) 
      {
      if(vrtx[i].lon < 0.0) sign=-1; 
      }
   if(sign==-1)
      {
      if(pt.lon > 180.) pt.lon=-360.;
      }
   else
      {
      if(pt.lon < 0.0) pt.lon+=360.;
      }
      
   
   for (i = 0, j = nvert-1; i < nvert; j = i++) 
      {
      if ( ((vrtx[i].lat>pt.lat) != (vrtx[j].lat>pt.lat)) &&
	    (pt.lon< (vrtx[j].lon-vrtx[i].lon) * (pt.lat-vrtx[i].lat) / 
            (vrtx[j].lat-vrtx[i].lat) + vrtx[i].lon) )
       c = !c;
      }
   return c;
   }

double distoline(struct geographic g,struct geographic t1, struct geographic t2,struct geographic *x,int *iplp)
   {
   struct cartesian h,a,b,p,d,p1,p2;
   double dis1,dis2,gdis;
   double t0,normb;

   h=sphtocart(geotosph(g));
   p1=sphtocart(geotosph(t1));
   p2=sphtocart(geotosph(t2));

   /* -- line x-a+tb -- */
   a=p1;
   b=add(p2,mult(p1,-1.));
   normb=norm(b);
   b=mult(b,1./normb);

   /* -- plane n.x=d -- */
   t0=scaprod(h,b)-scaprod(a,b);
   if(t0 >= 0.0 && t0 <= normb)
      {
      p=add(a,mult(b,t0));
      d=add(h,mult(p,-1.0));
      gdis=norm(d);
      *iplp=2;
      *x=sphtogeo(cartosph(p));
      }
   else
      {
      dis1=norm(add(h,mult(p1,-1.0)));
      dis2=norm(add(h,mult(p2,-1.0)));
      *iplp=3;
      if(dis1 <= dis2)
         {
         gdis=dis1;
         x=&t1;
         }
      else
         {
         gdis=dis2;
         x=&t2;
         }
      }
   return gdis;
   }

double distoplane(struct geographic g,struct geographic t1,struct geographic t2,struct geographic t3,struct geographic *x,int *iplp,double *strike,double *dip)
   {
   double gdis,gdis1,gdis2,gdis3;
   double ss1,ss2,ss3,dn,dp;
   struct geographic x1,x2,x3;
   struct cartesian r,s,s1,s2,s3,d,d1,d2,d3;
   struct cartesian h,north,n,p1,p2,p3,p;
   int iplp1,iplp2,iplp3;
   gdis=1.e24;
   /* -- convert to spherical coordinates */
   h=sphtocart(geotosph(g));
   p1=sphtocart(geotosph(t1));
   p2=sphtocart(geotosph(t2));
   p3=sphtocart(geotosph(t3));
   /* -- average location */
   r=mult(add(p1,add(p2,p3)),1./3.);
   /* -- setup three vectors spanning the plane */
   d1=add(p2,mult(p1,-1.0));
   d2=add(p3,mult(p2,-1.0));
   d3=add(p1,mult(p3,-1.0));
   /* -- normal to the plane and origin (n.x=d) */
   n=vecprod(d1,d2);
   n=mult(n,1./norm(n));
   /* -- while we're at it, also determine the strike/dip of the plane */
   /* -- nx, ny, nz is normal to the plane */
   /* -- get dip -- */
   *dip = acos(scaprod(r,n)/(norm(r)*norm(n)));
   *dip = *dip/CONV;
   if(*dip > 90.) *dip-=180;
   /* -- get strike -- */
   /* -- first get local north vector -- */
   if(r.z == 0.0)
      {
      north.z=1.0;
      north.x=0.0;
      north.y=0.0;
      }
   else
      {
      north.x=-r.x;
      north.y=-r.y;
      north.z=(r.x*r.x+r.y*r.y)/r.z;
      north=mult(north,r.z/fabs(r.z));
      }

   /* -- get the strike vector -- */
   /* -- sx, sy, sz is strike vector in global coords */
   s=vecprod(r,n);
   s=mult(s,1./norm(s));
   /* -- use inner product to get strike angle -- */
   *strike=acos(scaprod(north,s)/(norm(north)*norm(s)));
   *strike=*strike/CONV-180.;
   /* -- now get the old strik/dip/slip in the new coordinate system */
   /* -- n = updip, e = along strike, z = radial outward */
   /* -- orientation in local cartesian */

   /* -- back to our main problem */
   dn=scaprod(n,p1);
   if(dn < 0.0)
      {
      n=mult(n,-1.0);
      dn=-dn;
      }
   /* -- project point on plane */
   dp=scaprod(h,n);
   dp=dn-dp;
   p=add(h,mult(n,dp));
   /* -- determine whether point in triangle, all vector products should have */
   /* -- same sign */
   s1=vecprod(add(p,mult(p1,-1.)),d1);
   s2=vecprod(add(p,mult(p2,-1.)),d2);
   s3=vecprod(add(p,mult(p3,-1.)),d3);
   ss1=scaprod(s1,s2);
   ss2=scaprod(s2,s3);
   ss3=scaprod(s3,s1);

   if (ss1 > 0.0 && ss2 > 0.0 && ss3 > 0.0)
      {
      gdis=fabs(dp);
      *x=sphtogeo(cartosph(p));
      *iplp=1;
      }
   else
      {
      gdis1 = distoline(g,t1,t2,&x1,&iplp1);
      gdis2 = distoline(g,t2,t3,&x2,&iplp2);
      gdis3 = distoline(g,t3,t1,&x3,&iplp3);
      /*printf("%9.4f %9.4f %9.4f\n",gdis1,gdis2,gdis3);*/
      if(gdis2 < gdis1)
         {
         if(gdis2 < gdis3)
            {
            gdis=gdis2;
            *x=x2;
            *iplp=iplp2;
            }
         else
            {
            gdis=gdis3;
            *x=x3;
            *iplp=iplp3;
            }
         }
      else if(gdis1 < gdis3)
         {
         gdis=gdis1;
         *x=x1;
         *iplp=iplp1;
         }
      else
         {
         gdis=gdis3;
         *x=x3;
         *iplp=iplp3;
         }
      }
   return gdis;
   }

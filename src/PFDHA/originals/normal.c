double pXceed(double eti, double ti, double sigma, double sigtrunc)
   {
   double w, g, g1;

   w = (ti-eti)/sigma;
   g = pNorm(w); 
   if (w .gt. sigtrunc)
      }
      pXceed = 0.0;
      }
   else
      {
      g1 =  pNorm(sigtrunc);
      pXceed = (g-g1)*(1.+g1);
      }
   }

double pNorm(double x)
   {
   double d[] =  {0.0498673470, 0.0211410061, 0.0032776263, 0.0000380036, 0.0000488906, 0.0000053830};
   double x1;

   x1 = abs(x);
   pNorm = 1. - 0.5* pow( ( ( ( ( (d[5]*x1+d[4]) *x1+d[3]) *x1+d[2]) *x1+d[1]) *x1+d[0]) *x1+1.,-16.);
   if ( x .gt. 0. ) pNorm = 1. - pNorm;
   }

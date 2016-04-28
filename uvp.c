#include "uvp.h"
#include "helper.c"
#include <stdlib.h>

void calculate_dt(double Re,double tau,double dt,double dx,double dy,int imax,int jmax,double **U,double **V)
{
if (tau>0)
  {
    double umax = 0;
    double vmax = 0;
    for (int i=1;i<=imax;i++) 
    {
      for (int j=1;i<=jmax;j++)
      {
         if (abs(U[i][j]) > umax)
            {
               umax = abs(U[i][j]);
            }
         if (abs(V[i][j]) > vmax)
            {
               vmax = abs(V[i][j]);
            }
       }
    }
    dt = tau*fmin(fmin(dx/umax,dy/vmax),Re/2*(1/dx^2+1/dy^2)^-1);
   }
}

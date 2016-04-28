#include "uvp.h"
#include "helper.c"
#include <stdlib.h>

/* calculate the tine stepping based on the value of tau*/

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


/* Calculating right hand size of the pressure equation */

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{
for (int i=0;i<imax;i++)
   {
     for (int j=0;j<jmax;j++)
      {
         RS[i][j] = 1/dt*((F[i][j] - F[i-1][j])/dx + (G[i][j] - G[i][j-1])/dy);
      }
   }
}


/* updating the values of U and V */

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
for (int i=1;i<=imax-1;i++)
  {
    for (int j =1;j<=jmax;j++)
      {
        U[i][j] = F[i][j] - dt/dx*(P[i+1][j] - P[i][j]);
      }
   }

for (int i=1;i<=imax;i++)
  {
    for (int j =1;j<=jmax-1;j++)
      {
        V[i][j] = G[i][j] - dt/dy*(P[i][j+1] - P[i][j]);
      }
  }
}





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

void calculate_rs(double dt,double dx,double dy,  int imax,  int jmax,  double **F,  double **G,  double **RS)
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

void calculate_uv(double dt,  double dx,  double dy,  int imax,  int jmax,double **U,  double **V,  double **F, double **G,double **P)
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

/* calculating F and G*/

double du2dx(int i, int j)
  {
  return 1/(4*dx)*(((U[i][j]+U[i+1][j])^2-(U[i-1][j]+U[i][j])^2)+alpha/dy*(abs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j])-abs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j])));
  }

double duvdy(int i, int j)
  {
  return 1/(4*dy)*((V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))+alpha/dy*(abs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-abs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])));
  }

double d2udx2(int i, int j)
  {
  return (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx^2);
  }

double d2udy2(int i, int j)
  {
  return (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy^2);
  }

double d2vdx2(int i, int j)
  {
  return (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx^2);
  }

double d2vdy2(int i, int j)
  {
  return (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy^2);
  }

double dv2dy(int i, int j)
  {
  return (1/dy)*(((V[i][j] + V[i][j+1])/2)^2 - ((V[i][j-1] + V[i][j])/2)^2) + (alpha/dy)*(abs(V[i][j] + V[i][j+1])/2*(V[i][j] - V[i][j+1])/2 - abs(V[i][j-1] + V[i][j])/2*(V[i][j-1] - V[i][j])/2);
  }

double duvdx(int i, int j)
  {
  return (1/dx)*((U[i][j] + U[i][j+1])/2*(V[i][j] + V[i+1][j])/2 - (U[i-1][j] + U[i-1][j+1])/2*(V[i-1][j] + V[i][j])/2) + alpha/dx*(abs(U[i][j] + U[i][j+1])/2*(V[i][j] - V[i+1][j])/2 - abs(U[i-1][j] + U[i-1][j+1])/2*(V[i-1][j] - V[i][j])/2);
   }

// calculate_fg function
 
void calculate_fg(double Re,double GX,double GY,double alpha,double dt,double dx,  double dy,int imax,int jmax,  double **U,  double **V,  double **F,  double **G)
  {
   for (int i = 1; i <= imax-1; ++i)
  {
    for (int j = 1; j <= jmax; ++j)
    {
      F[i][j] = U + dt*(1/Re*(d2udx2(i,i)+d2udy2(i,j)-du2dx(i,j)-duvdy(i,j))+GX);
    }
  }
  for (int i =1;i<=imax;i++)
  {
    for (int j =1;j<=jmax-1;j++)
    {
      G[i][j] = V[i][j] + dt*(1/Re*(d2vdx2(i,j) + d2vdy2(i,j) - duvdx(i,j) - dv2dy(i,j) + GY);
    }
   }
  }



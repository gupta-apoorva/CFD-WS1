#include "boundary_val.h"
void boundaryvalues(int imax,int jmax,double **U,double **V, double** P, double** G, double** F)
{
	for(int i=1; i<=imax; i++)
	{
		U[i][0] = -U[i][1];
		U[i][jmax+1] = U[i][jmax];
                V[i][0] = -V[i][1];
		V[i][jmax+1] = -V[i][jmax];
                P[i][0] = P[i][1];
		P[i][jmax+1] = P[i][jmax];
                G[i][0] = V[i][0];
                G[i][jmax] = V[i][jmax]; 
                
	}
	for(int j=1; j<=jmax; j++)
	{
		U[0][j] = -U[1][j];
		U[imax+1][j] = -U[imax][j];
                V[0][j] = -V[1][j];
		V[imax+1][j] = -V[imax][j];
                P[0][j] = P[1][j];
		P[imax+1][j] = P[imax][j];
                F[0][j] = U[0][j];
                F[imax][j] = U[imax][j];
	
}
}

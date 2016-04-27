void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
)
{
	for(int i=1; i<=imax; i++)
	{
		U[i][0] = -U[i][1];
		U[i][jmax+i] = U[i][jmax];
	}
	for(int j=1; j<=jmax; j++)
	{
		U[0][j] = -U[1][j];
		U[imax+1][j] = -U[imax][j];
	}
	for(int i=1; i<=imax; i++)
	{
		V[i][0] = -V[i][1];
		V[i][jmax+i] = -V[i][jmax];
	}
	for(int j=1; j<=jmax; j++)
	{
		V[0][j] = -V[1][j];
		V[imax+1][j] = -V[imax][j];
	}




}

#include <mpi.h>
#include "mpi_residual.h"
#include "mpi_etc.h"
#include <bits/stdc++.h>
#define corner 5

using namespace std;

extern int rank, size, n_expanded;

double residual(double *a, double *aInv, int n) {

	int mem = n * n_expanded / size;

	MPI_Bcast(aInv, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (rank == 0)
	{
		print(a, n_expanded, n);
	}


	double *a_part;
	a_part = new double[mem];

	double *res_part;
	res_part = new double[mem];
	
	double *res_matrix;

	if (rank == 0)
	{	
		res_matrix = new double[n_expanded * n];
		memset(res_matrix, 0, n_expanded * n * sizeof(double));
	}
	
//	trans(aInv, n_expanded, n);

	MPI_Scatter(a, mem, MPI_DOUBLE,
    a_part, mem, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int j = 0; j < mem / n; j++)
		for (int i = 0; i < n; i++) {
			res_part[i + j * n] = 0;
			for (int k = 0; k < n; k++)
				res_part[i + j * n] += a[k + i * n] * aInv[j + k * n];
		}	

//	mult(a_part, aInv, res_part, mem / n, n);
	
//	print(res_part, n, m);
	
	MPI_Gather(res_part, mem, MPI_DOUBLE,
    res_matrix, mem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
	if (rank == 0)
	{
		print(res_matrix, n, n);
 		for (int i = 0; i < n; i++)
			res_matrix[i * n + i] -= 1;
		double res = 0;
		for (int i = 0; i < n; i++) {
			double preRes = 0;
			for (int j = 0; j < n; j++)
				preRes += fabs(res_matrix[i * n + j]);
			res = max(res, preRes);
		}
		delete[] res_matrix;
		return res;
	}
	return 0;
}

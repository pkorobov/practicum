#include <mpi.h>
#include "mpi_residual.h"
#include "mpi_etc.h"
#include <bits/stdc++.h>
#define corner 5

using namespace std;

extern int rank, size, n_expanded;

double residual(double *a, double *aInv, int n) {

	int mem = n * n_expanded / size;

	double *a_part = NULL;
	a_part = new double[mem];

	double *res_part = NULL;
	res_part = new double[mem];
	
	double *res_matrix = NULL;

	MPI_Bcast(aInv, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{	
		res_matrix = new double[n_expanded * n];
		memset(res_matrix, 0, n_expanded * n * sizeof(double));
	}

	trans(aInv, n_expanded, n);
//	trans(a, n_expanded, n);

	MPI_Scatter(a, mem, MPI_DOUBLE,
    a_part, mem, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = 0; i < mem / n; i++)
		for (int j = 0; j < n; j++) {
			res_part[i * n + j] = 0;
			for (int k = 0; k < n; k++)
				res_part[i * n + j] += a_part[i * n + k] * aInv[j * n + k];
		}	

//	mult(a_part, aInv, res_part, mem / n, n);
	
//	if (rank == 1)
//		print(res_part, mem / n, n);
	
	MPI_Gather(res_part, mem, MPI_DOUBLE,
    res_matrix, mem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
	if (rank == 0)
	{
//		print(res_matrix, n_expanded, n);
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

#include <mpi.h>
#include "mpi_residual.h"
#include "mpi_etc.h"
#include <bits/stdc++.h>
#define corner 5

using namespace std;

extern int rank, size, n_expanded;

double residual(double *a, double *aInv, int n) {

	int s = n_expanded / size;

	MPI_Bcast(aInv, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double *a_part;
	a_part = new double[n * s];

	double *res_part;
	res_part = new double[n * s];
	
	double *res_matrix;

	if (rank == 0)
	{	
		res_matrix = new double[n_expanded * n];
		memset(res_matrix, 0, n_expanded * n * sizeof(double));
	}
	
//	trans(aInv, n_expanded, n);

	MPI_Scatter(a, n * s, MPI_DOUBLE,
    a_part, n * s, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	trans(aInv, n, n);

	for (int j = 0; j < s; j++)
		for (int i = 0; i < n; i++) {
			res_part[i + j * n] = 0;
			for (int k = 0; k < n; k++)
				res_part[i + j * n] += a[k + i * n] * aInv[k + (j + s * rank) * n];
		}	
	
	MPI_Gather(res_part, n * s, MPI_DOUBLE,
    res_matrix, n * s, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
	if (rank == 0)
	{
//		print(res_matrix, n, n);
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

	delete[] a_part;
	delete[] res_part;

	return 0;
}

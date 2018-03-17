#include <mpi.h>
#include <bits/stdc++.h>
#include "mpi_etc.h"
#include "mpi_algo.h"
#include "mpi_residual.h"

#define corner 5
#define eps 1e-8

using namespace std;

int rank, size, n_expanded;

int main(int argc, char** argv)
{
	
	int n = 2000;
	double *a = NULL, *aInv = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	n_expanded = n % size == 0 ? n : n + size - n % size;
	a = new double[n * n_expanded];
	aInv = new double[n * n_expanded];
	memset(a, 0, n * n_expanded * sizeof(double));
	memset(aInv, 0, n * n_expanded * sizeof(double));
	
	generate(a, n, formula3);
	if (rank == 0){
		cout << "A" << endl;
//		print(a, n_expanded, n);
	}
	
	double t1, t2, dt;
	t1 = MPI_Wtime();

	if (solve(a, aInv, n) < 0)
		return -1;
	
	t2 = MPI_Wtime();
	dt = t2 - t1;
	if (rank == 0)
		cout << "Time: " << dt << " seconds" << endl;
	if (rank == 0){
		cout << "A^-1" << endl;
//		print(aInv, n, n);
	}	

	generate(a, n, formula3);

	double r = residual(a, aInv, n);
		if (rank == 0){
		cout << "Residual: " << r << endl;
	}

	delete[] a;
	delete[] aInv;
	MPI_Finalize();
	return 0;
}

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
	
	int n = 3;
	double *a = NULL, *aInv = NULL;
//	double a[] = {1, 2, 3, 4};
//	double aInv[] = {-2, 1, 1.5, -0.5};

//	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	n_expanded = n % size == 0 ? n : n + size - n % size;
	a = new double[n * n_expanded];
	aInv = new double[n * n_expanded];
	memset(a, 0, n * n_expanded * sizeof(double));
	memset(aInv, 0, n * n_expanded * sizeof(double));
	
	generate(a, n, formula2);
	generate(aInv, n, formula2);
	if (rank == 0){
		cout << "A" << endl;
		print(a, n, n);
	}
	solve(a, aInv, n);
	if (rank == 0){
		cout << "A^-1" << endl;
		print(aInv, n, n);
	}
	double r = residual(a, aInv, n);
	if (rank == 0){
		printf("%f\n", r);
	}
	delete[] a;
	delete[] aInv;
	MPI_Finalize();
	return 0;
}

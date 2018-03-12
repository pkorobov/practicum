#include <mpi.h>
#include <bits/stdc++.h>
#include "mpi_algo.h"
#define corner 5
#define eps 1e-8


using namespace std;

static double *Cos, *Sin;


int solve(double *a, double *aInv, int n) { 	
	
	
	double *q = NULL, *rInv = NULL, *a_part = NULL, *q_part = NULL;

	q = new double[n * n_expanded];
	rInv = new double[n * n];
	a_part = new double[n * n_expanded / size];
	q_part = new double[n * n_expanded / size];
	Cos = new double[n - 1];
	Sin = new double[n - 1];

   	memset(q, 0, n * n * sizeof(double));
	for (int i = 0; i < n; i++)
		q[i * n + i] = 1;		
	
	for (int k = 0; k < n - 1; k++) { 
	
		if (rank == 0) {				
		
			for (int l = k + 1; l < n; l++)	{	
				
				double x = a[k * n + k], y = a[k * n + l];
//				cout << x << " ! " << y << endl;
				double mod = sqrt(x * x + y * y);
				if (fabs(x) < eps && fabs(y) < eps) {
					continue;
				}

				Cos[l - k - 1] = x / mod;
				Sin[l - k - 1] = -y / mod;
				
//				cout << Cos[l - k - 1] << " " << Sin[l - k - 1] << endl;
				
//				a[k * n + k] = mod;
//				a[l * n + k] = 0;					
			}
		}

		MPI_Bcast(Cos, n - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(Sin, n - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		MPI_Scatter(a, n * n_expanded / size, MPI_DOUBLE, a_part, 
		n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Scatter(q, n * n_expanded / size, MPI_DOUBLE, q_part, 
		n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int m = 0; m < n_expanded / size; m++) {	

			for (int l = k + 1; l < n; l++)	{		
//			for (int m = from[thread_number]; m < to[thread_number]; m++) {	

				double temp_k, temp_l;
				if (m >= k) {
					temp_k = Cos[l - k - 1] * a_part[m * n + k] - Sin[l - k - 1] * a_part[m * n + l];
					temp_l = Sin[l - k - 1] * a_part[m * n + k] + Cos[l - k - 1] * a_part[m * n + l];
					a_part[m * n + k] = temp_k;
					a_part[m * n + l] = temp_l;
				}			
				temp_k = Cos[l - k - 1] * q_part[m * n + k] - Sin[l - k - 1] * q_part[m * n + l];
				temp_l = Sin[l - k - 1] * q_part[m * n + k] + Cos[l - k - 1] * q_part[m * n + l];
				q_part[m * n + k] = temp_k;
				q_part[m * n + l] = temp_l;
			}		
			
			if (rank == 0)
				print(a_part, n, n_expanded);
		}

		MPI_Gather(a_part, n * n_expanded / size, MPI_DOUBLE, a, 
		n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Gather(q_part, n * n_expanded / size, MPI_DOUBLE, q, 
		n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0 && fabs(a[k * n + k]) < eps) {
			printf("Error, degenerate matrix\n");
			return -1;
		}
	
	}
	
	if (rank == 0) {
		cout << "R" << endl;
		print(a, n, n);
		cout << "Q" << endl;
		print(q, n, n);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	for (int j = rank; j < n; j += size)
		for (int i = n - 1; i >= 0; i--) {	
			rInv[i * n + j] = (i == j ? 1 : 0);
			for (int k = i + 1; k < n; k++)
				rInv[i * n + j] -=  a[i * n + k] * rInv[k * n + j];
			rInv[i * n + j] /= a[i * n + i];
		}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0)
//		trans(q, n, n);
		mult(rInv, q, aInv, n, n);
	
	if (rank == 0) {
//		print(a, n, n);
//		print(q, n, n);
		delete[] rInv;
	}
	return 0;
}

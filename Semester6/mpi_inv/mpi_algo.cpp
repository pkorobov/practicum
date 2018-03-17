#include <mpi.h>
#include <bits/stdc++.h>
#include "mpi_algo.h"
#include "mpi_etc.h"
#define corner 5
#define eps 1e-8
#define  I(a, x, y) a[x + y * n]
#define  I_t(a, x, y) a[x * n + y]

using namespace std;

static double *Cos, *Sin;
extern int rank, size, n_expanded;


int solve(double *a, double *aInv, int n) { 	
	
	double *q = NULL, *a_part = NULL, *q_part = NULL;

	q = new double[n * n_expanded];
	a_part = new double[n * n_expanded / size];
	q_part = new double[n * n_expanded / size];
	Cos = new double[n - 1];
	Sin = new double[n - 1];

   	memset(q, 0, n * n * sizeof(double));
	for (int i = 0; i < n; i++)
		I(q, i, i) = 1;		
	
	for (int k = 0; k < n - 1; k++) { 
	
		if (rank == 0) {				
		
			for (int l = k + 1; l < n; l++)	{	
				
				double x = I(a, k, k), y = I(a, l, k);
//				cout << x << " ! " << y << endl;
				double mod = sqrt(x * x + y * y);
				if (fabs(x) < eps && fabs(y) < eps) {
					continue;
				}

				Cos[l - k - 1] = x / mod;
				Sin[l - k - 1] = -y / mod;
				
//				cout << Cos[l - k - 1] << " " << Sin[l - k - 1] << endl;
				
				I(a, k, k) = mod;
				I(a, l, k) = 0;					
			}
		}

		MPI_Bcast(Cos, n - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(Sin, n - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		MPI_Scatter(a, n * n_expanded / size, MPI_DOUBLE, a_part, 
		n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Scatter(q, n * n_expanded / size, MPI_DOUBLE, q_part, 
		n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		/*
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0) {
			printf("0:\n");
			print(a_part, n_expanded / size, n);
//			print(a_part, n, n_expanded / size);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 1) {
			printf("1:\n");
			print(a_part, n_expanded / size, n);
//			print(a_part, n, n_expanded / size);
		}
//		return -2;
*/
//		MPI_Barrier(MPI_COMM_WORLD);
		
		for (int m = 0; m < n_expanded / size; m++) {	
			for (int l = k + 1; l < n; l++)	{		
				double temp_k, temp_l;
				if (m + n_expanded / size * rank > k) {
					temp_k = Cos[l - k - 1] * I(a_part, k, m) - Sin[l - k - 1] * I(a_part, l, m);
					temp_l = Sin[l - k - 1] * I(a_part, k, m) + Cos[l - k - 1] * I(a_part, l, m);
					I(a_part, k, m) = temp_k;
					I(a_part, l, m) = temp_l;
				}			
				temp_k = Cos[l - k - 1] * I(q_part, k, m) - Sin[l - k - 1] * I(q_part, l, m);
				temp_l = Sin[l - k - 1] * I(q_part, k, m) + Cos[l - k - 1] * I(q_part, l, m);
				I(q_part, k, m) = temp_k;
				I(q_part, l, m) = temp_l;
/*				if (rank == 1) {
					printf("R>>>>>>>>>>>>>>\n");
					print(a_part, n_expanded / size, n);
					printf("-------------:\n");
	
					printf("Q>>>>>>>>>>>>>>\n");
					print(q_part, n_expanded / size, n);
					printf("-------------:\n");
				}
*/
//				print(a_part, n, n_expanded / size);
		}

		}


		MPI_Gather(a_part, n * n_expanded / size, MPI_DOUBLE, a, 
		n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Gather(q_part, n * n_expanded / size, MPI_DOUBLE, q, 
		n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				
        if (fabs(a[k * n + k]) < eps) {
//			print(a, n, n);
			printf("%d: Error, degenerate matrix\n", rank);
			return -1;
		}
	}

	MPI_Bcast(a, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(q, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

/*	
	if (rank == 0) {
		cout << "R" << endl;
		print(a, n, n);
		cout << "Q" << endl;
		print(q, n, n);
	}
*/
	double *aInv_part = NULL;
	aInv_part = new double[n * n_expanded / size];

	MPI_Scatter(aInv, n * n_expanded / size, MPI_DOUBLE, aInv_part, 
	n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);	
	
	int s = n_expanded / size;
	for (int j = 0; j < s && j + s * rank < n; j++)
		for (int i = n - 1; i >= 0; i--) {	
			if (i == n - 1)
			{
//				print(q, n, n);
//				printf("%f ------ %d ------- %d -------- %d\n", q[i + n * (j + s * rank)], j, rank, rank * s);
			}
			I(aInv_part, i, j) = q[i + n * (j + s * rank)]; // I(q, i, j + rank * s) всё ломает
//			printf(">>>>>>>\n%d\n-------\n", j + s * rank);
			for (int k = i + 1; k < n; k++)
				I(aInv_part, i, j) -=  I(a, i, k) * I(aInv_part, k, j);
			I(aInv_part, i, j) /= I(a, i, i);
		}

	MPI_Gather(aInv_part, n * n_expanded / size, MPI_DOUBLE, aInv, 
	n * n_expanded / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
	delete[] q;	
	return 0;
}

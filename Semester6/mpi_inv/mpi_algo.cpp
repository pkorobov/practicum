#include <mpi.h>
#include <bits/stdc++.h>
#include "mpi_algo.h"
#include "mpi_etc.h"
#define corner 5
#define eps 1e-8
#define  I(a, x, y) a[x + y * n]
#define  I_t(a, x, y) a[y * n + x]

using namespace std;

static double *Cos, *Sin;
extern int rank, size, n_expanded;


int solve(double *a, double *aInv, int n) { 	
	
	
	double *q = NULL, *rInv = NULL, *a_part = NULL, *aInv_part = NULL, *q_part = NULL;


	q = new double[n * n_expanded];
	aInv_part = new double[n * n_expanded / size];
	a_part = new double[n * n_expanded / size];
	q_part = new double[n * n_expanded / size];
	Cos = new double[n - 1];
	Sin = new double[n - 1];

	int s = n_expanded / size;

   	memset(q_part, 0, n * s * sizeof(double));
	for (int i = 0; i < s; i++)
		I(q_part, i + rank * s, i) = 1;		
//		q_part[i + rank * s + i * n] = 1;		
		
	double t1, t2, dt, ta, tb, tsum = 0;
	
	t1 = MPI_Wtime();
	
	ta = MPI_Wtime();

	MPI_Scatter(a, n * s, MPI_DOUBLE, a_part, 
	n * s, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	tb = MPI_Wtime();
	tsum += tb - ta;
	
	for (int k = 0; k < n - 1; k++) { 
		int root = -1;
		
		for (int t = 0; t < size; t++)
			if (t * s <= k && k < (t + 1) * s) 	
				root = t;

		if (rank == root) {				
			for (int l = k + 1; l < n; l++)	{	
				
				double x = I(a_part, k, k % s), y = I(a_part, l, k % s);
				double mod = sqrt(x * x + y * y);
				if (fabs(x) < eps && fabs(y) < eps) {
					continue;
				}

				Cos[l - k - 1] = x / mod;
				Sin[l - k - 1] = -y / mod;
				
				I(a_part, k, k % s) = mod;
				I(a_part, l, k % s) = 0;					
			}
		}
	
		MPI_Bcast(Cos, n - 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
		MPI_Bcast(Sin, n - 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

		for (int m = 0; m < s; m++) {	
			for (int l = k + 1; l < n; l++)	{		
				double temp_k, temp_l;
				if (m + s * rank > k) {
					temp_k = Cos[l - k - 1] * I(a_part, k, m) - Sin[l - k - 1] * I(a_part, l, m);
					temp_l = Sin[l - k - 1] * I(a_part, k, m) + Cos[l - k - 1] * I(a_part, l, m);
					I(a_part, k, m) = temp_k;
					I(a_part, l, m) = temp_l;
				}			
				temp_k = Cos[l - k - 1] * I(q_part, k, m) - Sin[l - k - 1] * I(q_part, l, m);
				temp_l = Sin[l - k - 1] * I(q_part, k, m) + Cos[l - k - 1] * I(q_part, l, m);
				I(q_part, k, m)  = temp_k;
				I(q_part, l, m)  = temp_l;
			}		
		}
	}

	t2 = MPI_Wtime();

	if (rank == 0)
	{
		dt = t2 - t1;
		cout << "Rotation time: " << dt << " seconds";
	}

	ta = MPI_Wtime();

	MPI_Allgather(a_part, n * s, MPI_DOUBLE, a, 
	n * s, MPI_DOUBLE, MPI_COMM_WORLD);

	MPI_Allgather(q_part, n * s, MPI_DOUBLE, q, 
	n * s, MPI_DOUBLE, MPI_COMM_WORLD);

	tb = MPI_Wtime();
	tsum += tb - ta;

	if (rank == 0) {
		cout << "R" << endl;
		print(a, n, n);
		cout << "Q" << endl;
		print(q, n, n);
	}

	t1 = MPI_Wtime();

	for (int j = 0; j < s && j + s * rank < n; j++)
		for (int i = n - 1; i >= 0; i--) {	
			I(aInv_part, i, j) = q[i + (j + s * rank) * n];//I(q, i, j);
			for (int k = i + 1; k < n; k++)
				I(aInv_part, i, j) -=  I(a, i, k) * I(aInv_part, k, j);
			I(aInv_part, i, j) /= I(a, i, i);
		}

	t2 = MPI_Wtime();
	if (rank == 0)
	{
		dt = t2 - t1;
		cout << "Reverse gauss time: " << dt << " seconds";
	}

	ta = MPI_Wtime();
		
	MPI_Gather(aInv_part, n * s, MPI_DOUBLE, aInv, 
	n * s, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	tb = MPI_Wtime();
	tsum += tb - ta;
	
	if (rank == 0)
		cout << "Mesages time: " << tsum << endl;
	delete[] q;
	delete[] q_part;
	delete[] a_part;
	delete[] aInv_part;
	delete[] Sin;
	delete[] Cos;
	return 0;
}

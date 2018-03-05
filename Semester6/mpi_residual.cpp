#include <mpi.h>
#include <bits/stdc++.h>
#define corner 5

using namespace std;

void multMatVec(double *m, double *vec, double *res, int n);
double residual(double *a, double *aInv, int n);
void mult(double *a, double *b, double *c, int n, int m);
void trans(double *a, int n, int m);
void print(double *a, int n, int m);
void generate(double *a, int n, double formula(int, int));
double formula1(int i, int j);
double formula2(int i, int j);
double formula3(int i, int j);

int rank, size, n_expanded;

int main(int argc, char** argv)
{
	
	int n = 9;
	double *a = NULL, *aInv = NULL;
//	double a[] = {1, 2, 3, 4};
//	double aInv[] = {-2, 1, 1.5, -0.5};

//	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	n_expanded = n + size - n % size;
	a = new double[n * n_expanded];
	aInv = new double[n * n_expanded];
	memset(a, 0, n * n_expanded * sizeof(double));
	memset(aInv, 0, n * n_expanded * sizeof(double));
	
	generate(a, n, formula2);
	generate(aInv, n, formula2);
	
	double r = residual(a, aInv, n);
	if (rank == 0)
		printf("%f\n", r);
	MPI_Finalize();
	return 0;
}


double residual(double *a, double *aInv, int n) {

	int mem = n * n_expanded / size;

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

	mult(a_part, aInv, res_part, mem / n, n);
	
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

void mult(double *a, double *b, double *c, int m, int n) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++) {
			c[i * n + j] = 0;
			for (int k = 0; k < n; k++)
				c[i * n + j] += a[i * n + k] * b[k * n + j];
		}	
}


void trans(double *a, int m, int n)
{
        for (int i = 0; i < m; i++)
                for (int j = i + 1; j < n; j++)
                {
                        double t = a[i * n + j];
                        a[i * n + j] = a[j * n + i];
                        a[j * n + i] = t;
                }
}

void print(double *a, int m, int n)
{
	printf("\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%-7.2f ", a[i * n + j]);
		printf("\n");	
	}
}

double formula1(int i, int j) {
	return 1.0 / (i + j + 1);
}

double formula2(int i, int j) {
	return fabs(i - j);
}

double formula3(int i, int j) {
	return max(i, j) + 1;
}

void generate(double *a, int n, double formula(int, int)) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i * n + j] = formula(i, j);
}


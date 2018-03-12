#include <mpi.h>
#include <bits/stdc++.h>
#define corner 5
#define eps 1e-8


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
int solve(double *a, double *aInv, int n);

int rank, size, n_expanded;
double *Cos, *Sin;

int main(int argc, char** argv)
{
	
	int n = 5;
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


double residual(double *a, double *aInv, int n) {

	int mem = n * n_expanded / size;

	if (rank == 0)
	{
//		print(a, n_expanded, n);
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

    delete[] a_part;
    delete[] res_part;
    
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
			printf("%-7.2f ", a[j * n + i]);
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


#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>

#define corner 5
#define eps 1e-16

using namespace std;

int solve(double *a, double *aInv, int n);
void generate(double *a, int n, double formula(int, int));
double formula1(int i, int j);
double formula2(int i, int j);
double formula3(int i, int j);
double residual(double *a, double *aInv, int n);
void mult(double *a, double *b, double *c, int n);
void trans(double *a, int n);
void print(double *a, int n);

int main(int argc, char** argv)
{  
	int n;
	double *a = NULL, *aInv = NULL;
	if (argc >= 2) {
		if (!strcmp(argv[1], "-formula")) {
			char *t;
			n = strtol(argv[3], &t, 10);
			a = new double[n * n];
			char f = argv[2][0];
			switch(f) {
				case '1' : generate(a, n, formula1); break;
				case '2' : generate(a, n, formula2); break;
				case '3' : generate(a, n, formula3); break;
//				case '4' : generate(a, n, formula4); break;
			}
		}
		else if (!strcmp(argv[1], "-file")) {
			char name[100];
			strncpy(name, argv[2], 100);
  			freopen(name, "r", stdin); 	
			cin >> n;
			cout << n;
			a = new double[n * n];
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					cin >> a[i * n + j];
		}
		else {
			cout << "Incorrect input" << endl;
			return 0;
		}
	}
	else if (argc == 1) {
		cin >> n;
		a = new double[n * n];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				cin >> a[i * n + j];
	}
	else {
		cout << "Incorrect input" << endl;
		return 0;
	}
	aInv = new double[n * n];
	printf("\nA\n");
	print(a, n);
	printf("\nA^-1\n");
	time_t start = clock();
	if (solve(a, aInv, n) == -1)
		return -2;
	time_t end = clock();
	double time = (end - start) / CLOCKS_PER_SEC; 
 	print(aInv, n);
	
	generate(a, n, formula2);
	printf("%e\n", residual(a, aInv, n));
	printf("Time: %dm, %ds:\n", (int)time / 60, (int)time % 60);
	delete[] a;
	delete[] aInv;
	return 0;
}

int solve(double *a, double *aInv, int n) {
	double *q, mod;
	q = new double[n * n];
    	memset(q, 0, n * n * sizeof(double));
	for (int i = 0; i < n; i++)
		q[i * n + i] = 1;
	for (int k = 0; k < n - 1; k++) { 
		for (int l = k + 1; l < n; l++)	{	
			double x = a[k * n + k], y = a[l * n + k];
			mod = sqrt(x * x + y * y);
			if (fabs(x) < eps && fabs(y) < eps) {
				continue;
			}
			double	Cos = x / mod;
			double	Sin = -y / mod;

			a[k * n + k] = mod;
			a[l * n + k] = 0;
			
			for (int m = 0; m < n; m++) {
				double temp_k, temp_l;
				if (m > k) {
					temp_k = Cos * a[k * n + m] - Sin * a[l * n + m];
					temp_l = Sin * a[k * n + m] + Cos * a[l * n + m];
					a[k * n + m] = temp_k;
					a[l * n + m] = temp_l;
				}			
				temp_k = Cos * q[k * n + m] - Sin * q[l * n + m];
				temp_l = Sin * q[k * n + m] + Cos * q[l * n + m];
				q[k * n + m] = temp_k;
				q[l * n + m] = temp_l;
			}
		}
		if (fabs(a[k * n + k]) < eps) {
			printf("Error, degenerate matrix\n");
			return -1;
		}
	}	
	print(a, n);

	for (int j = 0; j < n; j++)
		for (int i = n - 1; i >= 0; i--) {	
			aInv[i * n + j] = q[i * n + j];
			for (int k = i + 1; k < n; k++)
				aInv[i * n + j] -=  a[i * n + k] * aInv[k * n + j];
			aInv[i * n + j] /= a[i * n + i];
		}

	delete[] q;	
	return 0;
}

double residual(double *a, double *aInv, int n) {
	double *c;
	c = new double[n * n];
	memset(c, 0, n * n * sizeof(double));
	mult(a, aInv, c, n);
	for (int i = 0; i < n; i++)
		c[i * n + i] -= 1;
	double res = 0;
	for (int i = 0; i < n; i++) {
		double preRes = 0;
		for (int j = 0; j < n; j++)
			preRes += fabs(c[i * n + j]);
		res = max(res, preRes);
	}
	delete[] c;
	return res;
}

void print(double *a, int n)
{
	printf("\n");
	for (int i = 0; i < min(n, corner); i++) {
		for (int j = 0; j < min(n, corner); j++)
			printf("%-7.2f ", a[i * n + j]);
		if (n > corner)
			printf("%-7s %-7.2f ", "...", a[i * n + n - 1]);
		printf("\n");	
	}
	if (n > corner)
	{
		for (int j = 0; j <= corner + 1; j++)
			printf("%-7s ", "...");
		printf("\n");
		for (int j = 0; j < corner; j++)
			printf("%-7.2f ", a[(n - 1) * n + j]);
		printf("%-7s ", "...");
		printf("%-7.2f ", a[(n - 1) * n + n - 1]);
	}
	printf("\n");
}

void trans(double *a, int n)
{
        for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                {
                        double t = a[i * n + j];
                        a[i * n + j] = a[j * n + i];
                        a[j * n + i] = t;
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

void mult(double *a, double *b, double *c, int n) {
	trans(b, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			c[i * n + j] = 0;
			for (int k = 0; k < n; k++)
				c[i * n + j] += a[i * n + k] * b[j * n + k];
		}	
}

void generate(double *a, int n, double formula(int, int)) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i * n + j] = formula(i, j);
}


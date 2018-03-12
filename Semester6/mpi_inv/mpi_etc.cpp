#include <mpi.h>
#include <bits/stdc++.h>
#include "mpi_etc.h"
#define corner 5
#define  I(a, x, y) a[x * n + y]
#define  I_t(a, x, y) a[y * n + x]

using namespace std;

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

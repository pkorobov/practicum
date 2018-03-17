#include <mpi.h>
#include <bits/stdc++.h>
#include "mpi_etc.h"
#define corner 5
#define  I(a, x, y) a[x * n + y]
#define  I_t(a, x, y) a[y * n + x]

using namespace std;

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
	for (int i = 0; i < min(m, corner); i++) {
		for (int j = 0; j < min(n, corner); j++)
			printf("%-7.2f ", a[j * n + i]);
		if (n > corner)
			printf("%-7s %-7.2f ", "...", a[i + (n - 1) * n]);
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

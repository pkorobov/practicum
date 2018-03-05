#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include "others.h"
#define corner 5
#define maxn_threads 8

using namespace std;

struct thread_data_res
{
	double *a, *b, *c;
	int nthreads, thread_number, n;
};

double residual(double *a, double *aInv, int n, int nthreads) {
	
	thread_data_res arg[maxn_threads];
	pthread_t id[maxn_threads];

	double *c;
	c = new double[n * n];
	memset(c, 0, n * n * sizeof(double));

	for (int s = 0; s < nthreads; s++) {
		arg[s].a = a;
		arg[s].b = aInv;
		arg[s].c = c;
	
		arg[s].n = n;
		arg[s].thread_number = s;
		arg[s].nthreads = nthreads;
	}


	trans(aInv, n);

	for (int s = 0; s < nthreads; s++) {
		int rc = pthread_create(&id[s], NULL, mult, &arg[s]);		
		if(rc) {
			cout << "Break at creating thread! Error code: " << rc << endl;	
			return -1;
		}
	}

	for (int s = 0; s < nthreads; s++)
		pthread_join(id[s], NULL);

//	mult(a, aInv, c, n, 0, n);
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

int input(int argc, char **argv, double* &a, int &n, int &nthreads) {
	char *t;
	if (argc >= 3) {
		if (!strcmp(argv[1], "-formula")) {
			n = strtol(argv[3], &t, 10);
			nthreads = strtol(argv[4], &t, 10);
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
			nthreads = strtol(argv[3], &t, 10);
  			if(!freopen(name, "r", stdin)) {
				cin >> n;
				cout << n;
				a = new double[n * n];
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
						cin >> a[i * n + j];
			}
			else
				cout << "Error while reading file" << endl;
		}
		else {
			cout << "Incorrect input" << endl;
			return -1;
		}
	}
	else if (argc == 2) {
		cin >> n;

		nthreads = strtol(argv[1], &t, 10);

		a = new double[n * n];

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				cin >> a[i * n + j];
	}
	else {
		cout << "Incorrect input" << endl;
		return -1;
	}
	return 0;
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
		for (int j = i + 1; j < n; j++)	{
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

void mult(double *a, double *b, double *c, int n, int from, int to) {
	for (int i = from; i < n; i++)
		for (int j = 0; j < n; j++) {
			c[i * n + j] = 0;
			for (int k = 0; k < n; k++)
				c[i * n + j] += a[i * n + k] * b[j * n + k];
		}	
}

void *mult(void *data) {
	
	thread_data_res *unpacked = (thread_data_res *)data;

	double *a = unpacked->a;
	double *b = unpacked->b;
	double *c = unpacked->c;
	int from = unpacked->thread_number;
	int period = unpacked->nthreads;
	int n = unpacked->n;
	
	for (int i = from; i < n; i += period)
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


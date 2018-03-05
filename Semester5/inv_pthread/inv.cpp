#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include "others.h"

#define eps 1e-16
#define maxn_threads 8

using namespace std;


static pthread_barrier_t barrier;
//static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
//static pthread_cond_t condvar = PTHREAD_COND_INITIALIZER;

static int degenerate_flag = 0;
double *Cos, *Sin, mod;
double *rInv;


struct thread_data
{
	int n;
	int thread_number;
	int nthreads;
	double *a, *aInv, *q;
	int *from, *to;
};

void synchronize(int total_threads);
void *solve(void *data);

int main(int argc, char** argv)
{  
	int n, nthreads, rc;
	double *a = NULL, *aInv = NULL, *q = NULL;
	int *from = NULL, *to = NULL;
	thread_data arg[maxn_threads];
	pthread_t id[maxn_threads];

	input(argc, argv, a, n, nthreads);

	pthread_barrier_init(&barrier, NULL, nthreads);
		
	cout << nthreads << endl;
		
	printf("\nA\n");
	print(a, n);

	aInv = new double[n * n];
	q = new double[n * n];
	from = new int[nthreads];
	to = new int[nthreads];

	memset(q, 0, n * n * sizeof(double));
	for (int i = 0; i < n; i++)
		q[i * n + i] = 1;
	
	for (int i = 0; i < nthreads; i++)
	{
		if (i == 0)
			from[0] = 0;
		else
			from[i] = to[i-1];
		to[i] = from[i] + n / nthreads + (int)(i < n % nthreads);
	}
	
	for (int s = 0; s < nthreads; s++) {
		arg[s].a = a;
		arg[s].aInv = aInv;
		arg[s].q = q;
		arg[s].from = from;
		arg[s].to = to;

		arg[s].n = n;
		arg[s].thread_number = s;
		arg[s].nthreads = nthreads;
	}
	timespec t1, t2;	
	clock_gettime(CLOCK_REALTIME, &t1);
	time_t start = clock();	
	for (int s = 0; s < nthreads; s++) {
		rc = pthread_create(&id[s], NULL, solve, &arg[s]);		
		if(rc) {
			cout << "Break at creating thread! Error code: " << rc << endl;	
			return -1;
		}
	}

	for (int s = 0; s < nthreads; s++)
		pthread_join(id[s], NULL);

	clock_gettime(CLOCK_REALTIME, &t2);
	printf("\nA^-1\n");
 	print(aInv, n);

	input(argc, argv, a, n, nthreads);

	double time = (t2.tv_sec - t1.tv_sec)+(t2.tv_nsec - t1.tv_nsec)*1e-9;

	printf("%e\n", residual(a, aInv, n, nthreads));
	cout << "Time: " << time << " seconds" << endl;
	delete[] q;	
	delete[] a;
	delete[] aInv;
	return 0;
}

void *solve(void *data) {
	
	thread_data *unpacked = (thread_data *)data;
	
	double *a = unpacked->a;
	double *q = unpacked->q;
	double *aInv = unpacked->aInv;
	int *from = unpacked->from;
	int *to = unpacked->to;
	int n = unpacked->n;
	int thread_number = unpacked->thread_number;
	int nthreads = unpacked->nthreads;

	if (thread_number == 0){
		Cos = new double[n - 1];
		Sin = new double[n - 1];
		rInv = new double[n * n];
	}
	
	for (int k = 0; k < n - 1; k++) { 
	
		if (thread_number == 0) {				
		
			for (int l = k + 1; l < n; l++)	{	
				
				double x = a[k * n + k], y = a[l * n + k];
				mod = sqrt(x * x + y * y);
				if (fabs(x) < eps && fabs(y) < eps) {
					continue;
				}

				Cos[l - k - 1] = x / mod;
				Sin[l - k - 1] = -y / mod;

				a[k * n + k] = mod;
				a[l * n + k] = 0;					
			}
		}
/*		
		pthread_mutex_lock(&mutex);
	    cout << "Thread " << thread_number << " is waiting at point 1..." << endl;
		pthread_mutex_unlock(&mutex);
*/
		pthread_barrier_wait(&barrier);
/*
		pthread_mutex_lock(&mutex);
	    cout << "Thread " << thread_number << " is going at point 1!" << endl;
		pthread_mutex_unlock(&mutex);
*/	
		for (int l = k + 1; l < n; l++)	{		
			for (int m = from[thread_number]; m < to[thread_number]; m++) {	
//			for (int m = thread_number; m < n; m += nthreads) {	
				double temp_k, temp_l;
				if (m > k) {
					temp_k = Cos[l - k - 1] * a[k * n + m] - Sin[l - k - 1] * a[l * n + m];
					temp_l = Sin[l - k - 1] * a[k * n + m] + Cos[l - k - 1] * a[l * n + m];
					a[k * n + m] = temp_k;
					a[l * n + m] = temp_l;
				}			
				temp_k = Cos[l - k - 1] * q[k * n + m] - Sin[l - k - 1] * q[l * n + m];
				temp_l = Sin[l - k - 1] * q[k * n + m] + Cos[l - k - 1] * q[l * n + m];
				q[k * n + m] = temp_k;
				q[l * n + m] = temp_l;
			}		
		}
/*
		pthread_mutex_lock(&mutex);
	    cout << "Thread " << thread_number << " is waiting at point 2..." << endl;
		pthread_mutex_unlock(&mutex);
*/
		pthread_barrier_wait(&barrier);
/*
		pthread_mutex_lock(&mutex);
	    cout << "Thread " << thread_number << " is going at point 2!" << endl;
		pthread_mutex_unlock(&mutex);
*/
		if (fabs(a[k * n + k]) < eps) {
			printf("Error, degenerate matrix\n");
			degenerate_flag = 1;
			pthread_exit(NULL);
		}
	
	}


	for (int j = thread_number; j < n; j +=  nthreads)
		for (int i = n - 1; i >= 0; i--) {	
			rInv[i * n + j] = (i == j ? 1 : 0);
			for (int k = i + 1; k < n; k++)
				rInv[i * n + j] -=  a[i * n + k] * rInv[k * n + j];
			rInv[i * n + j] /= a[i * n + i];
		}
	if (thread_number == 0)
		trans(q, n);

	pthread_barrier_wait(&barrier);
	mult(rInv, q, aInv, n, from[thread_number], to[thread_number]);
	pthread_barrier_wait(&barrier);
	
	if (thread_number == 0) {
		print(a, n);
		print(q, n);
		delete[] rInv;
	}
	return NULL;
}

void synchronize(int total_threads){
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads) {
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} 
	else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads) {
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} 
	else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);
}

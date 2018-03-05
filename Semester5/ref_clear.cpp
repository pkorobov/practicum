#include <bits/stdc++.h>

#define corner 5
#define eps 1e-6
using namespace std;

pthread_barrier_t barrier;

int solve(double *a, int n) {
	for (int k = 0; k < n - 1; k++) {
		double s = 0;
		for (int l = k + 1; l < n; l++) 
			s += a[l * n + k] * a[l * n + k];
		double columnNorm = sqrt(a[k * n + k] * a[k * n + k] + s);
		double *x;
		x = new double[n - k];
		for (int l = 0; l < n - k; l++) 
			x[l] = a[(k + l) * n + k] - (l == 0 ? columnNorm : 0);
		double xNorm = sqrt(x[0] * x[0] + s);
		for (int l = 0; l < n - k; l++) 
			x[l] /= xNorm;
		a[k * n + k] = columnNorm;
		cout << columnNorm << endl;
		for (int l = k + 1; l < n; l++)
			a[l * n + k] = 0;
		for (int j = k + 1; j < n; j++) {
			double dot = 0;
			for (int i = k; i < n; i++)
				dot += x[i - k] * a[i * n + j];
			for (int i = k; i < n; i++) {
				a[i * n + j] -= 2 * x[i - k] * dot; 
			}
		}
//		double u_ij = (i == j ? 1 : 0) - x[i] * x[j]		
	}
}
void print(double *a, int n) {
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

int main(int argc, char **argv) {
	int n;
	double *a;
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
	if (argc == 2) {
		char *t;
		n = strtol(argv[1], &t, 10);
		a = new double[n * n];
		generate(a, n, formula2);
	}
	print(a, n);
	cout << endl;
	solve(a, n);
	print(a, n);
	return 0;
}


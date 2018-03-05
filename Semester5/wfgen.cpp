#include <bits/stdc++.h>

using namespace std;
#define delta 1e-16

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

void print(double *a, int n) {
	cout << "{";		
	for (int i = 0; i < n - 1; i++) {
		
		cout << "{" << a[i * n] << ", ";		
		for (int j = 1; j < n - 1; j++)
			cout << a[i * n + j] << ", ";
		cout << a[i * n + (n - 1)] << "}, " << endl;
	}
	
	cout << "{" << a[(n - 1) * n] << ", ";		
	for (int j = 1; j < n - 1; j++)
		cout << a[(n - 1) * n + j] << ", ";
	cout << a[(n - 1) * n + (n - 1)] << "}}" << endl;
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
	cout << "Eigenvalues[";
	print(a, n);
	cout << "]" << endl;;
	return 0;
}

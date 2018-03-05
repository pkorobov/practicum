#include <bits/stdc++.h>

#define corner 5
#define eps 1e-7
#define delta 1e-12

using namespace std;
void print(double *a, int n);
void printVector(double *a, int n);

int eigenValues(double *ans, double *diag, double *lowerDiag, double *upperDiag, double aNorm, int n);
int QRDecomposition(double *diag, double *lowerDiag, double *upperDiag, double *qCos, double *qSin, int n);
int reduction(double *a, int n);
void generate(double *a, int n, double formula(int, int));
void generate(double *a, int n, double formula(int, int, int));
double formula1(int i, int j);
double formula2(int i, int j);
double formula3(int i, int j);
double formula4(int i, int j, int n);
double norm(double *a, int n);

int main(int argc, char **argv) {
	int n;
	double *a, *lowerDiag, *diag, *upperDiag, *ans;

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
				case '4' : generate(a, n, formula4); break;
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
	double euclideanNorm = 0;
	for(int i = 0; i < n; i++)
		euclideanNorm += a[i] * a[i];
	euclideanNorm = sqrt(euclideanNorm);
	double trace = 0, trace2 = 0;
	for (int i = 0; i < n; i++)
		trace += a[i * n + i];
	double aNorm = norm(a, n);
	print(a, n);
	cout << endl;
//Приведение к трехдиагональной матрице
	reduction(a, n);
	print(a, n);
	ans = new double[n];
	diag = new double[n];
	lowerDiag = new double[n - 1];
	upperDiag = new double[n - 1];
	
	for (int i = 0; i < n; i++)
		diag[i] = (fabs(a[i * n + i]) > delta ? a[i * n + i] : 0);
	for (int i = 0; i < n - 1; i++)
		lowerDiag[i] = (fabs(a[i * n + i + 1]) > delta ? a[i * n + i + 1] : 0);	
	for (int i = 0; i < n - 1; i++)
		upperDiag[i] = (fabs(a[i * n + i + 1]) > delta ? a[i * n + i + 1] : 0);
	
	delete[] a;
	
//Собственные значения
	eigenValues(ans, diag, lowerDiag, upperDiag, aNorm, n);
	for (int i = 0; i < n; i++)
		trace2 += ans[i];
	printVector(ans, n);
	double det = 1;
	for (int i = 0; i < n; i++)
		det *= ans[i];
	cout << "Trace of input: " << trace << endl; 
	cout << "Trace of result: " << trace2 << endl; 
	cout << "Euclidean norm of input: " << euclideanNorm << endl;
	cout << "Det of result: " << det << endl;

	delete[] diag;
	delete[] upperDiag;
	delete[] lowerDiag;
	delete[] ans;
	
	return 0;
}


int QRDecomposition(double *diag, double *lowerDiag, double *upperDiag, double *qCos, double *qSin, int n) {
	for(int k = 0; k < n - 1; k++) {
		
		if (fabs(diag[k]) < delta) diag[k] = 0; 
		if (fabs(lowerDiag[k]) < delta) lowerDiag[k] = 0;
		if (fabs(diag[k + 1]) < delta) diag[k + 1] = 0; 
		if (fabs(upperDiag[k]) < delta) upperDiag[k] = 0;
		if (k < n - 2 && fabs(upperDiag[k + 1]) < delta) upperDiag[k + 1] = 0;

		double x = diag[k], y = lowerDiag[k];
		
        double mod = sqrt(x * x + y * y);
       	if (mod < delta) {
			qCos[k] = 1;
			qSin[k] = 0;
			continue;
		}
		
		qCos[k] = x / mod;
		qSin[k] = -y / mod;
		diag[k] = mod;
		lowerDiag[k] = 0;	
		
//		if (fabs(qCos[k]) < delta) qCos[k] = 0; 
//		if (fabs(qSin[k]) < delta) qSin[k] = 0;
		
		double temp1 = qCos[k] * upperDiag[k] - qSin[k] * diag[k + 1];
		double temp2 = qSin[k] * upperDiag[k] + qCos[k] * diag[k + 1];
        upperDiag[k] = temp1;
        diag[k + 1] = temp2;

		if (k < n - 2) {
			temp2 = qCos[k] * upperDiag[k + 1];
			upperDiag[k + 1] = temp2;
        }
        
	}
	return 0;
}

int eigenValues(double *ans, double *diag, double *lowerDiag, double *upperDiag, double aNorm, int n) {
	double *qCos, *qSin;
	qCos = new double[n - 1];
	qSin = new double[n - 1];
	int k = 0;
	memset(ans, 0, n - 1);
	int step = 0;
	while (n - k > 2) {
		step++;
/*		printVector(diag, n);
		printVector(lowerDiag, n - 1);
*/
		double s_k = diag[n - k - 1];
		
		for(int l = 0; l < n - k; l++)
			diag[l] -= s_k;

		QRDecomposition(diag, lowerDiag, upperDiag, qCos, qSin, n - k);
// A' = RQ
		
		for(int l = 0; l < n - k - 1; l++) {
			
			double temp1 = qCos[l] * diag[l] - qSin[l] * upperDiag[l];
			double temp2 = qSin[l] * diag[l] + qCos[l] * upperDiag[l];
			diag[l] = temp1;
       		upperDiag[l] = temp2;
			
			temp1 = qCos[l] * lowerDiag[l] - qSin[l] * diag[l + 1];
			temp2 = qSin[l] * lowerDiag[l] + qCos[l] * diag[l + 1];
			lowerDiag[l] = temp1;
          	diag[l + 1] = temp2;

		}	
		for(int l = 0; l < n - k - 1; l++)
			upperDiag[l] = lowerDiag[l];
		
		for(int l = 0; l < n - k; l++)
			diag[l] += s_k;
		
		if (fabs(lowerDiag[n - k - 2]) < eps * aNorm) {
			ans[n - k - 1] = diag[n - k - 1];
			k++;
		}
	}
	cout << "Iterations in general: " << step << endl;
	double b = -(diag[0] + diag[1]);
	double c = diag[0] * diag[1] - lowerDiag[0] * lowerDiag[0];
	double discriminant = b * b - 4 * c;
	if (discriminant < eps)	
		discriminant = 0;
	ans[0] = (-b - sqrt(discriminant))/2;
	ans[1] = (-b + sqrt(discriminant))/2;

	for (int i = 0; i < n; i++)
		if (fabs(ans[i]) < eps)	
			ans[i] = 0;
	delete[] qSin;
	delete[] qCos;
	return 0;	
}

int reduction(double *a, int n) {
	double *x, *y, *z;
	x = new double[n];
	y = new double[n];
	z = new double[n];

	for (int k = 0; k < n - 2; k++) {
		double s = 0;
		for (int l = k + 2; l < n; l++) 
			s += a[l * n + k] * a[l * n + k];
		double columnNorm = (sqrt(a[(k + 1) * n + k] * a[(k + 1) * n + k] + s) > delta ? sqrt(a[(k + 1) * n + k] * a[(k + 1) * n + k] + s) : 0);
		
		for (int l = 0; l < n - k - 1; l++) 
			x[l] = a[(k + l + 1) * n + k] - (l == 0 ? columnNorm : 0);
		double xNorm = (sqrt(x[0] * x[0] + s) > delta ? sqrt(x[0] * x[0] + s) : 0);
		
		if (xNorm < delta)
			continue;
		
		for (int l = 0; l < n - k - 1; l++) 
			x[l] /= xNorm;
		
		double dot = 0;

		a[(k + 1) * n + k] = columnNorm;
		a[k * n + (k + 1)] = columnNorm;
		for (int l = k + 2; l < n; l++) {
			a[l * n + k] = 0;
			a[k * n + l] = 0;
		}

/*
Счет к-го столбца руками
  		for (int i = k; i < n - 1; i++)
			dot += x[i - k] * a[(i + 1) * n + k];

		for (int i = k; i < n - 1; i++)
			a[(i + 1) * n + k] -= 2 * x[i - k] * dot; 
	
		for (int j = k; j < n - 1; j++)
			a[k * n + (j + 1)] -= 2 * x[j - k] * dot;
*/

		memset(y, 0, (n - k - 1) * sizeof(double));
		for (int m = 0; m < n - k - 1; m++)
			for (int l = 0; l < n - k - 1; l++) 
				y[m] += a[(k + 1 + m) * n + (k + 1 + l)] * x[l];

		dot = 0;
		for (int l = 0; l < n - k - 1; l++)
			dot += x[l] * y[l];
		dot *= 2;

		for (int l = 0; l < n - k - 1; l++)
			z[l] = 2 * y[l] - dot * x[l];
		
		for (int m = 0; m < n - k - 1; m++)
			for (int l = m; l < n - k - 1; l++) 
				a[(k + 1 + l) * n + (k + 1 + m)] -= z[l] * x[m] + x[l] * z[m];

		for (int m = 0; m < n - k - 1; m++)
			for (int l = 0; l < m; l++) 
				a[(k + 1 + l) * n + (k + 1 + m)] = a[(k + 1 + m) * n + (k + 1 + l)];
		
	}
	delete[] x;
	delete[] y;
	delete[] z;
	return 0;
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

double formula4(int i, int j, int n) {
	if (i == j && i != n - 1)
		return 1;
	else if (i == n - 1)
		return j + 1;
	else if (j == n - 1)
		return i + 1;
	return 0;
}

void generate(double *a, int n, double formula(int, int)) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i * n + j] = formula(i, j);
}

void generate(double *a, int n, double formula(int, int, int)) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i * n + j] = formula(i, j, n);
}

double norm(double *a, int n) {
	double res = 0;
	for (int i = 0; i < n; i++) {
		double preRes = 0;
		for (int j = 0; j < n; j++)
			preRes += fabs(a[i * n + j]);
		res = max(res, preRes);
	}
	return res;
}

void print(double *a, int n) {
	printf("\n");
	for (int i = 0; i < min(n, corner); i++) {
		for (int j = 0; j < min(n, corner); j++)
			printf("%-9.4f ", a[i * n + j]);
		if (n > corner)
			printf("%-9.4s %-9.4f ", "...", a[i * n + n - 1]);
		printf("\n");	
	}
	if (n > corner)
	{
		for (int j = 0; j <= corner + 1; j++)
			printf("%-9.4s ", "...");
		printf("\n");
		for (int j = 0; j < corner; j++)
			printf("%-9.4f ", a[(n - 1) * n+ j]);
		printf("%-9s ", "...");
		printf("%-9.4f ", a[(n - 1) * n + n - 1]);
	}
	printf("\n");
}

void printVector(double *a, int n){
	printf("\n");
	for (int i = 0; i < min(n, corner); i++)
			printf("%-9.4f ", a[i]);
	if (n > corner * 2)
		printf("%-9.4s ", "...");
	for (int i = max(n - corner, corner); i < n; i++)
		printf("%-9.4f ", a[i]);		
	printf("\n");
}

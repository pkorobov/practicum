#ifndef OTHERS
#define OTHERS
void generate(double *a, int n, double formula(int, int));
double formula1(int i, int j);
double formula2(int i, int j);
double formula3(int i, int j);
double residual(double *a, double *aInv, int n, int nthreads);
void mult(double *a, double *b, double *c, int n, int from, int to);
void trans(double *a, int n);
void print(double *a, int n);
int input(int argc, char** argv, double* &a, int &n, int &nthreads);
void *mult(void *data);

#endif


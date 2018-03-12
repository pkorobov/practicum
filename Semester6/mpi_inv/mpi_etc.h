#ifndef ETC
#define ETC

double residual(double *a, double *aInv, int n);
void mult(double *a, double *b, double *c, int n, int m);
void trans(double *a, int n, int m);
void print(double *a, int n, int m);
void generate(double *a, int n, double formula(int, int));
double formula1(int i, int j);
double formula2(int i, int j);
double formula3(int i, int j); 

#endif
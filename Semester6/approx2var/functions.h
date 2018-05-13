#ifndef FUNCTIONS_H
#define FUNCTIONS_H

int initValues(bool fromFile, double *x, double *values, double *der);

double f_0 (double x, double y);
double f_1 (double x, double y);
double f_2 (double x, double y);
double f_3 (double x, double y);
double f_4 (double x, double y);

double dfx_0 (double x, double y);
double dfx_1 (double x, double y);
double dfx_2 (double x, double y);
double dfy_2 (double x, double y);
double dfx_3 (double x, double y);
double dfy_3 (double x, double y);
double dfx_4 (double x, double y);
double dfy_4 (double x, double y);

double dfy_0 (double x, double y);
double dfy_1 (double x, double y);

void mult(double *a, double *b, double *c, int n);
void trans(double *a, int n);

#endif // FUNCTIONS_H

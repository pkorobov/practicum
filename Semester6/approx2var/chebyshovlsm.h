#ifndef CHEBYSHOVLSM_H
#define CHEBYSHOVLSM_H
#include <QtMath>
#include "abstractmethod.h"

class ChebyshovLSM : public abstractMethod
{
public:
    ChebyshovLSM();
    ChebyshovLSM(double a_, double b_, int n_, double *x_, double *values_, int N_,
                 double c_ = -1, double d_ = 1, int m_ = 5, double *y_ = nullptr);
    ~ChebyshovLSM();
    double (*f) (double);
    double method_compute(double x);
    double method_compute2v(double x, double y);

    void method_init();
    void method_init2v(double *x, double *y, double *values, int n, int m, int N);
    void calc_coefficients(double a, double b, double *x_, double *values, int n, int N, double *alpha);


private:
    int N;
    double *alpha, *betta, *gamma;
    double f_(double x);
    double U(int n, double x);
    double T(int n, double x, double a, double b);
    double IntA(int i, double x);
    double IntB(int i, double x);
};

#endif  // CHEBYSHOVLSM_H

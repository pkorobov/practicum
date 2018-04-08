#ifndef CHEBYSHOVLSM_H
#define CHEBYSHOVLSM_H
#include <QtMath>

class ChebyshovLSM
{
public:
    double a, b;
    double *x;
    int n, N;
    ChebyshovLSM();
    ChebyshovLSM(double (*f_)(double), double a_ = -1, double b_ = 1, int n_ = 5, int N_ = 10);
    ~ChebyshovLSM();
    double (*f) (double);
    double method_compute(double x);
    void method_init();
//    void change_func(double (*f_)(double));

private:
    double *alpha;
    double f_(double x);
    double U(int n, double x);
    double T(int n, double x);
    double IntA(int i, double x);
    double IntB(int i, double x);
};

#endif  // CHEBYSHOVLSM_H

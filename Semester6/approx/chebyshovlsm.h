#ifndef CHEBYSHOVLSM_H
#define CHEBYSHOVLSM_H
#include <QtMath>
#include "abstractmethod.h"

class ChebyshovLSM : public abstractMethod
{
public:
    int N;
    ChebyshovLSM();
    ChebyshovLSM(double a_, double b_, int n_, double *x_, double *values_, int N_);
    ~ChebyshovLSM();
    double (*f) (double);
    double method_compute(double x);
    void method_init();

private:
    double *alpha;
    double f_(double x);
    double U(int n, double x);
    double T(int n, double x);
    double IntA(int i, double x);
    double IntB(int i, double x);
};

#endif  // CHEBYSHOVLSM_H

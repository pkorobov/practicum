#ifndef CUBICSPLINES_H
#define CUBICSPLINES_H
#include "abstractmethod.h"

class cubicSplines : public abstractMethod
{
public:
    cubicSplines();
    cubicSplines(double a_, double b_, int n_, double *x_, double *values_,
                 double c_ = -1, double d_ = 1, int m_ = 5, double *y_ = nullptr);
    ~cubicSplines();

    void method_init(double *derivatives);
    double method_compute(double x);
    double method_init2v(double *x, double *y, double *values, int n, int m, int N);
    void calc_coefficients(double a, double b, double *x_, double *values, int n, int N, double *alpha);

private:
    double *state;
};

#endif // CUBICSPLINES_H

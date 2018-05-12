#ifndef CUBICSPLINES_H
#define CUBICSPLINES_H
#include "abstractmethod.h"
#include "QGLWidget"

class cubicSplines : public abstractMethod
{
public:
    cubicSplines();
    cubicSplines(double a_, double b_, int n_, double *x_, double *derx_, double *values_,
                 double c_ = -1, double d_ = 1, int m_ = 5, double *y_ = nullptr, double *dery_ = nullptr);
    ~cubicSplines();

    void method_init(double *derivatives);
    double method_compute(double x);
    double method_init2v(double *x, double *y, double *values, int n, int m, double *derx, double *dery);
    double method_compute2v(double x1, double x2);
    void calc_coefficients(double a, double b, double *x_, double *values, int n, double *derivatives, double *d);
private:
    double *derx, *dery;
    double *state;
    double gamma[200][200][16];
};

#endif // CUBICSPLINES_H

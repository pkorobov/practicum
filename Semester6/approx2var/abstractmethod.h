#ifndef ABSTRACTMETHOD_H
#define ABSTRACTMETHOD_H


class abstractMethod
{
protected:
    double a, b, c, d;
    double *x, *y, *values;
    int n, m;
public:
    abstractMethod();
    abstractMethod(double a_, double b_, int n_, double *x_, double *values_,
                   double c_ = -1, double d_ = 1, int m_ = 5, double *y = nullptr);
//    virtual double method_init(double *values) = 0;
//    virtual double method_compute(double x) = 0;
    virtual ~abstractMethod();
};

#endif // ABSTRACTMETHOD_H

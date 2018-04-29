#ifndef ABSTRACTMETHOD_H
#define ABSTRACTMETHOD_H


class abstractMethod
{
protected:
    double a, b;
    double *x, *values;
    int n;
public:
    abstractMethod();
    abstractMethod(double a_, double b_, int n_, double *x_, double *values_);
//    virtual double method_init(double *values) = 0;
    virtual double method_compute(double x) = 0;
    virtual ~abstractMethod();
};

#endif // ABSTRACTMETHOD_H

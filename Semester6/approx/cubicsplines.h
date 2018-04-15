#ifndef CUBICSPLINES_H
#define CUBICSPLINES_H
#include "abstractmethod.h"

class cubicSplines : public abstractMethod
{
public:
    cubicSplines();
    cubicSplines(double a_, double b_, int n_, double *x_, double *values_);
    ~cubicSplines();

    double method_init(double *derivatives);
    double method_compute(double x);
private:
    double *state;
};

#endif // CUBICSPLINES_H

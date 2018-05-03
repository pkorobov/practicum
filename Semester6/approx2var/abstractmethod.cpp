#include "abstractmethod.h"

abstractMethod::abstractMethod()
{
    ;
}

abstractMethod::abstractMethod(double a_, double b_, int n_, double *x_, double *values_,
double c_, double d_, int m_, double *y_)
{
    n = n_;
    a = a_;
    b = b_;
    m = m_;
    c = c_;
    d = d_;
    y = y_;
    values = values_;
}

abstractMethod::~abstractMethod()
{

}

#include "abstractmethod.h"

abstractMethod::abstractMethod()
{
    ;
}

abstractMethod::abstractMethod(double a_, double b_, int n_, double *x_, double *values_)
{
    n = n_;
    a = a_;
    b = b_;
    x = x_;
    values = values_;
}

abstractMethod::~abstractMethod()
{

}

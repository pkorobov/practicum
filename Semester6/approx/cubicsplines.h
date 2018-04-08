#ifndef CUBICSPLINES_H
#define CUBICSPLINES_H


class cubicSplines
{
public:
    cubicSplines();
    cubicSplines(double (*f_)(double), double a_ = -1.0, double b_ = 1.0, int n_ = 10);
    ~cubicSplines();
    double a, b;
    double *x;
    int n;
    double (*f) (double);

    double method_init(double *values, double *derivatives);
    double method_compute(double x);
private:
    double *state;
};

#endif // CUBICSPLINES_H

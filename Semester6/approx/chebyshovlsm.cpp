#include "chebyshovlsm.h"
#include <QtMath>
#include <QDebug>
#include <QTime>
#include "window.h"


//N - число узлов
//n - степень многочлена Чебышева
double f_default(double x)
{
//    return qSin(x);
//    return qSin(x) / qSqrt(1 + x * x);
      return qCos(x);
}

ChebyshovLSM::ChebyshovLSM()
{
    a = -10;
    b = 10;
    n = 8;
    N = 5;
    f = f_default;
    alpha = new double[N];
    x = new double[n];
    memset(alpha, 0, N * sizeof(double));
    memset(x, 0, n * sizeof(double));
}

ChebyshovLSM::ChebyshovLSM(double a_, double b_, int n_, double *x_, double *values_, int N_)
    : abstractMethod(a_, b_, n_, x_, values_)
{
    N = N_;
    alpha = new double[N];
    memset(alpha, 0, N * sizeof(double));
}

ChebyshovLSM::~ChebyshovLSM()
{
    delete[] alpha;
}
/*
double ChebyshovLSM::f_(double x)
{
    double y = ((a + b) + (b - a) * x) * 0.5;
    return f(y);
}
*/
//Многочлен чебышева второго рода
double ChebyshovLSM::U(int n, double x)
{
    if (n == 0)
        return 1;
    else if (n == 1)
        return 2 * x;
    else
        return 2 * x * U(n - 1, x) - U(n - 2, x);
}

double ChebyshovLSM::T(int i, double x)
{
    double z = 2 * (2 * x - (a + b)) / (b - a);
    if (i == 0)
        return 1;
    else if (i == 1)
        return z / 2;
    else
        return z * T(i - 1, x) - T(i - 2, x);
}

//Интеграл T_i(x) / (sqrt(1 - x^2)
double ChebyshovLSM::IntA(int i, double x)
{
    if (i == 0)
        return -qAcos(x);
    else
        return -1.0 / i * qSqrt(1 - x * x) * U(i - 1, x);
}

//Интеграл T_i(x) * x / (sqrt(1 - x^2)
double ChebyshovLSM::IntB(int i, double x)
{
    if (i == 0)
        return -qSqrt(1 - x * x) * U(0, x);
    else if (i == 1)
        return -0.5 * (qAcos(x) + 0.5 * qSqrt(1 - x * x) * U(1, x));
    else
        return -0.5 * qSqrt(1 - x * x) * (1.0 / (i - 1) * U(i - 2, x) + 1.0 / (i + 1) * U(i, x));
}

void ChebyshovLSM::method_init()
{
    double *c, *d, *u;
    double *x;

    x = new double[n];
    u = new double[n * N];
    c = new double[n * N];
    d = new double[n * N];

    for (int i = 0; i < n; i++)
        x[i] = (2 * ChebyshovLSM::x[i] - (b + a)) / (b - a);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j < n - 1) {
                double a_ij = IntA(i, x[j + 1]) - IntA(i, x[j]);
                double b_ij = IntB(i, x[j + 1]) - IntB(i, x[j]);
                c[i * n + j] = (a_ij * x[j + 1] - b_ij) / (x[j + 1] - x[j]);
                d[i * n + j] = (b_ij - a_ij * x[j]) / (x[j + 1] - x[j]);
            }

            if (j == 0)
                u[i * n + j] = c[i * n + j];
            else if (j == n - 1)
                u[i * n + j] = d[i * n + j - 1];
            else
                u[i * n + j] = c[i * n + j] + d[i * n + j - 1];
        }
    }

    for (int i = 0; i < N; i++)
    {
        alpha[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            alpha[i] += u[i * n + j] * values[j];
        }
        alpha[i] *=  2 / M_PI;
        if (i == 0) alpha[i] /= 2;
    }
    delete[] x;
    delete[] u;
    delete[] c;
    delete[] d;
}

double ChebyshovLSM::method_compute(double x)
{
    double value = 0;
    for (int i = 0; i < N; i++)
        value += alpha[i] * T(i, x);
    return value;
}

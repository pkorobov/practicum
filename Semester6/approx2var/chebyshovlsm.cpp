#include "chebyshovlsm.h"
#include <QtMath>
#include <QDebug>
#include <QTime>
#include "window.h"
#include <algorithm>
#include <cmath>

using namespace std;

//N - число узлов
//n - степень многочлена Чебышева

ChebyshovLSM::ChebyshovLSM()
{

}

ChebyshovLSM::ChebyshovLSM(double a_, double b_, int n_, double *x_, double *values_, int N_,
                           double c_, double d_, int m_, double *y_)
    : abstractMethod(a_, b_, n_, x_, values_, c_, d_, m_, y_)
{
    N = N_;
    alpha = new double[N];
    memset(alpha, 0, N * sizeof(double));
}

ChebyshovLSM::~ChebyshovLSM()
{
    delete[] alpha;
    delete[] gamma;
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

void ChebyshovLSM::calc_coefficients(double a, double b, double *x_, double *values, int n, int N, double *alpha)
{
    double *c, *d, *u;
    double *x;

    x = new double[n];
    u = new double[n * N];
    c = new double[n * N];
    d = new double[n * N];

    for (int i = 0; i < n; i++)
        x[i] = (2 * x_[i] - (b + a)) / (b - a);
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

void ChebyshovLSM::method_init2v(double *x, double *y, double *values, int n, int m, int N)
{
    int L = max(max(m, n), N);
    gamma = new double[L * L];
    double a = x[0], b = x[n - 1];
    double c = y[0], d = y[m - 1];
    for (int i = 0; i < n; i++)
    {
        double column[m], res[N];
        for (int k = 0; k < m; k++)
            column[k] = values[i + k * L];
        calc_coefficients(a, b, x, column, m, N, res);
        for (int k = 0; k < N; k++)
            gamma[i + k * L] = res[k];
    }
    for (int j = 0; j < N; j++)
    {
        double string[n];
        for (int k = 0; k < n; k++)
            string[k] = gamma[k + j * L];
        calc_coefficients(c, d, y, string, n, N, gamma + j * L);
    }
}

double ChebyshovLSM::method_compute(double x)
{
    double value = 0;
    for (int i = 0; i < N; i++)
        value += alpha[i] * T(i, x);
    return value;
}

double ChebyshovLSM::method_compute2v(double x, double y)
{
    int L = max(max(m, n), N);
    double value = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            value += gamma[i + j * L] * T(i, x) * T(j, y);
    return value;
}

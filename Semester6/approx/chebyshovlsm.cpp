#include "chebyshovlsm.h"
#include <QtMath>
#include <QDebug>
#include <QTime>
#include "window.h"


//N - число узлов
//n - степень многочлена Чебышева
double f_default(double x)
{
      return qCos(x);
}

ChebyshovLSM::ChebyshovLSM()
{
    a = -10;
    b = 10;
    n = 5;
    N = 8;
    f = f_default;
    state = new double[n];
    x = new double[N + 1];
    memset(state, 0, n * sizeof(double));
    memset(x, 0, (N + 1) * sizeof(double));
}

ChebyshovLSM::ChebyshovLSM(double (*f_)(double), double a_, double b_, int n_, int N_)
{
    f = f_;
    a = a_;
    b = b_;
    n = n_;
    N = N_;
    state = new double[n];
    x = new double[N + 1];
    memset(state, 0, n * sizeof(double));
    memset(x, 0, (N + 1) * sizeof(double));
}

/*void ChebyshovLSM::change_func(double (*f_)(double))
{
    f = f_;
}
*/
ChebyshovLSM::~ChebyshovLSM()
{
    delete[] state;
}

double ChebyshovLSM::f_(double x)
{
    double y = ((a + b) + (b - a) * x) * 0.5;
    return f(y);
}

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

    u = new double[n * (N + 1)];
    c = new double[n * (N + 1)];
    d = new double[n * (N + 1)];


    qsrand(QTime::currentTime().msec());

    double delta = (b - a) / N;
    for (int i = 0; i <= N; i++)
        x[i] = a + delta * i;

    for (int i = 0; i <= N; i++)
    {
        x[i] = (2 * x[i] - (b + a)) / (b - a);
//        x[i] = -1 + (double) qrand() / RAND_MAX * 2;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= N; j++)
        {
            if (j < N) {
                double a_ij = IntA(i, x[j + 1]) - IntA(i, x[j]);
                double b_ij = IntB(i, x[j + 1]) - IntB(i, x[j]);
                c[i * (N + 1) + j] = (a_ij * x[j + 1] - b_ij) / (x[j + 1] - x[j]);
                d[i * (N + 1) + j] = (b_ij - a_ij * x[j]) / (x[j + 1] - x[j]);
            }

            if (j == 0)
                u[i * (N + 1) + j] = c[i * (N + 1) + j];
            else if (j == N)
                u[i * (N + 1) + j] = d[i * (N + 1) + j - 1];
            else
                u[i * (N + 1) + j] = c[i * (N + 1) + j] + d[i * (N + 1) + j - 1];
        }
    }

    for (int i = 0; i < n; i++)
    {
        state[i] = 0.0;
        for (int j = 0; j <= N; j++)
        {
            state[i] += u[i * (N + 1) + j] * f_(x[j]);
        }
        state[i] *=  2 / M_PI;
        if (i == 0) state[i] /= 2;
//        qDebug() << i << " " << state[i];
    }
}

double ChebyshovLSM::method_compute(double x)
{
    double value = 0;
    for (int i = 0; i < n; i++)
        value += state[i] * T(i, x);
    return value;
}

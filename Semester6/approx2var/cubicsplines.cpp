#include "cubicsplines.h"
#include <QTime>
#include <QDebug>
#include <iostream>


cubicSplines::cubicSplines()
{
    state = new double[4 * n];
    memset(state, 0, 4 * n * sizeof(double));
}

cubicSplines::~cubicSplines()
{
    delete[] state;
}

cubicSplines::cubicSplines(double a_, double b_, int n_, double *x_, double *values_,
                           double c_, double d_, int m_, double *y_)
    : abstractMethod(a_, b_, n_, x_, values_, c_, d_, m_)
{
    state = new double[4 * n];
    memset(state, 0, 4 * n * sizeof(double));
}

void cubicSplines::method_init(double *derivatives)
{
    double *X;
    double *d;

    d = new double[n];
    X = new double[n * (n + 1)];
    memset(X, 0, n * (n + 1) * sizeof(double));

    X[0 * n + 0] = 1;
    X[n * n + 0] = derivatives[0];

    X[(n - 1) * n + (n - 1)] = 1;
    X[n * n + (n - 1)] = derivatives[1];
    for (int i = 1; i < n - 1; i++)
    {
        X[i * n + (i - 1)] = x[i + 1] - x[i];
        X[i * n + i] = 2 * (x[i + 1] - x[i - 1]);
        X[i * n + (i + 1)] = x[i] - x[i - 1];
        X[n * n + i] = 3 * (values[i] - values[i - 1]) / (x[i] - x[i - 1]) * (x[i + 1] - x[i]) +
                       3 * (values[i + 1] - values[i]) / (x[i + 1] - x[i]) * (x[i] - x[i - 1]);
    }
/*    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j < n; j++)
            qDebug() << X[i * n + j] << " " << (X + n * n + j);
        qDebug() << endl;
    }
*/    for (int k = 0; k < n; k++)
        for (int i = k + 1; i < n; i++)
        {
            double t = X[i * n + k] / X[k * n + k];
//            qDebug() << "[[[[]]]" << t;
            for (int j = k + 1; j <= n; j++)
                X[i * n + j] -= t * X[k * n + j];
            X[i * n + k] = 0;

            X[n * n + i] -= t * X[n * n + k];
        }

    for (int i = n - 1; i >= 0; i--)
    {
        d[i] = X[n * n + i];
        for (int j = i + 1; j < n; j++)
            d[i] -= d[j] * X[i * n + j];
        d[i] /= X[i * n + i];
    }

/*    for (int i = 0; i < n; i++)
        qDebug() << d[i];
*/
    for (int i = 0; i < n - 1; i++)
    {
        state[0 * n + i] = values[i]; //values[i];
        state[1 * n + i] = d[i];
        state[2 * n + i] = (3 * (values[i + 1] - values[i]) / (x[i + 1] - x[i]) - 2 * d[i] - d[i + 1]) / (x[i + 1] - x[i]);
        state[3 * n + i] = (d[i] + d[i + 1] - 2 * (values[i + 1] - values[i]) / (x[i + 1] - x[i])) / (x[i + 1] - x[i]) / (x[i + 1] - x[i]);
    }
/*    for (int j = 0; j < 4; j++)
    {
        for (int i = 0; i < n; i++)
            qDebug() << state[j * n + i];
        qDebug() << "!!!";
    }
*/
}

/*double cubicSplines::method_init2v(double *x, double *y, double *values, int n, int m, int N)
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
        double string[n], res[N];
        for (int k = 0; k < n; k++)
            string[k] = gamma[k + j * L];
        calc_coefficients(c, d, y, string, n, N, gamma + j * L);
    }
}
*/
double cubicSplines::method_compute(double y)
{
    for (int i = 0; i < n - 1; i++)
    {
        if (x[i] <= y && y <= x[i + 1])
            return state[0 * n + i] + state[1 * n + i] * (y - x[i]) +
            state[2 * n + i] * (y - x[i]) * (y - x[i]) + state[3 * n + i] * (y - x[i]) * (y - x[i]) * (y - x[i]);
    }
    return 0;
}

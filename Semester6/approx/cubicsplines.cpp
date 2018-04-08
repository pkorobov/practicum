#include "cubicsplines.h"
#include <QTime>

cubicSplines::cubicSplines()
{
    state = new double[4 * n];
    memset(state, 0, 4 * n * sizeof(double));
}

cubicSplines::~cubicSplines()
{
    delete[] state;
}


cubicSplines::cubicSplines(double (*f_)(double), double a_, double b_, int n_)
{
    f = f_;
    a = a_;
    b = b_;
    n = n_;
    state = new double[4 * n];
    x = new double[n];
    memset(state, 0, 4 * n * sizeof(double));
}

double cubicSplines::method_init(double *values, double *derivatives)
{
    double *X;
    double *d;

    d = new double[n];
    X = new double[n * (n + 1)];
    memset(X, 0, n * n * sizeof(double));

    double delta = (b - a) / (n - 1);
    for (int i = 0; i <= n - 1; i++)
        x[i] = a + delta * i;

    X[0 + 0 * n] = 1;
    X[0 + n * n] = derivatives[0];

    X[n - 1 + n * (n - 1)] = 1;
    X[n + n * (n - 1)] = derivatives[n];
    for (int i = 1; i < n - 1; i++)
    {
        X[i + n * (i - 1)] = x[i + 1] - x[i];
        X[i + n * i] = 2 * (x[i + 1] - x[i - 1]);
        X[i + n * (i + 1)] = x[i] - x[i - 1];
        X[i + n * n] = 3 * (f(x[i]) - f(x[i - 1])) / (x[i] - x[i - 1]) * (x[i + 1] - x[i]) +
                       3 * (f(x[i + 1]) - f(x[i])) / (x[i + 1] - x[i]) * (x[i] - x[i - 1]);
    }
    for (int k = 0; k < n; k++)
        for (int i = k + 1; k < n; k++)
        {
            double t = X[i + k * n] / X[k + k * n];
            X[i + k * n] = 0;
            for (int j = k + 1; j < n; j++)
                X[i + j * n] -= t * X[k + j * n];
        }

    for (int i = n - 1; i >= 0; i--)
    {
        d[i] = X[i + n * n];
        for (int j = i + 1; j < n; j++)
            d[i] -= d[j] * X[i + j * n];
        d[i] /= X[i + i * n];
    }
    for (int i = 0; i < n - 1; i++)
    {
        state[i + 0 * n] = f(x[i]);//values[i];
        state[i + 1 * n] = derivatives[i];
        state[i + 2 * n] = 3 * (f(x[i + 1]) - f(x[i])) / (x[i + 1] - x[i]) - 2 * d[i] - d[i + 1];
        state[i + 3 * n] = (d[i] + d[i + 1] - 2 * (f(x[i + 1]) - f(x[i])) / (x[i + 1] - x[i])) / (x[i + 1] - x[i]) / (x[i + 1] - x[i]);
    }
}

double cubicSplines::method_compute(double y)
{

    for (int i = 0; i < n - 1; i++)
    {
        if (x[i] <= y && y <= x[i + 1])
            return state[i] + state[i + n] * (y - x[i]) +
            state[i + 2 * n] * (y - x[i]) * (y - x[i]) + state[i + 3 * n] * (y - x[i]) * (y - x[i]) * (y - x[i]);
    }
}

#include "cubicsplines.h"
#include <QTime>
#include <QDebug>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "functions.h"

using namespace std;

cubicSplines::cubicSplines()
{
    state = new double[4 * n];
    memset(state, 0, 4 * n * sizeof(double));
}

cubicSplines::~cubicSplines()
{
    delete[] state;
}

cubicSplines::cubicSplines(double a_, double b_, int n_, double *x_, double *derx_, double *values_,
                           double c_, double d_, int m_, double *y_, double *dery_)
    : abstractMethod(a_, b_, n_, x_, values_, c_, d_, m_, y_)
{
    derx = derx_;
    dery = dery_;
    state = new double[4 * n];
    memset(state, 0, 4 * n * sizeof(double));
}

void cubicSplines::calc_coefficients(double *x, double *values, int n, double * derivatives, double *d)
{
    double *X;

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

    for (int k = 0; k < n; k++)
        for (int i = k + 1; i < n; i++)
        {
            double t = X[i * n + k] / X[k * n + k];
            for (int j = k + 1; j < n; j++)
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
}

void cubicSplines::method_init2v(double *x, double *y, double *values, int n, int m, double *derx, double *dery)
{
    double *Fx, *Fy, *Fxy;
    int L = max(n, m);
    Fx = new double[L * L];
    Fy = new double[L * L];
    Fxy = new double[L * L];
    double der[2] = {0, 0};

    for (int j = 0; j < m; j++)
    {
        double string[n], res[n];

        der[0] = derx[0 * m + j];
        der[1] = derx[1 * m + j];

        for (int i = 0; i < n; i++)
            string[i] = values[i + j * n];
        calc_coefficients(x, string, n, der, res);
//Решение при фиксированныех у
        for (int i = 0; i < n; i++)
            Fx[i + j * n] = res[i];
    }

    for (int i = 0; i < n; i++)
    {
        der[0] = dery[0 * n + i];
        der[1] = dery[1 * n + i];

        double col[m], res[m];
        for (int j = 0; j < m; j++)
            col[j] = values[i + j * n];

        calc_coefficients(y, col, m, der, res);
//Решение при фиксированныех х
        for (int j = 0; j < m; j++)
            Fy[i + j * n] = res[j];
    }


    for (int i = 0; i < n; i++)
    {
        der[0] = dery[0 * n + i];
        der[1] = dery[1 * n + i];

        double col[m], res[m];
        for (int j = 0; j < m; j++)
            col[j] = Fx[i + j * n];
        calc_coefficients(y, col, m, der, res);
        for (int j = 0; j < m; j++)
            Fxy[i + j * n] = res[j];
    }

    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < m - 1; j++)
        {

            double h = x[i + 1] - x[i];
            double A1[16] = {1, 0, 0, 0,
                        0, 1, 0, 0,
                        -3/h/h, -2/h, 3/h/h, -1/h,
                          2/h/h/h, 1/h/h, -2/h/h/h, 1/h/h};
             h = y[j + 1] - y[j];

             double A2[16] = {1, 0, 0, 0,
                         0, 1, 0, 0,
                         -3/h/h, -2/h, 3/h/h, -1/h,
                           2/h/h/h, 1/h/h, -2/h/h/h, 1/h/h};             

             double Fij[16] = {values[i + j * n], Fy[i + j * n], values[i + (j + 1) * n], Fy[i + (j + 1) * n],
                           Fx[i + j * n], Fxy[i + j * n], Fx[i + (j + 1) * n], Fxy[i + (j + 1) * n],
                           values[i + 1 + j * n], Fy[i + 1 + j * n], values[i + 1 + (j + 1) * n], Fy[i + 1 + (j + 1) * n],
                           Fx[i + 1 + j * n], Fxy[i + 1 + j * n], Fx[i + 1 + (j + 1) * n], Fxy[i + 1 + (j + 1) * n]};

             double res[16];
/*
             qDebug() << "Fij" << i << " " << j;
             for (int k = 0; k < 4; k++)
                 qDebug() << Fij[0 + k * 4] << Fij[1 + k * 4] << Fij[2 + k * 4] << Fij[3 + k * 4];
*/

             mult(A1, Fij, res, 4);
             trans(A2, 4);
             mult(res, A2, gamma[i][j], 4);
             trans(gamma[i][j], 4);
        }
        delete[] Fx;
        delete[] Fy;
        delete[] Fxy;
}

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

double cubicSplines::method_compute2v(double x1, double x2)
{
    x1 = min(x1, b);
    x1 = max(x1, a);

    x2 = min(x2, d);
    x2 = max(x2, c);

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < m - 1; j++)
            if (x[i] <= x1 && x1 <= x[i + 1]
             && y[j] <= x2 && x2 <= y[j + 1])
             {

                double value = gamma[i][j][0 + 0 * 4] + gamma[i][j][1 + 0 * 4] * (x1 - x[i]) + gamma[i][j][2 + 0 * 4] * (x1 - x[i]) * (x1 - x[i]) +
                       gamma[i][j][3 + 0 * 4] * (x1 - x[i]) * (x1 - x[i]) * (x1 - x[i]) +
                       gamma[i][j][0 + 1 * 4] * (x2 - y[j]) + gamma[i][j][0 + 2 * 4] * (x2 - y[j]) * (x2 - y[j]) +
                       gamma[i][j][0 + 3 * 4] * (x2 - y[j]) * (x2 - y[j]) * (x2 - y[j]) +
                       gamma[i][j][1 + 1 * 4] * (x1 - x[i]) * (x2 - y[j]) + gamma[i][j][1 + 2 * 4] * (x1 - x[i]) * (x2 - y[j]) * (x2 - y[j]) +
                       gamma[i][j][1 + 3 * 4] * (x1 - x[i]) * (x2 - y[j]) * (x2 - y[j]) * (x2 - y[j]) +
                       gamma[i][j][2 + 1 * 4] * (x1 - x[i]) * (x1 - x[i]) * (x2 - y[j]) +
                       gamma[i][j][2 + 2 * 4] * (x1 - x[i]) * (x1 - x[i]) * (x2 - y[j]) * (x2 - y[j]) +
                       gamma[i][j][2 + 3 * 4] * (x1 - x[i]) * (x1 - x[i]) * (x2 - y[j]) * (x2 - y[j]) * (x2 - y[j]) +
                       gamma[i][j][3 + 1 * 4] * (x1 - x[i]) * (x1 - x[i]) * (x1 - x[i]) * (x2 - y[j]) +
                       gamma[i][j][3 + 2 * 4] * (x1 - x[i]) * (x1 - x[i]) * (x1 - x[i]) * (x2 - y[j]) * (x2 - y[j]) +
                       gamma[i][j][3 + 3 * 4] * (x1 - x[i]) * (x1 - x[i]) * (x1 - x[i]) * (x2 - y[j]) * (x2 - y[j]) * (x2 - y[j]);
//                qDebug() << i << j << value;
                return value;
            }
    }
    return 0;
}

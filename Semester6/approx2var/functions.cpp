#include "functions.h"
#include <cmath>
#include <algorithm>

void mult(double *a, double *b, double *c, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            c[i * n + j] = 0;
            for (int k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
        }
}

void trans(double *a, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
        {
            std::swap(a[i + j * n], a[j + i * n]);
        }
}

double f_0 (double x, double y)
{
  return x + 2 * y + 1;
}

double dfx_0 (double x, double y)
{
  return 1 + x * 0 + y * 0;
}


double dfy_0 (double x, double y)
{
  return 2 + x * 0 + y * 0;
}


double f_1 (double x, double y)
{
  return x * x * x + y * y;
}

double dfx_1 (double x, double y)
{
  return 3 * x * x + y * 0;
}

double dfy_1 (double x, double y)
{
  return 2 * y + x * 0;
}

double f_2 (double x, double y)
{
    return tanh(x) * y;
}

double dfx_2 (double x, double y)
{
  return y / cosh(x) / cosh(x);
}

double dfy_2 (double x, double y)
{
  return tanh(x) + y * 0;
}

double f_3 (double x, double y)
{
  return exp(2 * x) * sin(5 * y) / (3 + x + y * y);
}

double dfx_3 (double x, double y)
{
  return (exp(2 * x) * (5 + 2 * x + 2 * y * y) * sin(5 * y)) / (3 + x + y * y) / (3 + x + y * y);
}

double dfy_3 (double x, double y)
{
  return (exp(2 * x) * (5 * (3 + x + y * y) * cos(5 * y) - 2 * y * sin(5 * y))) / (3 + x + y * y) / (3 + x + y * y);
}


double f_4 (double x, double y)
{
  return fabs(x) + fabs(y);
}

double dfx_4 (double x, double y)
{
    if (x > 0)
      return 1 + y * 0;
    else
      return -1;
}

double dfy_4 (double x, double y)
{
    return 0 * x * y;
}

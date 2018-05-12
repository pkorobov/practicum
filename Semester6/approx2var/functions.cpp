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
  return 1;
}


double dfy_0 (double x, double y)
{
  return 2;
}


double f_1 (double x, double y)
{
  return x * x * x + y * y;
}

double dfx_1 (double x, double y)
{
  return 3 * x * x;
}

double dfy_1 (double x, double y)
{
  return 2 * y;
}

double f_2 (double x, double y)
{
    return tanh(x) * y;
}

double dfx_2 (double x, double y)
{
  return 0;
}

double dfy_2 (double x, double y)
{
  return 0;
}

double f_3 (double x, double y)
{
  return exp(2 * x) * sin(5 * y) / (3 + x + y * y);
}

double dfx_3 (double x, double y)
{
  return ((1 + x * x) * cos(x) - x * sin(x)) / pow(1 + x * x, 1.5);
}

double dfy_3 (double x, double y)
{
  return ((1 + x * x) * cos(x) - x * sin(x)) / pow(1 + x * x, 1.5);
}


double f_4 (double x)
{
  return abs(x);
}

double df_4 (double x)
{
    if (x > 0)
      return 1;
    else
      return -1;
}

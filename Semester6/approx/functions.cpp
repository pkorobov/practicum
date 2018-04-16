#include "functions.h"
#include "window.h"
#include <cmath>
#include <QCoreApplication>
#include <QDebug>

double f_0 (double x)
{
  return x;
}

double df_0 (double x)
{

  return x * 0 + 1;
}

double f_1 (double x)
{
  return x * x * x;
}

double df_1 (double x)
{
  return 2 * x * x;
}

double f_2 (double x)
{
  return sin(x) / sqrt(1 + x * x);
}

double df_2 (double x)
{
  return ((1 + x * x) * cos(x) - x * sin(x)) / pow(1 + x * x, 1.5);
}

double f_3 (double x)
{
  return exp(-x*x);
}

double df_3 (double x)
{
  return -2 * x * exp(-x*x);
}

double f_4 (double x)
{
  return qAbs(x);
}

double df_4 (double x)
{
    if (x > 0)
      return 1;
    else
      return -1;
}

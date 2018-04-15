#include "functions.h"
#include "window.h"
#include <QtMath>
#include <QCoreApplication>
#include <QDebug>

double f_0 (double x)
{
  return x;
}

double df_0 (double x)
{
  return 1;
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
  return qSin(x) / qSqrt(1 + x * x);
}

double df_2 (double x)
{
  return ((1 + x * x) *qCos(x) - x * qSin(x)) / qPow(1 + x * x, 1.5);
}

double f_3 (double x)
{
  return qExp(-x*x);
}

double df_3 (double x)
{
  return -2 * x * qExp(-x*x);
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

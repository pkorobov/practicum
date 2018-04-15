#include <QPainter>
#include <QPen>
#include <QDebug>
#include <QtMath>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "window.h"
#include "chebyshovlsm.h"
#include "cubicsplines.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define DEFAULT_STEPS 100

using namespace std;

double f_0 (double x)
{
  return x;
}

double f_1 (double x)
{
  return x * x * x;
}

double f_2 (double x)
{
  return qSin(x) / qSqrt(1 + x * x);
}

double f_3 (double x)
{
  return qExp(-x*x);
}

double f_4 (double x)
{
  return qAbs(x);
}

Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;

//  steps = DEFAULT_STEPS;
  x = new double[1000];

  double delta_x = (b - a) / (n - 1);
  for (int i = 0; i < n; i++)
      x[i] = a + i * delta_x;

  steps = width();

  func_id = 0;

  change_func();
}

QSize Window::minimumSizeHint () const
{
  return QSize(100, 100);
}

QSize Window::sizeHint () const
{
  return QSize(1000, 1000);
}

int Window::parse_command_line(int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
  {
//    QString file = argv[2];
    fromFile = true;
    ifstream myfile;
    myfile.open(argv[2]);
    myfile >> n;
    for (int i = 0; i < n; i++)
    {
        myfile >> x[i];
        myfile >> values[i];
    }
    myfile.close();
  }
  char c;
  if (sscanf(argv[1], "%lf", &a) != 1
    || sscanf(argv[2], "%lf", &b) != 1
    || b - a < 1.e-6
    || argc > 3 && sscanf(argv[3], "%d", &n) != 1
    || n <= 0)
    return -2;
  else
  {
      double delta_x = (b - a) / (n - 1);
      for (int i = 0; i < n; i++)
         x[i] = a + i * delta_x;
  }
  qDebug() << a;
  qDebug() << b;
  return 0;
}

/// change current function for drawing
void Window::change_func()
{
  func_id = (func_id + 1) % 5;

  switch (func_id)
  {
      case 0:
        f_name = "f (x) = x";
        f = f_0;

        break;
      case 1:
        f_name = "f (x) = x * x * x";
        f = f_1;
        break;
      case 2:
        f_name = "f (x) = Sin(x) / Sqrt(1 + x*x)";
        f = f_2;
        break;
      case 3:
        f_name = "f (x) = Exp(-x*x)";
        f = f_3;
        break;
      case 4:
        f_name = "f (x) = |x|";
        f = f_4;
        break;
    }
    update ();
}

/// render graph
void Window::paintEvent(QPaintEvent * /* event */)
{
  QPainter painter (this);
  double x1, x2, y1, y2;
  double max_y, min_y;
  double *values;

  values = new double[n];

  for (int i = 0; i < n; i++)
       values[i] = f(x[i]);

  ChebyshovLSM algo1(a, b, n, x, values, 5);
  cubicSplines algo2(a, b, n, x, values);
  double arr[n + 1];
  for (int i = 0; i <= n; i++)
      arr[i] = 1;
  arr[0] = -1;
  arr[1] = -1;

  algo1.method_init();
  algo2.method_init(arr);

  max_y = min_y = 0;
  double delta_x = (b - a) / steps;
  double delta_y;

  if (fromFile)
      for (int i = 0; i < n; i++)
      {
          y1 = values[i];
          if (y1 < min_y)
            min_y = y1;
          if (y1 > max_y)
            max_y = y1;
      }
  else
      for (x1 = a; x1 - b < 1.e-6; x1 += delta_x)
      {
          y1 = f(x1);
          if (y1 < min_y)
            min_y = y1;
          if (y1 > max_y)
            max_y = y1;
      }
  for (x1 = a; x1 - b < 1.e-6; x1 += delta_x)
    {
      y2 = algo1.method_compute(x1);
      double y3 = algo2.method_compute(x1);
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;

      if (y2 < min_y)
        min_y = y2;
      if (y2 > max_y)
       max_y = y2;

      if (y3 < min_y)
       min_y = y3;
      if (y3 > max_y)
       max_y = y3;
    }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;

  painter.setPen(QPen(Qt::black, 0));
  // save current Coordinate System
  painter.save();

  // make Coordinate Transformations
  painter.translate(0.5 * width (), 0.5 * height ());
  painter.scale(width() / (b - a), -height() / (max_y - min_y));
  painter.translate(-0.5 * (a + b), -0.5 * (min_y + max_y));

  painter.drawLine(a - 1, 0, b + 1, 0);
  painter.drawLine(0, min_y - 1, 0, max_y + 1);

  painter.setPen(QPen(Qt::red, 0));
  qDebug() << min_y << max_y;
  // draw approximated line for graph
  if (fromFile)
  {
      x1 = x[0];
      y1 = values[0];
      for (int i = 1; i < n - 1; i++)
      {
          x2 = x[i];
          y2 = values[i];
          painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
          x1 = x2;
          y1 = y2;
      }
      x2 = x[n - 1];
      y2 = values[n - 1];
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
  }
  else
  {
      x1 = a;
      y1 = f(x1);
      for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
      {
          y2 = f(x2);
          painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

          x1 = x2, y1 = y2;
      }
      x2 = b;
      y2 = f(x2);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
  }

  painter.setPen(QPen(Qt::green, 0));

  delta_x = (b - a) / steps;

  x1 = a;
  y1 = algo1.method_compute(x1);
  for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
  {
      y2 = algo1.method_compute(x2);
      painter.drawLine (QPointF ((x1) , y1), QPointF (x2, y2));
      x1 = x2, y1 = y2;
  }
  x2 = b;
  y2 = algo1.method_compute(x2);
  painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));  // draw axis

  painter.setPen(QPen(Qt::blue, 0));

  x1 = a;
  y1 = algo2.method_compute(x1);
  for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
  {
      y2 = algo2.method_compute(x2);
   //  qDebug() << x2 << ", " << y2;
      painter.drawLine(QPointF ((x1) , y1), QPointF (x2, y2));
      x1 = x2, y1 = y2;
  }
  x2 = b;
  y2 = algo2.method_compute(x2);


  painter.drawLine(QPointF (x1, y1), QPointF (x2, y2));

  // restore previously saved Coordinate System
  painter.restore();

  // render function name

  painter.setPen("blue");
  painter.drawText(0, 20, f_name);
//  delete[] values;
}

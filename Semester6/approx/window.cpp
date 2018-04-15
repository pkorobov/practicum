#include <QPainter>
#include <QPen>
#include <QDebug>
#include <QtMath>
#include <QFile>
#include <QTextStream>
#include <QSpinBox>
#include <QPushButton>
#include <QMainWindow>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "window.h"
#include "chebyshovlsm.h"
#include "cubicsplines.h"
#include "functions.h"


#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define DEFAULT_STEPS 100

using namespace std;

Window::Window(QWidget *parent)
  : QWidget(parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;

  func_id = 0;
  change_func();

  steps = width();

  QPushButton *double_button = new QPushButton(this);
  double_button->move(120, 5);
  double_button->setText("Points x 2");

  QPushButton *halve_button = new QPushButton(this);
  halve_button->move(215, 5);
  halve_button->setText("Points / 2");

  QSpinBox *delta_box = new QSpinBox(this);
  delta_box->move(310, 5);

  //  delta->setText("Delta +/-");

  connect(double_button, SIGNAL(clicked()), this, SLOT(doublePoints()));
  connect(halve_button, SIGNAL(clicked()), this, SLOT(halvePoints()));
  connect(delta_box, SIGNAL(valueChanged(int)), this, SLOT(addDelta(int)));

  delta_box->setMinimum(-1);
  delta_box->setMaximum(n - 1);
}

QSize Window::minimumSizeHint() const
{
  return QSize(100, 100);
}

QSize Window::sizeHint() const
{
  return QSize(800, 800);
}

void Window::addDelta(int i)
{
    deltaPoint = i;
    update();
}

void Window::doublePoints()
{
    if (fromFile)
        return;

    n *= 2;
    update();
}

void Window::halvePoints()
{
    if (fromFile)
        return;

    if (n > 1)
        n /= 2;
    update();
}

int Window::parse_command_line(int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
  {
//    QString file = argv[2];
    fromFile = true;
    change_func();
    filename = argv[1];
    return 0;
  }
  char c;
  if (sscanf(argv[1], "%lf", &a) != 1
    || sscanf(argv[2], "%lf", &b) != 1
    || b - a < 1.e-6
    || argc > 3 && sscanf(argv[3], "%d", &n) != 1
    || n <= 0)
    return -2;
  return 0;
}

/// change current function for drawing
void Window::change_func()
{
  if (fromFile)
  {
       f_name = "Data from file";
       update();
       return;
  }
  func_id = (func_id + 1) % 5;

  switch (func_id)
  {
      case 0:
        f_name = "f (x) = x";
        f = f_0;
        df = df_0;
        break;
      case 1:
        f_name = "f (x) = x * x * x";
        f = f_1;
        df = df_1;
        break;
      case 2:
        f_name = "f (x) = Sin(x) / Sqrt(1 + x*x)";
        f = f_2;
        df = df_2;
        break;
      case 3:
        f_name = "f (x) = Exp(-x*x)";
        f = f_3;
        df = df_3;
        break;
      case 4:
        f_name = "f (x) = |x|";
        f = f_4;
        df = df_4;
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
  x = new double[n];
  values = new double[n];
  initValues(fromFile);

  if (deltaPoint > -1)
      values[deltaPoint] += 1;

  ChebyshovLSM algo1(a, b, n, x, values, 12);
  cubicSplines algo2(a, b, n, x, values);

  algo1.method_init();
  algo2.method_init(der);

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
      for (x1 = a; x1 < b - 1.e-6; x1 += delta_x)
      {
          y1 = f(x1);
          if (y1 < min_y)
            min_y = y1;
          if (y1 > max_y)
            max_y = y1;
      }
  for (x1 = a; x1 < b - 1.e-6; x1 += delta_x)
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
          x1 = x2, y1 = y2;
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
  for (x2 = x1 + delta_x; x2 < b - 1.e-6; x2 += delta_x)
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

  delta_x = 0.01;
  for (x2 = x1 + delta_x; x2 < b - 1.e-6; x2 += delta_x)
  {
      y2 = algo2.method_compute(x2);
      painter.drawLine(QPointF((x1) , y1), QPointF(x2, y2));
      x1 = x2, y1 = y2;
  }
  x2 = b;
  y2 = algo2.method_compute(x2);

  painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

  // restore previously saved Coordinate System
  painter.restore();

  // render function name

  painter.setPen("blue");
  painter.drawText(0, 20, f_name);
//  delete[] values;
}

int Window::initValues(bool fromFile)
{
    if (!fromFile)
    {
        double delta_x = (b - a) / (n - 1);
        for (int i = 0; i < n; i++)
        {
            x[i] = a + i * delta_x;
            values[i] = f(x[i]);
        }
        der[0] = df(x[0]);
        der[1] = df(x[n - 1]);
    }
    else
    {
        QFile fileIn(filename);
        if(fileIn.open(QIODevice::ReadOnly))
        {
            QTextStream fin(&fileIn);
            fin >> n;
            for (int i = 0; i < n; i++)
            {
                fin >> x[i];
                fin >> values[i];
            }
            fin  >> der[0];
            fin >> der[1];

            a = x[0];
            b = x[n - 1];
            fileIn.close();
            ifstream myfile;
        }
    }
}

#include <QtGui>
#include "window.h"
#include <QDebug>
#include <cmath>
#include <algorithm>
#include "chebyshovlsm.h"
#include "functions.h"

using namespace std;

double f(double x, double y)
{
    return  sin(x * y) / sqrt(1 + x * x * y * y);
}

Window::  Window(QWidget* pwgt/*= 0*/) : QGLWidget(pwgt), m_xRotate(-90), m_zRotate(0)
{
//    a = c = -1;
  //  b = d = 1;
    a = -1, b = 1, c = -1, d = 1;
    n = 10, m = 10, N = 5;
    func_id = 0;
}

int Window::parse_command_line(int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
  {
//    QString file = argv[2];
    fromFile = true;
    //change_func();
    filename = argv[1];
    return 0;
  }
  if (sscanf(argv[1], "%lf", &a) != 1
    || sscanf(argv[2], "%lf", &b) != 1
    || sscanf(argv[3], "%lf", &c) != 1
    || sscanf(argv[4], "%lf", &d) != 1
    || b - a < 1.0e-6 || d - c < 1.0e-6
    || (argc > 6 && (sscanf(argv[5], "%d", &n) != 1 || sscanf(argv[6], "%d", &m) != 1))
    || n <= 0 || m <= 0
    || (argc > 7 && sscanf(argv[6], "%d", &N) != 1)
    || N <= 0)
    return -2;
  return 0;
}

/*virtual*/void   Window::initializeGL()
{
    qglClearColor(Qt::gray);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    int size1 = n;
    int size2 = m;
    x = new double[size1];
    y = new double[size2];
    values = new double[size1 * size2];

    double delta_x = (b - a) / (n - 1);
    double delta_y = (d - c) / (m - 1);

    ChebyshovLSM algo(a, b, n, NULL, NULL, N);

    change_func();

    for (int i = 0; i < n; i++)
        x[i] = a + i * delta_x;
    for (int j = 0; j < m; j++)
        y[j] = c + j * delta_y;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
        {
            values[i + j * n] = f(x[i], y[j]);
        }

    algo.method_init2v(x, y, values, n, m, N);


    delta_x = (b - a) / (size1 - 1);
    delta_y = (d - c) / (size2 - 1);

    for (int i = 0; i < size1; i++)
        x[i] = a + i * delta_x;
    for (int j = 0; j < size2; j++)
        y[j] = c + j * delta_y;
    for (int i = 0; i < size1; i++)
        for (int j = 0; j < size2; j++)
        {
            values[i + j * size1] = algo.method_compute2v(x[i], y[j]);
        }
    algo.method_init2v(x, y, values, n, m, 10);
    graph = drawGraph(f);
    graph2 = drawGraph(x, y, values, size1, size2);
}

/*virtual*/void   Window::resizeGL(int nWidth, int nHeight)
{
    glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(min(a, c), max(b, d), min(a, c), max(b, d), 1.0, 10.0);
}

/*virtual*/void   Window::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    drawAxis();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -3.0);

    glRotatef(m_xRotate, 1.0, 0.0, 0.0);
    glRotatef(m_zRotate, 0.0, 0.0, 1.0);

    calculate();

    graph = drawGraph(f);
    graph2 = drawGraph(x, y, values, n, m);

    glCallList(graph);
    glCallList(graph2);
}

/*virtual*/void   Window::mousePressEvent(QMouseEvent* pe)
{
    m_ptPosition = pe->pos();
}

/*virtual*/void   Window::mouseMoveEvent(QMouseEvent* pe)
{
    m_xRotate += 180 * (GLfloat)(pe->y() - m_ptPosition.y()) / height();
    m_zRotate += 180 * (GLfloat)(pe->x() - m_ptPosition.x()) / width();
    updateGL();

    m_ptPosition = pe->pos();
}
void   Window::drawAxis()
{
    glLineWidth(1.0f);

    qglColor(Qt::red);
    glBegin(GL_LINES);

        glVertex3f( 1.0f,  0.0f,  0.0f);
        glVertex3f(-1.0f,  0.0f,  0.0f);

        qglColor(Qt::black);

        glVertex3f( 0.0f,  1.0f,  0.0f);
        glVertex3f( 0.0f, -1.0f,  0.0f);

        glVertex3f( 0.0f,  0.0f,  1.0f);
        glVertex3f( 0.0f,  0.0f, -1.0f);
    glEnd();
}

GLuint Window::drawGraph(double(*f)(double, double))
{
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
            qglColor(Qt::red);
            double delta_x = (b - a) / (n - 1);
            double delta_y = (d - c) / (m - 1);

            for (int i = 0; i < n - 1; i++)
            {
                glBegin(GL_QUADS);
                for (int j = 0; j < m - 1; j++)
                {
                    double x1 = a + delta_x * i;
                    double x2 = a + delta_x * (i + 1);
                    double y1 = c + delta_y * j;
                    double y2 = c + delta_y * (j + 1);
                    glVertex3f(x1, y1, f(x1, y1));
                    glVertex3f(x1, y2, f(x1, y2));
                    glVertex3f(x2, y2, f(x2, y2));
                    glVertex3f(x2, y1, f(x2, y1));
                }
                glEnd();
            }
            glLineWidth(1.5f);
            qglColor(Qt::black);
            for (int i = 0; i < n; i++)
            {
                glBegin(GL_LINE_STRIP);
                for (int j = 0; j < m; j++)
                {
                    double x = a + delta_x * i;
                    double y = c + delta_y * j;
                    glVertex3f(x, y, f(x, y));
               }
               glEnd();

            }
            for (int j = 0; j < m; j++)
            {
                glBegin(GL_LINE_STRIP);
                for (int i = 0; i < n; i++)
                {
                    double x = a + delta_x * i;
                    double y = c + delta_y * j;
                    glVertex3f(x, y, f(x, y));
                }
                glEnd();
            }
    glEndList();

    return list;
}

GLuint Window::drawGraph(double *x, double *y, double *values, int n, int m)
{
    a = x[0];
    b = x[n - 1];
    c = y[0];
    d = y[n - 1];

    GLuint list = glGenLists(1);

    glNewList(list, GL_COMPILE);
            qglColor(Qt::green);


            for (int i = 0; i < n - 1; i++)
            {
                glBegin(GL_QUADS);
                for (int j = 0; j < m - 1; j++)
                {
                    glVertex3f(x[i], y[j], values[i + j * n]);
                    glVertex3f(x[i], y[j + 1], values[i + (j + 1) * n]);
                    glVertex3f(x[i + 1], y[j + 1], values[i + 1 + (j + 1) * n]);
                    glVertex3f(x[i + 1], y[j], values[i + 1 + j * n]);
                }
                glEnd();
            }
            glLineWidth(1.5f);
            qglColor(Qt::black);
            for (int i = 0; i < n; i++)
            {
                glBegin(GL_LINE_STRIP);
                for (int j = 0; j < m; j++)
                {
                    glVertex3f(x[i], y[j], values[i + j * n]);
               }
               glEnd();

            }
            for (int j = 0; j < m; j++)
            {
                glBegin(GL_LINE_STRIP);
                for (int i = 0; i < n; i++)
                {
                    glVertex3f(x[i], y[j], values[i + j * n]);
                }
                glEnd();
            }
    glEndList();

    return list;
}

int Window::initValues()
{
    if (!fromFile)
    {
        double delta_x = (b - a) / (n - 1);
        double delta_y = (d - c) / (m - 1);
        double *x, *y, *values;
        x = new double[n];
        y = new double[m];
        values = new double[n * m];
        for (int i = 0; i < n - 1; i++)
            x[i] = a + i * delta_x;
        for (int j = 0; j < m - 1; j++)
            y[j] = c + j * delta_y;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                values[i + j * n] = f(x[i], y[j]);
            }

//        der[0] = df(x[0]);
//        der[1] = df(x[n - 1]);
    }
    else
    {
        QFile fileIn(filename);
        if(fileIn.open(QIODevice::ReadOnly))
        {
            QTextStream fin(&fileIn);
            fin >> n >> N;
            for (int i = 0; i < n; i++)
            {
/*                fin >> x[i];
                fin >> values[i];a
*/            }
/*            fin  >> der[0];
            fin >> der[1];

            a = x[0];
            b = x[n - 1];
*/            fileIn.close();
        }
    }
    return 0;
}

void Window::change_func()
{
//  if (fromFile)
  {
//       f_name = "Data from file";
//       update();
//       return;
  }
  func_id = (func_id + 1) % 3;

  switch (func_id)
  {
      case 0:
//        f_name = "f (x) = x";
        f = f_0;
//        df = df_0;
        break;
      case 1:
//        f_name = "f (x) = x * x * x";
        f = f_1;
//        df = df_1;
        break;
      case 2:
//        f_name = "f (x) = Sin(x) / Sqrt(1 + x*x)";
        f = f_2;
//        df = df_2;
        break;
    }
    update ();
}

void Window::calculate()
{
    double delta_x = (b - a) / (n - 1);
    double delta_y = (d - c) / (m - 1);

    ChebyshovLSM algo(a, b, n, NULL, NULL, N, c, d, m);

    for (int i = 0; i < n; i++)
        x[i] = a + i * delta_x;
    for (int j = 0; j < m; j++)
        y[j] = c + j * delta_y;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
        {
            values[i + j * n] = f(x[i], y[j]);
        }

    algo.method_init2v(x, y, values, n, m, N);

    double size1 = n, size2 = m;

    delta_x = (b - a) / (size1 - 1);
    delta_y = (d - c) / (size2 - 1);

    for (int i = 0; i < size1; i++)
        x[i] = a + i * delta_x;
    for (int j = 0; j < size2; j++)
        y[j] = c + j * delta_y;
    for (int i = 0; i < size1; i++)
        for (int j = 0; j < size2; j++)
        {
            values[i + j * n] = algo.method_compute2v(x[i], y[j]);
        }
}

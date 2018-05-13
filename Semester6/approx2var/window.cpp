#include <QtGui>
#include "window.h"
#include <QDebug>
#include <cmath>
#include <algorithm>
#include <QPushButton>
#include <QLabel>
#include <QFont>

using namespace std;

Window::Window(QWidget* pwgt) : QGLWidget(pwgt), m_xRotate(-90), m_zRotate(0)
{
    a = -1, b = 1, c = -1, d = 1;
    n = 10, m = 10, N = 5;
    func_id = 2;
    win_number = 0;
    win_name = "Original function";
    plot_original = false;
    scale = 1.0f;

    int s = 200;
    x = new double[s];
    y = new double[s];
    values = new double[s * s];

    derx = new double[2 * s];
    dery = new double[2 * s];
}

Window::~Window()
{
    delete[] x;
    delete[] y;
    delete[] values;
    delete[] derx;
    delete[] dery;
}

int Window::parse_command_line(int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (sscanf(argv[1], "%lf", &a) != 1
    || sscanf(argv[2], "%lf", &b) != 1
    || sscanf(argv[3], "%lf", &c) != 1
    || sscanf(argv[4], "%lf", &d) != 1
    || b - a < 1.0e-6 || d - c < 1.0e-6
    || (argc > 6 && (sscanf(argv[5], "%d", &n) != 1 || sscanf(argv[6], "%d", &m) != 1))
    || n <= 0 || m <= 0
    || (argc > 7 && sscanf(argv[7], "%d", &N) != 1)
    || N <= 0)
    return -2;
  return 0;
}

void Window::initializeGL()
{
    qglClearColor(Qt::gray);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    change_func();
}

void Window::resizeGL(int nWidth, int nHeight)
{
    glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double t = max(abs(a), abs(b));
    t = max(t, abs(c));
    t = max(t, abs(d));
    t += 2.0;
    glOrtho(-t, t, -t, t, -10.0, 10.0);
//    glFrustum(-t, t, -t, t, 1.0, 100.0);
}

void Window::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    drawAxis();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
   // glTranslatef(0.0, 0.0, -scale);

    glScalef(scale, scale, scale);
    glRotatef(m_xRotate, 1.0, 0.0, 0.0);
    glRotatef(m_zRotate, 0.0, 0.0, 1.0);
    for (int i = 0; i < 3; i++)
        if (win_number == i)
            glCallList(graph[i]);

    if (win_number == 3)
        glCallList(graphErr[0]);
    if (win_number == 4)
        glCallList(graphErr[1]);

    if (win_number > 0 && plot_original)
        glCallList(graph[0]);

    QFont font("Decorative", 8, QFont::Normal);

    qglColor(Qt::black);

    renderText(0, 20, f_name, font);
    if (win_number == 1 || win_number == 3)
        renderText(0, 40, "Grid size: " + QString::number(n) + "x" + QString::number(m) + ", [" + QString::number(a) + ", " +
        QString::number(b) + "] x [" + QString::number(c) + ", " + QString::number(d) + "]" + " , N = " + QString::number(N), font);
    else
        renderText(0, 40, "Grid size: " + QString::number(n) + "x" + QString::number(m) + ", [" + QString::number(a) + " ," +
        QString::number(b) + "] x [" + QString::number(c) + ", " + QString::number(d) + "]", font);

    renderText(0, 60, win_name, font);
    if (win_number == 3)
        renderText(0, 80, "Max absolute error: " + QString::number(maxErr1), font);
    if (win_number == 4)
        renderText(0, 80, "Max absolute error: " + QString::number(maxErr2), font);
}

void Window::mousePressEvent(QMouseEvent* pe)
{
    m_ptPosition = pe->pos();
}

void Window::mouseMoveEvent(QMouseEvent* pe)
{
    m_xRotate += 180 * (GLfloat)(pe->y() - m_ptPosition.y()) / height();
    m_zRotate += 180 * (GLfloat)(pe->x() - m_ptPosition.x()) / width();
    updateGL();

    m_ptPosition = pe->pos();
}

void Window::wheelEvent(QWheelEvent* pe)
{
    if ((pe->delta())>0)
        scale*=1.1;
    else if ((pe->delta())<0)
        scale/=1.1;

    updateGL();
}

void Window::drawAxis()
{
    glLineWidth(1.0f);

    qglColor(Qt::red);
    glBegin(GL_LINES);

        glVertex3f(a - 1,  0.0f,  0.0f);
        glVertex3f(b + 1,  0.0f,  0.0f);

        qglColor(Qt::blue);

        glVertex3f( 0.0f,  d + 1,  0.0f);
        glVertex3f( 0.0f, c - 1,  0.0f);

        qglColor(Qt::green);

        glVertex3f( 0.0f,  0.0f,  1.0);
        glVertex3f( 0.0f,  0.0f, -1.0);
    glEnd();
}


void Window::change_func()
{
  func_id = (func_id + 1) % 5;

  switch (func_id)
  {
      case 0:
        f_name = "f(x, y) = x + 2 * y + 1";
        f = f_0;
        dfx = dfx_0;
        dfy = dfy_0;
        break;
      case 1:
        f_name = "f(x, y) = x * x * x + y * y";
        f = f_1;
        dfx = dfx_1;
        dfy = dfy_1;
        break;
      case 2:
        f_name = "f(x, y) = tanh(x) * y";
        f = f_2;
        dfx = dfx_2;
        dfy = dfy_2;
        break;
      case 3:
        f_name = "f(x, y) = exp(2 * x) * sin(5 * y) / (3 + x + y * y)";
        f = f_3;
        dfx = dfx_3;
        dfy = dfy_3;
        break;
      case 4:
          f_name = "f(x, y) = |x| + |y|";
          f = f_4;
          dfx = dfx_4;
          dfy = dfy_4;
    }

    calculate();
    updateGL();
}

void Window::doubleN()
{
    if (N * 2 < 50)
    {
        N *= 2;
        calculate();
        updateGL();
    }
}

void Window::halveN()
{
    if (N / 2 > 0)
    {
        N /= 2;
        calculate();
        updateGL();
    }
}

void Window::doublePoints()
{
    if (n * 2 < 200 && m * 2 < 200)
    {
        m *= 2;
        n *= 2;
        calculate();
        updateGL();
    }
}

void Window::halvePoints()
{
    if (n / 2 > 1 && m / 2 > 1)
    {
        n /= 2;
        m /= 2;
        calculate();
        updateGL();
    }
}

void Window::closer()
{
    scale /= 1.1;
    updateGL();
}

void Window::further()
{
    scale *= 1.1;
    updateGL();
}

void Window::change_win()
{
    win_number = (win_number + 1) % 5;
    switch(win_number)
    {
        case 0:
            win_name = "Original function";
            break;
        case 1:
            win_name = "Chebyshov least square method";
            break;
        case 2:
            win_name = "Cubic splines method";
            break;
        case 3:
            win_name = "Error 1";
            break;
        case 4:
            win_name = "Error 2";
            break;
    }
    updateGL();
}

void Window::enable_original()
{
    plot_original = (plot_original + 1) % 2;
    updateGL();
}

GLuint Window::drawGraph(double(*f)(double, double), QColor color)
{
    GLuint list = glGenLists(1);

    int grid_split1 = (b - a) / 0.05;
    int grid_split2 = (d - c) / 0.05;
    int split = width();

    glNewList(list, GL_COMPILE);
            qglColor(color);
            double delta_x = (b - a) / (split - 1);
            double delta_y = (d - c) / (split - 1);

            for (int i = 0; i < split - 1; i++)
            {
                glBegin(GL_QUADS);
                for (int j = 0; j < split - 1; j++)
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
            for (int i = 0; i < split; i += grid_split1)
            {
                glBegin(GL_LINE_STRIP);
                for (int j = 0; j < split; j++)
                {
                    double x = a + delta_x * i;
                    double y = c + delta_y * j;
                    glVertex3f(x, y, f(x, y));
               }
               glEnd();

            }
            for (int j = 0; j < split; j += grid_split2)
            {
                glBegin(GL_LINE_STRIP);
                for (int i = 0; i < split; i++)
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


GLuint Window::drawGraph(abstractMethod & algo, QColor color)
{
    GLuint list = glGenLists(1);

    int grid_split1 = (b - a) / 0.05;
    int grid_split2 = (d - c) / 0.05;
    int split = width();

    glNewList(list, GL_COMPILE);
            qglColor(color);
            double delta_x = (b - a) / (split - 1);
            double delta_y = (d - c) / (split - 1);

            for (int i = 0; i < split - 1; i++)
            {
                glBegin(GL_QUADS);
                for (int j = 0; j < split - 1; j++)
                {
                    double x1 = a + delta_x * i;
                    double x2 = a + delta_x * (i + 1);
                    double y1 = c + delta_y * j;
                    double y2 = c + delta_y * (j + 1);
                    glVertex3f(x1, y1, algo.method_compute2v(x1, y1));
                    glVertex3f(x1, y2, algo.method_compute2v(x1, y2));
                    glVertex3f(x2, y2, algo.method_compute2v(x2, y2));
                    glVertex3f(x2, y1, algo.method_compute2v(x2, y1));
                }
                glEnd();
            }
            glLineWidth(1.5f);
            qglColor(Qt::black);
            for (int i = 0; i < split; i += grid_split1)
            {
                glBegin(GL_LINE_STRIP);
                for (int j = 0; j < split; j++)
                {
                    double x = a + delta_x * i;
                    double y = c + delta_y * j;
                    glVertex3f(x, y, algo.method_compute2v(x, y));
               }
               glEnd();

            }
            for (int j = 0; j < split; j += grid_split2)
            {
                glBegin(GL_LINE_STRIP);
                for (int i = 0; i < split; i++)
                {
                    double x = a + delta_x * i;
                    double y = c + delta_y * j;
                    glVertex3f(x, y, algo.method_compute2v(x, y));
                }
                glEnd();
            }
    glEndList();

    return list;
}

void Window::calculate()
{
    double delta_x = (b - a) / (n - 1);
    double delta_y = (d - c) / (m - 1);

    for (int i = 0; i < n; i++)
        x[i] = a + i * delta_x;
    for (int j = 0; j < m; j++)
        y[j] = c + j * delta_y;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
        {
            values[i + j * n] = f(x[i], y[j]);
        }

    for (int j = 0; j < m; j++)
    {
        derx[0 * m + j] = dfx(x[0], y[j]);
        derx[1 * m + j] = dfx(x[n - 1], y[j]);
    }

    for (int i = 0; i < n; i++)
    {
        dery[0 * n + i] = dfy(x[i], y[0]);
        dery[1 * n + i] = dfy(x[i], y[n - 1]);
    }

    ChebyshovLSM algo1(a, b, n, NULL, NULL, N, c, d, m);
    cubicSplines algo2(a, b, n, x, NULL, NULL, c, d, m, y);

    algo1.method_init2v(x, y, values, n, m, N);
    algo2.method_init2v(x, y, values, n, m, derx, dery);


    delta_x = (b - a) / (n - 1);
    delta_y = (d - c) / (m - 1);

    for (int i = 0; i < n; i++)
        x[i] = a + i * delta_x;
    for (int j = 0; j < m; j++)
        y[j] = c + j * delta_y;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
        {
            values[i + j * n] = f(x[i], y[j]);
        }
    graph[0] = drawGraph(f, QColor(200, 200, 50));
    graph[1] = drawGraph(algo1, QColor(100, 200, 50));
    graph[2] = drawGraph(algo2, QColor(100, 50, 100));

    graphErr[0] = drawError(f, algo1, QColor(100, 200, 50), maxErr1);
    graphErr[1] = drawError(f, algo2, QColor(100, 50, 100), maxErr2);
}

GLuint Window::drawError(double(*f)(double, double), abstractMethod & algo, QColor color, double & maxErr)
{
    maxErr = 0;

    GLuint list = glGenLists(1);

    int grid_split1 = (b - a) / 0.05;
    int grid_split2 = (d - c) / 0.05;
    int split = width();

    glNewList(list, GL_COMPILE);
            qglColor(color);
            double delta_x = (b - a) / (split - 1);
            double delta_y = (d - c) / (split - 1);

            for (int i = 0; i < split - 1; i++)
            {
                glBegin(GL_QUADS);
                for (int j = 0; j < split - 1; j++)
                {
                    double x1 = a + delta_x * i;
                    double x2 = a + delta_x * (i + 1);
                    double y1 = c + delta_y * j;
                    double y2 = c + delta_y * (j + 1);
                    glVertex3f(x1, y1, abs(f(x1, y1) - algo.method_compute2v(x1, y1)));
                    glVertex3f(x1, y2, abs(f(x1, y2) - algo.method_compute2v(x1, y2)));
                    glVertex3f(x2, y2, abs(f(x2, y2) - algo.method_compute2v(x2, y2)));
                    glVertex3f(x2, y1, abs(f(x2, y1) - algo.method_compute2v(x2, y1)));

                    maxErr = max(maxErr, abs(f(x1, y1) - algo.method_compute2v(x1, y1)));
                }
                glEnd();
            }
            glLineWidth(1.5f);
            qglColor(Qt::black);
            for (int i = 0; i < split; i += grid_split1)
            {
                glBegin(GL_LINE_STRIP);
                for (int j = 0; j < split; j++)
                {
                    double x = a + delta_x * i;
                    double y = c + delta_y * j;
                    glVertex3f(x, y, abs(f(x, y) - algo.method_compute2v(x, y)));
               }
               glEnd();

            }
            for (int j = 0; j < split; j += grid_split2)
            {
                glBegin(GL_LINE_STRIP);
                for (int i = 0; i < split; i++)
                {
                    double x = a + delta_x * i;
                    double y = c + delta_y * j;
                    glVertex3f(x, y, abs(f(x, y) - algo.method_compute2v(x, y)));
                }
                glEnd();
            }
    glEndList();

    return list;
}


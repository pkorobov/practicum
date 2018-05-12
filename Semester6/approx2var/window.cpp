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
    n = 10, m = 5, N = 5;
    func_id = 2;
    win_number = 0;
    win_name = "Original function";
    plot_original = false;
    fromFile = false;
    scale = 3.0f;

    int s = 200;
    x = new double[s];
    y = new double[s];
    values = new double[s * s];

    derx = new double[2 * s];
    dery = new double[2 * s];

    valuesAlgo1 = new double[s * s];
    valuesAlgo2 = new double[s * s];
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
    glFrustum(min(a, c), max(b, d), min(a, c), max(b, d), 1.0, 10.0);
}

void Window::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    drawAxis();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -scale);

    glRotatef(m_xRotate, 1.0, 0.0, 0.0);
    glRotatef(m_zRotate, 0.0, 0.0, 1.0);
    for (int i = 0; i < 3; i++)
        if (win_number == i)
            glCallList(graph[i]);
    if (win_number > 0 && plot_original)
        glCallList(graph[0]);

    QFont font("Decorative", 8, QFont::Normal);

    QString str1, str2;
    qglColor(Qt::black);

    renderText(0, 20, f_name, font);
    renderText(0, 40, "Grid size: " + str1.setNum(n) + "x" + str2.setNum(m), font);
    renderText(0, 60, win_name, font);
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
//  if (fromFile)
  {
//       f_name = "Data from file";
//       update();
//       return;
  }
  func_id = (func_id + 1) % 4;

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
    }

    calculate();
    updateGL();
}

void Window::doublePoints()
{
    if (fromFile)
        return;
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
    if (fromFile)
        return;

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
    scale -= 0.1;
    updateGL();
}

void Window::further()
{
    scale += 0.1;
    updateGL();
}

void Window::change_win()
{
    win_number = (win_number + 1) % 3;
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
    }
    updateGL();
}

void Window::enable_original()
{
    plot_original = (plot_original + 1) % 2;
    updateGL();
    qDebug() << plot_original;
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
            valuesAlgo1[i + j * n] = algo1.method_compute2v(x[i], y[j]);
            valuesAlgo2[i + j * n] = algo2.method_compute2v(x[i], y[j]);
        }
    graph[0] = drawGraph(f, QColor(200, 200, 50));
    graph[1] = drawGraph(algo1, QColor(100, 200, 50));
    graph[2] = drawGraph(algo2, QColor(100, 50, 100));

    /*    graph[1] = drawGraph(x, y, valuesAlgo1, n, m, QColor(100, 200, 50));
    graph[2] = drawGraph(x, y, valuesAlgo2, n, m, QColor(100, 50, 100));
*/
}

int Window::initValues()
{
    if (!fromFile)
    {
        double delta_x = (b - a) / (n - 1);
        double delta_y = (d - c) / (m - 1);
        double *x, *y;
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
                fin >> valuesAlgo1[i];a
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



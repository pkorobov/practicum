#pragma once

#include <qgl.h>
#include "chebyshovlsm.h"
#include "cubicsplines.h"
#include "functions.h"
#include <QString>

// ======================================================================
class   Window : public QGLWidget {
    Q_OBJECT
private:
    int n, m, N, size1, size2;
    double scale;
    double a, b, c, d;
    double *x, *y, *values, *valuesAlgo1, *valuesAlgo2, *derx, *dery;
    bool fromFile, plot_original;
    double(*f)(double, double);
    double(*dfx)(double, double);
    double(*dfy)(double, double);

    QString filename;
    QString f_name;
    QString win_name;
    GLuint  graph[3];
    GLfloat m_xRotate;
    GLfloat m_zRotate;
    QPoint  m_ptPosition;
    int initValues();
    void calculate();
    int func_id;
    int win_number;

protected:
    virtual void initializeGL();
    virtual void resizeGL(int nWidth, int nHeight);
    virtual void paintGL();
    virtual void mousePressEvent(QMouseEvent* pe);
    virtual void mouseMoveEvent (QMouseEvent* pe);
    void drawAxis();
    GLuint drawGraph(double(*f)(double, double), QColor color);
    GLuint drawGraph(abstractMethod & algo, QColor color);
    GLuint drawGraph(double *x, double *y, double *values, int n, int m, QColor color);

public:
    int parse_command_line(int argc, char *argv[]);
    Window(QWidget* pwgt = 0);
    virtual ~Window(){}

public slots:
    void enable_original();
    void change_win();
    void change_func();
    void doublePoints();
    void halvePoints();
    void closer();
    void further();
};


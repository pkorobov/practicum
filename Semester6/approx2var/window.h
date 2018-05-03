#pragma once

#include <QGLWidget>

// ======================================================================
class   Window : public QGLWidget {
    Q_OBJECT
private:
    int n, m, N;
    double a, b, c, d;
    double *x, *y, *values;
    bool fromFile;
    double(*f)(double, double);
    QString filename;
    GLuint  graph, graph2;
    GLfloat m_xRotate;
    GLfloat m_zRotate;
    QPoint  m_ptPosition;
    int initValues();
    void calculate();
    int func_id;

protected:
    virtual void initializeGL();
    virtual void resizeGL(int nWidth, int nHeight);
    virtual void paintGL();
    virtual void mousePressEvent(QMouseEvent* pe);
    virtual void mouseMoveEvent (QMouseEvent* pe);
    void drawAxis();
    GLuint drawGraph(double(*f)(double, double));
    GLuint drawGraph(double *x, double *y, double *values, int n, int m);

public:
    int parse_command_line(int argc, char *argv[]);
    Window(QWidget* pwgt = 0);
    virtual ~Window(){}

public slots:
    void change_func();
};


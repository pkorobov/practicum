// ======================================================================
//    Window.cpp
// ======================================================================
//                   This file is a part of the book 
//             "Qt 5.3 Professional programming with C++"
// ======================================================================
//  Copyright (c) 2014 by Max Schlee
//
//  Email : Max.Schlee@neonway.com
//  Blog  : http://www.maxschlee.com
//
//  Social Networks
//  ---------------
//  FaceBook : http://www.facebook.com/mschlee
//  Twitter  : http://twitter.com/Max_Schlee
//  2Look.me : http://2look.me/NW100003
//  Xing     : http://www.xing.com/profile/Max_Schlee
//  vk.com   : https://vk.com/max.schlee
// ======================================================================

#include <QtGui>
#include "window.h"
#include <QDebug>
#include <cmath>
#include <algorithm>

using namespace std;

double f(double x, double y)
{
    return x * x + y * y;
}

// ----------------------------------------------------------------------
Window::  Window(QWidget* pwgt/*= 0*/) : QGLWidget(pwgt), m_xRotate(0), m_yRotate(0)
{
    a = -1, b = 1, c = -1, d = 1;
    n = 10, m = 20;
}

// ----------------------------------------------------------------------
/*virtual*/void   Window::initializeGL()
{
    qglClearColor(Qt::gray);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    graph = drawGraph(f);
}

// ----------------------------------------------------------------------
/*virtual*/void   Window::resizeGL(int nWidth, int nHeight)
{
    glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(min(a, c), max(b, d), -1.0, 1.0, 1.0, 10.0);
}

// ----------------------------------------------------------------------
/*virtual*/void   Window::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    drawAxis();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -3.0);

    glRotatef(m_xRotate, 1.0, 0.0, 0.0);
    glRotatef(m_yRotate, 0.0, 0.0, 1.0);

    glCallList(graph);
}

// ----------------------------------------------------------------------
/*virtual*/void   Window::mousePressEvent(QMouseEvent* pe)
{
    m_ptPosition = pe->pos();
}

// ----------------------------------------------------------------------
/*virtual*/void   Window::mouseMoveEvent(QMouseEvent* pe)
{
    m_xRotate += 180 * (GLfloat)(pe->y() - m_ptPosition.y()) / height();
    m_yRotate += 180 * (GLfloat)(pe->x() - m_ptPosition.x()) / width();
    updateGL();

    m_ptPosition = pe->pos();
}
void   Window::drawAxis()
{
    glLineWidth(1.0f);

    qglColor(Qt::black);
    glBegin(GL_LINES);
        glVertex3f( 1.0f,  0.0f,  0.0f);
        glVertex3f(-1.0f,  0.0f,  0.0f);

        glVertex3f( 0.0f,  1.0f,  0.0f);
        glVertex3f( 0.0f, -1.0f,  0.0f);

        glVertex3f( 0.0f,  0.0f,  1.0f);
        glVertex3f( 0.0f,  0.0f, -1.0f);
    glEnd();
}
// ----------------------------------------------------------------------
GLuint Window::drawGraph(double(*f)(double, double))
{
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
            qglColor(Qt::green);
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

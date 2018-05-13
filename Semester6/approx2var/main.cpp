
#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>
#include "QDebug"

#include "window.h"

int main (int argc, char *argv[])
{
    QApplication app(argc, argv);
    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar = new QMenuBar(window);
    Window *graph_area = new Window(window);
    QAction *action;



    if (graph_area->parse_command_line(argc, argv))
    {
          QMessageBox::warning(0, "Wrong input arguments!", "Wrong input arguments!");
          return -1;
    }
    window->resize(400, 400);
//    graph_area->resize(400, 400);
    tool_bar->setMaximumHeight(50);

    action = tool_bar->addAction("Change function", graph_area, SLOT(change_func()));
    action->setShortcut(QString("Ctrl+C"));

    action = tool_bar->addAction("Exit", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    action = tool_bar->addAction("Change window", graph_area, SLOT(change_win()));
    action->setShortcut(QString("1"));

    action = tool_bar->addAction("Add delta", graph_area, SLOT(enable_delta()));
    action->setShortcut(QString("Ctrl+D"));

    action = tool_bar->addAction("Plot original", graph_area, SLOT(enable_original()));
    action->setShortcut(QString("0"));

    action = tool_bar->addAction("Halve points", graph_area, SLOT(halvePoints()));
    action->setShortcut(QString("Ctrl+O"));

    action = tool_bar->addAction("Double points", graph_area, SLOT(doublePoints()));
    action->setShortcut(QString("Ctrl+P"));

    action = tool_bar->addAction("Scale increase", graph_area, SLOT(further()));
    action->setShortcut(QString("Ctrl+U"));

    action = tool_bar->addAction("Scale decrease", graph_area, SLOT(closer()));
    action->setShortcut(QString("Ctrl+I"));

    action = tool_bar->addAction("Double N", graph_area, SLOT(doubleN()));
    action->setShortcut(QString("Ctrl+N"));

    action = tool_bar->addAction("Halve N", graph_area, SLOT(halveN()));
    action->setShortcut(QString("Ctrl+M"));

    window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Graph");

    window->show();

    qDebug() << "Ctrl+X - exit, Ctrl+C - change function, 1 - change window, Ctrl+P - double points";
    qDebug() << "Ctrl+O - halve points, Ctrl+U - scale increase, Ctrl+I - scale decrease";

    return app.exec();
}

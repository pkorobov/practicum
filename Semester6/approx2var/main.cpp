
#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>

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

    action = tool_bar->addAction("Change function", graph_area, SLOT(change_func()));
    action->setShortcut(QString("Ctrl+C"));

    action = tool_bar->addAction("Exit", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    tool_bar->setMaximumHeight(30);

    window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Graph");

    window->show();
    return app.exec();
}

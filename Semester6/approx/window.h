#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  const char *f_name;
  double a;
  double b;
  double *x, *values;
  double der[2];
  QString filename;

  bool fromFile = false;
  int n, steps;

public:
  double(*f) (double);
  double(*df) (double);
  int initValues(bool fromFile);
  Window(QWidget *parent);

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  int parse_command_line(int argc, char *argv[]);

public slots:
  void change_func();

protected:
  void paintEvent (QPaintEvent *event);
};

#endif

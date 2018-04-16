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
  int deltaPoint = -1;
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
signals:
  void nChanged(int n);
public slots:
  void change_func();
  void doublePoints();
  void halvePoints();
  void addDelta(int i);
protected:
  void paintEvent (QPaintEvent *event);
};

#endif

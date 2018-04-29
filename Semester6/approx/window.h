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
  int deltaPoint;
  QString filename;

  bool fromFile;
  int n, N, steps;

public:
  double(*f) (double);
  double(*df) (double);
  int initValues(bool fromFile);
  Window(QWidget *parent);

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  int parse_command_line(int argc, char *argv[]);
signals:
  void nChanged(QString);
  void error(QString);

public slots:
  void change_func();
  void doublePoints();
  void halvePoints();
  void addDelta(int);

protected:
  void paintEvent (QPaintEvent *event);
};

#endif

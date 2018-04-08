
QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = approx
TEMPLATE = app

DEFINES += QT_DEPRECATED_WARNINGS

SOURCES += \
        main.cpp \
        window.cpp \
    chebyshovlsm.cpp \
    cubicsplines.cpp

HEADERS += \
        window.h \
    chebyshovlsm.h \
    cubicsplines.h

FORMS += \
        window.ui

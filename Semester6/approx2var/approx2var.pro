TEMPLATE     = app
QT          += widgets opengl
SOURCES += \
        main.cpp \
        window.cpp \
    chebyshovlsm.cpp \
    cubicsplines.cpp \
    abstractmethod.cpp \
    functions.cpp

HEADERS += \
        window.h \
    chebyshovlsm.h \
    abstractmethod.h \
    cubicsplines.h \
    functions.h

TARGET = approx2var


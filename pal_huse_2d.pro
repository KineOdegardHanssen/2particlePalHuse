TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    ph_system.cpp \
    ph_evolve.cpp \
    ph_running.cpp

HEADERS += \
    ph_system.h \
    ph_evolve.h \
    ph_running.h

LIBS += -larmadillo -llapack -lblas
QMAKE_CXXFLAGS += -Wall -std=c++0x


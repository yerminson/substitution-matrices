#-------------------------------------------------
#
# Project created by QtCreator 2015-05-22T00:00:47
#
#-------------------------------------------------

QT       += core gui


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = BlosumAlgorithm
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    blosum.cpp \
    block.cpp \
    configurations.cpp

HEADERS  += mainwindow.h \
    blosum.h \
    block.h \
    configurations.h

FORMS    += mainwindow.ui

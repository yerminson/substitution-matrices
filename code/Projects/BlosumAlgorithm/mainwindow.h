#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
  void on_fileButton_clicked();

  void on_calculateMatrixButton_clicked();

  void on_noClusterRadioButton_clicked();

  void on_existingClusterRadioButton_clicked();

  void on_existingWeightsRadioButton_clicked();

  void on_PositionBasesWeigthsRadioButton_clicked();

  void on_clusterPercentageRadioButton_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H

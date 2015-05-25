#include "mainwindow.h"
#include "ui_mainwindow.h"



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::on_fileButton_clicked()
{
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                  "/home"
                                                  );
  ui->fileNameLabel->setText(fileName);
}

void MainWindow::on_calculateMatrixButton_clicked()
{
  QMessageBox validationBox;
  QString fileName = ui->fileNameLabel->text();
  if(fileName.isEmpty()){
    validationBox.setText("Please select a database file to create the BLOSUM matrix");
    validationBox.setIcon(QMessageBox::Critical);
    validationBox.exec();
  }
  int fileNameSize = fileName.size();
  char* fileNameValue = new char[fileNameSize+1];
  strcpy(fileNameValue,fileName.toStdString().c_str());
  qDebug() << fileNameValue;
  QString minimunStrengh = ui->minimunStrenghInputText->toPlainText();
  if(minimunStrengh.isEmpty()){
    validationBox.setText("Please write a value for the Minimun Strengh");
    validationBox.setIcon(QMessageBox::Critical);
    validationBox.exec();
  }
  bool ok = true;
  int minimunStrenghValue = minimunStrengh.toInt(&ok);
  qDebug() << minimunStrenghValue;
  if(!ok){
    validationBox.setText("Please write a numeric vaule for the Minimun Strengh");
    validationBox.setIcon(QMessageBox::Critical);
    validationBox.exec();
  }

  QString maximunStrengh = ui->maximunStrenghInputText->toPlainText();
  if(maximunStrengh.isEmpty()){
    validationBox.setText("Please write a value for the Maximun Strengh");
    validationBox.setIcon(QMessageBox::Critical);
    validationBox.exec();
  }
  ok = true;
  int maximunStrenghValue = maximunStrengh.toInt(&ok);
  qDebug() << maximunStrenghValue;
  if(!ok){
    validationBox.setText("Please write a numeric vaule for the Maximun Strengh");
    validationBox.setIcon(QMessageBox::Critical);
    validationBox.exec();
  }

  if(maximunStrenghValue < minimunStrenghValue){
    validationBox.setText("The Maximun Strengh must be greater than Minimun Strengh");
    validationBox.setIcon(QMessageBox::Critical);
    validationBox.exec();

  }

  QString scale = ui->scaleInputText->toPlainText();
  if(scale.isEmpty()){
    validationBox.setText("Please write a value for the Scale");
    validationBox.setIcon(QMessageBox::Critical);
    validationBox.exec();
  }
  ok = true;
  int scaleValue = scale.toInt(&ok);
  qDebug() << scaleValue;
  if(!ok){
    validationBox.setText("Please write a numeric vaule for the scale");
    validationBox.setIcon(QMessageBox::Critical);
    validationBox.exec();
  }
  if(ok && scaleValue <= 0 || scaleValue > 100){
    validationBox.setText("Remember the scale must be a value between 1 and 100");
    validationBox.setIcon(QMessageBox::Warning);
    validationBox.exec();
  }
  bool noClustering = ui->noClusterRadioButton->isChecked();
  bool existingCluster = ui->existingClusterRadioButton->isChecked();
  bool positionBasedWeights = ui->PositionBasesWeigthsRadioButton->isChecked();
  bool clusterPercentage = ui->clusterPercentageRadioButton->isChecked();
  int Cluster = -1;
  double PBParameterValue = 0.0;
  if(noClustering){
    Cluster = -2;
  }else if(existingCluster){
    Cluster = -3;
  }else if(positionBasedWeights){
    Cluster = -4;
    QString PBParameter = ui->positionBasedWeightsInputText->toPlainText();
    if(PBParameter.isEmpty()){
      validationBox.setText("Please write a value for the PBParameter (n)");
      validationBox.setIcon(QMessageBox::Critical);
      validationBox.exec();
    }
    ok = true;
    PBParameterValue = PBParameter.toDouble(&ok);
    qDebug() << PBParameterValue;
    if(!ok){
      validationBox.setText("Please write a numeric vaule for the PBParameter");
      validationBox.setIcon(QMessageBox::Critical);
      validationBox.exec();
    }
    if(ok && PBParameterValue <= 0.0 || PBParameterValue > 9.0){
      validationBox.setText("Remember the scale must be a value between 1 and 9");
      validationBox.setIcon(QMessageBox::Warning);
      validationBox.exec();
    }
  }else if(clusterPercentage){
    Cluster = -1;
    QString ClusterText = ui->clusterPercentageInputText->toPlainText();
    if(ClusterText.isEmpty()){
      validationBox.setText("Please write a value for the Cluster");
      validationBox.setIcon(QMessageBox::Critical);
      validationBox.exec();
    }
    ok = true;
    Cluster = ClusterText.toDouble(&ok);
    qDebug() << Cluster;
    if(!ok){
      validationBox.setText("Please write a numeric value for the Cluster");
      validationBox.setIcon(QMessageBox::Critical);
      validationBox.exec();
    }
    if(ok && Cluster <= 0 || Cluster > 100){
      validationBox.setText("Remember the Cluster must be a value between 1 and 100");
      validationBox.setIcon(QMessageBox::Warning);
      validationBox.exec();
    }
  }


}

void MainWindow::on_noClusterRadioButton_clicked()
{
    ui->clusterPercentageInputText->setDisabled(true);
    ui->positionBasedWeightsInputText->setDisabled(true);
    ui->clusterPercentageInputText->setText("");
}

void MainWindow::on_existingClusterRadioButton_clicked()
{
     ui->clusterPercentageInputText->setDisabled(true);
     ui->positionBasedWeightsInputText->setDisabled(true);
     ui->clusterPercentageInputText->setText("");
     ui->positionBasedWeightsInputText->setText("");
}

void MainWindow::on_existingWeightsRadioButton_clicked()
{
     ui->clusterPercentageInputText->setDisabled(true);
     ui->positionBasedWeightsInputText->setDisabled(true);
     ui->clusterPercentageInputText->setText("");
     ui->positionBasedWeightsInputText->setText("");
}

void MainWindow::on_PositionBasesWeigthsRadioButton_clicked()
{
     ui->clusterPercentageInputText->setDisabled(true);
     ui->positionBasedWeightsInputText->setDisabled(false);
     ui->clusterPercentageInputText->setText("");

}

void MainWindow::on_clusterPercentageRadioButton_clicked()
{
     ui->clusterPercentageInputText->setDisabled(false);
     ui->positionBasedWeightsInputText->setDisabled(true);
     ui->positionBasedWeightsInputText->setText("");
}

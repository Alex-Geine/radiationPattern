#include "mainwindow.h"

#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->widget->set_array_range(11,11);

    ui->label->setStyleSheet("QLabel { background-color : rgb(192,192,192) ;}");
}

MainWindow::~MainWindow()
{
    delete ui;
}

QVector<double> linspace(double left, double right, size_t size)
{
    QVector<double> x(size);
    for (int i = 0; i < size; ++i) {
        x[i] = left + i*(right-left)/ double(size);
    }
    return x;
}

void MainWindow::on_pushButton_2_clicked()
{
    double l = 1;
    double R = ui->DoubleSpinBox_2->value();
    double k = ui->DoubleSpinBox->value();
    const int N = 400;

    double* ampl = g_getDiagram(ui->widget->array, l, k, R, N);
    for (int i = 0; i < N*N; ++i)
         ampl[i] = std::log10(ampl[i] + 1.);

    auto image = create_QImage(ampl,N,N);
    ui->label->setPixmap(QPixmap::fromImage(image));
}


void MainWindow::on_pushButton_clicked()
{
    ui->pushButton_2->setEnabled(true);
}


void MainWindow::on_pushButton_4_clicked()
{
    ui->widget->array_set_empty();
    ui->widget->repaint();
}


void MainWindow::on_pushButton_3_clicked()
{
    ui->widget->array_set_all();
    ui->widget->repaint();
}


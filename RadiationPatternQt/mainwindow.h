#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <RadiationPattern.h>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    bool isLog = false;

private slots:
    void on_pushButton_2_clicked();

    void on_pushButton_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_3_clicked();


    QImage create_QImage(double* data, unsigned w, unsigned h)
    {
/*
        for (auto i = 0; i < w * h; ++i)
            std::cout << data[i] << " ";
        std::cout << std::endl;

        std::cout << "After log: " << std::endl;
        for (auto i = 0; i < w * h; ++i)
            std::cout << data[i] << " ";
        std::cout << std::endl;
*/
        uint8_t* image = new uint8_t[w*h];
        double max = data[0];
        for (auto i = 0; i < w * h; ++i)
             if (max < data[i])
                 max = data[i];

        //std::cout << "max: " << max << std::endl;

        //std::cout << "Image: " << std::endl;
        for (auto i = 0; i < w*h; ++i)
        {
            image[i] = data[i] / max * 255.;
        //    std::cout << (uint32_t)image[i] << " ";
        }
        //std::cout << std::endl;
        QImage imageOut(w,h,QImage::Format_Grayscale8);
        for (auto i = 0; i < w; ++i)
            for (auto j =0 ;j < h; ++j)
            {
                uint8_t pixel = image[i * h + j];
                auto c = qRgb(pixel,pixel,pixel);
                imageOut.setPixel(i,j,c);
            }

        return imageOut;
    }

};
#endif // MAINWINDOW_H

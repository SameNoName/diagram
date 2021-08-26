#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <cmath>
#include <iostream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

double theory(double x, double m)
{
    if (abs(x-1) < 0.001){
        double d = (m+1)/2;
        return d;
    }
    else return ((m+1)*pow(x,m+3) - (m+2)*pow(x,m+2) + x)/((1-pow(x,m+2))*(1-x));
}

double inverseGen(double lambda)
{
    double r = (rand()%1000);
    r /= 1000;
    double rslt = -1/lambda*log(r);
    if (rslt > 8) rslt = 8;
    return rslt;
}

double theoryErlang(double t, double lambda, double r)
{
    const double e = exp(1);
    double answer = lambda * pow(lambda*t, r-1) * pow(e, -lambda*t);
    for (int i = 2; i <= r-1; ++i)
    {
        answer /= i;
    }
    return answer;
}

double generatorErlang(double lambda, double r)
{
    double sum = 0;
    for (int j = 0; j < r; ++j)
    {
        double rnd = (rand()%1000);
        rnd /= 1000;
        sum += log(rnd);
    }
    double rslt = -1/(lambda*r)*sum;
    if (rslt > 64) rslt = 64;
    return rslt;
}

void MainWindow::on_pushButton_clicked()
{
    srand(time(NULL));
    QVector<double> x, y, y1, y2;
    int a = ui->aValue->value();
    int b = ui->bValue->value();
    int m = ui->mValue->value();
    double step = ui->stepValue->value();
    int l = ui->lValue->value();
    double r = ui->rValue->value();

    y1.push_back(0);
    y2.push_back(0);
    for (double lambda = step; lambda <= b*1.01; lambda += step)
    {
        QVector<double> timeVector, nCurrent;
        int counter = 0;
        double currentTime = 0;
        double nPlusTime;//время, через которое количество заявок в системе увеличится на одну
        double nMinusTime;//время, через которое количество заявок в системе уменьшится на одну
        int n = 0;
        QVector<double> nVector;
        while (currentTime < l*5)
        {
            timeVector.push_back(currentTime);
            nCurrent.push_back(n);
            ++counter;
            if (n == 0)
            {
                nPlusTime = 0;
                nMinusTime = 0;
                currentTime += inverseGen(lambda);
                //currentTime += generatorErlang(lambda, r);
                n=1;
            }
            else
            {
                if (abs(nPlusTime) < 0.001)
                {
                    nPlusTime = inverseGen(lambda);
                    //nPlusTime += generatorErlang(lambda, r);
                }
                if (abs(nMinusTime) < 0.001)
                {
                    nMinusTime = inverseGen(1);
                    //nMinusTime += generatorErlang(1, r);
                }
                if (nPlusTime > nMinusTime)
                {
                    currentTime += nMinusTime;
                    nPlusTime -= nMinusTime;
                    nMinusTime = 0;
                    --n;
                    if (n < 0) n = 0;
                }
                else if (nPlusTime < nMinusTime)
                {
                    currentTime += nPlusTime;
                    nMinusTime -= nPlusTime;
                    nPlusTime = 0;
                    ++n;
                    if (n > m + 1) n = m + 1;
                }
                else if (nPlusTime == nMinusTime)
                {
                    currentTime += nPlusTime;
                    nPlusTime = 0;
                    nMinusTime = 0;
                }
            }
        }

        int five = 0;
        for (int i = 0; i < timeVector.size(); ++i)
        {
            if (timeVector[i] > five)
            {
                if (i == 0)
                {
                    nVector.push_back(nCurrent[0]);
                }
                else
                {
                    for (int j = 1; j <= i; ++j)
                    {
                        if (i - j == 0)
                        {
                            nVector.push_back(nCurrent[0]);
                            break;
                        }
                        if (timeVector[i-j] <= five+0.001)
                        {
                            nVector.push_back(nCurrent[i - j]);
                            break;
                        }
                    }
                }
                five += 5;
            }
        }

        double avgN = 0;
        for (int i = 0; i < nVector.size(); ++i)
        {
            avgN += nVector[i];
        }
        avgN /= nVector.size();
        y1.push_back(avgN);

        double disp = 0;
        for (int i = 0; i < nVector.size(); ++i)
        {
            disp += (nVector[i] - avgN) * (nVector[i] - avgN);
        }
        disp /= (nVector.size() + 1);

        double s = sqrt(disp);
        y2.push_back(s);
    }

    for (double X=a; X<=b*1.01; X+=step)
    {
        x.push_back(X);
        y.push_back(theory(X,m));
        //y.push_back(theoryErlang(X, m, r));
    }

    double maxY = m+1;
    for (int i = 0; i<y1.size(); ++i)
    {
       if (y1[i] > maxY) maxY = y1[i];
       if (y2[i] > maxY) maxY = y2[i];
    }

    ui->widget->clearGraphs();
    ui->widget->legend->setVisible(true);
    ui->widget->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);

    ui->widget->addGraph();
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->graph(0)->setData(x, y);
    ui->widget->graph(0)->setName("Теоретический график");

    ui->widget->addGraph();
    ui->widget->graph(1)->setPen(QPen(Qt::blue));
    ui->widget->graph(1)->setLineStyle(QCPGraph::lsNone);
    ui->widget->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->widget->graph(1)->setData(x, y1);
    ui->widget->graph(1)->setName("График, построенный методом МК");

    ui->widget->addGraph();
    ui->widget->graph(2)->setPen(QPen(Qt::green));
    ui->widget->graph(2)->setData(x, y2);
    ui->widget->graph(2)->setName("Стандартное отклонение");

    ui->widget->xAxis->setLabel("Lambda");
    ui->widget->yAxis->setLabel("N");

    ui->widget->xAxis->setRange(a-((b-a)/10),b+((b-a)/10));
    ui->widget->yAxis->setRange(0, 1.1* maxY);
    ui->widget->replot();
}

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <fmt/base.h>
#include <fmt/format.h>
#include <boost/numeric/odeint.hpp>
#include <QDebug>

using namespace boost::numeric::odeint;

using state_type = std::vector<double>;

// Definição do sistema de equações diferenciais
struct harmonic_oscillator {
    void operator()(const state_type& x, state_type& dxdt, double t) {
        dxdt[0] = x[1];              // dx/dt = v
        dxdt[1] = -x[0];             // dv/dt = -x
    }
};

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

void MainWindow::on_pushButton_clicked()
{
    static int i = 0;
    auto s = fmt::format("Push {}", i++);
    ui->pushButton->setText(s.c_str());

    // Definindo o estado inicial
    state_type x(2);
    x[0] = 1.0; // Posição inicial
    x[1] = 0.0; // Velocidade inicial

    // Definindo o stepper
    typedef controlled_runge_kutta<runge_kutta_dopri5<state_type>> stepper_type;
    stepper_type stepper;

    double t = 0.0;    // Tempo inicial
    double t_end = 10.0; // Tempo final
    double dt = 1;   // Tamanho inicial do passo

    while (t < t_end) {
        //Tries one step of step size dt.
        // If the step was successful, success is returned, the resulting state is written to x, the new time is stored in t and dt now contains a new (possibly larger) step-size for the next step.
        // If the error was too big, rejected is returned and the results are neglected - x and t are unchanged and dt now contains a reduced step-size to be used for the next try.

        auto res = stepper.try_step(harmonic_oscillator(), x, t, dt);

        if (res == success) {
            qDebug() << "t: " << t << ", x: " << x[0] << ", v: " << x[1] << "passo: " << dt;
        }
        else {

            qDebug() << "fail - new dt: " << dt;
        }
    }
}


#include <iostream>
#include <fstream>

#include <ddp-actuator-solver/ddpsolver.hh>
#include "dctemp.hh"

#include <time.h>
#include <sys/time.h>


using namespace std;
using namespace Eigen;

#define STATE_NB 5
#define COMMAND_NB 1

int main()
{
    DDPSolver<double,STATE_NB,COMMAND_NB>::stateVec_t xinit,xDes,x;
    DDPSolver<double,STATE_NB,COMMAND_NB>::commandVec_t u;

    xinit << 0.0,0.0,25.0,0.0,25.0;
    xDes << 0.5,0.0,0.0,0.0,0.0;

    double t_end = 3*60;
    double dt=1e-3;
    unsigned int T = (unsigned int)(t_end/dt);
    double err=0.0;
    double err_int = 0.0;
    double Kp = 1.0e-2;
    double Ki = 1.0*1e-4;

    DCTemp model(dt);


    ofstream fichier1("results1.csv",ios::out | ios::trunc);
    if(fichier1)
    {
        fichier1 << "q,qdot,T,tau_ext,T_a,u" << endl;
        fichier1 << T << "," << STATE_NB << "," << COMMAND_NB << endl;
        err = 300.0 - x(1, 0);
        err_int += err;
        u << Kp*err + Ki*err_int;
        if(u(0,0)>5.0)  u << 5.0;
        if(u(0,0)<-5.0) u << -5.0;
        x = xinit;
        fichier1 << x(0, 0) << "," << x(1, 0) << "," << x(2, 0) << ","
                 << x(3, 0) << "," << x(4, 0) << "," << u(0, 0) << endl;
        for (unsigned int i = 1; i < T; i++)
        {
            //x(3,0) = 0.1*sin(0.5*M_PI*(x(0,0)/180.0 ));
            err = 300.0 - x(1, 0);
            err_int += err;
            u << Kp*err + Ki*err_int;
            if(u(0,0)>5.0)  u << 5.0;
            if(u(0,0)<-5.0) u << -5.0;
            x = model.computeNextState(dt, x, u);
            fichier1 << x(0, 0) << "," << x(1, 0) << "," << x(2, 0) << ","
                     << x(3, 0) << "," << x(4, 0) << "," << u(0, 0) << endl;
        }
        fichier1.close();
    }
    else
        cerr << "erreur ouverte fichier" << endl;



    return 0;

}

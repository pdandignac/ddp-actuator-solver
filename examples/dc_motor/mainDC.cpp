#include <iostream>
#include <fstream>

#include <ddp-actuator-solver/ddpsolver.hh>
#include "dcmotor.hh"
#include "cost.hh"

#include <time.h>
#include <sys/time.h>


using namespace std;
using namespace Eigen;

#define STATE_NB 3
#define COMMAND_NB 1

int main()
{
    struct timeval tbegin,tend;
    double texec=0.0;
    DDPSolver<double,STATE_NB,COMMAND_NB>::stateVec_t xinit,xDes,x;
    DDPSolver<double,STATE_NB,COMMAND_NB>::commandVec_t u;

    xinit << 0.0,0.0,0.0;
    xDes << 0.0,150.0,0.0;

    double t_end = 0.05;//3*60;
    double dt=1e-5;
    unsigned int T = (unsigned int)(t_end/dt);
    unsigned int iterMax = 100;
    double stopCrit = 10;
    DDPSolver<double,STATE_NB,COMMAND_NB>::stateVecTab_t xList;
    DDPSolver<double,STATE_NB,COMMAND_NB>::commandVecTab_t uList;
    DDPSolver<double,STATE_NB,COMMAND_NB>::traj lastTraj;

    DCMotor model(dt);
    DCMotor* noisyModel=NULL;
    CostDC cost;
    DDPSolver<double,STATE_NB,COMMAND_NB> solver(model,cost,DISABLE_FULLDDP,ENABLE_QPBOX);
    solver.FirstInitSolver(xinit,xDes,T,dt,iterMax,stopCrit);

    int N = 100;
    gettimeofday(&tbegin,NULL);
    solver.solveTrajectory();
    gettimeofday(&tend,NULL);

    lastTraj = solver.getLastSolvedTrajectory();
    xList = lastTraj.xList;
    uList = lastTraj.uList;
    unsigned int iter = lastTraj.iter;


    texec=((double)(1000*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000)))/1000.;
    texec /= N;

    cout << endl;
    cout << "temps d'execution total du solveur ";
    cout << texec << endl;
    cout << "temps d'execution par pas de temps ";
    cout << texec/T << endl;
    cout << "Nombre d'itérations : " << iter << endl;


    ofstream fichier1("results1.csv",ios::out | ios::trunc);
    if(fichier1)
    {
        fichier1 << "q,qdot,tau,i" << endl;
        fichier1 << T << "," << STATE_NB << "," << COMMAND_NB << endl;
        u << uList[0];
        x = xinit;
        fichier1 << x(0, 0) << "," << x(1, 0) << "," << x(2, 0) << ","
                 << u(0, 0) << endl;
        for (unsigned int i = 1; i < T; i++)
        {
            u << uList[i];
            x = model.computeNextState(dt, x, u);
            fichier1 << x(0, 0) << "," << x(1, 0) << "," << x(2, 0) << ","
                     << u(0, 0) << endl;
        }
        fichier1.close();
    }
    else
        cerr << "erreur ouverte fichier" << endl;



    return 0;


}

#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>

#include <ddp-actuator-solver/ddpsolver.hh>
#include "romeosimpleactuator.hh"
#include "costfunctionromeoactuator.hh"


using namespace std;
using namespace Eigen;

#define STATE_NB 4
#define COMMAND_NB 1

int main()
{
    struct timeval tbegin,tend;
    double texec=0.0;
    DDPSolver<double,STATE_NB,COMMAND_NB>::stateVec_t xinit,xDes,x;
    DDPSolver<double,STATE_NB,COMMAND_NB>::commandVec_t u;

    xinit << 0.0,0.0,0.0,0.0;
    xDes << 1.0,0.0,100.0,0.0;

    unsigned int T = 3000;
    double dt=1e-3;
    unsigned int iterMax = 25;
    double stopCrit = 1000;//1e-5;
    DDPSolver<double,STATE_NB,COMMAND_NB>::stateVecTab_t xList;
    DDPSolver<double,STATE_NB,COMMAND_NB>::commandVecTab_t uList;
    DDPSolver<double,STATE_NB,COMMAND_NB>::traj lastTraj;

    RomeoSimpleActuator romeoActuatorModel(dt);
    RomeoSimpleActuator* romeoNoisyModel=NULL;
    CostFunctionRomeoActuator costRomeoActuator;
    DDPSolver<double,STATE_NB,COMMAND_NB> testSolverRomeoActuator(romeoActuatorModel,costRomeoActuator,DISABLE_FULLDDP,DISABLE_QPBOX);
    testSolverRomeoActuator.FirstInitSolver(xinit,xDes,T,dt,iterMax,stopCrit);

    int N = 100;
    gettimeofday(&tbegin,NULL);
    testSolverRomeoActuator.solveTrajectory();
    gettimeofday(&tend,NULL);

    lastTraj = testSolverRomeoActuator.getLastSolvedTrajectory();
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
        fichier1 << "q,q_dot,theta,theta_dot,u" << endl;
        fichier1 << T << "," << STATE_NB << "," << COMMAND_NB << endl;
        x = xinit;
        fichier1 << x(0, 0) << "," << x(1, 0) << "," << x(2, 0) << ","
                 << x(3, 0) << "," << uList[0] << endl;
        for (unsigned int i = 1; i < T; i++)
        {
            x = romeoActuatorModel.computeNextState(dt, x, uList[i - 1]);
            fichier1 << x(0, 0) << "," << x(1, 0) << "," << x(2, 0) << ","
                     << x(3, 0) << "," << uList[i - 1] << endl;
        }
        fichier1 << xList[T](0, 0) << "," << xList[T](1, 0) << "," << xList[T](2, 0) << "," << xList[T](3, 0) << ","
                 << uList[T - 1](0, 0) << endl;
        fichier1.close();
    }
    else
        cerr << "erreur ouverte fichier" << endl;

    ofstream fichier2("results2.csv",ios::out | ios::trunc);
    if(fichier2)
    {
        fichier2 << "tau,tauDot,q,qDot,u" << endl;
        fichier2 << T << ',' << 0 << endl;
        for(int j=0;j<0;j++)
        {
            romeoNoisyModel = new RomeoSimpleActuator(dt,1);
            x = xinit;
            for (unsigned int i = 1; i < T; i++)
            {
                x = romeoNoisyModel->computeNextState(dt, x, uList[i - 1]);
                fichier2 << x(0, 0) << "," << x(1, 0) << "," << x(2, 0) << ","
                         << x(3, 0) << "," << uList[i - 1] << endl;
            }
            fichier2 << xList[T](0, 0) << "," << xList[T](1, 0) << "," << xList[T](2, 0) << "," << xList[T](3, 0) << ","
                     << uList[T - 1](0, 0) << endl;
            delete romeoNoisyModel;
        }
        fichier2.close();
    }
    else
        cerr << "erreur ouverte fichier" << endl;


    return 0;

}

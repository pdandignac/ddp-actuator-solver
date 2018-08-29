#include <math.h>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <sys/time.h>

#include <iostream>

#include "dctemp.hh"

/*
 * x0 -> actuator position
 * x1 -> actuator speed
 * x2 -> motor temperature
 * x3 -> external torque
 * x4 -> ambiant temperature
 */

DCTemp::DCTemp(double& mydt,bool noiseOnParameters)
{
    stateNb=5;
    commandNb=1;
    dt = mydt;

    if(!noiseOnParameters)
    {
        J = 1e-3*1e-1;
        K_M=2.5*60.3e-2;
        f_VL=1e-3*3e1;
        R_TA=0.3*1.0e1;
        R_TH=1.84;
        C_TH=1e-2*1.77e3;
    }
    else
    {
        J = 1e-1;
        K_M=60.3e-1;
        f_VL=3e1;
        R_TA=1.0e0;
        R_TH=1.84;
        C_TH=1.77e3;
    }

    Id.setIdentity();


    fu.setZero();
    fx.setZero();

    fxx[0].setZero();
    fxx[1].setZero();
    fxx[2].setZero();
    fxx[3].setZero();
    fxx[4].setZero();

    fxu[0].setZero();
    fxu[0].setZero();

    fuu[0].setZero();
    fux[0].setZero();
    fxu[0].setZero();

    QxxCont.setZero();
    QuuCont.setZero();
    QuxCont.setZero();

    lowerCommandBounds << -1.0;
    upperCommandBounds << 1.0;
}

DCTemp::stateVec_t DCTemp::computeDeriv(double& , const stateVec_t& X, const commandVec_t &U)
{
    dX[0] = X[1];
    dX[1] = (K_M/J)*U[0] - (f_VL/J)*X[1] - (1.0/J)*X[3];
    dX[2] = (R_TA/C_TH)*U[0]*U[0] - (1.0/(R_TH*C_TH))*(X[2]-X[4]);
    dX[3] = 0.0;
    dX[4] = 0.0;
    //std::cout << dX.transpose() << std::endl;
    return dX;
}

DCTemp::stateVec_t DCTemp::computeNextState(double& dt, const stateVec_t& X,const commandVec_t& U)
{
    k1 = computeDeriv(dt,X,U);
    k2 = computeDeriv(dt,X+(dt/2)*k1,U);
    k3 = computeDeriv(dt,X+(dt/2)*k2,U);
    k4 = computeDeriv(dt,X+dt*k3,U);
    x_next = X + (dt/6)*(k1+2*k2+2*k3+k4);
    return x_next;
}

void DCTemp::computeAllModelDeriv(double& dt, const stateVec_t& X,const commandVec_t& U)
{
    double dh = 1e-7;
    stateVec_t Xp,Xm;
    commandVec_t Up,Um;
    Xp = X;
    Xm = X;
    Up = U;
    Um = U;
    for(unsigned int i=0;i<stateNb;i++)
    {
        Xp[i] += dh/2;
        Xm[i] -= dh/2;
        fx.col(i) = (computeNextState(dt, Xp, U) - computeNextState(dt, Xm, U))/dh;
        Xp = X;
        Xm = X;
    }
    for(unsigned int i=0;i<commandNb;i++)
    {
        Up[i] += dh/2;
        Um[i] -= dh/2;
        fu.col(i) = (computeNextState(dt, X, Up) - computeNextState(dt, X, Um))/dh;
        Up = U;
        Um = U;
    }
}

DCTemp::stateMat_t DCTemp::computeTensorContxx(const stateVec_t& )
{
    return QxxCont;
}

DCTemp::commandMat_t DCTemp::computeTensorContuu(const stateVec_t& )
{
    return QuuCont;
}

DCTemp::commandR_stateC_t DCTemp::computeTensorContux(const stateVec_t& )
{
    return QuxCont;
}

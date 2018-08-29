#include <math.h>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <sys/time.h>

#include <iostream>

#include "dcmotor.hh"

/*
 * x0 -> actuator position
 * x1 -> actuator speed
 * x2 -> external torque
 */

DCMotor::DCMotor(double& mydt,bool noiseOnParameters)
{
    stateNb=3;
    commandNb=1;
    dt = mydt;

    if(!noiseOnParameters)
    {
        J = 1e-5;
        K_M = 60.3e-3;
        f = 3e-3;
    }
    else
    {
        J = 1e-5;
        K_M = 60.3e-3;
        f = 3e-3;
    }

    Id.setIdentity();
    A <<    0.0,1.0,0.0,
            0.0,-f/J,-1.0/J,
            0.0,0.0,0.0;
    B << 0.0,(K_M/J),0.0;

    fx.setZero();
    fx = (dt*A+Id);

    fu.setZero();
    fu = dt*B;

    fxx[0].setZero();
    fxx[1].setZero();
    fxx[2].setZero();

    fxu[0].setZero();
    fxu[0].setZero();

    fuu[0].setZero();
    fux[0].setZero();
    fxu[0].setZero();

    QxxCont.setZero();
    QuuCont.setZero();
    QuxCont.setZero();

    lowerCommandBounds << -10.0;
    upperCommandBounds << 10.0;
}

DCMotor::stateVec_t DCMotor::computeDeriv(double& , const stateVec_t& X, const commandVec_t &U)
{
    dX[0] = X[1];
    dX[1] = (K_M/J)*U[0] - (1.0/J)*X[2];
    dX[2] = 0.0;
    //std::cout << dX.transpose() << std::endl;
    return dX;
}

DCMotor::stateVec_t DCMotor::computeNextState(double& dt, const stateVec_t& X,const commandVec_t& U)
{
    /*k1 = computeDeriv(dt,X,U);
    k2 = computeDeriv(dt,X+(dt/2.0)*k1,U);
    k3 = computeDeriv(dt,X+(dt/2.0)*k2,U);
    k4 = computeDeriv(dt,X+dt*k3,U);
    x_next = X + (dt/6.0)*(k1+2.0*k2+2.0*k3+k4);*/
    x_next = fx*X + fu*U;
    return x_next;
}

void DCMotor::computeAllModelDeriv(double& dt, const stateVec_t& X,const commandVec_t& U)
{
    double dh = 1e-10;
    stateVec_t Xp,Xm;
    commandVec_t Up,Um;
    Xp = X;
    Xm = X;
    Up = U;
    Um = U;
    /*for(unsigned int i=0;i<stateNb;i++)
    {
        Xp[i] += dh/2.0;
        Xm[i] -= dh/2.0;
        fx.col(i) = (computeNextState(dt, Xp, U) - computeNextState(dt, Xm, U))/dh;
        Xp = X;
        Xm = X;
    }*/
        /*std::cout << fx << std::endl;
        std::cout << "-" << std::endl;
        std::cout << fu << std::endl; 
        std::cout << std::endl;*/
}

DCMotor::stateMat_t DCMotor::computeTensorContxx(const stateVec_t& )
{
    return QxxCont;
}

DCMotor::commandMat_t DCMotor::computeTensorContuu(const stateVec_t& )
{
    return QuuCont;
}

DCMotor::commandR_stateC_t DCMotor::computeTensorContux(const stateVec_t& )
{
    return QuxCont;
}

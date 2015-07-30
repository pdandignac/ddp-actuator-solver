#include "romeosimpleactuator.h"
#include <math.h>

#define pi M_PI

RomeoSimpleActuator::RomeoSimpleActuator()
{
    Id.setIdentity();

    A.setZero();
    A(0,1) = 1.0;
    A(2,3) = 1.0;
    A(1,0) = -((k/Jl)+(k/(Jm*R*R)));
    A(1,1) = -(fvm/Jm);
    A(1,3) = -((fvm*k)/Jm);
    A(3,0) = 1.0/Jl;

    A13atan = (2.0*Jm*R/(pi*Jl))*Cf0;
    A33atan = (2.0/(pi*Jl))*Cf0;

    B <<  0.0,
          k/(R*Jm),
          0.0,
          0.0;

    fxBase << 1.0,      1.0,      0.0,      0.0,
              -(k/Jl)-(k/(Jm*R*R)),      -(fvm/Jm),      0.0,      -(fvm*k)/Jm,
              0.0,      0.0,      1.0,      1.0,
              1.0/Jl,      0.0,      0.0,      0.0;
    fx.setZero();

    fxx[0].setZero();
    fxx[1].setZero();
    fxx[2].setZero();
    fxx[3].setZero();

    fxu[0].setZero();
    fxu[0].setZero();
    fuBase << 0.0,
              k/(R*Jm),
              0.0,
              0.0;
    fu.setZero();
    fuu[0].setZero();
    fux[0].setZero();
    fxu[0].setZero();

    QxxCont.setZero();
    QuuCont.setZero();
    QuxCont.setZero();
}


stateVec_t RomeoSimpleActuator::computeNextState(double& dt, const stateVec_t& X,const commandVec_t& U)
{
    stateMat_t Ad = (A*dt+Id);
    stateR_commandC_t Bd = dt*B;
    stateVec_t result = Ad*X + Bd*U;
    result(1,0)+=dt*A13atan*atan(a*X(3,0));
    result(3,0)+=dt*A33atan*atan(a*X(3,0));

    return result;
}

void RomeoSimpleActuator::computeAllModelDeriv(double& dt, const stateVec_t& X,const commandVec_t& U)
{
    fx = fxBase;

    fx(0,1) *= dt;
    fx(1,0) *= dt;
    fx(1,1) *= dt;
    fx(1,1) += 1.0;
    fx(1,3) *= dt;
    fx(1,3) += ((2*dt*Jm*R)/(pi*Jl))*Cf0*(a/(1.0+(a*a*X(3,0)*X(3,0))));
    fx(2,3) *= dt;
    fx(3,0) *= dt;
    fx(3,3) = 1.0-((2*dt*Cf0)/(pi*Jl))*(a/(1.0+(a*a*X(3,0)*X(3,0))));

    fu = dt*fuBase;

    fxx[3](1,3) = -((2*dt*Jm*R)/(pi*Jl))*Cf0*((2*a*a*a*X(3,0))/((1+(a*a*X(3,0)*X(3,0)))*(1+(a*a*X(3,0)*X(3,0)))));
    fxx[3](3,3) = +((2*dt*Cf0)/(pi*Jl))*((2*a*a*a*X(3,0))/((1+(a*a*X(3,0)*X(3,0)))*(1+(a*a*X(3,0)*X(3,0)))));
}

stateMat_t RomeoSimpleActuator::computeTensorContxx(const stateVec_t& nextVx)
{
    QxxCont = nextVx[3]*fxx[3];
    return QxxCont;
}

commandMat_t RomeoSimpleActuator::computeTensorContuu(const stateVec_t& nextVx)
{
    return QuuCont;
}

commandR_stateC_t RomeoSimpleActuator::computeTensorContux(const stateVec_t& nextVx)
{
    return QuxCont;
}

/// accessors ///
unsigned int RomeoSimpleActuator::getStateNb()
{
    return stateNb;
}

unsigned int RomeoSimpleActuator::getCommandNb()
{
    return commandNb;
}

stateMat_t& RomeoSimpleActuator::getfx()
{
    return fx;
}

stateTens_t& RomeoSimpleActuator::getfxx()
{
    return fxx;
}

stateR_commandC_t& RomeoSimpleActuator::getfu()
{
    return fu;
}

stateR_commandC_commandD_t& RomeoSimpleActuator::getfuu()
{
    return fuu;
}

stateR_stateC_commandD_t& RomeoSimpleActuator::getfxu()
{
    return fxu;
}

stateR_commandC_stateD_t& RomeoSimpleActuator::getfux()
{
    return fux;
}
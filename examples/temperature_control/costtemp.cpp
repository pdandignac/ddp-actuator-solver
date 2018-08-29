#include "costtemp.hh"
#include <iostream>

CostTemp::CostTemp()
{
    Q <<    0.0,    0.0,    0.0,    0.0,    0.0,
            0.0,    1.0e-3, 0.0,    0.0,    0.0,
            0.0,    0.0,    1.0e30, 0.0,    0.0,
            0.0,    0.0,    0.0,    0.0,    0.0,
            0.0,    0.0,    0.0,    0.0,    0.0;
    R << 1.0e-8;
    a = 2.0;
    Tmax = 100.0;

    lxx = Q;
    luu = R;
    lux << 0.0,0.0,0.0,0.0,0.0;
    lxu << 0.0,0.0,0.0,0.0,0.0;
    lx.setZero();
}

void CostTemp::computeAllCostDeriv(const stateVec_t& X,const stateVec_t& Xdes, const commandVec_t& U)
{
    lx = Q*(X-Xdes);
    //std::cout << lx(1,0) << "-" << lx(2,0) << std::endl;
    //std::cout << (X-Xdes).transpose() << std::endl;
    //lx(2,0) += exp((1/a)*(X[2]-Tmax));
    //lxx(2,2) = Q(2,2) + (1/a)*exp((1/a)*(X[2]-Tmax));
    lu = R*U;
}

void CostTemp::computeFinalCostDeriv(const stateVec_t& X,const stateVec_t& Xdes)
{
    lx = 1.0*Q*(X-Xdes);
    lxx = 1.0*Q;
}

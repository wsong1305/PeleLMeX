#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real P_mean    = 101325.0_rt;
    amrex::Real v_in      = 28.0_rt;
    amrex::Real T_in      = 900.0_rt;
    amrex::Real inflow_diam = 900.0_rt;

};
#endif

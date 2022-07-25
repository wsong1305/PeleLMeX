#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean   = 101325.0_rt;
   amrex::Real T_mean   = 300.0_rt;
   amrex::Real Vcoflow  = 0.0_rt;

   amrex::Real ps_halfangle = 30.0_rt;
   amrex::Real ps_r0        = 0.0005_rt;
   amrex::Real ps_rmsvel    = 5.0_rt;

   amrex::Real part_temp    = 0.0;
   amrex::Real mass_flow_rate = 0.0;

   amrex::Real jet_start_time = 0.0;
   amrex::Real jet_end_time = 1.0;

   amrex::Real floating_injection_mass = 0.0;
   amrex::Real part_dia_min = 1.0e-7;
};
#endif

#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{

  amrex::Real P_mean = 101325.0;
  amrex::Real splitx = 0.01;
  amrex::Real midtanh = 0.001;
  amrex::Real widthtanh = 0.001;
  amrex::Real Zst = 0.055;
  amrex::Real T_in = 300.0;
  amrex::Real U_b = 5;
  amrex::Real U_c = 0.5;
  amrex::Real x0 = 0.005;
  amrex::Real D = 0.0032;
  amrex::Real H = 0.0001;
  amrex::Real pertmag = 0.001;
  amrex::Real standoff = 0.01;
  amrex::Real amplification = 100;

  amrex::Real Do_swirler = 0.0;
  amrex::Real Di_swirler = 0.0;
  amrex::Real D_jet = 0.0;
  amrex::Real phi_jet = 0.0;
  amrex::Real swirl_angle = 0.0;

  amrex::Real U_jet = 0.0;

  int bathID = -1;
  int fuelID = -1;
  int oxidID = -1;
    
  bool ignition = false;
  amrex::Real T0 = 900.;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
  // Jet parameters
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_jet = {{0.0}};
  amrex::Real jet_start_time = 0.;
  amrex::Real jet_end_time = 10000.;
  amrex::Real jet_dia = 1.E-5;
  amrex::Real jet_vel = 50.;
  // This is the spread of the jet
  amrex::Real spread_angle = 20.;
  
  amrex::Real part_weibull_k = 3.;
  amrex::Real part_mean_dia = 1.E-6;
  amrex::Real part_stdev_dia = 0.;
  amrex::Real mass_flow_rate = 2.3E-3;
  amrex::Real part_temp = 300.;
  // This is the polar angle of the jet from the Z axis
  // 0 is perpendicular to the lower jet boundary
  amrex::Real jet_angle = 0.;
  // These are bounds on the azimuthal angle
  amrex::Real jet_azi_lo = 0.;
  amrex::Real jet_azi_hi = 360.;
  // Center of the jet
  amrex::RealVect jet_cent;
};
#endif

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
    amrex::Real Tb_init = 300.0;
    amrex::Real U_b = 5;
    amrex::Real U_c = 0.5;
    amrex::Real x0 = 0.005;
    amrex::Real D = 0.0032;
    amrex::Real H = 0.0001;
    amrex::Real pertmag = 0.001;
    amrex::Real standoff = 0.01;
    amrex::Real vel_fluct = 0;

    amrex::Real Do_swirler = 0.0;
    amrex::Real Di_swirler = 0.0;
    amrex::Real D_jet = 0.0;

    amrex::Real phi_in = 0.0;        // equivalence ratio at the start of the simulation
    amrex::Real phi_target = 0.0;    // equiv ratio desired after successful ignition
    amrex::Real t_start = 0.0;       // time to start changing the equiv ratio to target
    amrex::Real dt_phi_change = 0.0; // time interval to equiv ratio to target

    bool ignition = false;

    amrex::Real swirl_angle = 0.0;

    amrex::Real jet_angle = 0.0;
    amrex::Real U_jet = 0.0;
    amrex::Real jet_radius = 0.0;
    amrex::Real jet_loc_x = 0.0;
    amrex::Real jet_loc_y = 0.0;
    amrex::Real jet_loc_z = 0.0;

    int bathID = -1;
    int fuelID = -1;
    int oxidID = -1;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> Ys = {{0.0}};

    amrex::Real jet_vel = 50.;
    amrex::Real jet_dia = 1.E-5;
    amrex::Real part_mean_dia = 5.E-6;
    amrex::Real part_stdev_dia = 0.;
    amrex::Real mass_flow_rate = 2.3E-3;
    amrex::Real part_temp = 300.;
    amrex::Real jet_start_time = 0.;
    amrex::Real jet_end_time = 10000.;
    amrex::Real spray_angle = 0.42;
    amrex::Real jet_dx_mod = 50.;
    amrex::Real gas_jet_vel = 50.;
    amrex::Real gas_jet_dia = 1.E-5;
 //   amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_jet = {{0.0}};
 //   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> jet_cents;
};
#endif

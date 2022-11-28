#ifndef _EBUSERDEFINED_H_
#define _EBUSERDEFINED_H_

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(const Geometry& geom,
              const int required_coarsening_level,
              const int max_coarsening_level)
{
    // ParmParse your geometry parameters
    // --> fixed geometry of the challenge problem
    //Geometric parameters. Default values
    amrex::ParmParse pp("eb2");
    Real d_air = 0.005;
    Real d_fuel = 0.001;
    Real th_fuel = 0.0002;
    Real th_inlet = 0.002;
    pp.query("d_air", d_air);
    pp.query("d_fuel", d_fuel);
    pp.query("th_fuel", th_fuel);
    pp.query("th_inlet", th_inlet);

    amrex::EB2::CylinderIF air_pipe(0.5*d_air,th_inlet*2, 2, {AMREX_D_DECL(0.0, 0.0, 0.0)}, true);
    amrex::EB2::CylinderIF fuel_pipe_int(0.5*d_fuel-th_fuel,th_inlet*2, 2, {AMREX_D_DECL(0.0, 0.0, 0.0)}, true);
    amrex::EB2::CylinderIF fuel_pipe_out(0.5*d_fuel+th_fuel,th_inlet*2, 2, {AMREX_D_DECL(0.0, 0.0, 0.0)}, false);

    amrex::EB2::PlaneIF bottom_plane({AMREX_D_DECL(0.0, 0.0, th_inlet)},{AMREX_D_DECL(0.0,0.0,1.0)},false);

    auto fuel_pipe = amrex::EB2::makeIntersection(fuel_pipe_int,fuel_pipe_out);
    auto pipes = amrex::EB2::makeIntersection(air_pipe, bottom_plane);
    auto pipes_2 = amrex::EB2::makeUnion(pipes,fuel_pipe);

    auto gshop = amrex::EB2::makeShop(pipes_2);
    amrex::EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
}
#endif
#endif

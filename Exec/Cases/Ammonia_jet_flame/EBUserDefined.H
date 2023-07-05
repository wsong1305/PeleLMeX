#ifndef _EBUSERDEFINED_H_
#define _EBUSERDEFINED_H_

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(const Geometry& geom,
              const int required_coarsening_level,
              const int max_coarsening_level)
{
    // ParmParse your geometry parameters
    // --> fixed geometry of the challenge problem

    const amrex::Real *problo = geom.ProbLo();
    const amrex::Real *probhi = geom.ProbHi();
    amrex::Real dy = geom.CellSize()[1]; 
    amrex::Real dx = geom.CellSize()[0]; 

    amrex::Print() << "** Initialising flame holder geometry... \n";
    
    amrex::ParmParse pp("EB");
    
    amrex::Real in_diam = 0.0;
    pp.query("in_diam",in_diam);
    
    amrex::Real out_diam = 0.0;
    pp.query("out_diam",out_diam);

    amrex::Real plane_nozzle = 0.0;
    pp.query("nozzle_plane_height" ,plane_nozzle);
    
    amrex::Real plane_cavity = 0.0;
    pp.query("cavity_plane_height" ,plane_cavity);

    amrex::Real flame_holder_diam = 0.0;
    pp.query("flame_holder_diam",flame_holder_diam);

    EB2::PlaneIF out_vert_right({AMREX_D_DECL(out_diam/2.0,0.,0.)},
                                   {AMREX_D_DECL(-1.,0.,0.)});

    EB2::PlaneIF in_vert_right({AMREX_D_DECL(in_diam/2.0,0.,0.)},
                                   {AMREX_D_DECL(1.,0.,0.)});
    
    EB2::PlaneIF nozzle_plane_height({AMREX_D_DECL(0.,plane_nozzle,0.)},
                                  {AMREX_D_DECL(0.,-1.,0.)});
    
    EB2::PlaneIF cavity_plane_height({AMREX_D_DECL(0.,plane_cavity,0.)},
                                  {AMREX_D_DECL(0.,-1.,0.)});

    EB2::PlaneIF flame_holder_diam_plane({AMREX_D_DECL(flame_holder_diam/2.0,0.,0.)},
                                    {AMREX_D_DECL(1.,0.,0.)});

    
    auto outer_geom = EB2::makeIntersection(out_vert_right,nozzle_plane_height);
    auto nozzle     = EB2::makeIntersection(outer_geom,in_vert_right);
    auto geometry   = EB2::makeIntersection(nozzle,cavity_plane_height);

    // auto inner_geom2 = EB2::makeIntersection(inner_geom,in_vert_right);

    auto polys = EB2::makeUnion(geometry,flame_holder_diam_plane);

    auto pr     = EB2::lathe(polys);
    auto shop   = EB2::makeShop(pr);

    // EB2::Build(shop, geom, max_coarsening_level, max_coarsening_level);
    EB2::Build(shop, geom, required_coarsening_level,max_coarsening_level);

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBState(const amrex::Real xEBface[AMREX_SPACEDIM],
           amrex::Real s_ext[NVAR],
           const amrex::Real /*time*/,
           amrex::GeometryData const& geomdata,
           ProbParm const& prob_parm)
{
  // Provide the state on the EB centroid provided
  // by xEBface
  // Not the entire has to be specified. For instance, if using isothermal
  // boundary, only s_ext[TEMP] will be used.
  // Aborting here (message will not show on GPUs)
  // Abort("Default setEBState called. Copy and edit this file for your needs");
  
  amrex::Real r = sqrt(pow((0.0-xEBface[0]),2)+pow((0.0-xEBface[1]),2));

   s_ext[TEMP] = prob_parm.T_in;

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBType(const amrex::Real * /*xEBface[AMREX_SPACEDIM]*/,
          amrex::Real & EBflagType,
          amrex::GeometryData const& /*geomdata*/,
          ProbParm const& /*prob_parm*/)
{
  // Provide a factor applied on the EB face diffusion
  // coefficient: 0.0 will turn off the diffusion (homogeneous Neumann BC)
  // while 1.0 will let it unchanged. Intermediate values can be specified
  // for partially Dirichlet BCs
  EBflagType = 1.0;
  // Abort("Default setEBType called. Copy and edit this file for your needs");
}

#endif
#endif

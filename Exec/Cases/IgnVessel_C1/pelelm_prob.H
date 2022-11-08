#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>
#include <interpolator.H>
#include <PeleLM.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Phi(ProbParm const& prob_parm,
               const amrex::Real phi_in,
               amrex::Real Y[])
{
   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real Xt[NUM_SPECIES] = {0.0};
   amrex::Real a = 0.0;
   // if ( prob_parm.fuelID == NXC7H16_ID ) {
   //   a = 11.0;
   // }
   Xt[prob_parm.oxidID] = 1.0 / ( 1.0 + phi_in / a + 0.79 / 0.21 );
   Xt[prob_parm.fuelID] = phi_in * Xt[prob_parm.oxidID] / a;
   Xt[prob_parm.bathID] = 1.0 - Xt[prob_parm.oxidID] - Xt[prob_parm.fuelID];

   eos.X2Y(Xt,Y);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Ksi(ProbParm const& prob_parm,
               const amrex::Real ksi_in,
               amrex::Real Y[])
{
   amrex::Real phi = ksi_in / ( 1.0 - amrex::min(ksi_in,0.999999) ) * (1.0 - prob_parm.Zst ) / prob_parm.Zst;

   set_Y_from_Phi(prob_parm,phi,Y);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k,
                 int is_incompressible,
                 amrex::Array4<amrex::Real> const& state,
                 amrex::Array4<amrex::Real> const& aux,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(,
                 const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1]);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real pert = 0.0;
    amrex::Real y1, y2, z1, z2;
    amrex::Real massfrac_eq[NUM_SPECIES] = {0.0};

    amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real molefrac_air[NUM_SPECIES] = {0.0};


    massfrac_air[prob_parm.oxidID] = 0.233;
    massfrac_air[prob_parm.bathID] = 0.767;

    massfrac[prob_parm.oxidID] = 0.164;
    massfrac[prob_parm.bathID] = 0.72;
    massfrac[CO2_ID] = 0.094;
    massfrac[H2O_ID] = 0.022;


    state(i,j,k,TEMP) = prob_parm.T_in;
    //for (int n = 0; n < NUM_SPECIES; ++n) {
    //    molefrac[n] = molefrac_air[n];
    //}

    //eos.X2Y(molefrac, massfrac);

    state(i,j,k,VELX) = 0;
    state(i,j,k,VELY) = 0;
    state(i,j,k,VELZ) = 0;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,TEMP), massfrac, state(i,j,k,RHOH));
    state(i,j,k,RHOH) *= 1.0e-4 * state(i,j,k,DENSITY);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
        state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx      = geomdata.CellSize();
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> thermo_state = {0.0};

  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};


  auto eos = pele::physics::PhysicsType::eos();

  //radial distance from the injector center
  amrex::Real r = sqrt(pow((0.0-x[0]),2)+pow((0.0-x[1]),2)); 

  //compute tanh function
  amrex::Real eta = (1.0 - std::tanh(2.0*(r-prob_parm.D_jet*0.5)/(0.3*prob_parm.D_jet*0.5)))*0.5;
  if(eta < 0.0){
    eta = 0.0;
  }

  pele::physics::PMF::pmf(pmf_data, r, r, thermo_state);

  s_ext[VELX] = 0.0;
  s_ext[VELY] = 0.0; 
  s_ext[VELZ] = thermo_state[3]; //eta*prob_parm.U_jet;

  //Using tanh function to specify radial temperature distribution
  // s_ext[TEMP] = (1.0-eta)*prob_parm.T_in+prob_parm.T_jet*eta;
  s_ext[TEMP] = thermo_state[1];

  for (int n = 0; n < NUM_SPECIES; n++){
    massfrac[n] = thermo_state[n+4];
  }

  amrex::Real rho_cgs, P_cgs, RhoH_temp;
  P_cgs = prob_parm.P_mean * 10.0;


  eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
  s_ext[DENSITY] = rho_cgs * 1.0e3;

  eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
  s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int /*i*/, int /*j*/, int /*k*/,
           amrex::Array4<amrex::Real> const& /*beta*/,
           amrex::GeometryData const& /*geomdata*/,
           amrex::Box const& /*dBox*/,
           const int  /*dir*/,
           const int  /*state_comp*/,
           const int  /*nComp*/)
{
}
#endif

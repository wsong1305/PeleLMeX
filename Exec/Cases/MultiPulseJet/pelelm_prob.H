#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k,
                 int is_incompressible,
                 amrex::Array4<amrex::Real> const& state,
                 amrex::Array4<amrex::Real> const& aux,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    state(i,j,k,TEMP) = prob_parm.T_ox;

    for (int n = 0; n < NUM_SPECIES; n++)
    {
        massfrac[n] = prob_parm.Y_ox[n];
    }

    state(i,j,k,VELX) = 0.0;
    state(i,j,k,VELY) = 0.0;
    state(i,j,k,VELZ) = 0.0;

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i,j,k,TEMP), massfrac, h_cgs);
    state(i,j,k,RHOH) = h_cgs * 1.0e-4 * state(i,j,k,DENSITY);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
time_func(amrex::Real time, amrex::Real start, amrex::Real duration, amrex::Real tau, amrex::Real dwell)
{
  amrex::Real eta_time;
  amrex::Real eta_t1;
  amrex::Real eta_t2;
  amrex::Real inj1_end;
  amrex::Real inj2_beg;
  amrex::Real inj2_end;

  inj1_end = start + duration;
  inj2_beg = duration + dwell;
  inj2_end = inj2_beg + duration;

  eta_t1 = 0.5 * (-std::tanh((time - inj1_end)/(duration*tau)));
  eta_t2 = 0.5 * (std::tanh((time - inj2_beg)/(duration*tau)) - std::tanh((time-inj2_end)/(duration*tau)));

  eta_time = eta_t1 + eta_t2;

  return eta_time;
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
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real radius = {0.0};
  amrex::Real enthalpy;
  int ierr;

  // if (sgn == 1 && idir == 2) {
  if (sgn == 1) {

    // Jet function
    radius = std::sqrt(std::pow(x[0]-prob_parm.center_xy[0],2) + std::pow(x[1]-prob_parm.center_xy[1],2));
    amrex::Real eta = 0.5*(1.0 - std::tanh(2.0*(radius-prob_parm.D*0.5)/(0.1*prob_parm.D*0.5)));
    eta *= time_func(time, prob_parm.inj_start, prob_parm.inj_dur, prob_parm.tau, prob_parm.dwell);


    s_ext[VELX] = 0.0;
    s_ext[VELY] = 0.0;
    s_ext[VELZ] = eta * prob_parm.v_in + (1.0-eta) * prob_parm.v_co;

    s_ext[TEMP] = eta * prob_parm.T_fu + (1.0-eta) * prob_parm.T_ox;
    

    for (int n = 0; n < NUM_SPECIES; n++)
    {
      massfrac[n] = eta * prob_parm.Y_fuel[n] + (1.0-eta) * prob_parm.Y_ox[n];
      // massfrac[n] = prob_parm.Y_ox[n];
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
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif

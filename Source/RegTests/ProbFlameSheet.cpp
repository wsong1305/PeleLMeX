#include "ProbFlameSheet.H"
#include "AMReX_ParmParse.H"
#include "PeleLM_Index.H"
#include "PMF.H"

using namespace amrex;

void
FlameSheet::initParams(void **a_params, void **a_params_d, amrex::Real &a_pInit)
{
   m_FSparams = new FSParams{};

   Print() << " Initializing FlameSheet case \n";

   amrex::ParmParse pp("prob");
   pp.query("P_mean",   m_FSparams->P_mean);
   pp.query("standoff", m_FSparams->standoff);
   pp.query("pertmag",  m_FSparams->pertmag);
   pmf_data.initialize();

   a_pInit = m_FSparams->P_mean;

   m_FSparams_d = (FSParams*)The_Arena()->alloc(sizeof(FSParams));
   Gpu::copy(Gpu::hostToDevice, m_FSparams, m_FSparams+1, m_FSparams_d);
   
   *a_params = static_cast<void*>(m_FSparams);
   *a_params_d = static_cast<void*>(m_FSparams_d);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
FlameSheet::bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int /*m_nAux*/,
  amrex::Vector<amrex::Real> &s_ext,
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  void* a_params,
  pele::physics::PMF::PmfData::DataContainer const *a_pmf_data)
{
  FSParams* params = static_cast<FSParams*>(a_params);
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  auto eos = pele::physics::PhysicsType::eos();
  if (sgn == 1) {
    pele::physics::PMF::pmf(a_pmf_data,prob_lo[idir], prob_lo[idir], pmf_vals);

    s_ext[VELX] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
    s_ext[VELY] = pmf_vals[1]*1e-2;
#elif ( AMREX_SPACEDIM == 3 )
    s_ext[VELY] = 0.0;
    s_ext[VELZ] = pmf_vals[1]*1e-2;
#endif

    s_ext[TEMP] = pmf_vals[0];
    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    eos.X2Y(molefrac, massfrac);

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = params->P_mean * 10.0;

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
FlameSheet::initdata_k(int i, int j, int k,
                       amrex::Array4<amrex::Real> const& state,
                       amrex::GeometryData const& geomdata,
                       void* a_params,
                       pele::physics::PMF::PmfData::DataContainer const *a_pmf_data)
{
     
    FSParams* params = static_cast<FSParams*>(a_params);
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
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert = 0.0;

    if (params->pertmag > 0.0) 
    {

#if ( AMREX_SPACEDIM == 2 )
       pert = params->pertmag * 
             (1.0 * std::sin(2 * Pi * 4 * x / Lx) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / Lx) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / Lx) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / Lx) +
              0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / Lx));
    }

    amrex::Real y1 = (y - params->standoff - 0.5*dx[1] + pert)*100;
    amrex::Real y2 = (y - params->standoff + 0.5*dx[1] + pert)*100;

#elif ( AMREX_SPACEDIM == 3 )

       pert = params->pertmag *
             (1.0 * std::sin(2 * Pi * 4 * x / Lx) *
                std::sin(2 * Pi * 5 * y / Ly) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / Lx) *
                std::sin(2 * Pi * 4 * (y - 0.0053765) / Ly) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / Lx) *
                std::sin(2 * Pi * 3 * (y - 0.02137) / Ly) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / Lx) *
                std::sin(2 * Pi * 6 * (y - 0.018) / Ly) +
              0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / Lx));
    }

    amrex::Real y1 = (z - params->standoff - 0.5*dx[2] + pert)*100;
    amrex::Real y2 = (z - params->standoff + 0.5*dx[2] + pert)*100;
#endif


    pele::physics::PMF::pmf(a_pmf_data,y1, y2, pmf_vals);

    state(i,j,k,TEMP) = pmf_vals[0];;

    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    eos.X2Y(molefrac, massfrac);

    state(i,j,k,VELX) = 0.0;
#if ( AMREX_SPACEDIM == 2 )
    state(i,j,k,VELY) = pmf_vals[1]*1e-2;
#elif ( AMREX_SPACEDIM == 3 )
    state(i,j,k,VELY) = 0.0;
    state(i,j,k,VELZ) = pmf_vals[1]*1e-2;
#endif

    amrex::Real P_cgs = params->P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i,j,k,TEMP), massfrac, h_cgs);
    state(i,j,k,RHOH) = h_cgs * 1.0e-4 * state(i,j,k,DENSITY);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
       state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

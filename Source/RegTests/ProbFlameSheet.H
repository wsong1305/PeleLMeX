#ifndef PH_FLAMESHEET_H
#define PH_FLAMESHEET_H

#include "ProblemHelper.H"

class FlameSheet : public pele::physics::ProblemHelper::Register<FlameSheet>
{
public:
  static const std::string identifier() { return "FlameSheet"; }

  // Override virtuals
  void initParams(void **a_params, void **a_params_d, amrex::Real &a_pInit) override;

  void bcnormal(const amrex::Real* x,
                const int m_nAux,
                amrex::Vector<amrex::Real> &s_ext,
                const int idir,
                const int sgn,
                const amrex::Real time,
                amrex::GeometryData const& geomdata,
                void* a_params,
                pele::physics::PMF::PmfData::DataContainer const * pmf_data) override;

  void initdata_k(int i, int j, int k,
                  amrex::Array4<amrex::Real> const& state,
                  amrex::GeometryData const& geomdata,
                  void* a_params,
                  pele::physics::PMF::PmfData::DataContainer const * pmf_data) override;

  // FlameSheet params
  struct FSParams
  {
     amrex::Real P_mean   = 101325.0;
     amrex::Real standoff = 0.0;
     amrex::Real pertmag  = 0.0004;
  };

  FSParams *m_FSparams;
  FSParams *m_FSparams_d;
};
#endif

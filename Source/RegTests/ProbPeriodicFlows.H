#ifndef PH_PERIODICFLOWS_H
#define PH_PERIODICFLOWS_H

#include "ProblemHelper.H"

class PeriodicFlows : public pele::physics::ProblemHelper::Register<PeriodicFlows>
{
public:
  static const std::string identifier() { return "PeriodicFlows"; }

  // Override virtuals
  void initParams(void **a_params, void **a_params_d, amrex::Real &a_pInit) override;

  void initdata_k(int i, int j, int k,
                  amrex::Array4<amrex::Real> const& state,
                  amrex::GeometryData const& geomdata,
                  void* a_params,
                  pele::physics::PMF::PmfData::DataContainer const *a_pmf_data) override;

  // PeriodicFlows params
  struct PFParams
  {
     // Shared params
     int probType = 0;  
     amrex::Real T_mean = 298.0;
     amrex::Real P_mean = 101325.0;
     amrex::Real meanFlowMag = 0.0;
     int  meanFlowDir = 1;

     // CoVo params
     amrex::Real rvort = 0.07;
     amrex::Real xvort = 0.5;
     amrex::Real yvort = 0.5;
     amrex::Real forcevort = 6.0;

     // CoGau & DiffGau params
     amrex::Real rgauss = 0.1;
     amrex::Real xgauss = 0.5;
     amrex::Real ygauss = 0.5;
     amrex::Real ampgauss = 0.1;
     int gauss_type = 0;

     amrex::Real gaussTime = 0.001;
     amrex::Real gaussDiff = 0.000001;
  };

  PFParams *m_PFparams;
  PFParams *m_PFparams_d;
};
#endif

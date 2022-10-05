#include "ProbPeriodicFlows.H"
#include "AMReX_ParmParse.H"
#include "PeleLM_Index.H"

using namespace amrex;

void
PeriodicFlows::initParams(void **a_params, void **a_params_d, amrex::Real &a_pInit)
{
   m_PFparams = new PFParams{};

   Print() << " Initializing PeriodicFlows case \n";
   
   amrex::ParmParse pp("prob");
   std::string type;
   pp.query("type", type);
   pp.query("T_mean", m_PFparams->T_mean);
   pp.query("P_mean", m_PFparams->P_mean);

   if ( type == "ConvectedVortex" ) {
      m_PFparams->probType = 0;
      pp.query("rvort", m_PFparams->rvort);
      pp.query("xvort", m_PFparams->xvort);
      pp.query("yvort", m_PFparams->yvort);
      pp.query("forcevort", m_PFparams->forcevort);
      pp.query("meanFlowDir", m_PFparams->meanFlowDir);
      pp.query("meanFlowMag", m_PFparams->meanFlowMag);
   } else if ( type == "ConvectedGaussian" ) {
      m_PFparams->probType = 1;
      pp.query("gaussian_rad", m_PFparams->rgauss);
      pp.query("gaussian_x0", m_PFparams->xgauss);
      pp.query("gaussian_y0", m_PFparams->ygauss);
      pp.query("gaussian_ampl", m_PFparams->ampgauss);
      std::string gtype;
      pp.query("gaussian_type", gtype);
      if ( gtype == "Spec" ) {
         m_PFparams->gauss_type = 0;
      } else if ( gtype == "Temp" ) { 
         m_PFparams->gauss_type = 1;
      } else {
         amrex::Print() << " Unknown prob.gaussian_type ! Should be Spec or Temp \n";
         amrex::Abort();
      }
      pp.query("meanFlowDir", m_PFparams->meanFlowDir);
      pp.query("meanFlowMag", m_PFparams->meanFlowMag);
   } else if ( type == "ConvectedTanH" ) {
      m_PFparams->probType = 3;
      pp.query("tanh_rad", m_PFparams->rgauss);
      pp.query("tanh_x0", m_PFparams->xgauss);
      pp.query("tanh_y0", m_PFparams->ygauss);
      pp.query("tanh_ampl", m_PFparams->ampgauss);
      std::string gtype;
      pp.query("tanh_type", gtype);
      if ( gtype == "Spec" ) {
         m_PFparams->gauss_type = 0;
      } else if ( gtype == "Temp" ) { 
         m_PFparams->gauss_type = 1;
      } else {
         amrex::Print() << " Unknown prob.tanh_type ! Should be Spec or Temp \n";
         amrex::Abort();
      }
      pp.query("meanFlowDir", m_PFparams->meanFlowDir);
      pp.query("meanFlowMag", m_PFparams->meanFlowMag);
   } else if ( type == "DiffusedGaussian" ) {
      m_PFparams->probType = 2;
      pp.query("gaussian_time", m_PFparams->gaussTime);
      pp.query("gaussian_diffusivity", m_PFparams->gaussDiff);
      pp.query("gaussian_x0", m_PFparams->xgauss);
      pp.query("gaussian_y0", m_PFparams->ygauss);
      pp.query("gaussian_ampl", m_PFparams->ampgauss);
      std::string gtype;
      pp.query("gaussian_type", gtype);
      if ( gtype == "Spec" ) {
         m_PFparams->gauss_type = 0;
      } else if ( gtype == "Temp" ) {
         m_PFparams->gauss_type = 1;
      } else {
         amrex::Print() << " Unknown prob.gaussian_type ! Should be Spec or Temp \n";
         amrex::Abort();
      }
   } else { 
       amrex::Print() << " Unknown prob.type ! Should be ConvectedVortex, ConvectedGaussian or DiffusedGaussian \n";
       amrex::Abort();
   }
   
   a_pInit = m_PFparams->P_mean;

   m_PFparams_d = (PFParams*)The_Arena()->alloc(sizeof(PFParams));
   Gpu::copy(Gpu::hostToDevice, m_PFparams, m_PFparams+1, m_PFparams_d);
   
   *a_params = static_cast<void*>(m_PFparams);
   *a_params_d = static_cast<void*>(m_PFparams_d);
}

void
PeriodicFlows::initdata_k(int i, int j, int k,
                        amrex::Array4<amrex::Real> const& state,
                        amrex::GeometryData const& geomdata,
                        void* a_params,
                        pele::physics::PMF::PmfData::DataContainer const* /*a_pmf_data*/)
{
    PFParams* params = static_cast<PFParams*>(a_params);
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    AMREX_D_TERM(const amrex::Real x_c = prob_lo[0] + 0.5*Lx;,
                 const amrex::Real y_c = prob_lo[1] + 0.5*Ly;,
                 const amrex::Real z_c = prob_lo[2] + 0.5*Lz);

    auto eos = pele::physics::PhysicsType::eos();
    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    // Species IDs
#ifdef O2_ID
    int sp1_ID = O2_ID;
#endif
#ifdef N2a_ID
    int sp1_ID = N2a_ID;
#endif
    int sp2_ID = N2_ID;

    massfrac[sp1_ID] = 0.233;
    massfrac[sp2_ID] = 0.767;

    if ( params->probType == 0 ) { // CoVo
       const amrex::Real deltax = x - params->xvort;
       const amrex::Real deltay = y - params->yvort;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = params->rvort * params->rvort;
       const amrex::Real u_vort = -params->forcevort*deltay/r_sq * exp(-d_sq/r_sq/2.);
       const amrex::Real v_vort = params->forcevort*deltax/r_sq * exp(-d_sq/r_sq/2.);
       const amrex::Real w_vort = 0.;

       switch(params->meanFlowDir) {
         case 1  :
            AMREX_D_TERM(state(i,j,k,VELX) = params->meanFlowMag + u_vort;,
                         state(i,j,k,VELY) = v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -1  :
            AMREX_D_TERM(state(i,j,k,VELX) = -params->meanFlowMag + u_vort;,
                         state(i,j,k,VELY) = v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case 2  :
            AMREX_D_TERM(state(i,j,k,VELX) = u_vort;,
                         state(i,j,k,VELY) = params->meanFlowMag + v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -2  :
            AMREX_D_TERM(state(i,j,k,VELX) = u_vort;,
                         state(i,j,k,VELY) = -params->meanFlowMag + v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case 3  :
            AMREX_D_TERM(state(i,j,k,VELX) = params->meanFlowMag + u_vort;,
                         state(i,j,k,VELY) = params->meanFlowMag + v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -3  :
            AMREX_D_TERM(state(i,j,k,VELX) = -params->meanFlowMag + u_vort;,
                         state(i,j,k,VELY) = -params->meanFlowMag + v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
       }

       state(i,j,k,TEMP) = params->T_mean;

    } else if ( params->probType == 1 ) { // CoGau

       switch(params->meanFlowDir) {
         case 1  :
            AMREX_D_TERM(state(i,j,k,VELX) = params->meanFlowMag;,
                         state(i,j,k,VELY) = 0.0;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case -1  :                        
            AMREX_D_TERM(state(i,j,k,VELX) = -params->meanFlowMag;,
                         state(i,j,k,VELY) = 0.0;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case 2  :                         
            AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                         state(i,j,k,VELY) = params->meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case -2  :                        
            AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                         state(i,j,k,VELY) = -params->meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case 3  :                         
            AMREX_D_TERM(state(i,j,k,VELX) = params->meanFlowMag;,
                         state(i,j,k,VELY) = params->meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case -3  :                        
            AMREX_D_TERM(state(i,j,k,VELX) = -params->meanFlowMag;,
                         state(i,j,k,VELY) = -params->meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;
       }

       state(i,j,k,TEMP) = params->T_mean;
       const amrex::Real deltax = x - params->xgauss;
       const amrex::Real deltay = y - params->ygauss;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = params->rgauss * params->rgauss;
       if ( params->gauss_type == 0 ) { // Spec
          massfrac[sp1_ID] += massfrac[sp1_ID] * params->ampgauss * std::exp(-d_sq/r_sq);
          massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
       } else if ( params->gauss_type == 1 ) { // Temp
          state(i,j,k,TEMP) = params->T_mean * ( 1.0 + params->ampgauss * std::exp(-d_sq/r_sq) );
       } else {
          amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
       }
    } else if ( params->probType == 3 ) { // CoTanH

       switch(params->meanFlowDir) {
         case 1  :
            AMREX_D_TERM(state(i,j,k,VELX) = params->meanFlowMag;,
                         state(i,j,k,VELY) = 0.0;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case -1  :                        
            AMREX_D_TERM(state(i,j,k,VELX) = -params->meanFlowMag;,
                         state(i,j,k,VELY) = 0.0;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case 2  :                         
            AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                         state(i,j,k,VELY) = params->meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case -2  :                        
            AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                         state(i,j,k,VELY) = -params->meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case 3  :                         
            AMREX_D_TERM(state(i,j,k,VELX) = params->meanFlowMag;,
                         state(i,j,k,VELY) = params->meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;                         
         case -3  :                        
            AMREX_D_TERM(state(i,j,k,VELX)  = -params->meanFlowMag;,
                         state(i,j,k,VELY)  = -params->meanFlowMag;,
                         state(i,j,k,VELZ)  = 0.0);
            break;
       }

       state(i,j,k,TEMP) = params->T_mean;
       const amrex::Real deltax = x - params->xgauss;
       const amrex::Real deltay = y - params->ygauss;
       const amrex::Real rad = std::sqrt(deltax*deltax + deltay*deltay);
       const amrex::Real eta = 0.5*(1.0 - std::tanh(2.0*(rad-params->rgauss)/(0.05*params->rgauss)));
       if ( params->gauss_type == 0 ) { // Spec
          massfrac[sp1_ID] += massfrac[sp1_ID] * params->ampgauss * eta;
          massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
       } else if ( params->gauss_type == 1 ) { // Temp
          state(i,j,k,TEMP) = params->T_mean * ( 1.0 + params->ampgauss * eta );
       } else {
          amrex::Abort("Unknown tanh_type: should be either Temp or Spec");
       }
    } else if ( params->probType == 2 ) { // DiffGau

       AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                    state(i,j,k,VELY) = 0.0;,
                    state(i,j,k,VELZ) = 0.0);

       state(i,j,k,TEMP) = params->T_mean;
       // Set up the Gaussian as the solution of the diffusion of a delta dirac initial distribution
       // after time gaussTime of a diffusion process with constant diffusion coeff. gaussDiff. Use ampgauss
       // to avoid undershoot/overshoot of the 2-species mixture.
       const amrex::Real deltax = x - params->xgauss;
       const amrex::Real deltay = y - params->ygauss;
       amrex::Real denom_inv = params->ampgauss / std::sqrt(4*Pi*params->gaussTime*params->gaussDiff);     
       const amrex::Real d_sq = deltax*deltax;
       const amrex::Real r_sq = 4.0*params->gaussTime*params->gaussDiff;
       if ( params->gauss_type == 0 ) { // Spec
          massfrac[sp1_ID] = denom_inv * std::exp(-d_sq/r_sq);
          massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
       } else if ( params->gauss_type == 1 ) { //Temp
          // TODO: find a better way, one with an analytical solution
          state(i,j,k,TEMP) = params->T_mean * ( 1.0 + params->ampgauss * std::exp(-d_sq/r_sq) );
       } else {
          amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
       }
    }

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

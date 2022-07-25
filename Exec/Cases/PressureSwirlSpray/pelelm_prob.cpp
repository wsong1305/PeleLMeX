#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   // Gas params
   pp.query("P_mean", PeleLM::prob_parm->P_mean);
   pp.query("T_mean", PeleLM::prob_parm->T_mean);
   pp.query("Vcoflow", PeleLM::prob_parm->Vcoflow);

   // Spray params
   pp.query("mass_flow_rate",PeleLM::prob_parm->mass_flow_rate);
   pp.query("spray_temp",PeleLM::prob_parm->part_temp);
   pp.query("spray_start_time",PeleLM::prob_parm->jet_start_time);
   pp.query("spray_end_time",PeleLM::prob_parm->jet_end_time);
   pp.query("pressure_swirl_theta",PeleLM::prob_parm->ps_halfangle);
   pp.query("pressure_swirl_r0",PeleLM::prob_parm->ps_r0);
   pp.query("pressure_swirl_rmsvel",PeleLM::prob_parm->ps_rmsvel);
}

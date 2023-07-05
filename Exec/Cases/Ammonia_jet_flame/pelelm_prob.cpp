#include <PeleLM.H>
#include <pelelm_prob.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   pp.query("P_mean",   prob_parm->P_mean);
   pp.query("v_in",  prob_parm->v_in);
   pp.query("T_in",  prob_parm->T_in);

   amrex::ParmParse ppEB("EB");
   ppEB.query("in_diam",prob_parm->inflow_diam);
   
}

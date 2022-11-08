#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuContainers.H>

void PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

   pp.query("P_mean"      , PeleLM::prob_parm->P_mean);
   pp.query("T_in"        , PeleLM::prob_parm->T_in);
   pp.query("U_jet"       , PeleLM::prob_parm->U_jet);
   pp.query("D_jet"       , PeleLM::prob_parm->D_jet);
   pp.query("T_jet"       , PeleLM::prob_parm->T_jet);
   
   PeleLM::prob_parm->bathID = N2_ID;  
   PeleLM::prob_parm->fuelID = POSF11498_ID;  
   PeleLM::prob_parm->oxidID = O2_ID; 
   
   PeleLM::pmf_data.initialize();
 

}

#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

   pp.query("P_mean"       , PeleLM::prob_parm->P_mean);
   pp.query("Zst"          ,    PeleLM::prob_parm->Zst);
   pp.query("T_in"         ,   PeleLM::prob_parm->T_in);
   pp.query("U_b"          ,    PeleLM::prob_parm->U_b);
   pp.query("standoff"     , PeleLM::prob_parm->standoff);


   PeleLM::prob_parm->bathID = N2_ID;  
   PeleLM::prob_parm->fuelID = NXC7H16_ID;  
   //PeleLM::prob_parm->fuelID = NC7H16_ID;  
   PeleLM::prob_parm->oxidID = O2_ID; 

   PeleLM::pmf_data.initialize();
}

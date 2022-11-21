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
   pp.query("pertmag"      ,  PeleLM::prob_parm->pertmag);
   pp.query("amplification",  PeleLM::prob_parm->amplification);
   pp.query("Do_swirler"   ,  PeleLM::prob_parm->Do_swirler);
   pp.query("Di_swirler"   ,  PeleLM::prob_parm->Di_swirler);
   pp.query("D_jet"        ,  PeleLM::prob_parm->D_jet);
   pp.query("U_jet"        ,  PeleLM::prob_parm->U_jet);
   pp.query("phi_jet"      ,  PeleLM::prob_parm->phi_jet);
   pp.query("swirl_angle"  ,  PeleLM::prob_parm->swirl_angle);
   pp.query("ignition"     ,  PeleLM::prob_parm->ignition);
   std::vector<amrex::Real> ign_loc(3,0.);
   pp.queryarr("ign_region", ign_loc);
   for (int n =0; n < 3; n++){
	   PeleLM::prob_parm->ign_region[n] = ign_loc[n];
   }

   PeleLM::prob_parm->bathID = N2_ID;  
   PeleLM::prob_parm->fuelID = NXC7H16_ID;  
   //PeleLM::prob_parm->fuelID = NC7H16_ID;  
   PeleLM::prob_parm->oxidID = O2_ID; 
}

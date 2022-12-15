#include <PeleLM.H>
#include <pelelm_prob.H>


void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   pp.query("P_mean"       ,  PeleLM::prob_parm->P_mean);
   pp.query("Zst"          ,  PeleLM::prob_parm->Zst);
   pp.query("T_in"         ,  PeleLM::prob_parm->T_in);
   pp.query("U_b"          ,  PeleLM::prob_parm->U_b);

   amrex::Real phi_in;
   pp.query("phi_in"       ,  phi_in);
   pp.query("phi_target"   ,  PeleLM::prob_parm->phi_target);
   pp.query("t_start"      ,  PeleLM::prob_parm->t_start);
   pp.query("dt_phi_change",  PeleLM::prob_parm->dt_phi_change);

   pp.query("ignition"      ,  PeleLM::prob_parm->ignition);

   std::vector<amrex::Real> region(AMREX_SPACEDIM);
   pp.getarr("ign_region", region);

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
    PeleLM::prob_parm->ign_region[idir] = region[idir];


   PeleLM::prob_parm->phi_in = phi_in;

   auto problo = geom[0].ProbLo();
   auto probhi = geom[0].ProbHi();

   PeleLM::prob_parm->splitx = 0.5 * (problo[0] + probhi[0]);
   PeleLM::prob_parm->midtanh = 0.6 * (problo[0] + probhi[0]);
   PeleLM::prob_parm->widthtanh = 0.05 * (problo[0] + probhi[0]); 

   PeleLM::prob_parm->bathID = N2_ID;  
   PeleLM::prob_parm->fuelID = POSF10325_ID;  
   PeleLM::prob_parm->oxidID = O2_ID; 

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real Xt[NUM_SPECIES] = {0.0};
   amrex::Real massfrac[NUM_SPECIES] = {0.0};
   amrex::Real a = 16.5;
   

   // Xt[O2_ID] = 1.0 / ( 1.0 + phi_in / a + 0.79 / 0.21 );
   // Xt[POSF10325_ID] = phi_in * Xt[O2_ID] / a;
   // Xt[N2_ID] = 1.0 - Xt[O2_ID] - Xt[POSF10325_ID];

   Xt[O2_ID] = 0.21;
   Xt[N2_ID] = 0.79;

   eos.X2Y(Xt,massfrac);


  for (int n = 0; n < NUM_SPECIES; n++)
    (PeleLM::prob_parm->Ys)[n] = massfrac[n];

  PeleLM::pmf_data.initialize();

}

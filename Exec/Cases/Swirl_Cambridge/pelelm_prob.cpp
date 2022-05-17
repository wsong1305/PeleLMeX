#include <PeleLM.H>
#include <pelelm_prob.H>

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
   pp.query("jet_angle"    ,  PeleLM::prob_parm->jet_angle);
   pp.query("phi_jet"      ,  PeleLM::prob_parm->phi_jet);
   pp.query("swirl_angle"  ,  PeleLM::prob_parm->swirl_angle);

   //Spray stuff
   pp.query("jet_dx_mod", PeleLM::prob_parm->jet_dx_mod);
   pp.get("jet_dia", PeleLM::prob_parm->jet_dia);
   amrex::Vector<amrex::Real> jet_loc = {{0.0, 0.0, 0.0}};
   pp.queryarr("jet_location",jet_loc);
   pp.get("part_mean_dia", PeleLM::prob_parm->part_mean_dia);
   pp.query("part_stdev_dia", PeleLM::prob_parm->part_stdev_dia);
   pp.get("part_temp", PeleLM::prob_parm->part_temp);
   pp.query("mass_flow_rate", PeleLM::prob_parm->mass_flow_rate);
   pp.get("spray_angle_deg", PeleLM::prob_parm->spray_angle);
   pp.query("gas_jet_dia", PeleLM::prob_parm->gas_jet_dia);
   pp.query("gas_jet_vel", PeleLM::prob_parm->gas_jet_vel);
   pp.query("jet_vel", PeleLM::prob_parm->jet_vel);
   std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
   in_Y_jet[0] = 1.;
   pp.queryarr("jet_mass_fracs", in_Y_jet);
   amrex::Real sumY = 0.;
   for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
     PeleLM::prob_parm->Y_jet[spf] = in_Y_jet[spf];
     sumY += in_Y_jet[spf];
   }
   amrex::Print() << " Spray injection at: x = " << jet_loc[0] << " y = " << jet_loc[1] << " z = " << jet_loc[2] << std::endl;
   // Convert to radians
   PeleLM::prob_parm->spray_angle *= M_PI / 180.;
   amrex::Real dom_len = geom[0].ProbHi(0) - geom[0].ProbLo(0);
   amrex::Real yloc = jet_loc[0];
   amrex::Real xloc = jet_loc[1];
   amrex::Real zloc = jet_loc[2];
   AMREX_D_TERM(PeleLM::prob_parm->jet_cents[0] = xloc;,
                PeleLM::prob_parm->jet_cents[1] = yloc;,
                PeleLM::prob_parm->jet_cents[2] = zloc;)

   // PeleLM::pmf_data.initialize();

   auto problo = geom[0].ProbLo();
   auto probhi = geom[0].ProbHi();

   PeleLM::prob_parm->splitx = 0.5 * (problo[0] + probhi[0]);
   PeleLM::prob_parm->midtanh = 0.6 * (problo[0] + probhi[0]);
   PeleLM::prob_parm->widthtanh = 0.05 * (problo[0] + probhi[0]); 

   PeleLM::prob_parm->bathID = N2_ID;  
   PeleLM::prob_parm->fuelID = NC7H16_ID;  
   PeleLM::prob_parm->oxidID = O2_ID; 
}

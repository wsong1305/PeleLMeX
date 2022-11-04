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
   
   std::string infile;
   pp.query("datafile"    , infile);

  //check if the file exists
  std::ifstream iffile(infile);
  const std::string memfile =
    static_cast<std::stringstream const&>(std::stringstream() << iffile.rdbuf())
      .str();
  if (!amrex::FileSystem::Exists(infile)) {
    amrex::Abort("injection_file does not exist");
  }

  //count the lines for memory allocation
  iffile.close();
  std::istringstream iss(memfile);
  std::string firstline, remaininglines;
  std::getline(iss, firstline);
  int line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }

  amrex::Print() << "Found " << line_count << " lines in " << infile << std::endl;
  PeleLM::prob_parm->nlines = line_count;

  int nRadialValues = 1;
  // std::array<amrex::Vector<amrex::Real>, ncolumns> data;
  amrex::Vector<amrex::Real> r;
  amrex::Vector<amrex::Real> scalars;

  r      .resize(line_count);
  scalars.resize(line_count);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);

  for (unsigned int i = 0; i < line_count; ++i) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> r[i];
    for (int col = 0; col < nRadialValues; ++col) {
      sinput >> scalars[i+col+1];
    }
  }

  amrex::Print() << "Finished reading data from " << infile << std::endl;

  PeleLM::prob_parm->radial_loc    = (amrex::Real*) amrex::The_Arena()->alloc(line_count*sizeof(amrex::Real));
  PeleLM::prob_parm->mixfrac_value = (amrex::Real*) amrex::The_Arena()->alloc(line_count*nRadialValues*sizeof(amrex::Real));

  amrex::Gpu::copy(amrex::Gpu::hostToDevice,       r.begin(),       r.end(), PeleLM::prob_parm->radial_loc);
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, scalars.begin(), scalars.end(), PeleLM::prob_parm->mixfrac_value);


}

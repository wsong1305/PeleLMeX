#include <PeleLMDeriveFunc.H>
#include <PeleLM_Index.H>
#include <PelePhysics.H>
#include <mechanism.H>
#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

//
// User-defined derived variables list
//
Vector<std::string> pelelm_setuserderives()
{
  //Vector<std::string> var_names({"derUserDefine_null"});
  return {"derUserDefine_null"}; //var_names;
}

//
// User-defined derived definition
//
void pelelm_deruserdef (PeleLM* /*a_pelelm*/, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                        const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                        const Geometry& /*geom*/, Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    // Checks
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);

    // Array4s
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);

    /* Implement your own function in the lambda below
       Example of derived are available in Source/PeleLMDeriveFunc.cpp
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = in_dat(i,j,k,DENSITY) * in_dat(i,j,k,VELX);
    });
    */

    Abort("Using derUserDefine derived requires providing a definition in local DeriveUserDefined.cpp");
}

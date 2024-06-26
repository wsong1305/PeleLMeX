#ifndef PELELM_USERKEYS_H
#define PELELM_USERKEYS_H
#include <AMReX_Array.H>

#include <string>
#include <map>

/**
 * \brief struct holding PeleLMeX physical BC options
 */
struct BoundaryCondition
{
  BoundaryCondition() = default;
  enum {
    BCInterior = 0,
    BCInflow,
    BCOutflow,
    BCSymmetry,
    BCSlipWallAdiab,
    BCNoSlipWallAdiab,
    BCSlipWallIsotherm,
    BCNoSlipWallIsotherm
  };
  const std::map<const std::string, int> str2int = {
    {"interior", BCInterior},
    {"inflow", BCInflow},
    {"outflow", BCOutflow},
    {"symmetry", BCSymmetry},
    {"slipwalladiab", BCSlipWallAdiab},
    {"noslipwalladiab", BCNoSlipWallAdiab},
    {"slipwallisotherm", BCSlipWallIsotherm},
    {"noslipwallisotherm", BCNoSlipWallIsotherm}};
  const amrex::Array<std::string, 2> searchKey{"lo_bc", "hi_bc"};
};
const BoundaryCondition boundarycondition;

/**
  * \brief struct holding PeleLMeX NS solver options
           default is LowMachNumber
  */
struct NSSolver
{
  NSSolver() = default;
  enum { LowMachNumber = 0, Incompressible };
  const std::map<std::string, int> str2int = {
    {"lowmachnumber", LowMachNumber},
    {"incompressible", Incompressible},
    {"default", LowMachNumber}};
  const amrex::Array<std::string, 1> searchKey{"ns_solver"};
};
const NSSolver nssolver;

/**
  * \brief struct holding PeleLMeX interpolator options
           default is PiecewiseLinearConserv
  */
struct Interpolator
{
  Interpolator() = default;
  enum {
    PiecewiseLinearConserv = 0,
    PiecewiseLinearConservMinMax,
    PiecewiseConstant
  };
  const std::map<const std::string, int> str2int = {
    {"pwlinearconserv", PiecewiseLinearConserv},
    {"pwlinearconservminmax", PiecewiseLinearConservMinMax},
    {"pwconstant", PiecewiseConstant},
    {"default", PiecewiseLinearConserv}};
  const amrex::Array<std::string, 1> searchKey{"regrid_interp_method"};
};
const Interpolator interpolator;

/**
  * \brief struct holding PeleLMeX LES SGS model options
           default is None
  */
struct LESModel
{
  LESModel() = default;
  enum { None = 0, Smagorinsky, WALE, Sigma };
  const std::map<const std::string, int> str2int = {
    {"none", None},
    {"smagorinsky", Smagorinsky},
    {"wale", WALE},
    {"sigma", Sigma},
    {"default", None}};
  const amrex::Array<std::string, 1> searchKey{"les_model"};
};
const LESModel lesmodel;

/**
  * \brief struct holding PeleLMeX advection scheme options
           default is Godunov_PLM
  */
struct AdvectionScheme
{
  AdvectionScheme() = default;
  enum {
    Godunov_PLM = 0,
    Godunov_PPM,
    Godunov_PPM_WENOZ,
    Godunov_PPM_NOLIM,
    Godunov_BDS
  };
  const std::map<const std::string, int> str2int = {
    {"godunov_plm", Godunov_PLM},
    {"godunov_ppm", Godunov_PPM},
    {"godunov_ppm_wenoz", Godunov_PPM_WENOZ},
    {"godunov_ppm_nolim", Godunov_PPM_NOLIM},
    {"godunov_bds", Godunov_BDS},
    {"default", Godunov_PLM}};
  const amrex::Array<std::string, 1> searchKey{"advection_scheme"};
};
const AdvectionScheme advscheme;

/**
  * \brief struct holding PeleLMeX load balancing cost options
           default is Ncell
  */
struct LoadBalanceCost
{
  LoadBalanceCost() = default;
  enum {
    Ncell = 0,
    ChemFunctCallAvg,
    ChemFunctCallMax,
    ChemFunctCallSum,
    UserDefinedDerivedAvg,
    UserDefinedDerivedSum,
    Timers
  };
  const std::map<const std::string, int> str2int = {
    {"ncell", Ncell},
    {"chemfunctcall_avg", ChemFunctCallAvg},
    {"chemfunctcall_max", ChemFunctCallMax},
    {"chemfunctcall_sum", ChemFunctCallSum},
    {"userdefined_avg", UserDefinedDerivedAvg},
    {"userdefined_sum", UserDefinedDerivedSum},
    {"default", Ncell}};
  const amrex::Array<std::string, 2> searchKey{
    "load_balancing_cost_estimate", "chem_load_balancing_cost_estimate"};
};
const LoadBalanceCost lbcost;

/**
  * \brief struct holding PeleLMeX load balancing method options
           default is SFC
  */
struct LoadBalanceMethod
{
  LoadBalanceMethod() = default;
  enum { SFC = 0, Knapsack, RoundRobin };
  const std::map<const std::string, int> str2int = {
    {"sfc", SFC},
    {"knapsack", Knapsack},
    {"roundrobin", RoundRobin},
    {"default", SFC}};
  const amrex::Array<std::string, 2> searchKey{
    "load_balancing_method", "chem_load_balancing_method"};
};
const LoadBalanceMethod lbmethod;
#endif

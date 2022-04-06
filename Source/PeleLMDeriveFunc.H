#ifndef PELELM_DERIVEFUNC_H_
#define PELELM_DERIVEFUNC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>
#ifdef PELE_USE_EFIELD
#include <PeleLMEFDeriveFunc.H>
#endif

class PeleLM;

// Single level derived
void pelelm_dermassfrac (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_dertemp (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                     const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                     const amrex::Geometry& geomdata,
                     amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derkineticenergy (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                              const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                              const amrex::Geometry& geomdata,
                              amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_deravgpress (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_dermgvort (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                       const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                       const amrex::Geometry& geomdata,
                       amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_dermixfrac (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                        const amrex::Geometry& geomdata,
                        amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derprogvar (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                        const amrex::Geometry& geomdata,
                        amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

// Multi level derived
// Surface Density Function
void pelelm_derSDF(PeleLM* a_pelelm,
                   amrex::Vector<amrex::MultiFab *> dervec, int dcomp, int ncomp,
                   amrex::Vector<const amrex::MultiFab *> statevec, amrex::Vector<const amrex::MultiFab *> pressvec,
                   const amrex::Vector<amrex::Geometry>  &a_geoms, amrex::Real time, const amrex::Vector<amrex::BCRec> &a_bcrecs);
#endif

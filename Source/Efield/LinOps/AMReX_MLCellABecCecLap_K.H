#ifndef AMREX_MLCELLABECCECLAP_K_H
#define AMREX_MLCELLABECCECLAP_K_H
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLCellABecCecLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLCellABecCecLap_2D_K.H>
#else
#include <AMReX_MLCellABecCecLap_3D_K.H>
#endif

#endif

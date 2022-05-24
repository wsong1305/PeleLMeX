PELEMP_HOME=/lustre/scratch/bsouzas/PeleMP/
PELELMEX_HOME=/lustre/scratch/bsouzas/PeleLMeX/
AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex
PELE_PHYSICS_HOME  ?= ${PELELMEX_HOME}/Submodules/PelePhysics
AMREX_HYDRO_HOME   ?= ${PELELMEX_HOME}/Submodules/AMReX-Hydro

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = TRUE
TINY_PROFILE    = TRUE
USE_EB          = TRUE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE

# PelePhysics
Chemistry_Model = heptane_fc
Eos_Model       = Fuego
Transport_Model = Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX

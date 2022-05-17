TOP = ../../../../PeleLMeX
AMREX_HOME         ?= ${TOP}/Submodules/amrex
PELELMEX_HOME      ?= ${TOP}
PELE_PHYSICS_HOME  ?= ${TOP}/Submodules/PelePhysics
AMREX_HYDRO_HOME   ?= ${TOP}/Submodules/AMReX-Hydro
PELEMP_HOME=/lustre/scratch/bsouzas/PeleMP

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE
USE_EB	        = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE
USE_TURBINFLOW  = TRUE

# PelePhysics
Chemistry_Model = heptane_3sp
Eos_Model		= Fuego
Transport_Model = Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX

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
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = TRUE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE
USE_TURBINFLOW  = FALSE

# PelePhysics
Chemistry_Model = ARLMech-HC-F24
Eos_Model	= Fuego
Transport_Model = Simple

USE_EB          = TRUE
CEXE_headers    += EBUserDefined.H

USE_PARTICLES = FALSE
SPRAY_FUEL_NUM = 1

ifeq ($(USE_CUDA),TRUE)
DEFINES += -DAMREX_GPUS_PER_NODE=2
endif

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX

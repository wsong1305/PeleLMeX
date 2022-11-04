PELEMP_HOME=/lustre/scratch/bsouzas/PeleMP/

TOP = ../../../..
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
# AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex_hypre
AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex
PELE_PHYSICS_HOME  ?= ${PELELMEX_HOME}/Submodules/PelePhysics
AMREX_HYDRO_HOME   ?= ${PELELMEX_HOME}/Submodules/AMReX-Hydro
# HYPRE_HOME         ?= ${PELELMEX_HOME}/Submodules/hypre
HYPRE_HOME          = /lustre/scratch/bsouzas/PeleLMeX/Submodules/hypre/install

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE
USE_EB          = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE
USE_HYPRE       = FALSE

# PelePhysics
Chemistry_Model = c1
# Chemistry_Model = heptane_fc
# Chemistry_Model = air
Eos_Model       = Fuego
Transport_Model = Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX

ifeq ($(USE_CUDA),TRUE)
DEFINES += -DAMREX_GPUS_PER_NODE=2
endif


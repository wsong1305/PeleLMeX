
TOP = //lustre/scratch/bsouzas/
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
# AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex_hypre
AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex
PELE_PHYSICS_HOME  ?= ${PELELMEX_HOME}/Submodules/PelePhysics
AMREX_HYDRO_HOME   ?= ${PELELMEX_HOME}/Submodules/AMReX-Hydro
HYPRE_HOME         = /pscratch/sd/e/esclapez/hypre_install

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE
USE_EB          = TRUE

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
# Chemistry_Model = heptane_fc
Chemistry_Model = air
Eos_Model       = Fuego
Transport_Model = Simple


include $(PELELMEX_HOME)/Exec/Make.PeleLMeX

ifeq ($(USE_CUDA),TRUE)
DEFINES += -DAMREX_GPUS_PER_NODE=4
endif

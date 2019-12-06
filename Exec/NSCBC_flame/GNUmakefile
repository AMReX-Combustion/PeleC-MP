PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE

DIM        = 3

COMP	   = gcc

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_REACT  = TRUE

USE_ASCENT = TRUE

# Include support for participating or tracer particles
USE_PARTICLES = TRUE

SPRAY_COMPONENTS = 8

DEFINES += -DSPRAY_COMPONENTS=$(SPRAY_COMPONENTS)

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
Chemistry_Model := dodecane_lu

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib
#Transport_dir := Constant

Bpack   := ./Make.package 
Blocs   := .

include ${PELEC_HOME}/Exec/Make.PeleC

ifeq ($(USE_PARTICLES), TRUE)
DEFINES += -DPARTICLES
endif
ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif

ifeq ($(USE_ASCENT), TRUE)
DEFINES +=-DUSE_ASCENT 
ASCENT_INSTALL_DIR=/gpfs/alpine/proj-shared/csc340/larsen/stats_no_cuda_shared/ascent/install
include $(ASCENT_INSTALL_DIR)/share/ascent/ascent_config.mk
include $(AMREX_HOME)/Src/Extern/Conduit/Make.package
ASCENT_LINK_FLAGS = $(subst -pthread,, $(ASCENT_MPI_LIB_FLAGS))
INCLUDE_LOCATIONS += $(ASCENT_INCLUDE_FLAGS)
INCLUDE_LOCATIONS += /ccs/home/arientm/ECP/repoversion/amrex/Src/Extern/Conduit
LIBRARIES += $(ASCENT_LINK_FLAGS) 
endif
# To add lblas
LIBRARY_LOCATIONS +=/gpfs/alpine/proj-shared/csc340/arientm/sw4_uberenv/uberenv_libs/spack/opt/spack/linux-rhel7-ppc64le/gcc-5.4.0/netlib-lapack-3.8.0-2zkztxddylrbg4yua4bgm7afwcsgyh34/lib64/
override XTRALIBS += -lblas

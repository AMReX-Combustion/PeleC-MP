PRECISION  = DOUBLE
PROFILE    = FALSE
USE_EB     = FALSE
DEBUG      = FALSE

DIM        = 2

COMP	   = gcc

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_CUDA   = TRUE

# Include particles (currently spray) capability/sources
USE_PARTICLES = TRUE

#The belos spec is Hacky max. 
SPRAY_COMPONENTS = 8

DEFINES += -DSPRAY_COMPONENTS=$(SPRAY_COMPONENTS)

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego 

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego 

Chemistry_Model := decane_3sp

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib
#Transport_dir := Constant

Bpack   := ./Make.package

include ${PELEC_HOME}/Exec/Make.PeleC

ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif



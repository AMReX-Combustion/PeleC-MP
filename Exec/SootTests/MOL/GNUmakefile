PRECISION  = DOUBLE
PROFILE    = FALSE

DIM        = 2

COMP	   = gcc

USE_REACT = TRUE

USE_MPI    = FALSE
USE_OMP    = FALSE
USE_GRAV   = FALSE

USE_EB = TRUE

HYP_TYPE = MOL

USE_SOOT = TRUE

#DEBUG = TRUE
DEBUG = FALSE

# If this is changes, must run a make clean and rerun make
#NUM_SOOT_MOMENTS = 3
NUM_SOOT_MOMENTS = 6

DEFINES += -DNUM_SOOT_MOMENTS=$(NUM_SOOT_MOMENTS)

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
Chemistry_Model := SootReaction

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib

Bpack   := ./Make.package 
Blocs   := .

include ${PELEC_HOME}/Exec/Make.PeleC

ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif

#ifndef SPRAY_KERNELS_H_
#define SPRAY_KERNELS_H_


#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void passive_move(SprayParticleContainer::ParticleType& p,
                  amrex::Array4<amrex::Real const> const& state_arr,
                  const int rho_idx,
                  const int xmom_idx,
                  const int ymom_idx,
                  const int zmom_idx,
                  const amrex::Real* dx,
                  const amrex::Real* plo,
                  const amrex::Real dt)
{


}


#endif

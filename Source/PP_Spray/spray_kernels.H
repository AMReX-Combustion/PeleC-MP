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

    amrex::Real x = (p.pos(0) - plo[0]) / dx[0];
    amrex::Real y = (p.pos(1) - plo[1]) / dx[1];
    amrex::Real z = (p.pos(2) - plo[2]) / dx[2];

    int i = std::floor(x);
    int j = std::floor(y);
    int k = std::floor(z);

    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;

    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi; 

    amrex::Real coef_lll = wx_lo * wy_lo * wz_lo;
    amrex::Real coef_hll = wx_hi * wy_lo * wz_lo;
    amrex::Real coef_lhl = wx_lo * wy_hi * wz_lo;
    amrex::Real coef_hhl = wx_hi * wy_hi * wz_lo;
    amrex::Real coef_llh = wx_lo * wy_lo * wz_hi;
    amrex::Real coef_hlh = wx_hi * wy_lo * wz_hi;
    amrex::Real coef_lhh = wx_lo * wy_hi * wz_hi;
    amrex::Real coef_hhh = wx_hi * wy_hi * wz_hi;

    amrex::Real u_p = 
                coef_lll*state_arr(i-1, j-1, k-1, xmom_idx)/state_arr(i-1, j-1, k-1, rho_idx) + 
                coef_llh*state_arr(i-1, j-1, k  , xmom_idx)/state_arr(i-1, j-1, k  , rho_idx) + 
                coef_lhl*state_arr(i-1, j,   k-1, xmom_idx)/state_arr(i-1, j,   k-1, rho_idx) + 
                coef_lhh*state_arr(i-1, j,   k  , xmom_idx)/state_arr(i-1, j,   k  , rho_idx) + 
                coef_hll*state_arr(i,   j-1, k-1, xmom_idx)/state_arr(i,   j-1, k-1, rho_idx) + 
                coef_hlh*state_arr(i,   j-1, k  , xmom_idx)/state_arr(i,   j-1, k  , rho_idx) + 
                coef_hhl*state_arr(i,   j,   k-1, xmom_idx)/state_arr(i,   j,   k-1, rho_idx) + 
                coef_hhh*state_arr(i,   j,   k  , xmom_idx)/state_arr(i,   j,   k  , rho_idx) ;
}


#endif

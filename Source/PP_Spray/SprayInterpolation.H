#ifndef _SPRAYINTERPOLATION_H_
#define _SPRAYINTERPOLATION_H_

// Slightly modified from MFIX code

/****************************************************************
 Functions for interpolation on non-EB mesh
 ***************************************************************/

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void trilinear_interp (const amrex::RealVect& pos,
                       const amrex::RealVect& plo,
                       const amrex::RealVect& dxi,
                       amrex::IntVect* indx_array,
                       amrex::Real* weights)
{
    const amrex::RealVect lx = (pos - plo)*dxi + 0.5;
    const amrex::IntVect ijk = lx.floor();

    AMREX_D_TERM(int i = ijk[0];, int j = ijk[1];, int k = ijk[2];);


    const amrex::RealVect sx_hi = lx - ijk;
    const amrex::RealVect sx_lo = 1 - sx_hi;
    amrex::GpuArray<amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>,2> ssv;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      ssv[dir][0] = sx_lo[dir];
      ssv[dir][1] = sx_hi[dir];
    }
    int cc = 0;
#if AMREX_SPACEDIM == 3
    for(int kk(-1); kk<1; kk++){
#endif
#if AMREX_SPACEDIM > 1
      for (int jj(-1); jj<1; jj++){
#endif
        for (int ii(-1); ii<1; ii++){
          AMREX_D_TERM(indx_array[cc][0] = i + ii;,
                       indx_array[cc][1] = j + jj;,
                       indx_array[cc][2] = k + kk;);
          weights[cc] = AMREX_D_TERM(ssv[0][ii+1],*ssv[1][jj+1],*ssv[2][kk+1]);
          cc++;
        }
#if AMREX_SPACEDIM > 1
      }
#endif
#if AMREX_SPACEDIM == 3
    }
#endif
}

#ifdef AMREX_USE_EB

/****************************************************************
 Functions for the Newtons solver
 ***************************************************************/

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real f (const int dir,
               const amrex::RealVect& pos,
               const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
               const amrex::Real& xi,
               const amrex::Real& eta,
               const amrex::Real& zeta)
{
    amrex::Real a0 = (nodes[0][dir] - pos[dir]);
    amrex::Real a1 = (nodes[1][dir] - nodes[0][dir]);
    amrex::Real a2 = (nodes[2][dir] - nodes[0][dir]);
    amrex::Real a3 = (nodes[4][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] - nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

    return a0 + a1*xi + a2*eta + a3*zeta + a4*xi*eta + a5*xi*zeta + a6*eta*zeta + a7*xi*eta*zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dfdxi (const int dir,
                   const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
                   const amrex::Real& xi,
                   const amrex::Real& eta,
                   const amrex::Real& zeta)
{
    amrex::Real a1 = (nodes[1][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] - nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

    return a1 + a4*eta + a5*zeta + a7*eta*zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dfdeta (const int dir,
                    const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
                    const amrex::Real& xi,
                    const amrex::Real& eta,
                    const amrex::Real& zeta)
{

    amrex::Real a2 = (nodes[2][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] - nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

    return a2 + a4*xi + a6*zeta + a7*xi*zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dfdzeta (const int dir,
                     const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
                     const amrex::Real& xi,
                     const amrex::Real& eta,
                     const amrex::Real& zeta)
{

    amrex::Real a3 = (nodes[4][dir] - nodes[0][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] - nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

    return a3 + a5*xi + a6*eta + a7*xi*eta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_interp_mapping (const amrex::RealVect& pos,
                         const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
                         amrex::Real& xi,
                         amrex::Real& eta,
                         amrex::Real& zeta)
{

#if(0)
    amrex::Print() << 0
                   << "  xi: " << xi << "  "
                   << " eta: " << eta << "  "
                   << "zeta: " << zeta << std::endl;
#endif

    int lc(0);
    amrex::Real err(1.0);

    while(err > 1.0e-3 && lc < 10) {

        amrex::Real f0 = f(0, pos, nodes, xi, eta, zeta);
        amrex::Real f1 = f(1, pos, nodes, xi, eta, zeta);
        amrex::Real f2 = f(2, pos, nodes, xi, eta, zeta);

        amrex::Real df0dxi   = dfdxi(0, nodes, xi, eta, zeta);
        amrex::Real df0deta  = dfdeta(0, nodes, xi, eta, zeta);
        amrex::Real df0dzeta = dfdzeta(0, nodes, xi, eta, zeta);

        amrex::Real df1dxi   = dfdxi(1, nodes, xi, eta, zeta);
        amrex::Real df1deta  = dfdeta(1, nodes, xi, eta, zeta);
        amrex::Real df1dzeta = dfdzeta(1, nodes, xi, eta, zeta);

        amrex::Real df2dxi   = dfdxi(2, nodes, xi, eta, zeta);
        amrex::Real df2deta  = dfdeta(2, nodes, xi, eta, zeta);
        amrex::Real df2dzeta = dfdzeta(2, nodes, xi, eta, zeta);

        amrex::Real detJ = df0dxi   * (df1deta*df2dzeta - df1dzeta*df2deta)
                         - df0deta  * (df1dxi*df2dzeta  - df1dzeta*df2dxi)
                         + df0dzeta * (df1dxi*df2deta   - df1deta*df2dxi);

        amrex::Real detJ_xi = f0       * (df1deta*df2dzeta - df1dzeta*df2deta)
                            - df0deta  * (f1*df2dzeta      - df1dzeta*f2)
                            + df0dzeta * (f1*df2deta       - df1deta*f2);

        amrex::Real detJ_eta = df0dxi   * (f1*df2dzeta     - df1dzeta*f2)
                             - f0       * (df1dxi*df2dzeta - df1dzeta*df2dxi)
                             + df0dzeta * (df1dxi*f2       - f1*df2dxi);

        amrex::Real detJ_zeta = df0dxi   * (df1deta*f2     - f1*df2deta)
                              - df0deta  * (df1dxi*f2      - f1*df2dxi)
                              + f0       * (df1dxi*df2deta - df1deta*df2dxi);

#if(0)
        amrex::Print() << std::endl << std::endl;
        amrex::Print() << "detJ      " << detJ << std::endl;
        amrex::Print() << "detJ_xi   " << detJ_xi << std::endl;
        amrex::Print() << "detJ_eta  " << detJ_eta << std::endl;
        amrex::Print() << "detJ_zeta " << detJ_zeta << std::endl;
        amrex::Print() << std::endl << std::endl;
#endif
        amrex::Real new_xi   = xi   - detJ_xi   / detJ;
        amrex::Real new_eta  = eta  - detJ_eta  / detJ;
        amrex::Real new_zeta = zeta - detJ_zeta / detJ;

        lc += 1;
        err = amrex::max(amrex::Math::abs(xi-new_xi), amrex::Math::abs(eta-new_eta), amrex::Math::abs(zeta-new_zeta));

        xi   = new_xi;
        eta  = new_eta;
        zeta = new_zeta;

#if(0)
    amrex::Print() << lc
                   << "  xi: " << xi << "  "
                   << " eta: " << eta << "  "
                   << "zeta: " << zeta << "  "
                   << " err: " << err << std::endl;
#endif
    }

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fe_interp(const amrex::RealVect& pos,
               const int ip, const int jp, const int kp,
               const amrex::RealVect& dx,
               const amrex::RealVect& dxi,
               const amrex::RealVect& plo,
               amrex::Array4<const EBCellFlag> const& flags,
               amrex::Array4<const amrex::Real> const& ccent,
               amrex::Array4<const amrex::Real> const& bcent,
               amrex::Array4<const amrex::Real> const& apx,
               amrex::Array4<const amrex::Real> const& apy,
               amrex::Array4<const amrex::Real> const& apz,
               amrex::Array4<const amrex::Real> const& vfrac,
               amrex::IntVect* indx_array,
               amrex::Real* weights)
{

    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

    const amrex::Real cdist_x = pos[0] - (ip + 0.5 + ccent(ip,jp,kp,0))*dx[0] - plo[0];
    const amrex::Real cdist_y = pos[1] - (jp + 0.5 + ccent(ip,jp,kp,1))*dx[1] - plo[1];
    const amrex::Real cdist_z = pos[2] - (kp + 0.5 + ccent(ip,jp,kp,2))*dx[2] - plo[2];
    const amrex::Real vtol = 0.05;

    // Distance between particle and cell centoid.
    const amrex::Real cdist =
      std::sqrt(cdist_x*cdist_x + cdist_y*cdist_y + cdist_z*cdist_z);

    // Before doing anything fancy, just check that the particle isn't overlapping
    // the cell center. If it is, use the cell value.
    if(cdist < tolerance){
      indx_array[0] = {ip,jp,kp};
      weights[0] = 1.;
      for (int aindx = 1; aindx < 8; ++aindx) {
        indx_array[aindx] = {ip,jp,kp};
        weights[aindx] = 0.;
      }
      return;
    }

    amrex::Real par_dot_EB(2.);
    amrex::Real cent_dot_EB(1.);

    // If the particle is inside a cut-cell, verify that it is on the correct
    // side of the EB before trying to interpolate.
    if(flags(ip,jp,kp).isSingleValued()){

        const Real adx = apx(ip+1,jp  ,kp  ) - apx(ip  ,jp  ,kp  );
        const Real ady = apy(ip  ,jp+1,kp  ) - apy(ip  ,jp  ,kp  );
        const Real adz = apz(ip  ,jp  ,kp+1) - apz(ip  ,jp  ,kp  );

        const Real apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);

        // Temporary sanity check.
        AMREX_ASSERT_WITH_MESSAGE( apnorm > std::numeric_limits<amrex::Real>::epsilon(),
            "Invalid EB normal. Not sure how to interpolate!");
        // EB normal
        const amrex::Real apnorminv = 1. / apnorm;
        const amrex::RealVect normal = {adx*apnorminv, ady*apnorminv, adz*apnorminv};

        // Projection of vector pointing from EB centroid to particle onto EB normal
        par_dot_EB = (pos[0] - (ip + 0.5 + bcent(ip,jp,kp,0))*dx[0] - plo[0])*normal[0] +
                     (pos[1] - (jp + 0.5 + bcent(ip,jp,kp,1))*dx[1] - plo[1])*normal[1] +
                     (pos[2] - (kp + 0.5 + bcent(ip,jp,kp,2))*dx[2] - plo[2])*normal[2];

        // Projection of vector pointing from EB centroid to cell centroid onto EB normal
        cent_dot_EB = (ccent(ip,jp,kp,0) - bcent(ip,jp,kp,0))*dx[0]*normal[0] +
                      (ccent(ip,jp,kp,1) - bcent(ip,jp,kp,1))*dx[1]*normal[1] +
                      (ccent(ip,jp,kp,2) - bcent(ip,jp,kp,2))*dx[2]*normal[2];
        // Temporary sanity check
        AMREX_ASSERT_WITH_MESSAGE( cent_dot_EB > tolerance,
            "cent_dot_EB < tolerance ... this makes no sense!");
    }
    // Particle position relative to cell center [-0.5, 0.5]
    amrex::Real gx = (pos[0]-plo[0])*dxi[0] - (ip + 0.5);
    amrex::Real gy = (pos[1]-plo[1])*dxi[1] - (jp + 0.5);
    amrex::Real gz = (pos[2]-plo[2])*dxi[2] - (kp + 0.5);

    // Use the centroid location of the cell containing the particle
    // to guess the interpolation stencil.
    const int i = (gx < ccent(ip,jp,kp,0)) ? ip : ip + 1;
    const int j = (gy < ccent(ip,jp,kp,1)) ? jp : jp + 1;
    const int k = (gz < ccent(ip,jp,kp,2)) ? kp : kp + 1;

    const int di = i - ip; // 0 or 1
    const int dj = j - jp; // 0 or 1
    const int dk = k - kp; // 0 or 1

    // Count the number of non-connected cells in the stencil
    int aindx = 0;
    int covered = 0;
    for(int kk(-1); kk<1; kk++){
      for(int jj(-1); jj<1; jj++){
        for(int ii(-1); ii<1; ii++){
          if(not flags(ip,jp,kp).isConnected(di+ii,dj+jj,dk+kk))
            covered += 1;
        }
      }
    }
    // A negative value implies that the particle is 'behind' the EB (and
    // therefore inside the wall). Although this shouldn't occur often, it could
    // for fast moving particles.  Set the fluid velocity to zero and use the
    // cell value for remaining items.
    if((flags(ip,jp,kp).isSingleValued() and (par_dot_EB <= tolerance)) ||
       flags(ip,jp,kp).isCovered()){
        amrex::Abort("SprayInterpolation.H::fe_interp(): Particle has penetrated an EB wall");
    // The particle is near the EB. It is either
    // 1) between the cell centroid and the EB, or
    // 2) one or more cells in the 2^3 stencil is covered.
    // Either way, interpolating to the particle is not straight forward so
    // (for now) we do a 1D interpolation to the EB.

    } else if((covered > 0) or (par_dot_EB < cent_dot_EB)){
        for (int aindx = 0; aindx < 8; ++aindx) {
          indx_array[aindx] = {ip,jp,kp};
          weights[aindx] = 0.;
        }
        weights[0] = 1.;
    } else {
        amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 8> nodes;

        // Here we pack the local 2x2x2 stencil into a local array. Note that
        // the node ordering is not consistent with the previous version.
        // Nodes 2 and 3 have been swapped as well as nodes 6 and 7. This was
        // to allow for a more compact for-loop filling.

        int lc(0);

        for(int kk(-1); kk<1; kk++){
            for(int jj(-1); jj<1; jj++){
                for(int ii(-1); ii<1; ii++){
                    nodes[lc][0] = (i + ii + 0.5 + ccent(i+ii,j+jj,k+kk,0))*dx[0];
                    nodes[lc][1] = (j + jj + 0.5 + ccent(i+ii,j+jj,k+kk,1))*dx[1];
                    nodes[lc][2] = (k + kk + 0.5 + ccent(i+ii,j+jj,k+kk,2))*dx[2];
                    lc+=1;
                }
            }
        }

        amrex::Real xi   = (pos[0] - nodes[0][0] - plo[0])*dxi[0];
        amrex::Real eta  = (pos[1] - nodes[0][1] - plo[1])*dxi[1];
        amrex::Real zeta = (pos[2] - nodes[0][2] - plo[2])*dxi[2];

        get_interp_mapping(pos-plo, nodes, xi, eta, zeta);
        indx_array[0] = {i-1,j-1,k-1};
        indx_array[1] = {i,  j-1,k-1};
        indx_array[2] = {i,  j,  k-1};
        indx_array[3] = {i-1,j,  k-1};
        indx_array[4] = {i-1,j-1,k};
        indx_array[5] = {i,  j-1,k};
        indx_array[6] = {i,  j,  k};
        indx_array[7] = {i-1,j,  k};
        weights[0] = 1 - xi - eta - zeta + xi*eta + xi*zeta + eta*zeta - xi*eta*zeta;
        weights[1] = xi*(1 - eta - zeta + eta*zeta);
        weights[2] = xi*eta - xi*eta*zeta;
        weights[3] = eta*(1 - xi - zeta + xi*zeta);
        weights[4] = zeta*(1 - xi - eta + xi*eta);
        weights[5] = xi*zeta - xi*eta*zeta;
        weights[6] = xi*eta*zeta;
        weights[7] = eta*zeta - xi*eta*zeta;
        // Rudimentary attempt to prevent source terms from being added
        // to small volume cells to avoid having to redistribute later
        // Will probably have to add a redistribute function
        amrex::Real rw = 0.;
        for (int aindx = 0; aindx < 8; ++aindx) {
          if (vfrac(indx_array[aindx]) < vtol) {
            weights[aindx] = 0.;
          }
          rw += weights[aindx];
        }
        rw = 1./rw;
        for (int aindx = 0; aindx < 8; ++aindx) {
          weights[aindx] *= rw;
        }
    }
    return;
}

#endif // AMREX_USE_EB
#endif // _SPRAYINTERPOLATION_H_



#ifndef _SOOTMODEL_DERIVE_H_
#define _SOOTMODEL_DERIVE_H_
#include <AMReX_FArrayBox.H>
#include <SootModel.H>

void soot_largeparticledata (const amrex::Box& bx, amrex::FArrayBox& slfab, int dcomp, int ncomp,
			     const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
			     amrex::Real time, const int* bcrec, int level);

void soot_genvars (const amrex::Box& bx, amrex::FArrayBox& slfab, int dcomp, int ncomp,
		   const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
		   amrex::Real time, const int* bcrec, int level);

#endif

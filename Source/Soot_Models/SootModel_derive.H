#ifndef _SOOTMODEL_DERIVE_H_
#define _SOOTMODEL_DERIVE_H_
#include <AMReX_FArrayBox.H>
#include <SootModel.H>

void soot_largepartnumdens (const amrex::Box& bx, amrex::FArrayBox& nlfab, int dcomp, int ncomp,
			    const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
			    amrex::Real time, const int* bcrec, int level);

void soot_largepartmeanvol (const amrex::Box& bx, amrex::FArrayBox& vlfab, int dcomp, int ncomp,
			    const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
			    amrex::Real time, const int* bcrec, int level);

void soot_largepartsurfarea (const amrex::Box& bx, amrex::FArrayBox& slfab, int dcomp, int ncomp,
			     const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
			     amrex::Real time, const int* bcrec, int level);

#endif

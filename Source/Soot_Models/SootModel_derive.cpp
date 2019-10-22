
#include <PeleC.H>
#include <SootModel.H>
#include "SootModel_derive.H"

using namespace amrex;

void soot_largepartnumdens (const amrex::Box& bx, amrex::FArrayBox& nlfab, int dcomp, int ncomp,
			    const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
			    amrex::Real time, const int* bcrec, int level)
{
  auto const dat = datafab.array();
  auto       nl  = nlfab.array();
  const auto lo = lbound(bx);
  const auto hi = ubound(bx);
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
	Real M00 = dat(i,j,k,PeleC::FirstSootVar);
	Real N0 = dat(i,j,k,PeleC::FirstSootVar + PeleC::NumSootVars - 1);
	nl(i,j,k,dcomp) = M00 - N0;
      }
    }
  }
}

void soot_largepartmeanvol (const amrex::Box& bx, amrex::FArrayBox& vlfab, int dcomp, int ncomp,
			    const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
			    amrex::Real time, const int* bcrec, int level)
{
  auto const dat = datafab.array();
  auto       vl  = vlfab.array();
  const auto lo = lbound(bx);
  const auto hi = ubound(bx);
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {	
	Real M00 = dat(i,j,k,PeleC::FirstSootVar);
	Real M10 = dat(i,j,k,PeleC::FirstSootVar + 1);
	Real N0 = dat(i,j,k,PeleC::FirstSootVar + PeleC::NumSootVars - 1);
	vl(i,j,k,dcomp) = (M10 - N0*SootModel::V0)/(M00 - N0 + 1.E-30);
      }
    }
  }
}

void soot_largepartsurfarea (const amrex::Box& bx, amrex::FArrayBox& slfab, int dcomp, int ncomp,
			     const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
			     amrex::Real time, const int* bcrec, int level)
{
  auto const dat = datafab.array();
  auto       sl  = slfab.array();
  const auto lo = lbound(bx);
  const auto hi = ubound(bx);
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {	
	Real M00 = dat(i,j,k,PeleC::FirstSootVar);
	Real M01 = dat(i,j,k,PeleC::FirstSootVar + 2);
	Real N0 = dat(i,j,k,PeleC::FirstSootVar + PeleC::NumSootVars - 1);
	sl(i,j,k,dcomp) = (M01 - N0*SootModel::S0)/(M00 - N0 + 1.E-30);
      }
    }
  }
}


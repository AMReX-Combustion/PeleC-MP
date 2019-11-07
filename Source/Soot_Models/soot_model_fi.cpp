
#include <AMReX_Print.H>
#include <SootModel.H>

using namespace amrex;
extern "C"
{
  void fi_init_soot_vars (Real* vals)
  {
    Vector<Real> moments_vec(NUM_SOOT_MOMENTS + 1);
    SootModel::initialSmallMomVals(moments_vec);
    for (int mom = 0; mom != NUM_SOOT_MOMENTS + 1; ++mom)
      {
	vals[mom] = moments_vec[mom];
      }
  }
}

#include <PeleLM.H>
#include <pelelm_prob.H>

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("ref_p", PeleLM::prob_parm->P_mean);
  pp.query("standoff", PeleLM::prob_parm->standoff);
  std::string pmf_datafile;
  pp.query("pmf_datafile", pmf_datafile);
  int pmf_do_average = 1;
  PMF::read_pmf(pmf_datafile, pmf_do_average);
  amrex::Real moments[NUM_SOOT_MOMENTS + 1] = {0.0};
  if (PeleLM::do_soot_solve) {
    SootData* const sd = PeleLM::soot_model->getSootData();
    sd->initialSmallMomVals(moments);
  }
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; ++n) {
    PeleLM::prob_parm->soot_vals[n] = moments[n];
  }
}
}

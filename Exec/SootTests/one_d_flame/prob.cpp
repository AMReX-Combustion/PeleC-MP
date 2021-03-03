#include "prob.H"

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

void
read_pmf(const std::string& myfile)
{
  std::string firstline;
  std::string secondline;
  std::string remaininglines;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  amrex::Vector<std::string> pmf_names;
  pmf_names.resize(variable_count);
  amrex::Vector<std::string> spec_names;
  EOS::speciesNames(spec_names);
  int tempCol = -1;
  int rhoCol = -1;
  int velCol = -1;
  int xCol = -1;
  pos1 = 0;
  // pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    if (i == 0) xCol = i;
    else if (pmf_names[i] == "U" || pmf_names[i] == "u") velCol = i - 1;
    else if (pmf_names[i] == "T" || pmf_names[i] == "temp") tempCol = i - 1;
    else if (pmf_names[i] == "rho") rhoCol = i - 1;
    else if (pmf_names[i] != spec_names[i - 4]) amrex::Abort("Variables do not match");
    pos1 = pos2 + 1;
  }

  if (rhoCol < 0 || tempCol < 0 || velCol < 0) amrex::Abort("rho, T, and U were not all found");
  int specCol = std::max(tempCol, std::max(rhoCol, velCol)) + 1;
  PeleC::prob_parm_device->tempCol = tempCol;
  PeleC::prob_parm_device->rhoCol = rhoCol;
  PeleC::prob_parm_device->velCol = velCol;
  PeleC::prob_parm_device->specCol = specCol;
  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  PeleC::prob_parm_device->pmf_N = line_count;
  PeleC::prob_parm_device->pmf_M = variable_count - 1;
  PeleC::prob_parm_host->h_pmf_X.resize(PeleC::prob_parm_device->pmf_N);
  PeleC::prob_parm_host->pmf_X.resize(PeleC::prob_parm_device->pmf_N);
  PeleC::prob_parm_host->h_pmf_Y.resize(
    PeleC::prob_parm_device->pmf_N * PeleC::prob_parm_device->pmf_M);
  PeleC::prob_parm_host->pmf_Y.resize(
    PeleC::prob_parm_device->pmf_N * PeleC::prob_parm_device->pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  for (unsigned int i = 0; i < PeleC::prob_parm_device->pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->h_pmf_X[i];
    for (unsigned int j = 0; j < PeleC::prob_parm_device->pmf_M; j++) {
      sinput >>
        PeleC::prob_parm_host->h_pmf_Y[j * PeleC::prob_parm_device->pmf_N + i];
    }
  }

  // Renormalize the mass fractions so they sum to 1
  const int N = PeleC::prob_parm_device->pmf_N;
  for (unsigned int i = 0; i < PeleC::prob_parm_device->pmf_N; i++) {
    amrex::Real sumY = 0.;
    for (int n = 0; n < NUM_SPECIES; n++) {
      const int col = specCol + n;
      amrex::Real Yval = PeleC::prob_parm_host->h_pmf_Y[N * col + i];
      Yval = std::min(1., std::max(Yval, 0.));
      sumY += Yval;
    }
    sumY = 1. / sumY;
    for (int n = 0; n < NUM_SPECIES; n++) {
      const int col = specCol + n;
      amrex::Real Yval = PeleC::prob_parm_host->h_pmf_Y[N * col + i];
      Yval = std::min(1., std::max(Yval, 0.));
      PeleC::prob_parm_host->h_pmf_Y[N * col + i] = Yval*sumY;
    }
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_X.begin(),
    PeleC::prob_parm_host->h_pmf_X.end(), PeleC::prob_parm_host->pmf_X.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_Y.begin(),
    PeleC::prob_parm_host->h_pmf_Y.end(), PeleC::prob_parm_host->pmf_Y.begin());
  PeleC::prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
}

void
init_bc()
{
  amrex::Real vt;
  amrex::Real ek;
  amrex::Real T;
  amrex::Real rho;
  amrex::Real e;
  amrex::Real massfrac[NUM_SPECIES];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  const int tempCol = PeleC::prob_parm_device->tempCol;
  const int velCol = PeleC::prob_parm_device->velCol;
  const int rhoCol = PeleC::prob_parm_device->rhoCol;
  const int specCol = PeleC::prob_parm_device->specCol;
  const amrex::Real yl = 0.0;
  const amrex::Real yr = 0.0;
  pmf(yl, yr, pmf_vals, *PeleC::prob_parm_device);
  amrex::Real mysum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = amrex::max<amrex::Real>(0.0, pmf_vals[specCol + n]);
    mysum += massfrac[n];
  }
  massfrac[N2_ID] = 1.0 - (mysum - massfrac[N2_ID]);
  T = pmf_vals[tempCol];
  PeleC::prob_parm_device->vn_in = pmf_vals[velCol];
  const amrex::Real p = PeleC::prob_parm_device->pamb;

  EOS::PYT2RE(p, massfrac, T, rho, e);

  vt = PeleC::prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  PeleC::prob_parm_device->fuel_state[URHO] = rho;
  PeleC::prob_parm_device->fuel_state[UMX] = rho * vt;
  PeleC::prob_parm_device->fuel_state[UMY] = 0.0;
  PeleC::prob_parm_device->fuel_state[UMZ] = 0.0;
  PeleC::prob_parm_device->fuel_state[UEINT] = rho * e;
  PeleC::prob_parm_device->fuel_state[UEDEN] = rho * (e + ek);
  PeleC::prob_parm_device->fuel_state[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    PeleC::prob_parm_device->fuel_state[UFS + n - 1] = rho * massfrac[n];
  }
}

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  std::string pmf_datafile;

  amrex::ParmParse pp("prob");
  pp.get("pamb", PeleC::prob_parm_device->pamb);
  pp.query("phi_in", PeleC::prob_parm_device->phi_in);
  pp.query("T_in", PeleC::prob_parm_device->T_in);
  pp.query("vn_in", PeleC::prob_parm_device->vn_in);
  pp.get("pmf_datafile", pmf_datafile);

  AMREX_D_TERM(
  PeleC::prob_parm_device->L[0] = probhi[0] - problo[0];,
  PeleC::prob_parm_device->L[1] = probhi[1] - problo[1];,
  PeleC::prob_parm_device->L[2] = probhi[2] - problo[2];);

  read_pmf(pmf_datafile);

  init_bc();
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
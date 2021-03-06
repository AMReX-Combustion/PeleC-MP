#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  amrex::Real Stmod = 5.;
  amrex::Real rhoRatio = 1000.;
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("reynolds", PeleC::h_prob_parm_device->reynolds);
  pp.query("mach", PeleC::h_prob_parm_device->mach);
  pp.query("convecting", PeleC::h_prob_parm_device->convecting);
  pp.query("init_redist", PeleC::prob_parm_host->numRedist);
  pp.query("ref_p", PeleC::h_prob_parm_device->p0);
  pp.query("ref_T", PeleC::h_prob_parm_device->T0);
  pp.query("st_mod", Stmod);
  pp.query("num_particles", PeleC::prob_parm_host->partNum);
  pp.query("density_ratio", rhoRatio);

  // Define the length scale
  PeleC::h_prob_parm_device->L = probhi[0] - problo[0];

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p0, massfrac, PeleC::h_prob_parm_device->T0,
    PeleC::h_prob_parm_device->rho0, eint);
  eos.RTY2Cs(
    PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac, cs);
  eos.TY2Cp(PeleC::h_prob_parm_device->T0, massfrac, cp);

  amrex::Real refL = PeleC::h_prob_parm_device->L;
  PeleC::h_prob_parm_device->v0 = PeleC::h_prob_parm_device->mach * cs;
  pele::physics::transport::TransParm trans_parm;

  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  const amrex::Real mu = PeleC::h_prob_parm_device->rho0 *
                         PeleC::h_prob_parm_device->v0 * refL /
                         PeleC::h_prob_parm_device->reynolds;
  trans_parm.const_viscosity = mu;
  trans_parm.const_conductivity = mu * cp / PeleC::h_prob_parm_device->prandtl;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif

  const amrex::Real St_num = Stmod / (8. * M_PI);
  amrex::Real refU = PeleC::h_prob_parm_device->v0;
  amrex::Real partRho;
  amrex::ParmParse ppp("particles");
  ppp.get("fuel_rho", partRho);
  amrex::Real refRho = PeleC::h_prob_parm_device->rho0;
  if (std::abs(rhoRatio - partRho / refRho) > 10.) {
    amrex::Print() << "Restart solution with particles.fuel_rho = "
                   << refRho * rhoRatio << std::endl;
    amrex::Abort();
  }
  // Time scale for Eulerian phase
  amrex::Real tau_g = refL / refU;
  amrex::Real tau_d = St_num * tau_g;
  amrex::Real dia = std::sqrt(18. * mu * tau_d / partRho);
  amrex::Real Re_d = dia * refU * PeleC::h_prob_parm_device->rho0 / mu;
  amrex::Real error = 1000.;
  const amrex::Real tol = 1.E-6;
  int k = 0;
  const int maxIter = 500;
  // Iteratively solve for the diameter of the spray droplets
  while (error > tol) {
    amrex::Real oldDia = dia;
    Re_d = dia * refU * PeleC::h_prob_parm_device->rho0 / mu;
    amrex::Real C_D = 24. / Re_d;
    if (Re_d > 1.)
      C_D *= (1. + std::pow(Re_d, 2. / 3.) / 6.);
    dia = 0.75 * refRho * C_D * refU * tau_d / partRho;
    error = std::abs(oldDia - dia) / dia;
    k++;
    if (k > maxIter) {
      amrex::Abort("Failed to converge a particle diameter");
    }
  }
  PeleC::prob_parm_host->partDia = dia;
  PeleC::prob_parm_host->partTemp = PeleC::h_prob_parm_device->T0;
  // Output IC
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "number of particles: "
                      << PeleC::prob_parm_host->partNum[0];
    for (int dir = 1; dir != AMREX_SPACEDIM; ++dir)
      amrex::Print(ofs) << ", " << PeleC::prob_parm_host->partNum[dir];
    amrex::Print(ofs) << std::endl;
    amrex::Print(ofs) << "rho0: " << PeleC::h_prob_parm_device->rho0 << std::endl;
    amrex::Print(ofs) << "cs: " << cs << std::endl;
    amrex::Print(ofs) << "U: " << PeleC::h_prob_parm_device->v0 << std::endl;
    amrex::Print(ofs) << "mu: " << trans_parm.const_viscosity << std::endl;
    amrex::Print(ofs) << "Re: " << PeleC::h_prob_parm_device->reynolds
                      << std::endl;
    amrex::Print(ofs) << "Stokes number: " << Stmod << "*Stc " << std::endl;
    amrex::Print(ofs) << "particle diameter: " << PeleC::prob_parm_host->partDia
                      << std::endl;
    amrex::Print(ofs) << "tau_d: " << tau_d << std::endl;
    amrex::Print(ofs) << "Re_d: " << Re_d << std::endl;
    amrex::Print(ofs) << "final time (72 tau_g): " << 72. * tau_g << std::endl;
    ofs.close();
  }
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


// AMReX include statements
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

// PeleC include statements
#include <PeleC.H>
#include <PeleC_F.H>

// PeleC-Physics include statements
#include <chemistry_file.H>

// PeleC-MP include statements
#include <SootModel.H>

using namespace amrex;

// This file will contains functions and initializations relevant to 
// the soot surface reactions

void
SootModel::initializeReactData()
{
  if (m_sootVerbosity)
    {
      Print() << "SootModel::initReactData(): Initializing reaction"
	      << " and species data for soot model" << std::endl;
    }
  // Retrieve the universal gas constant
  Real ruc, pa; // Unused, only needed for chemistry inputs
  CKRP(&Rgas, &ruc, &pa);
  // Number of gas phase species associated with the surface reactions
  const int ngs = GasSpecIndx::numGasSpecs;
  // Relevant species names for the surface chemistry model
  // TODO: Currently must correspond to GasSpecIndx enum in SootModel.H
  m_gasSpecNames = {"H2", "H", "OH", "H2O", "CO", "C2H2", "O2", m_PAHname};
  // Molecular weight of relevant surface reaction species
  m_gasMW.resize(ngs);
  // Number of accounted for PAH particles (hard-coded)
  const int numPAH = 3;
  // Names of PAH species
  // TODO: Currently can only handle naphthalene(C10H8), phenathrene(C14H10), 
  // and pyrene(C16H10) and can only handle 1 PAH inception species
  std::string PAH_names[] = {"A2", "A3", "A4"};
  // Corresponding sticking coefficients
  const Real PAH_gammas[] = {0.001, 0.0150, 0.025};
  // Average number of C atoms per dimer
  const Real PAH_numC[] = {10., 14., 16.};

  // Determine which PAH inception species is being used
  bool corrPAH = false;
  for (int cpah = 0; cpah != numPAH; ++cpah)
    {
      if (PAH_names[cpah] == m_PAHname)
	{
	  m_gammaStick = PAH_gammas[cpah];
	  m_dimerVol = 2.*PAH_numC[cpah];
	  corrPAH = true;
	}
    }
  m_nuclVol = 2.*m_dimerVol;
  m_nuclSurf = std::pow(m_nuclVol, 2./3.);
  if (!corrPAH)
    {
      Abort(m_PAHname + " not recognized as PAH inception species");
    }
  // Retrieve the molecular weights
  Vector<Real> species_mm(PeleC::NumSpec);
  get_mw(species_mm.dataPtr());
  Vector<std::string> species_names(PeleC::NumSpec);
  // Loop over all species
  for (int i = 0; i != PeleC::NumSpec; ++i)
    {
      int len = 20;
      Vector<int> int_spec_names(len);
      get_spec_names(int_spec_names.dataPtr(), &i, &len);
      char char_spec_names[len+1];
      for (int j = 0; j != len; ++j)
	{
	  char_spec_names[j] = int_spec_names[j];
	}
      char_spec_names[len] = '\0';
      species_names[i] = std::string(char_spec_names);
      // Check if species is the PAH inceptor
      if (species_names[i] == m_PAHname)
	{
	  m_PAHindx = i;
	}
      // Check if species matches others for surface reactions
      for (int sootSpec = 0; sootSpec != ngs; ++sootSpec)
	{
	  if (species_names[i] == m_gasSpecNames[sootSpec])
	    {
	      m_gasSpecRefs[sootSpec] = i;
	      m_gasMW[sootSpec] = species_mm[i];
	    }
	}
    }
  // Return error if no PAH inception species is specified in mechanism
  if (m_PAHindx == -1)
    {
      Abort("PAH inception species was not found in PelePhysics mechanism");
    }
  // Return error if not all soot species are present
  for (int sootSpec = 0; sootSpec != ngs; ++sootSpec)
    {
      if (m_gasSpecRefs[sootSpec] == -1)
	{
	  Abort("Species " + m_gasSpecNames[sootSpec] + 
		" must be present in PelePhysics mechanism");
	}
    }
  BL_ASSERT(m_gasSpecNames.size() == ngs);
  // Surface reaction information
  const int nsr = 7; // TODO: Should be variable
  m_numSurfReacts = nsr;
  A_f.assign(nsr, 0.);
  n_f.assign(nsr, 0.);
  ER_f.assign(nsr, 0.);
  A_b.assign(nsr, 0.);
  n_b.assign(nsr, 0.);
  ER_b.assign(nsr, 0.);
  // Number of reactants and products for each reaction
  rNum.assign(nsr, 0);
  pNum.assign(nsr, 0);
  // Corresponding species index for reactants and products
  // For example: nIndx_f[8][1] would be the GasSpeciesIndx of
  // the 2nd gaseous reactant in the 9th reaction
  nIndx_f.assign(nsr, {-1, -1, -1});
  nIndx_b.assign(nsr, {-1, -1, -1});
  // Represent the stoichiometric coefficient corresponding to
  // the nIndx.
  nu_f.assign(nsr, {0., 0., 0.});
  nu_b.assign(nsr, {0., 0., 0.});

  // Soot reaction indices
  sIndx_f.assign(nsr, -1);
  sIndx_b.assign(nsr, -1);
}

// Fill reaction data
// TODO: Currently hard-coded
void
SootModel::fillReactionData()
{  
  // Units are CGS
  // Demonstration of reaction indices using a fake reaction
  // Soot-H + 2OH + C2H2 <=> Soot-* + 4H + 2CO
  // For the i-th reaction
  // rNum[i] = 2; // 2 gas phase reactants
  // nIndx_f[i][0] = GasSpecIndx::indxOH;
  // nu_f[i][0] = 2.; // Since there are 2 moles of OH
  // nIndx_f[i][1] = GasSpecIndx::indxC2H4;
  // nu_f[i][1] = 1; // Since there is 1 mole of C2H2
  // sIndx_f[i] = SootIndx::indxSootH;

  // pNum[i] = 2; // 2 gas phase products
  // nIndx_b[i][0] = GasSpecIndx::indxH;
  // nu_b[i][0] = 4.; // Since there are 4 moles of H
  // nIndx_b[i][1] = GasSpecIndx::indxCO;
  // nu_b[i][1] = 2.; // Since there are 2 moles of CO
  // sIndx_b[i] = SootIndx::indxSootS;

  // 1. Soot-H + OH <=> Soot-* + H2O
  A_f[0] = 6.72E1;
  n_f[0] = 3.33;
  ER_f[0] = 6.09E10/Rgas;
  rNum[0] = 1;
  nIndx_f[0][0] = GasSpecIndx::indxOH;
  nu_f[0][0] = 1.;
  sIndx_f[0] = SootIndx::indxSootH;

  A_b[0] = 6.44E-1;
  n_b[0] = 3.79;
  ER_b[0] = 27.96E10/Rgas;
  pNum[0] = 1;
  nIndx_b[0][0] = GasSpecIndx::indxH2O;
  nu_b[0][0] = 1.;
  sIndx_b[0] = SootIndx::indxSootS;

  // 2. Soot-H + H <=> Soot-* + H2
  A_f[1] = 1.0E8;
  n_f[1] = 1.80;
  ER_f[1] = 68.42E10/Rgas;
  rNum[1] = 1;
  nIndx_f[1][0] = GasSpecIndx::indxH;
  nu_f[1][0] = 1.;
  sIndx_f[1] = SootIndx::indxSootH;

  A_b[1] = 8.68E4;
  n_b[1] = 2.36;
  ER_b[1] = 25.46E10/Rgas;
  pNum[1] = 1;
  nIndx_b[1][0] = GasSpecIndx::indxH2;
  nu_b[1][0] = 1.;
  sIndx_b[1] = SootIndx::indxSootS;

  // 3. Soot-H <=> Soot-* + H
  A_f[2] = 1.13E16;
  n_f[2] = -0.06;
  ER_f[2] = 476.05E10/Rgas;
  rNum[2] = 0;
  sIndx_f[2] = SootIndx::indxSootH;

  A_b[2] = 4.17E13;
  n_b[2] = 0.15;
  ER_b[2] = 0.;
  pNum[2] = 1;
  nIndx_b[2][0] = GasSpecIndx::indxH;
  nu_b[2][0] = 1.;
  sIndx_b[2] = SootIndx::indxSootS;

  // 4. Soot-* + C2H2 => Soot-H
  A_f[3] = 2.52E9;
  n_f[3] = 1.10;
  ER_f[3] = 17.13E10/Rgas;
  rNum[3] = 1;
  nIndx_f[3][0] = GasSpecIndx::indxC2H2;
  nu_f[3][0] = 1.;
  sIndx_f[3] = SootIndx::indxSootS;

  sIndx_b[3] = SootIndx::indxSootH;

  // 5. Soot-* + O2 => Soot-* + 2CO
  A_f[4] = 2.20E12;
  n_f[4] = 0.;
  ER_f[4] = 31.38E10/Rgas;
  rNum[4] = 1;
  nIndx_f[4][0] = GasSpecIndx::indxO2;
  nu_f[4][0] = 1.;
  sIndx_f[4] = SootIndx::indxSootS;

  pNum[4] = 1;
  nIndx_b[4][0] = GasSpecIndx::indxCO;
  nu_b[4][0] = 2.;
  sIndx_b[4] = SootIndx::indxSootS;

  // 6. Soot-H + OH => Soot-H + CO
  // TODO: This transforms the Arrhenius formulation to be
  // reaction probability, 8.94*sqrt(T)*probGamma*A
  Real probGamma = 0.13;
  // FIXME: Find out what the units for this are
  A_f[5] = 8.94*probGamma*avogadros*100.;
  n_f[5] = 0.5;
  ER_f[5] = 0.;
  rNum[5] = 1;
  nIndx_f[5][0] = GasSpecIndx::indxOH;
  nu_f[5][0] = 1.;
  sIndx_f[5] = SootIndx::indxSootH;

  A_b[5] = 0.;
  n_b[5] = 0.;
  ER_b[5] = 0.;
  pNum[5] = 1;
  nIndx_b[5][0] = GasSpecIndx::indxCO;
  nu_b[5][0] = 1.;
  sIndx_b[5] = SootIndx::indxSootH;

  // 7. A4 + A4 => DIMER
  // TODO: Makes use of Arrhenius form similar to last reaction
  A_f[6] = m_betaDimerFact*std::sqrt(m_colFact)*0.5*m_gammaStick;
  n_f[6] = 0.5;
  ER_f[6] = 0.;
  rNum[6] = 1;
  nIndx_f[6][0] = GasSpecIndx::indxPAH;
  nu_f[6][0] = 2.;
  sIndx_f[6] = -1; // No soot in reactants

  m_reactDataFilled = true;
}

// Compute the surface and gas phase chemistry rates
void
SootModel::chemicalSrc(const Real&         T,
		       const Vector<Real>& xi_n,
		       const Vector<Real>& moments,
		       const Vector<Real>& momFV,
		       Real&               k_sg,
		       Real&               k_ox,
		       Real&               k_o2,
		       Vector<Real>&       omega_src,
		       Real&               rho_src,
		       Real&               eng_src) const
{
  // Number of surface reactions
  const int nsr = m_numSurfReacts;
  // Number of gas species
  const int ngs = GasSpecIndx::numGasSpecs;
  // Number of soot types
  const int nss = SootIndx::numSootSpecs;
  // Forwad and back reaction rates
  Vector<Real> k_fwd(nsr);
  Vector<Real> k_bkwd(nsr);
  // Creation and destruction terms
  Vector<Real> w_fwd(nsr, 0.);
  Vector<Real> w_bkwd(nsr, 0.);
  const Real invT = 1./T;
  // Loop over reactions
  for (int i = 0; i != nsr; ++i)
    {
      k_fwd[i] = A_f[i]*std::pow(T, n_f[i])*std::exp(-ER_f[i]*invT);
      k_bkwd[i] = A_b[i]*std::pow(T, n_b[i])*std::exp(-ER_b[i]*invT);
      Real fwdM = 1.;
      for (int j = 0; j != rNum[i]; ++j)
	{
	  const int rIndx = nIndx_f[i][j];
	  const Real rExp = nu_f[i][j];
	  fwdM *= std::pow(xi_n[rIndx], rExp);
	}
      w_fwd[i] = k_fwd[i]*fwdM;
      Real bkwdM = 1.;
      for (int j = 0; j != pNum[i]; ++j)
	{
	  const int pIndx = nIndx_b[i][j];
	  const Real pExp = nu_b[i][j];
	  bkwdM *= std::pow(xi_n[pIndx], pExp);
	}
      w_bkwd[i] = k_bkwd[i]*bkwdM;
    }
  // Factor r for the quasi-steady state concentration of radical sites
  // TODO: This will depend on the surface reactions, currently hardcoded
  Real fSootStar = computeRadSiteConc(xi_n, k_fwd, k_bkwd);
  computeSurfRates(xi_n, w_fwd, w_bkwd, fSootStar, k_sg, k_ox, k_o2);
  // Determine the concentration of hydrogenated and radical soot surface sites
  // Quasi-steady state for surface radical sites on soot
  Vector<Real> C_Soot(nss);
  C_Soot[SootIndx::indxSootS] = fSootStar*m_SootDensityC*moments[2];
  C_Soot[SootIndx::indxSootH] = (1. - fSootStar)*m_SootDensityC*moments[2];

  // TODO: The last two reaction (6 and 7) are special 
  // they are not treated the same in this loop
  // TODO: Currently assumes soot is involved in every reaction on both sides
  for (int i = 0; i != nsr - 2; ++i)
    {
      int sootIndx = sIndx_f[i];
      BL_ASSERT(sootIndx >= 0);
      Real xi_soot = C_Soot[sootIndx];
      w_fwd[i] *= xi_soot;
      sootIndx = sIndx_b[i];
      BL_ASSERT(sootIndx >= 0);
      xi_soot = C_Soot[sootIndx];
      w_bkwd[i] *= xi_soot;
    }
  Real surf = S0*fracMom(0., 1., momFV);
  w_fwd[5] *= surf;
  // Loop over gas species and solve for omegas
  for (int i = 0; i != nsr; ++i)
    {
      // Creation/destruction rates for the i-th reaction
      Real rate = w_fwd[i] - w_bkwd[i];
      // Loop over reactant species and subtract nu'_(i,r)*rate
      for (int r = 0; r != rNum[i]; ++r)
	{
	  const int rIndx = nIndx_f[i][r];
	  Real nup = nu_f[i][r];
	  omega_src[rIndx] -= nup*rate;
	}
      // Loop over product species and add nu"_(i,r)*rate
      for (int p = 0; p != pNum[i]; ++p)
	{
	  const int pIndx = nIndx_b[i][p];
	  Real nupp = nu_b[i][p];
	  omega_src[pIndx] += nupp*rate;
	}
    }
  // Loop over each species to convert to proper units
  // and to compute continuity source term
  rho_src = 0.;
  eng_src = 0.; // Unused since energy already contains sensible and chemical
  for (int i = 0; i != ngs; ++i)
    {
      const Real cmm = m_gasMW[i];
      omega_src[i] *= cmm;
      rho_src += omega_src[i];
    }
}

// Return fSootStar
Real
SootModel::computeRadSiteConc(const Vector<Real>& xi_n,
			      const Vector<Real>& k_fwd,
			      const Vector<Real>& k_bkwd) const
{
  Real C_OH = xi_n[GasSpecIndx::indxOH];
  Real C_H = xi_n[GasSpecIndx::indxH];
  Real C_H2 = xi_n[GasSpecIndx::indxH2];
  Real C_H2O = xi_n[GasSpecIndx::indxH2O];
  Real C_C2H2 = xi_n[GasSpecIndx::indxC2H2];
  Real r1 = (k_fwd[0]*C_OH + k_fwd[1]*C_H + k_fwd[2]);
  Real r2 = (k_bkwd[0]*C_H2O + k_bkwd[1]*C_H2 + k_bkwd[2]*C_H 
	       + k_fwd[3]*C_C2H2);
  return r1/(r2 + r1);
}

// Compute the surface chemistry rates (1/s)
void
SootModel::computeSurfRates(const Vector<Real>& xi_n,
			    const Vector<Real>& w_fwd,
			    const Vector<Real>& w_bkwd,
			    const Real&         fSootStar,
			    Real&               k_sg,
			    Real&               k_ox,
			    Real&               k_o2) const
{
  k_sg = w_fwd[3]*fSootStar;
  k_ox = w_fwd[4]*fSootStar + w_fwd[5]*0.5/m_SootChi;
  k_o2 = w_fwd[4]*fSootStar;
}

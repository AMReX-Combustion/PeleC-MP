
// Standard library includes
#include <string>
#include <map>

// AMReX include statements
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_Derive.H>

// PeleC include statements
#include <PeleC.H>
#include <PeleC_F.H>

// PeleC-MP include statements
#include <SootModel.H>
#include <SootModel_derive.H>

using namespace amrex;

// Box definitions for adding plot outputs
//static Box same_box (const Box& b) { return b; }
// Read soot related inputs
void
SootModel::readSootParams()
{
  ParmParse pp("soot");
  pp.get("incept_pah", m_PAHname);
  pp.query("v", m_sootVerbosity);
  if (m_sootVerbosity)
    {
      Print() << "SootModel::readSootParams(): Reading input parameters"
	      << std::endl;
    }
  m_readSootParams = true;
}

// Define all member data
void
SootModel::defineMemberData()
{
  if (m_sootVerbosity)
    {
      Print() << "SootModel::defineMemberData(): Defining member data"
	      << std::endl;
    }
  // Assign the many constant values
  const int numSootMom = NUM_SOOT_MOMENTS;
  // The volume and surface area for the smallest soot particles
  m_V0 = m_SootMolarMass/(avogadros*m_SootDensity); // cm^3
  m_S0 = std::pow(36.*M_PI, 1./3.)*std::pow(m_V0, 2./3.); // cm^2

  // Soot density (unitless)
  m_SootDensityC = m_SootChi*m_S0;

  // Factors used throughout
  m_colFact = M_PI*Rgas/(2.*avogadros*m_SootDensity);
  m_colFact23 = std::pow(m_V0, 2./3.);
  m_colFact16 = std::pow(m_V0, 1./6.); // Units: cm^0.5
  m_colFactPi23 = std::pow(6./M_PI, 2./3.);

  // Compute V_nucl and V_dimer to fractional powers
  for (int i = 0; i != 9; ++i)
    {
      Real exponent = 2.*(Real)i - 3.;
      m_dimerExp6[i] = std::pow(m_dimerVol, exponent/6.);
    }
  for (int i = 0; i != 11; ++i)
    {
      Real exponent = (Real)i - 3.;
      m_nuclVolExp3[i] = std::pow(m_nuclVol, exponent/3.);
      m_nuclVolExp6[i] = std::pow(m_nuclVol, exponent/6.);
    }

  for (int i = 0; i != numSootMom; ++i)
    {
      // Used to convert moments to mol of C
      MomUnitConv[i] = std::pow(m_V0, MomOrderV[i])*
	std::pow(m_S0, MomOrderS[i])*avogadros;
      // Used for computing nucleation source term
      m_momFact[i] = std::pow(m_nuclVol, MomOrderV[i])*
	std::pow(m_nuclSurf, MomOrderS[i]);
    }
  // and to convert the weight of the delta function
  MomUnitConv[numSootMom] = avogadros;

  // Coagulation, oxidation, and fragmentation factors
  m_lambdaCoagFact = 
    1./(std::pow(6.*m_SootMolarMass/(M_PI*m_SootDensity*avogadros), 1./3.));
  for (int i = 0; i != numSootMom; ++i)
    {
      const Real expFact = MomOrderV[i] + 2./3.*MomOrderS[i];
      const Real factor = (std::pow(2., expFact - 1.) - 1.);
      m_ssfmCoagFact[i] = 
	std::pow(2., 2.5)*factor*std::pow(m_nuclVol, expFact + 1./6.);
      m_sscnCoagFact[i] = factor*std::pow(m_nuclVol, expFact);
      m_smallOxidFact[i] = std::pow(m_nuclVol, expFact - 1./3.);
      m_fragFact[i] = std::pow(2., 1. - MomOrderV[i] - MomOrderS[i]) - 1.;
    }
  m_ssfmCoagFact[numSootMom] = std::pow(2., 2.5)*std::pow(1./m_nuclVol, 0.5)*
    getNuclExp3(2);
  m_smallOxidFact[numSootMom] = getNuclExp3(-1);

  // Beta and dimer factors
  // cm^0.5 mol^-1
  const Real dnfact = 4.*std::sqrt(2.)*m_colFact16*m_colFactPi23*avogadros;
  m_betaDimerFact = dnfact*std::pow(0.5*m_dimerVol, 1./6.);
  m_betaNuclFact = 2.2*dnfact*getDimerExp6(1);

  /// Condensation factor
  m_condFact = std::sqrt((1./m_nuclVol) + (1./m_dimerVol))*
    std::pow((std::pow(m_nuclVol, 1./3.) + std::pow(m_dimerVol, 1./3.)), 2.);
}

// Add derive plot variables
void
SootModel::addSootDerivePlotVars(DeriveList&           derive_lst,
				 const DescriptorList& desc_lst)
{
  
  // Number density associated with the second mode (large particles)
  derive_lst.add("NL",IndexType::TheCellType(),1,soot_largepartnumdens,DeriveRec::TheSameBox);
  derive_lst.addComponent("NL",desc_lst,0,PeleC::Density,PeleC::NUM_STATE);

  // Mean volume associated with the second mode (large particles)
  derive_lst.add("VL",IndexType::TheCellType(),1,soot_largepartmeanvol,DeriveRec::TheSameBox);
  derive_lst.addComponent("VL",desc_lst,0,PeleC::Density,PeleC::NUM_STATE);

  // Mean surface area associated with the second mode (large particles)
  derive_lst.add("SL",IndexType::TheCellType(),1,soot_largepartsurfarea,DeriveRec::TheSameBox);
  derive_lst.addComponent("SL",desc_lst,0,PeleC::Density,PeleC::NUM_STATE);
}

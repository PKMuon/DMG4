#include "AnnihilationStepLimiter.hh"
#include "G4TransportationProcessType.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"

#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "DarkMatterParametersRegistry.hh"

AnnihilationStepLimiter::AnnihilationStepLimiter(DarkMatterAnnihilation* m_DarkMatterAnnihilation_in,const G4String& aName)
  : G4VDiscreteProcess(aName,fUserDefined),
    m_DarkMatterAnnihilation(m_DarkMatterAnnihilation_in),
    AnnihilationMaxStep(DBL_MAX),
    factor(1.)
{
  SetProcessSubType(2);
}


AnnihilationStepLimiter::~AnnihilationStepLimiter()
{}


G4bool AnnihilationStepLimiter::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle==*(G4Positron::Definition()));
}


G4double AnnihilationStepLimiter::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*){
  return DBL_MAX;
}

G4double 
  AnnihilationStepLimiter::PostStepGetPhysicalInteractionLength(
                                       const G4Track& aTrack,
                                       G4double, // previousStepSize
                                       G4ForceCondition* condition  )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  G4double DensityMat = aTrack.GetMaterial()->GetDensity() / (g / cm3);
  G4double etot = aTrack.GetTotalEnergy() / GeV; //this is the total energy of the positron at the beginning of the step

  if (m_DarkMatterAnnihilation->EmissionAllowed(etot, DensityMat)) {

    //G4 convention: put immediately units
    G4double Mres=m_DarkMatterAnnihilation->GetMA()*GeV;

    //Track energy in G4 units
    G4double E=aTrack.GetTotalEnergy();
    G4double Ekin=aTrack.GetKineticEnergy();

    //compute dEdX and delta0
    G4double dEdX=emCal.GetDEDX(Ekin,aTrack.GetParticleDefinition(),aTrack.GetMaterial());
    G4double delta0=1./dEdX;
    delta0 = delta0 / factor;

    //compute the step length (see:https://gitlab.cern.ch/P348/DMG4/-/work_items/20)
    G4double maxStepDen=dSigmadEoverSigma(E);
    AnnihilationMaxStep=delta0/maxStepDen;
    //std::cout << "Annihilation max step: " << AnnihilationMaxStep << std::endl;

    return AnnihilationMaxStep;
  }
  return DBL_MAX;
}


G4VParticleChange*
  AnnihilationStepLimiter::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&  )
// Do Nothing
//
{
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

G4double AnnihilationStepLimiter::GetMaxEloss(G4double E)
{
  //std::cout << "Max E loss: " << 1./dSigmadEoverSigma(E) << std::endl;
  return 1./dSigmadEoverSigma(E);
}


//E is the TOTAL POSITRON energy
G4double AnnihilationStepLimiter::dSigmadEoverSigma(G4double E)
{
  G4double Mres=m_DarkMatterAnnihilation->GetMA()*GeV;
  G4double Eres=(Mres*Mres-2*CLHEP::electron_mass_c2*CLHEP::electron_mass_c2)/(2*CLHEP::electron_mass_c2);
  G4double W=m_DarkMatterAnnihilation->Width()*GeV;
  G4double Wstar=W*(Mres/(2*CLHEP::electron_mass_c2));


  G4double ret=1./Wstar;
  if (fabs(E-Eres)>Wstar){
    ret=2*fabs(E-Eres)/((E-Eres)*(E-Eres)+Wstar*Wstar);
  }

  //G4cout<<"dS: "<<Mres/GeV<<" "<<Eres/GeV<<" "<<W/GeV<<" "<<Wstar/GeV<<" "<<E/GeV<<" "<<ret<<G4endl;
  return ret;
}

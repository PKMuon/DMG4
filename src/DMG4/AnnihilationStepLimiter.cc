#include "AnnihilationStepLimiter.hh"
#include "G4TransportationProcessType.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4SystemOfUnits.hh"

#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"

AnnihilationStepLimiter::AnnihilationStepLimiter(DarkMatterAnnihilation* m_DarkMatterAnnihilation_in,const G4String& aName)
  : G4StepLimiter(aName),
    m_DarkMatterAnnihilation(m_DarkMatterAnnihilation_in),
    AnnihilationMaxStep(DBL_MAX)
{
  factor=5;
}


AnnihilationStepLimiter::~AnnihilationStepLimiter()
{}


G4double 
  AnnihilationStepLimiter::PostStepGetPhysicalInteractionLength(
                                       const G4Track& aTrack,
                                       G4double, // previousStepSize
                                       G4ForceCondition* condition  )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  G4double DensityMat = aTrack.GetMaterial()->GetDensity() / (g / cm3);
  G4double ekin = aTrack.GetKineticEnergy() / GeV; //this is the energy of the positron at the beginning of the step

  G4cout<<"Calling AnnihilationStepLimiter::PostStepGetPhysicalInteractionLength "<<ekin<<G4endl;

  if (m_DarkMatterAnnihilation->EmissionAllowed(ekin, DensityMat)) {

    //G4 convention: put immediately units
    G4double Mres=m_DarkMatterAnnihilation->GetMA()*GeV;
    G4double Eres=(Mres*Mres)/CLHEP::electron_mass_c2;
    G4double W=m_DarkMatterAnnihilation->Width()*GeV;
    G4double Wstar=W*(Mres/(2*CLHEP::electron_mass_c2));

    //Track energy in G4 units
    G4double E=aTrack.GetKineticEnergy();

    //compute dEdX and delta0
    G4double dEdX=emCal.ComputeTotalDEDX(aTrack.GetKineticEnergy(),aTrack.GetParticleDefinition(),aTrack.GetMaterial());
    G4double delta0=Wstar/dEdX;
    delta0 = delta0 / factor;

    //compute the step length (see:https://gitlab.cern.ch/P348/DMG4/-/work_items/20)
    G4double maxStepDen=dSigmadEoverSigma(E);
    AnnihilationMaxStep=delta0/maxStepDen;

    G4cout<<"AnnihilationStepLimiter kinE: "<<E/GeV<<" Eres: "<<Eres/GeV<<" Wstar: "<<Wstar/GeV<<" delta0: "<<delta0/cm<<G4endl;

    return AnnihilationMaxStep;
  }
  G4cout<<"DEFAULT BIG STEP: "<<DBL_MAX<<G4endl;
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

G4double AnnihilationStepLimiter::GetMaxEloss(G4double E){
  G4double Mres=m_DarkMatterAnnihilation->GetMA()*GeV;
  G4double W=m_DarkMatterAnnihilation->Width()*GeV;
  G4double Wstar=W*(Mres/(2*CLHEP::electron_mass_c2));

  G4double maxEloss=Wstar/dSigmadEoverSigma(E);
  return maxEloss;
}

G4double AnnihilationStepLimiter::dSigmadEoverSigma(G4double E){
  G4double Mres=m_DarkMatterAnnihilation->GetMA()*GeV;
  G4double Eres=(Mres*Mres)/CLHEP::electron_mass_c2;
  G4double W=m_DarkMatterAnnihilation->Width()*GeV;
  G4double Wstar=W*(Mres/(2*CLHEP::electron_mass_c2));

  G4double ret=1;
  if (fabs(E-Eres)>Wstar){
    ret=Wstar*2*fabs(E-Eres)/((E-Eres)*(E-Eres)+Wstar*Wstar);
  }
  return ret;

}


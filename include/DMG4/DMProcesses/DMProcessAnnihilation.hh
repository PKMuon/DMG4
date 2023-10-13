#pragma once

#include <G4VDiscreteProcess.hh>

class DarkMatter;
class DarkMatterAnnihilation;
class G4ParticleDefinition;
class DarkMatterParametersFactory;
class AnnihilationStepLimiter;

class DMProcessAnnihilation : public G4VDiscreteProcess
{
  public:

    DMProcessAnnihilation(DarkMatterAnnihilation* DarkMatterPointer, G4ParticleDefinition* theDMParticlePtrIn, G4double BiasSigmaFactor);

    // Implements final state parameters when the process won.
    virtual G4VParticleChange* PostStepDoIt( const G4Track &, const G4Step & ) override;

    virtual G4double GetMeanFreePath( const G4Track & aTrack,
                                      G4double previousStepSize,
                                      G4ForceCondition * condition ) override;

    virtual G4bool IsApplicable(const G4ParticleDefinition &) override;

  private:

    DarkMatterAnnihilation* myDarkMatterAnnihilation;
    G4ParticleDefinition* theDMParticlePtr;
    G4double BiasSigmaFactor;

    DarkMatterParametersFactory* DMpar;
    G4int iBranchingType;
    G4double mChi,mChi1,mChi2;

    //This part is used to take into account the energy dependence of the cross section along the step
    G4double CrossSectionStepVal; //the cross section used by the GetMeanFreePath method. Should be larger than the cross section along the step
    G4double CrossSectionStepE;   //the energy used to compute the above
    G4double xi;

    //Part for the step limit
    AnnihilationStepLimiter *m_limiter;

    void SetStepLimiter();
};

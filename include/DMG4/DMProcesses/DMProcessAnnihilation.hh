#pragma once

#include <map>
#include <utility>

#include <G4VDiscreteProcess.hh>

#include <fstream>

class DarkMatter;
class DarkMatterAnnihilation;
class G4ParticleDefinition;
class DarkMatterParametersRegistry;
class AnnihilationStepLimiter;

class G4Element;


class DMProcessAnnihilation : public G4VDiscreteProcess
{
  public:

    DMProcessAnnihilation(DarkMatterAnnihilation* DarkMatterPointer, G4ParticleDefinition* theDMParticlePtrIn, G4double BiasSigmaFactor,AnnihilationStepLimiter *LimiterIn=0);

    // Implements final state parameters when the process won.
    virtual G4VParticleChange* PostStepDoIt( const G4Track &, const G4Step & ) override;

    virtual G4double GetMeanFreePath( const G4Track & aTrack,
                                      G4double previousStepSize,
                                      G4ForceCondition * condition ) override;

    virtual G4bool IsApplicable(const G4ParticleDefinition &) override;

    ~DMProcessAnnihilation();

  private:

    const G4Element *GetRandomElement(const G4ElementVector *elms);

    DarkMatterAnnihilation* myDarkMatterAnnihilation;
    G4ParticleDefinition* theDMParticlePtr;
    G4double BiasSigmaFactor;

    DarkMatterParametersRegistry* DMpar;
    G4int iBranchingType;
    G4double mChi,mChi1,mChi2;

    //This part is used to take into account the energy dependence of the cross section along the step
    G4double CrossSectionStepVal; //the cross section used by the GetMeanFreePath method. Should be larger than the cross section along the step
    G4double CrossSectionStepE;   //the energy used to compute the above
    G4double xi;

    //Part for the step limit
    AnnihilationStepLimiter *m_limiter;


    //Part for the atomic motion
    //Key: Z, payload: for each shell index, a vector of the kinetic energies of the electrons in that shell, extracted by MC, in GeV
    std::map<G4int, std::map<G4int,std::vector<G4double> > > shellElectronEnergies;

    //Key: Z, payload: for each shell index, the number of electrons in that shell
    std::map<G4int, std::map<G4int,G4int> > shellElectronZ;


    //The largest kinetic energy of all atomic electrons generated
    G4double maxShellElectronEnergy;

    G4double GetOneRandomEleEnergy(G4double B);
    std::vector<G4double> SimulateElectronEnergies(const G4Element *elm,G4int is);

};

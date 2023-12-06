#pragma once

#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4EmCalculator.hh"


class DarkMatterAnnihilation;

class AnnihilationStepLimiter : public G4VDiscreteProcess
{
  public:

     AnnihilationStepLimiter(DarkMatterAnnihilation*,const G4String& processName ="AnnihilationStepLimiter");

     virtual ~AnnihilationStepLimiter();

     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

     virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            );
                            

     virtual G4bool IsApplicable(const G4ParticleDefinition&);
     virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);



     void SetMaxStep(G4double MaxStepIn) {AnnihilationMaxStep = MaxStepIn;}
     double GetMaxStep() {return AnnihilationMaxStep;}
     G4double GetMaxEloss(G4double E);

     void SetFactor(double m_factor){
           factor=(m_factor >= 1 ? m_factor : 1.);
           G4cout<<"Annihilation step limiter, factor set to: "<<factor<<G4endl;
     }

  private:
  
  // hide assignment operator as private 
    AnnihilationStepLimiter(AnnihilationStepLimiter&);
    AnnihilationStepLimiter& operator=(const AnnihilationStepLimiter& right);

    G4double dSigmadEoverSigma(G4double E);


    G4double AnnihilationMaxStep;
    DarkMatterAnnihilation* m_DarkMatterAnnihilation;
    G4EmCalculator emCal;
    double factor;



};


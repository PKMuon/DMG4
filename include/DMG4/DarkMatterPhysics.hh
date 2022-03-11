#pragma once

#include <G4VPhysicsConstructor.hh>

class DarkMatter;
class DarkMatterPhysicsMessenger;

class DarkMatterPhysics : public G4VPhysicsConstructor {
  public:

    DarkMatterPhysics(DarkMatterPhysicsMessenger* Msg = nullptr);
    ~DarkMatterPhysics();
    bool DarkMatterPhysicsConfigure();

    //A.C. I introduced this method to allow to pass at run-time the A' mass (in GeV), the mChi/mA mass ratio, and the bias
    //The units for Amass must be GeV!!!
    DarkMatterPhysics(double Amass,double ratio,double alphaD,double Bias);
    bool DarkMatterPhysicsConfigureWithPars(double Amass,double ratio,double alphaD,double Bias);


    // Should call initial constructor of particle singletons
    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    DarkMatter* GetDarkMatterPointer() {return myDarkMatter;}
    G4double GetBiasSigmaFactor() {return BiasSigmaFactor;}

    // function for the messenger
    void SetBiasSigmaFactor0(G4double BiasSigmaFactor0) {fBiasSigmaFactor0 = BiasSigmaFactor0;}
    void SetEThresh(G4double EThresh) {fEThresh = EThresh;}
    void SetDMMass(G4double DMMass) {fDMMass = DMMass;}
    void SetDMParticle(G4String DMParticle) {fDMParticle = DMParticle;}
    void SetSigmaNorm(G4double SigmaNorm) {fSigmaNorm = SigmaNorm;}
    void SetANucl(G4double ANucl) {fANucl = ANucl;}
    void SetZNucl(G4double ZNucl) {fZNucl = ZNucl;}
    void SetDensity(G4double Density) {fDensity = Density;}
    void SetEpsilon(G4double Epsilon) {fEpsilon = Epsilon;}
    void SetDecay(G4int Decay) {fDecay = Decay;}

  private:
    DarkMatterPhysics(const DarkMatterPhysics &) = delete;
    DarkMatterPhysics & operator=(const DarkMatterPhysics &) = delete;
  private:
    DarkMatter* myDarkMatter;
    G4double BiasSigmaFactor;

    DarkMatterPhysicsMessenger* fMessenger; // pointer to the messenger
    G4double fEThresh;
    G4double fBiasSigmaFactor0;
    G4double fDMMass;
    G4String fDMParticle;
    G4double fSigmaNorm;
    G4double fANucl;
    G4double fZNucl;
    G4double fDensity;
    G4double fEpsilon;
    G4int    fDecay;
};

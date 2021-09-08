#pragma once

#include <G4VPhysicsConstructor.hh>

class DarkMatter;

class DarkMatterPhysics : public G4VPhysicsConstructor {
  public:
    DarkMatterPhysics();
    ~DarkMatterPhysics();
    bool DarkMatterPhysicsConfigure();
    // Should call initial constructor of particle singletons
    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    DarkMatter* GetDarkMatterPointer() {return myDarkMatter;}
    G4double GetBiasSigmaFactor() {return BiasSigmaFactor;}
  private:
    DarkMatterPhysics(const DarkMatterPhysics &) = delete;
    DarkMatterPhysics & operator=(const DarkMatterPhysics &) = delete;
  private:
    DarkMatter* myDarkMatter;
    G4double BiasSigmaFactor;
};

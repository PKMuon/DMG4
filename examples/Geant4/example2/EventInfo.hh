#pragma once

#include <G4ThreeVector.hh>
#include <G4Step.hh>
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>


struct ParticleInfo {
  ParticleInfo() : pID(0), E(0), p(0,0,0) {;}
  ParticleInfo(G4double pID_in, G4double E_in, G4ThreeVector p_in) : pID(pID_in), E(E_in), p(p_in) {;}

  G4double pID;
  G4double E; // Saved in GeV
  G4ThreeVector p;  // Saved in GeV/c

  void Reset() {
    // Reset all values
    pID = 0;
    E = 0;
    p.set(0,0,0);
  }

};

std::ofstream& operator<<(std::ofstream& stream, const ParticleInfo& particle);

struct DMEventInfo {
  G4ThreeVector productionVertex; 
  G4ThreeVector decayVertex; 
  ParticleInfo parent;
  std::vector<ParticleInfo> daughters;

  void AddParentInfo(const G4Step* aStep) {
    productionVertex = aStep->GetPreStepPoint()->GetPosition();
    parent.pID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    parent.E = aStep->GetTrack()->GetTotalEnergy()/GeV;
    parent.p = aStep->GetTrack()->GetMomentum()/GeV;
  }

  void AddDaughterInfo(const G4Step* aStep) {
    decayVertex = aStep->GetPreStepPoint()->GetPosition();
    for(const auto& daughter : *aStep->GetSecondaryInCurrentStep()) {
      daughters.emplace_back(
          daughter->GetParticleDefinition()->GetPDGEncoding(),
          daughter->GetTotalEnergy()/GeV,
          daughter->GetMomentum()/GeV
          );
    }
  }

  void DumpInfo(std::ofstream& outFile) {
    if(outFile.is_open()) {
      outFile << productionVertex.x() << " " << productionVertex.y() << " " << productionVertex.z() << std::endl;
      outFile << decayVertex.x() << " " << decayVertex.y() << " " << decayVertex.z() << std::endl;

      // Dump Parent info
      outFile << parent;

      outFile << daughters.size() << std::endl;

      // Dump Daughters info
      for (const auto &daughter : daughters) {
        outFile << daughter;
      }
    }
  }
  
  void Reset() {
    // Reset all values
    productionVertex.set(0,0,0);
    decayVertex.set(0,0,0);
    daughters.clear();
    parent.Reset();
  }
};

class EventInfo
{
  public:
    EventInfo();
    ~EventInfo();

    void AddParentInfo(const G4Step* aStep);
    void AddDaughterInfo(const G4Step* aStep);
    G4bool DumpInfo(std::ofstream& outFile);
    void Reset();

  private:
    std::map<G4int, DMEventInfo> dmEvents;
};

#include "EventInfo.hh"

#include "G4Step.hh"
#include "G4ios.hh"


std::ofstream& operator<<(std::ofstream& stream, const ParticleInfo& particle) {
  // Output format for file and std::cout
  stream << particle.pID << " " << particle.E << " " << particle.p.x() << " " << particle.p.y() << " " << particle.p.z() << std::endl;
  return stream;
}


EventInfo::EventInfo() {
  dmEvents.clear();
}

void EventInfo::AddParentInfo(const G4Step* aStep) {
  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4String particleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4cout << "INFO: EventInfo: Particle (" << particleName << ", " << trackID << ")  added to map!" << G4endl;
  // Add new DM event's parent info
  dmEvents[trackID].AddParentInfo(aStep);
}

void EventInfo::AddDaughterInfo(const G4Step* aStep) {
  G4int trackID = aStep->GetTrack()->GetTrackID();
  if(auto it = dmEvents.find(trackID); it != dmEvents.end()) {
    // Add new DM event's daughter info
    dmEvents[trackID].AddDaughterInfo(aStep);
  } else {
    G4cerr << "ERROR! EventInfo: Parent ID " << trackID << " not present in map!" << G4endl;
  }
}

G4bool EventInfo::DumpInfo(std::ofstream& outFile) {

  // Check file status
  if(outFile.is_open()) {
    outFile << "DMINFO " << dmEvents.size() << std::endl;

    for(auto& dmEvt : dmEvents) {
      dmEvt.second.DumpInfo(outFile);
    }
  } else {
    G4cerr << "ERROR! EventInfo: File not open!" << G4endl;
  }

  return outFile.good();

}

void EventInfo::Reset() {
  dmEvents.clear();
}

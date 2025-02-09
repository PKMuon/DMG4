#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4GenericMessenger;

using namespace std;


class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    ~DetectorConstruction();

  public:
  
    G4VPhysicalVolume* Construct();

    G4Material* GetGeneralAir() {return GeneralAir;}
    G4Material* GetGeneralAluminium() {return GeneralAluminium;}
    G4Material* GetGeneralIron() {return GeneralIron;}
    G4Material* GetGeneralLead() {return GeneralLead;}
    G4Material* GetGeneralVacuum() {return GeneralVacuum;}

    inline G4VPhysicalVolume* GetphysiWorld() {return physiWorld;}
    inline G4LogicalVolume* GetlogicWorld() {return logicWorld;}

    void SetAEmission(G4int AE) {AEmission = AE;}
    G4int GetAEmission() {return AEmission;}

  private:

    void DefineMaterials();

  private:

    G4Material* GeneralAir;
    G4Material* GeneralHelium;
    G4Material* GeneralAluminium;
    G4Material* GeneralIron;
    G4Material* GeneralCopper;
    G4Material* GeneralTungsten;
    G4Material* GeneralLead;
    G4Material* GeneralVacuum;

    G4Box*             solidWorld;
    G4LogicalVolume*   logicWorld;
    G4VPhysicalVolume* physiWorld;
    G4VPhysicalVolume* physiTube;

    G4Box*             solidBox;
    G4Tubs*            solidTube;
    G4LogicalVolume*   logicTube;
    G4LogicalVolume*   logicBox;

    G4int AEmission;
};

#endif

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

#include "G4ios.hh"


DetectorConstruction::DetectorConstruction()
{
}


DetectorConstruction::~DetectorConstruction()
{
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

// Elements: ----------------

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

// Materials: ----------------

  // Air
  density = 1.29*mg/cm3;
  //density = 0.001*mg/cm3;
  GeneralAir = new G4Material(name="GeneralAir", density, nel=2);
  GeneralAir->AddElement(elN, .7);
  GeneralAir->AddElement(elO, .3);

  // Aluminium
  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  GeneralAluminium = new G4Material(name="GeneralAluminium", z=13., a, density);

  // Iron
  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  GeneralIron = new G4Material(name="GeneralIron", z=26., a, density);

  // Lead
  a = 207.19*g/mole;
  density = 11.35*g/cm3;
    GeneralLead = new G4Material(name="GeneralLead", z=82., a, density);

  // Vacuum
  a = 1.*g/mole;
  density = 1.e-15*g/cm3;
  GeneralVacuum = new G4Material(name="GeneralVacuum",z=1.,a,density);

//--------- Sizes of the principal geometrical components (solids)  ---------

//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

  //------------------------------
  // World
  //------------------------------

  G4double WorldWidth = 200.*cm;
  G4double WorldHeight = 200.*cm;
  G4double WorldDepth = 4000.*cm;

  solidWorld= new G4Box("World",0.5*WorldWidth,0.5*WorldHeight,0.5*WorldDepth);
  logicWorld= new G4LogicalVolume( solidWorld, GeneralAir, "WorldLV", 0, 0, 0);

  //  Must place the World Physical volume unrotated at (0,0,0).
  //
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 "WorldPV",       // its name
                                 logicWorld,      // its logical volume
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume

  // Lead block:

  G4LogicalVolume*   logicECAL;
  G4VPhysicalVolume* physiECAL;
  G4double ECALWidth;
  G4double ECALHeight;
  G4double ECALStart;
  G4double ECALX;
  G4double ECALY;
  G4double ECALDepth;

  ECALStart=4600;
  ECALX=0.;
  ECALY=0.;
  ECALWidth=200.;
  ECALHeight=200.;
  ECALDepth=400.;

  solidBox = new G4Box("ECALbox",0.5*ECALWidth, 0.5*ECALHeight, 0.5*ECALDepth);

  logicECAL = new G4LogicalVolume(solidBox, GeneralLead, "ECAL_LV box", 0, 0, 0);

  physiECAL = new G4PVPlacement(0,                             // no rotation
       G4ThreeVector(ECALX, ECALY, ECALStart + 0.5*ECALDepth), // at
                                          "ECAL_DAMI",         // its name
                                           logicECAL,          // its logical volume
                                           physiWorld,         // its mother  volume
                                           false,              // no boolean operations
                                           0);                 // copy number

  return physiWorld;
}

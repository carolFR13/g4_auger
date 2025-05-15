#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction(){


}

DetectorConstruction::~DetectorConstruction(){

    
}

G4VPhysicalVolume *DetectorConstruction::Construct(){ // we are defining here our main function Construct

    G4bool checkOverlaps = true;

    G4NistManager *nist = G4NistManager::Instance();
    G4Material *worldMat = nist-> FindOrBuildMaterial("G4_AIR");
    G4Material*leadMat = nist -> FindOrBuildMaterial("G4_Pb");
    G4Material *detMat = nist -> FindOrBuildMaterial("G4_SODIUM_IODIDE");

    G4double xWorld = 5. * m;
    G4double yWorld = 5. * m;
    G4double zWorld = 5. * m;

    G4Box *solidWorld = new G4Box("solidWorld", 0.5 * xWorld, 0.5 * yWorld, 0.5* zWorld);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, checkOverlaps);

    return physWorld;


}
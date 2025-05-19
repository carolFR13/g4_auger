#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"


DetectorConstruction::DetectorConstruction(){


}

DetectorConstruction::~DetectorConstruction(){

    
}

G4VPhysicalVolume *DetectorConstruction::Construct(){ // we are defining here our main function Construct

    G4bool checkOverlaps = true;

    G4NistManager *nist = G4NistManager::Instance();
    G4Material *worldMat = nist-> FindOrBuildMaterial("G4_AIR");
    G4Material *waterMat = nist -> FindOrBuildMaterial("G4_WATER");

    G4double xWorld = 5. * m;
    G4double yWorld = 5. * m;
    G4double zWorld = 5. * m;

    G4Box *solidWorld = new G4Box("solidWorld", 0.5 * xWorld, 0.5 * zWorld, 0.5* yWorld);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, checkOverlaps);
    
    G4double radius = 1.8 * m;
    G4double height = 1.2 * m;
    G4Tubs *solidTank = new G4Tubs("solidTank", 0., radius, 0.5 * height, 0., 360. *deg);
    G4LogicalVolume *logicTank = new G4LogicalVolume(solidTank, waterMat, "logicTank");
    
    G4VisAttributes *tankVisAtt = new G4VisAttributes(G4Colour(0.2, 0.7, 1.0, 0.4));
    tankVisAtt->SetVisibility(true);
    tankVisAtt->SetForceSolid(true);
    tankVisAtt->SetForceAuxEdgeVisible(true);
    logicTank->SetVisAttributes(tankVisAtt);
    
    G4RotationMatrix *rotY = new G4RotationMatrix();
    rotY->rotateX(90. * deg);
    
    
    new G4PVPlacement(rotY, G4ThreeVector(0., 0., 0.), logicTank, "physTank", logicWorld, false, 0, checkOverlaps);
    

    return physWorld;


}
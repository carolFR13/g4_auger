#ifndef DETECTORCONSTRUCTION
#define DETECTORCONSTRUCTION

#include "G4VUserDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"


#include "G4SubtractionSolid.hh"

class DetectorConstruction : public G4VUserDetectorConstruction{

    public: 
        DetectorConstruction();
        virtual ~DetectorConstruction(); // virtual function cause its already defined in the parent class 

        virtual G4VPhysicalVolume* Construct();

    private:
        void LoadOpticalPropertiesFromFile(
            std::vector<G4double>& photonEnergies,
            std::vector<G4double>& refractiveIndices,
            std::vector<G4double>& absorptionLengths);
};

#endif
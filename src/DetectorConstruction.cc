#include "DetectorConstruction.hh"

#include <fstream>
#include <sstream>
#include <numeric>  // for std::iota
#include <algorithm> // for std::sort

DetectorConstruction::DetectorConstruction()
: flogicPMT(nullptr)
{}

DetectorConstruction::~DetectorConstruction(){}

void DetectorConstruction::LoadOpticalPropertiesFromFile(
    std::vector<G4double>& photonEnergies,
    std::vector<G4double>& refractiveIndices,
    std::vector<G4double>& absorptionLengths)
{
    std::ifstream file("../data/water_index.txt");
    if (!file.is_open()) {
        G4Exception("DetectorConstruction", "FileError", FatalException, "Could not open water_index.txt");
        return;
    }

    std::string line;
    std::map<G4double, G4double> n_map;
    std::map<G4double, G4double> k_map;
    bool readingN = false, readingK = false;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line.find("wl\tn") != std::string::npos) {
            readingN = true;
            readingK = false;
            continue;
        } else if (line.find("wl\tk") != std::string::npos) {
            readingN = false;
            readingK = true;
            continue;
        }

        std::istringstream iss(line);
        double wavelength_um, value;
        if (!(iss >> wavelength_um >> value)) continue;

        if (readingN) {
            n_map[wavelength_um] = value;
        } else if (readingK) {
            k_map[wavelength_um] = value;
        }
    }

    file.close();

    // Combine both maps using only wavelengths present in both
    for (const auto& [wavelength_um, n_val] : n_map) {
        if (k_map.find(wavelength_um) == k_map.end()) continue;
        double k_val = k_map[wavelength_um];

        // Convert wavelength in μm to photon energy in eV
        double wavelength_cm = wavelength_um * 1e-4;
        double energy = (1.2398419843320026 / wavelength_um) * eV;  // energy in eV
        photonEnergies.push_back(energy);
        refractiveIndices.push_back(n_val);

        // Compute absorption length from k: α = 4πk/λ → L = 1/α
        double alpha = (4.0 * pi * k_val) / wavelength_cm; // in cm^-1
        double absLength = 1.0 / alpha; // in cm
        absorptionLengths.push_back(absLength * cm); // convert to G4 units
    }

    // Sort by increasing energy (required by Geant4)
    std::vector<size_t> indices(photonEnergies.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return photonEnergies[a] < photonEnergies[b];
    });

    auto reorder = [&](std::vector<G4double>& vec) {
        std::vector<G4double> sorted;
        for (size_t idx : indices) sorted.push_back(vec[idx]);
        vec = std::move(sorted);
    };

    reorder(photonEnergies);
    reorder(refractiveIndices);
    reorder(absorptionLengths);
}

void DetectorConstruction::DefineMaterials() {
    auto* nist = G4NistManager::Instance();

    worldMat = nist->FindOrBuildMaterial("G4_AIR");
    waterMat = nist->FindOrBuildMaterial("G4_WATER");
    plasticMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    glassMat = nist->FindOrBuildMaterial("G4_GLASS_PLATE");

    // setting optical properties of the water to be pure water
    G4MaterialPropertiesTable* waterMPT = new G4MaterialPropertiesTable();

    std::vector<G4double> energies;
    std::vector<G4double> rindices;
    std::vector<G4double> absLengths;

    LoadOpticalPropertiesFromFile(energies, rindices, absLengths); // ya lo tienes hecho

    waterMPT->AddProperty("RINDEX", energies.data(), rindices.data(), energies.size());
    waterMPT->AddProperty("ABSLENGTH", energies.data(), absLengths.data(), energies.size());

    waterMat->SetMaterialPropertiesTable(waterMPT);

    // optical propierties for air (constant)
    G4MaterialPropertiesTable* airMPT = new G4MaterialPropertiesTable();
    G4double rAir[2] = {1.0, 1.0};
    G4double eAir[2] = {2.0 * eV, 4.0 * eV};
    airMPT->AddProperty("RINDEX", eAir, rAir, 2);
    worldMat->SetMaterialPropertiesTable(airMPT);
}




G4VPhysicalVolume *DetectorConstruction::Construct(){ // we are defining here our main function Construct

    G4bool checkOverlaps = true;

    DefineMaterials(); // we define the materials we are going to use

    // solid volumes definitions

    G4double xWorld = 5. * m;
    G4double yWorld = 5. * m;
    G4double zWorld = 5. * m;

    G4Box *solidWorld = new G4Box("solidWorld", 0.5 * xWorld, 0.5 * zWorld, 0.5* yWorld);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, checkOverlaps);
    

    // external metallic tank to hold the water
    G4double tankRadius = 1.8 * m;
    G4double tankHeight = 1.2 * m;
    G4double wallThickness = 0.02 * m;

    // Create a mother volume for the entire tank assembly
    G4Tubs* solidTankAssembly = new G4Tubs("solidTankAssembly",
                                          0., tankRadius +  wallThickness + 2. * mm,
                                          0.5 * (tankHeight + wallThickness + 2. * mm),
                                          0., 360.*deg);

    G4LogicalVolume* logicTankAssembly = new G4LogicalVolume(solidTankAssembly, worldMat, "logicTankAssembly");

    // we define the metallic shell as the result of a substraction
    G4Tubs* outerTube = new G4Tubs("outerTube",
                               0, (tankRadius + wallThickness),
                               0.5 * (tankHeight + wallThickness),  
                               0, 360. * deg);

    G4Tubs* innerTube = new G4Tubs("innerTube",
                               0, tankRadius,
                               0.5 * tankHeight , 
                               0.0, 360. * deg);

    G4SubtractionSolid* solidTank = new G4SubtractionSolid("solidTank", outerTube, innerTube);


    G4LogicalVolume* logicTank = new G4LogicalVolume(solidTank, plasticMat, "logicTank");

    // Rotation of 90 degrees to place the prism in the g4 axis definition
    auto rotation = new G4RotationMatrix();
    rotation->rotateX(90.*deg);

    // Place the entire assembly in the world
    new G4PVPlacement(rotation, 
                    G4ThreeVector(0., 0., 0.), 
                    logicTankAssembly, 
                    "physTankAssembly", 
                    logicWorld, 
                    false, 
                    0, 
                    checkOverlaps);

    // Place the steel tank shell inside the assembly
    new G4PVPlacement(0, 
                    G4ThreeVector(0., 0., 0.), 
                    logicTank,           
                    "physTank", 
                    logicTankAssembly,  
                    false, 
                    0, 
                    checkOverlaps);

    // internal water tub. height computed from the water volume

    // 12000 L water = 12 m^3 
    // V = pi * r^2 * h -> h = V / (pi * r^2)
    // h = 12 / (pi * 1.8^2) = 1.179 m

    G4double waterHeight = 1.179 * m;

    G4Tubs* solidWater = new G4Tubs("solidWater", 
                            0., tankRadius - 1.*mm,
                            0.5 * waterHeight, 
                            0., 360.*deg);

    G4LogicalVolume* logicWater = new G4LogicalVolume(solidWater, waterMat, "logicWater");

    G4double zOffset = - 0.5 * (tankHeight - waterHeight) ;
    
    new G4PVPlacement(0, 
                    G4ThreeVector(0., 0., zOffset), 
                    logicWater, 
                    "physWater", 
                    logicTankAssembly, 
                    false, 
                    0, 
                    checkOverlaps);


    // PMTs construction

    G4double pmtRadius = 11.5 * cm;
    G4double pmtHeight = 1 * cm;

    G4Tubs* solidPMT = new G4Tubs("solidPMT", 0., pmtRadius, 0.5 * pmtHeight, 0., 360.*deg);
    flogicPMT = new G4LogicalVolume(solidPMT, worldMat, "logicPMT");

    G4double placementRadius = 0.75 * m;  // distance from the center to each PMT
    G4double z_pos = 0.5 * (waterHeight + pmtHeight);

    std::vector<G4ThreeVector> PMT_positions = {
        G4ThreeVector( placementRadius, 0., z_pos),
        G4ThreeVector(-0.5 * placementRadius,  0.866 * placementRadius, z_pos),
        G4ThreeVector(-0.5 * placementRadius, -0.866 * placementRadius, z_pos)
        };

    for (size_t i = 0; i < PMT_positions.size(); ++i) {
        new G4PVPlacement(0,
                        PMT_positions[i],
                        flogicPMT,
                        "physPMT",
                        logicTankAssembly,  
                        false,
                        i,
                        checkOverlaps);      
        }


    auto* pmtSurface = new G4OpticalSurface("WaterToPMTSurface");
    pmtSurface->SetType(dielectric_metal);
    pmtSurface->SetFinish(polished);
    pmtSurface->SetModel(glisur);

    // Define arrays for optical properties
    const G4int NUM_ENTRIES = 2;
    G4double photonEnergy[NUM_ENTRIES] = {2.0*eV, 4.0*eV};
    G4double reflectivity[NUM_ENTRIES] = {0.9, 0.9};
    G4double efficiency[NUM_ENTRIES] = {0.8, 0.8};

    auto* surfaceMPT = new G4MaterialPropertiesTable();
    surfaceMPT->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, NUM_ENTRIES);
    surfaceMPT->AddProperty("EFFICIENCY", photonEnergy, efficiency, NUM_ENTRIES);
    
    pmtSurface->SetMaterialPropertiesTable(surfaceMPT);

    new G4LogicalSkinSurface("PMTSurface", flogicPMT, pmtSurface);
    
    G4VisAttributes *waterVisAtt = new G4VisAttributes(G4Colour(0.2, 0.7, 1.0, 0.3));
    waterVisAtt->SetVisibility(true);
    waterVisAtt->SetForceSolid(true);
    waterVisAtt->SetForceAuxEdgeVisible(true);
    logicWater->SetVisAttributes(waterVisAtt);

    G4VisAttributes* tankVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.6)); // greyish
    tankVisAtt->SetVisibility(true);
    tankVisAtt->SetForceSolid(true);
    tankVisAtt->SetForceAuxEdgeVisible(true);
    logicTank ->SetVisAttributes(tankVisAtt);

    // make the assembly invisible
    G4VisAttributes* assemblyVisAtt = new G4VisAttributes();
    assemblyVisAtt->SetVisibility(false);
    logicTankAssembly->SetVisAttributes(assemblyVisAtt);

    G4VisAttributes* pmtVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.9)); // pure white
    pmtVisAtt->SetVisibility(true);
    pmtVisAtt->SetForceSolid(true);
    pmtVisAtt->SetForceAuxEdgeVisible(true);
    flogicPMT->SetVisAttributes(pmtVisAtt);

    return physWorld;
}


void DetectorConstruction::ConstructSDandField(){

    G4cout << "Constructing sensitive detectors..." << G4endl;

    if(flogicPMT != nullptr) {
        G4SDManager *sdManager = G4SDManager::GetSDMpointer();
        SensitiveDetector *sensDet = new SensitiveDetector("SensitiveDetector");
        sdManager->AddNewDetector(sensDet);
        flogicPMT->SetSensitiveDetector(sensDet);
        G4cout << "Sensitive detector set for logicDetector" << G4endl;

    } else {
        G4cerr << "Error: logicPMT is null in ConstructSDandField" << G4endl;
    }
    
}

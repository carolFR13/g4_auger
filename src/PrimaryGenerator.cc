#include "PrimaryGenerator.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include <cmath>
#include "G4PhysicalConstants.hh"

PrimaryGenerator::PrimaryGenerator(){

    fParticleGun = new G4ParticleGun(1);


    // Particle type

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable -> FindParticle("mu-");

    fParticleGun -> SetParticleEnergy(4. * GeV);
    fParticleGun -> SetParticleDefinition(particle);

}


PrimaryGenerator::~PrimaryGenerator(){

    delete fParticleGun;
}

void PrimaryGenerator::GeneratePrimaries(G4Event *anEvent){

    // Particle position: circle of radius 2m in the xz plane at y = 1.5m
    // G4double radius = 2.0 * m;
    // G4double r = radius * std::sqrt(G4UniformRand());
    // G4double phi1 = 2. * pi * G4UniformRand();
    
    // G4double x = r * std::cos(phi1);
    // G4double z = r * std::sin(phi1);
    // G4double y = 1.5 * m;



    // Particle position: spherical casquette of radius 2m

    G4double r = 1.5 * m; // constant
    G4double theta = 0.5 * pi * G4UniformRand();  // from 0 to pi/2
    G4double phi = 2.0 * pi * G4UniformRand();  // from 0 to 2pi

    G4double x = r * std::sin(theta) * std::cos(phi);
    G4double z = r * std::sin(theta) * std::sin(phi);
    G4double y = r * std::cos(theta);

    G4ThreeVector pos(x, y, z);
    fParticleGun->SetParticlePosition(pos);


    // Particle direction
    G4double theta_mom;
    do {
        theta_mom = std::acos(G4UniformRand());
    } while (G4UniformRand() > std::pow(std::cos(theta_mom), 2));
    
    G4double phi_mom = 2. * pi * G4UniformRand();

    G4double px = std::sin(theta_mom) * std::cos(phi_mom);
    G4double pz = std::sin(theta_mom) * std::sin(phi_mom);
    G4double py = -std::cos(theta_mom);  // negativo: dirección hacia abajo
    G4ThreeVector mom(px,py,pz);
    fParticleGun->SetParticleMomentumDirection(mom);

    // Create vertex
    fParticleGun -> GeneratePrimaryVertex(anEvent);
    
}
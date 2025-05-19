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

    // Particle position

    G4double radius = 2.0 * m;
    G4double r = radius * std::sqrt(G4UniformRand());
    G4double phi1 = 2. * pi * G4UniformRand();
    
    G4double x = r * std::cos(phi1);
    G4double z = r * std::sin(phi1);
    G4double y = 1.5 * m;

    G4ThreeVector pos(x, y, z);
    fParticleGun->SetParticlePosition(pos);

    // Particle direction
    
    G4double theta;
    do {
        theta = std::acos(G4UniformRand());
    } while (G4UniformRand() > std::pow(std::cos(theta), 2));
    
    G4double phi2 = 2. * pi * G4UniformRand();

    G4double px = std::sin(theta) * std::cos(phi2); 
    G4double pz = std::sin(theta) * std::sin(phi2); 
    G4double py = -std::cos(theta);  // we make sure the particle is forwared in y negative direction

    G4ThreeVector mom(px,py,pz);
    fParticleGun->SetParticleMomentumDirection(mom);

    // Create vertex

    fParticleGun -> GeneratePrimaryVertex(anEvent);
    
}
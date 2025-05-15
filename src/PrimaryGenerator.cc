#include "PrimaryGenerator.hh"

PrimaryGenerator::PrimaryGenerator(){

    fParticleGun = new G4ParticleGun(1);

    // Particle position

    G4double x = 0. * m;
    G4double y = 0. * m;
    G4double z = 0. * m;

    G4ThreeVector pos(x, y, z);

    // Particle direction

    G4double px = 0. ; 
    G4double py = 0. ; 
    G4double pz = 1. ;  // we make sure the particle is forwared in z direction

    G4ThreeVector mom(px,py,pz);

    // Particle type

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable -> FindParticle("mu-");

    fParticleGun -> SetParticlePosition(pos);
    fParticleGun -> SetParticleMomentumDirection(mom);
    fParticleGun -> SetParticleEnergy(4. * MeV);
    fParticleGun -> SetParticleDefinition(particle);

}


PrimaryGenerator::~PrimaryGenerator(){

    delete fParticleGun;
}

void PrimaryGenerator::GeneratePrimaries(G4Event *anEvent){

    // Create vertex

    fParticleGun -> GeneratePrimaryVertex(anEvent);
    
}
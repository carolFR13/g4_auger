#include "PhysicsList.hh"
#include "G4OpticalPhysics.hh"

PhysicsList::PhysicsList(){
    // EM Physics 
    RegisterPhysics(new G4EmStandardPhysics());

    RegisterPhysics(new G4OpticalPhysics());

}

PhysicsList::~PhysicsList(){}
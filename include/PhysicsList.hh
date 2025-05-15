#ifndef PHYSICSLISTS_HH
#define PHYSICSLISTS_HH

#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"

class PhysicsList : public G4VModularPhysicsList{

    public: 
    PhysicsList();
    ~PhysicsList();
};


#endif
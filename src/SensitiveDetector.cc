#include "SensitiveDetector.hh"



SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name) {

    fTotalEnergyDeposited = 0;
}

SensitiveDetector::~SensitiveDetector(){}


void SensitiveDetector::Initialize(G4HCofThisEvent *){

    fTotalEnergyDeposited = 0;
    G4cout << "Starting new event" << G4endl;

}


void SensitiveDetector::EndOfEvent(G4HCofThisEvent *){

    G4cout << "Deposited energy: " << fTotalEnergyDeposited << G4endl;

}

G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory * ){

    G4Track* track = aStep->GetTrack();

    if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        G4ThreeVector pos = track->GetPosition();
        G4double time = track->GetGlobalTime();
        G4int copyID = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
        G4cout << "[PMT " << copyID << "] hit by photon!" << G4endl;
        G4cout << "Position: " << pos << ", Time: " << time / ns << " ns" << G4endl;

        // Terminar el fotón: simula que el PMT lo absorbe
        track->SetTrackStatus(fStopAndKill);

        return true;
    }
    
    return false;

}
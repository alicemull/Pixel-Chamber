//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2PrimaryGeneratorAction.cc
/// \brief Implementation of the B2PrimaryGeneratorAction class

#include "B2PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include <G4GeneralParticleSource.hh>
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Utility/memory.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::B2PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction()
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);
  // UNCOMMENT WHEN GPS IS NEEDED //
  //fParticleGun = new G4GeneralParticleSource();

  // default particle kinematic

  G4ParticleDefinition* particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("proton");

  fParticleGun->SetParticleDefinition(particleDefinition);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  fParticleGun->SetParticleEnergy(400.0*GeV);
  // UNCOMMENT WHEN GPS IS NEEDED //
  // REMEMBER TO DI IN .hh FILE TOO //
  //fParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));    
 //fParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(3.0 * GeV);
  //G4RandGauss::setTheSeed(1784);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::~B2PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore.

  G4double worldZHalfLength = 0;
  //G4double worldYHalfLength = 0;
  G4LogicalVolume* worldLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = NULL;
  if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  if ( worldBox ) worldZHalfLength = worldBox->GetZHalfLength();
  //if ( worldBox ) worldYHalfLength = worldBox->GetYHalfLength();
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  G4double PixXDim = 29.24 *micrometer;
  G4double PixYDim = 50 *micrometer;
  G4double PixZDim = 26.88 *micrometer;
  G4int ColNo = 1024; 
  G4int RawNo = 512;
  G4int NbOfSensors = 216;
  G4double Pi = 3.14159265358979323846;
  G4double Theta1 = 2*Pi*G4UniformRand();
  G4double r1 = G4UniformRand()*mm;
  G4double x_i = PixXDim*1024/2;
  //G4double y_i = r1*sin(Theta1);
  //G4double z_i = r1*cos(Theta1);
  G4double mean = 0.;
  G4double sigma_i = 0.4;
  
  
  
  //CLHEP::HepRandom::setTheSeed(17);
  G4double y_i =G4RandGauss::shoot(mean, sigma_i)*mm;
  G4double z_i = G4RandGauss::shoot(mean,sigma_i)*mm;
  G4double r = G4UniformRand()*0.2*mm;
  G4double Theta = 2*Pi*G4UniformRand();
  G4double x_f = x_i;
  //  G4double y_f = r*cos(Theta);
  // G4double z_f = r*sin(Theta);
  G4double px = x_f;
  // G4double py =abs(y_f-y_i);
  // G4double pz = abs(z_f - z_i);
  G4double sigma_f = 0.02;
  G4double py = G4RandGauss::shoot(y_i, sigma_f);
  G4double pz = G4RandGauss::shoot(z_i, sigma_f);
  G4double norm = sqrt(px*px+py*py+pz*pz);
  fParticleGun->SetParticlePosition(G4ThreeVector(-x_i, y_i, z_i));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px/norm,py/norm,pz/norm));
  //fParticleGun->SetParticlePosition(G4ThreeVector(0., worldYHalfLength, 0.));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

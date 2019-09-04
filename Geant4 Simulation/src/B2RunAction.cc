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
/// \file B2RunAction.cc
/// \brief Implementation of the B2RunAction class

#include "B2RunAction.hh"
#include "B4Analysis.hh"
#include "G4SteppingControl.hh"
#include "G4StepPoint.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2RunAction::B2RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1);
 
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2RunAction::~B2RunAction()
{
delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(false);

        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
 

  analysisManager->SetFileName("B2");
  analysisManager->OpenFile("B2");
  analysisManager->SetFirstNtupleId(1);
  analysisManager->CreateNtuple("Interest", "Values");


  analysisManager->CreateNtupleIColumn("PixelNb");
  analysisManager->CreateNtupleIColumn("TrackID");
  analysisManager->CreateNtupleIColumn("PartPDGCode");
  analysisManager->CreateNtupleIColumn("SensorNb");
  analysisManager->CreateNtupleDColumn("Edep"); 
  analysisManager->CreateNtupleIColumn("EventID"); 
  analysisManager->CreateNtupleDColumn("PixelCentreXPosition"); 
  analysisManager->CreateNtupleDColumn("PixelCentreYPosition"); 
  analysisManager->CreateNtupleDColumn("PixelCentreZPosition"); 
  analysisManager->CreateNtupleIColumn("BigSensorNb");
  analysisManager->CreateNtupleDColumn("Momentum_x");
  analysisManager->CreateNtupleDColumn("Momentum_y");
  analysisManager->CreateNtupleDColumn("Momentum_z");
  
  analysisManager->FinishNtuple();
  G4cout << "Ntuple-1 created" << G4endl;
  //analysisManager->OpenFile("B2");


  //analysisManager->SetFileName("B2_D");
 
  analysisManager->CreateNtuple("Interest_D", "Values_D");
  
  analysisManager->CreateNtupleIColumn("EventID_D"); 
  analysisManager->CreateNtupleIColumn("i");
  analysisManager->CreateNtupleIColumn("j");
  analysisManager->CreateNtupleIColumn("k");
  analysisManager->CreateNtupleIColumn("ParticlePDGCode_D");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("TrackID_1");
  analysisManager->CreateNtupleDColumn("TrackID_2");
  analysisManager->CreateNtupleDColumn("TrackID_3");
  analysisManager->CreateNtupleDColumn("Momentum_x_D");
  analysisManager->CreateNtupleDColumn("Momentum_y_D");
  analysisManager->CreateNtupleDColumn("Momentum_z_D");
 
  

  analysisManager->FinishNtuple();
  G4cout << "Ntuple-2 created" << G4endl;
  // analysisManager->OpenFile("B2_D");
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2RunAction::EndOfRunAction(const G4Run* )
{
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


	// save histograms & ntuple
	//
	analysisManager->Write();
	analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

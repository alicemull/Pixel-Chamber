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
/// \file B2MagneticField.hh
/// \brief Definition of the B2MagneticField class

#ifndef B2MagneticField_H
#define B2MagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

class G4GenericMessenger;

/// Magnetic field

class B2MagneticField : public G4MagneticField
{
public:
  B2MagneticField();
  virtual ~B2MagneticField();
    
  virtual void GetFieldValue(const G4double point[4],double* bField ) const;
    
  // void SetField(G4double bx, G4double by, G4double bz) {
  //   fBx = bx;
  //   fBy = by;
  //   fBz = bz;
  // }

  void SetField(G4double val)  { fBy = val; }
  void SetXField(G4double val) { fBx = val; }
  void SetYField(G4double val) { fBy = val; }
  void SetZField(G4double val) { fBz = val; }
 
  G4double GetFieldX() const { return fBx; }
  G4double GetFieldY() const { return fBy; }
  G4double GetFieldZ() const { return fBz; }
   
private:
  void DefineCommands();

  G4GenericMessenger* fMessenger;
  G4double fBx;
  G4double fBy;
  G4double fBz;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

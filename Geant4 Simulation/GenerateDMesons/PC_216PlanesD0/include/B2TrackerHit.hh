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
/// \file B2TrackerHit.hh
/// \brief Definition of the B2TrackerHit class

#ifndef B2TrackerHit_h
#define B2TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"


/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class B2TrackerHit : public G4VHit
{
public:
  B2TrackerHit();
  B2TrackerHit(const B2TrackerHit&);
  virtual ~B2TrackerHit();
  
  // operators
  const B2TrackerHit& operator=(const B2TrackerHit&);
  G4bool operator==(const B2TrackerHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  // methods from base class
  virtual void Draw();
  virtual void Print();

  
  // Set methods
  void SetTrackID  (G4int track)      { fTrackID = track; };
  void SetPixelNb  (G4int pix)        { fPixelNb = pix; };
  void SetBigSensorNb (G4int Bigsens) { fBigSensorNb = Bigsens; };
  void SetEdep     (G4double de)      { fEdep = de; };
  void SetETot     (G4double etot)    { fEtot = etot;};
  void SetPos      (G4ThreeVector xyz){ fPos = xyz; };
  void SetParticleName (G4String name){ fPartName = name;};
  void SetPDGCode (G4int PDGc)        { fPartPDG = PDGc;};
  void SetPixelPos (G4ThreeVector x_py_pz_p){ fPixelPos = x_py_pz_p; };
  void SetEventID  (G4int ev)         { fEventID = ev; };
  void SetPixelPosX (G4double x_p){ fPixelPosx = x_p; };
  void SetPixelPosY (G4double y_p){ fPixelPosy = y_p; };
  void SetPixelPosZ (G4double z_p){ fPixelPosz = z_p; };
  void SetVert(G4ThreeVector vxvyvz){ fVertex = vxvyvz;};
  void SetVertexPosX (double x_v){ fVertexPosx = x_v; };
  void SetVertexPosY (double y_v){ fVertexPosy = y_v; };
  void SetVertexPosZ (double z_v){ fVertexPosz = z_v; };  
  void SetSensorNb (G4int sens)       { fSensorNb = sens; };
  void SetMom(G4ThreeVector pxpypz) { fMomentum = pxpypz;};
  void SetMomentum_x (double px) {fMomentum_x = px;};
  void SetMomentum_y (double py) {fMomentum_y = py;};
  void SetMomentum_z (double pz) {fMomentum_z = pz;};
  
  

  // Get methods
  G4int GetTrackID() const     { return fTrackID; };
  G4int GetPixelNb() const     { return fPixelNb; };
  G4int GetSensorNb() const    { return fSensorNb; };
  G4double GetEdep() const     { return fEdep; };
  G4double GetETot() const     {return fEtot;};
  G4ThreeVector GetPos() const { return fPos; };
  G4String GetParticleName() const     { return fPartName;};
  G4int GetPDGCode() const     { return fPartPDG;};
  G4ThreeVector GetPixelPos() const { return fPixelPos; };
  G4int GetEventID() const     { return fEventID; };
  G4double GetPixelPosX() const { return fPixelPosx; };
  G4double GetPixelPosY() const { return fPixelPosy; };
  G4double GetPixelPosZ() const { return fPixelPosz; };
  G4ThreeVector GetVert() const { return fVertex; };
  G4double GetVertexPosX() const { return fVertexPosx; };
  G4double GetVertexPosY() const { return fVertexPosy; };
  G4double GetVertexPosZ() const { return fVertexPosz; };
  G4int GetBigSensorNb() const    { return fBigSensorNb; };
  G4ThreeVector GetMom() const { return fMomentum; };
  G4double GetMomentum_x() const { return fMomentum_x;};
  G4double GetMomentum_y() const { return fMomentum_y;};
  G4double GetMomentum_z() const { return fMomentum_z;};


  private:

  G4int         fTrackID;
  G4int         fPixelNb;
  G4int         fSensorNb;
  G4double      fEdep;
  G4double      fEtot;
  G4ThreeVector fPos;
  G4String      fPartName;
  G4int         fPartPDG;
  G4ThreeVector fPixelPos;
  G4int         fEventID;
  G4double      fPixelPosx;
  G4double      fPixelPosy;
  G4double      fPixelPosz;
  G4ThreeVector fVertex;
  G4double      fVertexPosx;
  G4double      fVertexPosy;
  G4double      fVertexPosz;
  G4int         fBigSensorNb;
  G4ThreeVector fMomentum;
  G4double fMomentum_x;
  G4double fMomentum_y;
  G4double fMomentum_z;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<B2TrackerHit> B2TrackerHitsCollection;

extern G4ThreadLocal G4Allocator<B2TrackerHit>* B2TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* B2TrackerHit::operator new(size_t)
{
  if(!B2TrackerHitAllocator)
      B2TrackerHitAllocator = new G4Allocator<B2TrackerHit>;
  return (void *) B2TrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void B2TrackerHit::operator delete(void *hit)
{
  B2TrackerHitAllocator->FreeSingle((B2TrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

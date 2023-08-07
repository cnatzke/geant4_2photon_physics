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
// $Id$
//
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4TwoPhotonEvaporation.hh                                         //
//  Author: C. Natzke                                                         //
//  Based on file: G4PhotonEvaporation.cc                                     //
//  Date:  2022-04-19                                                         //
//  Description: performs isomeric transition for excited states of           //
//               radioactive nuclei by emitting two gammas                    //
//               and returns daughter particles the rest frame                //
//               of the parent nucleus                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4TWOPHOTONEVAPORATION_HH
#define G4TWOPHOTONEVAPORATION_HH 1

#include "globals.hh"
#include "G4VEvaporationChannel.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4Threading.hh"
#include "Randomize.hh"

const G4int MAXDEPOINT2G = 10;
const G4int MAXGRDATA2G = 300;

class G4TwoPhotonTransition;

class G4TwoPhotonEvaporation : public G4VEvaporationChannel
{

public:
  explicit G4TwoPhotonEvaporation(G4TwoPhotonTransition *ptr = nullptr);

  virtual ~G4TwoPhotonEvaporation();

  virtual void Initialise() final;

  // two-photon emission
  virtual G4FragmentVector *EmittedFragments(G4Fragment *theNucleus) final;

  // compute emission probability for both continum and discrete cases
  // must be called before any method above
  virtual G4double GetEmissionProbability(G4Fragment *theNucleus) final;

  virtual G4double GetFinalLevelEnergy(G4int Z, G4int A, G4double energy) final;

  virtual G4double GetUpperLevelEnergy(G4int Z, G4int A) final;

  void SetTwoPhotonTransition(G4TwoPhotonTransition *);

  virtual void SetICM(G4bool);

  virtual void RDMForced(G4bool);

  inline void SetVerboseLevel(G4int verbose);

  inline void SetMultipoleMixingRatio(G4float multipoleMixing);

  inline void SetAngularRatio(G4float angularRatio);

  inline G4int GetVacantShellNumber() const;

  G4TwoPhotonEvaporation(const G4TwoPhotonEvaporation &right) = delete;
  const G4TwoPhotonEvaporation &operator=(const G4TwoPhotonEvaporation &right) = delete;

private:
  void InitialiseGRData();

  G4FragmentVector *GenerateGammas(G4Fragment *nucleus);

  inline void InitialiseLevelManager(G4int Z, G4int A);

  G4NuclearLevelData *fNuclearLevelData;
  const G4LevelManager *fLevelManager;
  G4TwoPhotonTransition *fTransition;

  G4int fVerbose;
  G4int theZ;
  G4int theA;
  G4int fPoints;
  G4int fCode;
  size_t fIndex;

  static G4float GREnergy2G[MAXGRDATA2G];
  static G4float GRWidth2G[MAXGRDATA2G];

  G4double fCummProbability[MAXDEPOINT2G];

  G4double fLevelEnergyMax;
  G4double fExcEnergy;
  G4double fProbability;
  G4double fStep;
  G4double fMaxLifeTime;
  G4double fRelativeBR;
  G4double fMultipoleMixing;
  G4double fAngularRatio;

  G4double Tolerance;

  G4bool fICM;
  G4bool fRDM;
  G4bool fSampleTime;
  G4bool fIsomerFlag;
  G4bool fTwoPhotonTransition;
  G4bool isInitialised;

#ifdef G4MULTITHREADED
  static G4Mutex PhotonEvaporationMutex;
#endif
};

inline void G4TwoPhotonEvaporation::SetVerboseLevel(G4int verbose)
{
  fVerbose = verbose;
}

inline void
G4TwoPhotonEvaporation::SetMultipoleMixingRatio(G4float multipoleMixing)
{
  fMultipoleMixing = multipoleMixing;
}

inline void
G4TwoPhotonEvaporation::SetAngularRatio(G4float angularRatio)
{
  fAngularRatio = angularRatio;
}

inline void
G4TwoPhotonEvaporation::SetRelativeBR(G4float branchingRatio)
{
  fRelativeBR = branchingRatio;
}

inline void
G4TwoPhotonEvaporation::InitialiseLevelManager(G4int Z, G4int A)
{
  if (Z != theZ || A != theA)
  {
    theZ = Z;
    theA = A;
    fIndex = 0;
    fLevelManager = fNuclearLevelData->GetLevelManager(theZ, theA);
    fLevelEnergyMax = fLevelManager ? fLevelManager->MaxLevelEnergy() : 0.0;
  }
}

inline G4int G4TwoPhotonEvaporation::GetVacantShellNumber() const
{
  return vShellNumber;
}

#endif

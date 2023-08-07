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
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      Triumf, Vancouver, Canada
//
//      File name:     G4TwoPhotonTransition
//
//      Author:        C. Natzke
//
//      Creation date: 2022-06-20
//
//      Modifications:
//
// -------------------------------------------------------------------

#ifndef G4TwoPhotonTransition_HH
#define G4TwoPhotonTransition_HH 1

#include "globals.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4PolarizationTransition.hh"
#include "G4FragmentVector.hh"
#include "Randomize.hh"

class G4TwoPhotonTransition
{
public:
  explicit G4TwoPhotonTransition();

  virtual ~G4TwoPhotonTransition();

  virtual G4FragmentVector *SampleTransition(G4Fragment *nucleus,
                                             G4double newExcEnergy,
                                             G4double multipoleRatio,
                                             G4double angularRatio,
                                             G4int shell,
                                             G4bool isGamma,
                                             G4bool isTwoPhoton);

  virtual void SampleEnergy(G4double totalTransEnergy);

  virtual void SampleDirection();

  inline G4ThreeVector GetPhoton0Direction()
  {
    return fDirectionPhoton0;
  }

  inline G4ThreeVector GetPhoton1Direction()
  {
    return fDirectionPhoton1;
  }

  inline G4double GetPhoton0Energy()
  {
    return fGamma0Energy;
  }

  inline G4double GetPhoton1Energy()
  {
    return fGamma1Energy;
  }

  inline void SetVerbose(G4int val)
  {
    fVerbose = val;
  };

private:
  void SetUpEnergySpectrumSampler(G4double transitionEnergy, G4float multipoleRatio);
  void SetUpAngularDistributionSampler(G4float alphaE1, G4float chi);
  void CreateRotationMatrix(const G4ThreeVector &vector1, const G4ThreeVector &vector2);
  G4ThreeVector RotateVector(const G4ThreeVector &vector);
  G4ThreeVector SphericalToCartesian(const G4double &theta, const G4double &phi);

  G4TwoPhotonTransition(const G4TwoPhotonTransition &right) = delete;
  const G4TwoPhotonTransition &operator=(const G4TwoPhotonTransition &right) = delete;
  G4bool operator==(const G4TwoPhotonTransition &right) const = delete;
  G4bool operator!=(const G4TwoPhotonTransition &right) const = delete;

  G4int fDim;
  G4double fGamma0Energy, fGamma1Energy;
  G4ThreeVector fDirection;
  G4ThreeVector fDirectionPhoton0, fDirectionPhoton1;
  std::vector<std::vector<G4double>> fRotationMatrix;
  G4RandGeneral *energySpectrumSampler;
  G4RandGeneral *angularDistributionSampler;

  G4int fVerbose;
  G4double fMultipoleRatio, fAngularRatio;

protected:
};

#endif

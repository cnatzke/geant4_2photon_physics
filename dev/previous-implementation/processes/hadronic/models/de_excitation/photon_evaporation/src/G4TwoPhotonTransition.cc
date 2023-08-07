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
//      TRIUMF, Vancouver, Canada
//
//      File name:     G4TwoPhotonTransition
//
//      Author C. Natzke 2022-06-20
//
// -------------------------------------------------------------------

#include "G4TwoPhotonTransition.hh"
#include "G4AtomicShells.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4LorentzVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4TwoPhotonTransition::G4TwoPhotonTransition()
    : fDim(3), fDirectionPhoton1(0., 0., 0.), fDirectionPhoton2(0., 0., 0.), fVerbose(0)
{
  fRotationMatrix.resize(fDim, std::vector<G4double>(fDim, 0.));
}

G4TwoPhotonTransition::~G4TwoPhotonTransition()
{
}

std::vector<G4Fragment *>
G4TwoPhotonTransition::SampleTransition(G4Fragment *nucleus,
                                        G4double newExcEnergy,
                                        G4double multipoleRatio,
                                        G4double angularRatio)
{
  fMultipoleRatio = multipoleRatio;
  fAngularRatio = angularRatio;
  G4Fragment *resultGamma1 = nullptr;
  G4Fragment *resultGamma2 = nullptr;
  std::vector<G4Fragment *> resultsVector;

  G4double totalTransEnergy = nucleus->GetExcitationEnergy() - newExcEnergy;
  if (fVerbose > 2)
  {
    G4cout << "G4TwoPhotonTransition::GenerateGamma - totalTransEnergy(MeV)= "
           << totalTransEnergy << "  Eexnew= " << newExcEnergy
           << G4endl;
  }
  if (totalTransEnergy <= 0.0)
  {
    G4Exception("G4TwoPhotonTransition::SampleTransition()", "HAD_TWOPHOTON_000", FatalException, "Total transition energy is less than 0.0");
  }

  // select secondaries
  G4ParticleDefinition *gamma1 = G4Gamma::Gamma();
  G4ParticleDefinition *gamma2 = G4Gamma::Gamma();

  // get energies of photons and emission directions
  SampleEnergy(totalTransEnergy);
  SampleDirection();

  // Do complete Lorentz computation
  G4LorentzVector lv = nucleus->GetMomentum();
  G4double mass = nucleus->GetGroundStateMass() + newExcEnergy;
  G4double ecm = lv.mag();
  G4ThreeVector bst = lv.boostVector();

  ecm = std::max(ecm, mass);
  G4double energy = 0.5 * ((ecm - mass) * (ecm + mass)) / ecm;
  G4double momPhoton1 = energy * (eGamma1 / totalTransEnergy);
  G4double momPhoton2 = energy - momPhoton1;

  G4LorentzVector photon1FourMom(momPhoton1 * fDirectionPhoton1.x(),
                                 momPhoton1 * fDirectionPhoton1.y(),
                                 momPhoton1 * fDirectionPhoton1.z(), eGamma1);

  G4LorentzVector photon2FourMom(momPhoton2 * fDirectionPhoton2.x(),
                                 momPhoton2 * fDirectionPhoton2.y(),
                                 momPhoton2 * fDirectionPhoton2.z(), eGamma2);

  // G4cout << momPhoton1 << ", " << momPhoton2 << ", " << eGamma1 << ", " << eGamma2 << G4endl;

  // residual
  energy = std::max(ecm - energy, mass);
  lv.set(-(momPhoton1 * fDirectionPhoton1.x() + momPhoton2 * fDirectionPhoton2.x()),
         -(momPhoton1 * fDirectionPhoton1.y() + momPhoton2 * fDirectionPhoton2.y()),
         -(momPhoton1 * fDirectionPhoton1.z() + momPhoton2 * fDirectionPhoton2.z()),
         energy);

  // Lab system transform for short lived level
  lv.boost(bst);

  // modified primary fragment
  nucleus->SetExcEnergyAndMomentum(newExcEnergy, lv);

  // gamma or e- are produced
  photon1FourMom.boost(bst);
  photon2FourMom.boost(bst);

  resultGamma1 = new G4Fragment(photon1FourMom, gamma1);
  resultGamma2 = new G4Fragment(photon2FourMom, gamma2);

  resultsVector.push_back(resultGamma1);
  resultsVector.push_back(resultGamma2);

  /*
  if (fVerbose > 2)
  {
    G4cout << "G4TwoPhotonTransition::SampleTransition : " << *resultsVector << G4endl;
    G4cout << "       Left nucleus: " << *nucleus << G4endl;
  }
  */

  return resultsVector;
}

void G4TwoPhotonTransition::SampleEnergy(G4double totalTransEnergy)
{
  // find continuous energy distribution of photons
  SetUpEnergySpectrumSampler(totalTransEnergy, fMultipoleRatio);

  if (energySpectrumSampler)
  {

    eGamma1 = totalTransEnergy * energySpectrumSampler->shoot(G4Random::getTheEngine());
    eGamma2 = totalTransEnergy - eGamma1; // keV

    if (fVerbose > 2)
    {
      G4cout << "G4TwoPhotonTransition::eGamma1: " << eGamma1 << " | "
             << "G4TwoPhotonTransition::eGamma2: " << eGamma2 << " | "
             << "G4TwoPhotonTransition::totalTransEnergy: " << eGamma1 + eGamma2
             << G4endl;
    }
  }
  else
  {
    G4Exception("G4TwoPhotonTransition::SampleTransition()", "HAD_TWOPHOTON_001", FatalException, "No initialized energy spectrum sampler");
  }

  energySpectrumSampler = NULL;
}

void G4TwoPhotonTransition::SampleDirection()
{
  std::vector<std::vector<G4double>> rotationMatrix;

  G4double alphaE1 = fAngularRatio;
  G4double chi = 1.0;

  SetUpAngularDistributionSampler(alphaE1, chi);

  if (angularDistributionSampler)
  {
    // the first photon is emitted in a random direction and sets the polarization axis
    G4double theta1 = pi * G4UniformRand();
    G4double phi1 = twopi * G4UniformRand();

    // sample the angle between the two photons
    G4double theta2 = pi * angularDistributionSampler->shoot(G4Random::getTheEngine());
    G4double phi2 = twopi * G4UniformRand();

    G4ThreeVector labAxis = G4ThreeVector(0, 0, 1.);
    G4ThreeVector polarizationAxis = SphericalToCartesian(theta1, phi1);
    G4ThreeVector emissionVector = SphericalToCartesian(theta2, phi2);

    // find rotation matrix that maps lab axis to the polarization axis
    CreateRotationMatrix(labAxis, polarizationAxis);

    // rotate emission vector to polarization axis
    G4ThreeVector rotatedEmissionVector = RotateVector(emissionVector);

    // Check the angle between vectors
    if (fVerbose > 2)
    {
      G4cout << "###### Angle SANITY CHECK ######" << G4endl;
      G4cout << "Theta: " << theta2 << " | Angle: " << polarizationAxis.angle(rotatedEmissionVector) << G4endl;
      G4cout << "#########################" << G4endl;
    }

    fDirectionPhoton1.set(polarizationAxis.getX(), polarizationAxis.getY(), polarizationAxis.getZ());
    fDirectionPhoton2.set(rotatedEmissionVector.getX(), rotatedEmissionVector.getY(), rotatedEmissionVector.getZ());
  }
  else
  {
    G4Exception("G4TwoPhotonTransition::SampleTransition()", "HAD_TWOPHOTON_002", FatalException, "No initialized angular distribution sampler");
  }

  angularDistributionSampler = NULL;
}

void G4TwoPhotonTransition::SetUpEnergySpectrumSampler(G4double transitionEnergy, G4float multipoleRatio)
{
  if (transitionEnergy > 0)
  {
    // Array to store spectrum pdf
    G4int npti = 100;
    G4double *pdf = new G4double[npti];

    G4double e;             // energy of one photon
    G4double percentDipole; // percentage of decay that is pure dipole
    G4double d, q, f;       // normalized pdf
    for (G4int ptn = 0; ptn < npti; ptn++)
    {
      // Sample energy range
      e = transitionEnergy * (G4double(ptn) + 0.5) / G4double(npti);

      // Build numberical pdf

      // Normalized pdf for pure dipole transition
      d = 140 * std::pow(e, 3) * pow((transitionEnergy - e), 3) / std::pow(transitionEnergy, 7);

      // Normalized pdf for pure quadrupole transition
      q = 2772 * std::pow(e, 5) * pow((transitionEnergy - e), 5) / std::pow(transitionEnergy, 11);

      // mixing ratio
      percentDipole = 1. / (1 + std::pow(multipoleRatio, 2));
      f = percentDipole * d + (1 - percentDipole) * q;

      pdf[ptn] = f;
    }
    energySpectrumSampler = new G4RandGeneral(pdf, npti);
    delete[] pdf;
  }
}

void G4TwoPhotonTransition::SetUpAngularDistributionSampler(G4float alphaE1, G4float chi)
{
  // Array to store spectrum pdf
  G4int npti = 100;
  G4double *pdf = new G4double[npti];

  G4double theta; // angle between photons
  G4double w, f;  // angular distribution function

  G4double norm = 2 / (3 * CLHEP::pi);                                            // normalizing factor
  G4double coeff = 4 * alphaE1 * chi / (std::pow(alphaE1, 2) + std::pow(chi, 2)); // interference term coefficient

  for (G4int ptn = 0; ptn <= npti; ptn++)
  {
    // Sample angular range
    theta = CLHEP::pi * (G4double(ptn) + 0.5) / G4double(npti);

    // Build numberical pdf
    // Normalized pdf for angle between photons
    // J. Kramp, D. Habs, R. Kroth, M. Music, J. Schirmer, D. Schwalm, and C. Broude, Nuclear Two-Photon Decay in 0+→0+ Transitions, Nuclear Physics, Section A 474, 412 (1987).
    w = 1 + coeff * std::cos(theta) + std::pow(std::cos(theta), 2);
    f = norm * w;

    pdf[ptn] = f;
  }
  angularDistributionSampler = new G4RandGeneral(pdf, npti);
  delete[] pdf;
}

void G4TwoPhotonTransition::CreateRotationMatrix(const G4ThreeVector &vector1, const G4ThreeVector &vector2)
{
  G4ThreeVector a, b, v;
  G4double c;
  G4double d[3][3];
  G4int i, j, k;

  // first normalize the vectors
  a = vector1.unit();
  b = vector2.unit();

  v = a.cross(b);
  c = a.dot(b);

  G4double iMat[3][3] = {{1., 0, 0},
                         {0, 1., 0},
                         {0, 0, 1.}};

  std::vector<std::vector<G4double>> rotationMatrix{{0., -v.z(), v.y()},
                                                    {v.z(), 0., -v.x()},
                                                    {-v.y(), v.x(), 0.0}};

  // set product matrix to zero
  for (i = 0; i < fDim; i++)
  {
    for (j = 0; j < fDim; j++)
    {
      d[i][j] = 0.;
    }
  }

  // matrix multiplication
  for (i = 0; i < fDim; i++)
  {
    for (j = 0; j < fDim; j++)
    {
      for (k = 0; k < fDim; k++)
      {
        d[i][j] += rotationMatrix[i][k] * rotationMatrix[k][j];
      }
    }
  }

  // building the rotation matrix
  for (i = 0; i < fDim; i++)
  {
    for (j = 0; j < fDim; j++)
    {
      rotationMatrix[i][j] = rotationMatrix[i][j] + iMat[i][j] + d[i][j] * (1.0 - c) / std::pow(v.mag(), 2);
    }
  }

  fRotationMatrix = rotationMatrix;

} // end FindRotationMatrix

G4ThreeVector G4TwoPhotonTransition::RotateVector(const G4ThreeVector &vector)
{
  G4int i, j;
  G4double tempVal;
  G4ThreeVector tempVec(0., 0., 0.);

  for (i = 0; i < fDim; i++)
  {
    tempVal = 0.;
    for (j = 0; j < fDim; j++)
    {
      tempVal += fRotationMatrix[i][j] * vector[j];
    }
    tempVec[i] = tempVal;
  }

  return tempVec;

} // end RotateVector

G4ThreeVector G4TwoPhotonTransition::SphericalToCartesian(const G4double &theta, const G4double &phi)
{
  G4ThreeVector res = G4ThreeVector(0., 0., 0.);

  // find trigometric values
  G4double cosTheta = std::cos(theta);
  G4double sinTheta = std::sin(theta);
  G4double cosPhi = std::cos(phi);
  G4double sinPhi = std::sin(phi);

  res.set(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
  return res;

} // end SphericalToCartesian
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
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4LorentzVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4TwoPhotonTransition::G4TwoPhotonTransition()
    : fDim(3), fDirection(0., 0., 0.), fDirectionPhoton0(0., 0., 0.), fDirectionPhoton1(0., 0., 0.), fVerbose(0)
{
  fRotationMatrix.resize(fDim, std::vector<G4double>(fDim, 0.));
}

G4TwoPhotonTransition::~G4TwoPhotonTransition()
{
}

G4FragmentVector *
G4TwoPhotonTransition::SampleTransition(G4Fragment *nucleus,
                                        G4double newExcEnergy,
                                        G4double multipoleRatio,
                                        G4double angularRatio,
                                        G4int shell,
                                        G4bool isGamma,
                                        G4bool isTwoPhoton)
{
  fMultipoleRatio = multipoleRatio;
  fAngularRatio = angularRatio;
  G4Fragment *resultGamma0 = nullptr;
  G4Fragment *resultGamma1 = nullptr;
  G4FragmentVector *results;
  G4double bond_energy = 0.0;

  if (!isGamma)
  {
    if (!isTwoPhoton)
    {
      if (0 <= shell)
      {
        G4int Z = nucleus->GetZ_asInt();
        if (Z <= 100)
        {
          G4int idx = (G4int)shell;
          idx = std::min(idx, G4AtomicShells::GetNumberOfShells(Z) - 1);
          bond_energy = G4AtomicShells::GetBindingEnergy(Z, idx);
        }
      }
    }
  }
  G4double eTrans = nucleus->GetExcitationEnergy() - newExcEnergy - bond_energy;

  if (fVerbose > 2)
  {
    G4cout << "G4TwoPhotonTransisition::GenerateDecayParticles - Etrans(MeV)= "
           << eTrans << "  Eexnew= " << newExcEnergy
           << " Ebond= " << bond_energy << G4endl;
  }
  if (eTrans <= 0.0)
  {
    eTrans += bond_energy;
    bond_energy = 0.0;
  }
  // Do complete Lorentz computation
  G4LorentzVector lv = nucleus->GetMomentum();
  G4ThreeVector bst = lv.boostVector();

  // select secondaries
  G4ParticleDefinition *part0;
  G4ParticleDefinition *part1;

  if (isGamma)
  {
    part0 = G4Gamma::Gamma();
  }
  else if (isTwoPhoton)
  {
    part0 = G4Gamma::Gamma();
    part1 = G4Gamma::Gamma();
  }
  else
  {
    part0 = G4Electron::Electron();
    G4int ne = std::max(nucleus->GetNumberOfElectrons() - 1, 0);
    nucleus->SetNumberOfElectrons(ne);
  }

  if (isTwoPhoton)
  {
    // sets direction and energy sharing of decay
    SampleEnergy(eTrans);
    SampleDirection();

    G4ThreeVector bst = lv.boostVector();

    // Updated momentum calculations
    // first photon 4-vector
    G4LorentzVector photon0FourMom(fGamma0Energy * fDirectionPhoton0.x(),
                                   fGamma0Energy * fDirectionPhoton0.y(),
                                   fGamma0Energy * fDirectionPhoton0.z(), fGamma0Energy);

    // second photon 4-vector
    G4LorentzVector photon1FourMom(fGamma1Energy * fDirectionPhoton1.x(),
                                   fGamma1Energy * fDirectionPhoton1.y(),
                                   fGamma1Energy * fDirectionPhoton1.z(), fGamma1Energy);

    // residual nucleus
    G4double resNucleusX = -fGamma0Energy * fDirectionPhoton0.x() - fGamma1Energy * fDirectionPhoton1.x();
    G4double resNucleusY = -fGamma0Energy * fDirectionPhoton0.y() - fGamma1Energy * fDirectionPhoton1.y();
    G4double resNucleusZ = -fGamma0Energy * fDirectionPhoton0.z() - fGamma1Energy * fDirectionPhoton1.z();

    lv.set(resNucleusX, resNucleusY, resNucleusZ, photon0Energy + photon1Energy);
    // Lab system transform for short lived level
    lv.boost(bst);

    // modified primary fragment
    nucleus->SetExcEnergyAndMomentum(newExcEnergy, lv);

    photon0FourMom.boost(bst);
    photon1FourMom.boost(bst);

    resultGamma0 = new G4Fragment(photon0FourMom, part0);
    resultGamma1 = new G4Fragment(photon1FourMom, part1);

    results->push_back(resultGamma0);
    results->push_back(resultGamma1);
  }
  else
  {
    fDirection = G4RandomDirection();
    G4double mass = nucleus->GetGroundStateMass() + newExcEnergy;
    G4double emass = part0->GetPDGMass();

    // 2-body decay in rest frame
    G4double ecm = lv.mag();
    G4ThreeVector bst = lv.boostVector();
    if (!isGamma)
    {
      ecm += (CLHEP::electron_mass_c2 - bond_energy);
    }

    // G4cout << "Ecm= " << ecm << " mass= " << mass << " emass= " << emass << G4endl;

    ecm = std::max(ecm, mass + emass);
    G4double energy = 0.5 * ((ecm - mass) * (ecm + mass) + emass * emass) / ecm;
    G4double mom = (emass > 0.0) ? std::sqrt((energy - emass) * (energy + emass))
                                 : energy;

    // emitted gamma or e-
    G4LorentzVector res4mom(mom * fDirection.x(),
                            mom * fDirection.y(),
                            mom * fDirection.z(), energy);
    // residual
    energy = std::max(ecm - energy, mass);
    lv.set(-mom * fDirection.x(), -mom * fDirection.y(), -mom * fDirection.z(), energy);

    // Lab system transform for short lived level
    lv.boost(bst);

    // modified primary fragment
    nucleus->SetExcEnergyAndMomentum(newExcEnergy, lv);

    // gamma or e- are produced
    res4mom.boost(bst);
    G4Fragment *result = new G4Fragment(res4mom, part0);
    results->push_back(result);

    // G4cout << " DeltaE= " << e0 - lv.e() - res4mom.e() + emass
    //	 << "   Emass= " << emass << G4endl;
    if (fVerbose > 2)
    {
      G4cout << "G4TwoPhotonTransition::SampleTransition : " << *result << G4endl;
      G4cout << "       Left nucleus: " << *nucleus << G4endl;
    }
  } // end !isTwoPhoton

  return results;
}

/*
if (fVerbose > 2)
{
  G4cout << "G4TwoPhotonTransition::SampleTransition() transitionEnergy= " << eTrans
         << G4endl;
}
if (fVerbose > 2)
{
  G4cout << "G4TwoPhotonTransition::SampleTransition : " << *results << G4endl;
  G4cout << "       Left nucleus: " << *nucleus << G4endl;
}
*/

void G4TwoPhotonTransition::SampleEnergy(G4double transitionEnergy)
{
  // find continuous energy distribution of photons
  SetUpEnergySpectrumSampler(transitionEnergy, fMultipoleRatio);

  if (energySpectrumSampler)
  {

    fGamma0Energy = transitionEnergy * energySpectrumSampler->shoot(G4Random::getTheEngine());
    fGamma1Energy = transitionEnergy - fGamma1Energy; // keV

    if (fVerbose > 2)
    {
      G4cout << "G4TwoPhotonTransition::fGamma0Energy= " << fGamma0Energy << " | "
             << "G4TwoPhotonTransition::fGamma1Energy= " << fGamma1Energy << " | "
             << "G4TwoPhotonTransition::sumEnergy= " << fGamma0Energy + fGamma1Energy
             << "G4TwoPhotonTransition::transitionEnergy= " << transitionEnergy
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
  // * chi could be tuneable, but for now it is hardcoded to simplify the math
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

    fDirectionPhoton0.set(polarizationAxis.getX(), polarizationAxis.getY(), polarizationAxis.getZ());
    fDirectionPhoton1.set(rotatedEmissionVector.getX(), rotatedEmissionVector.getY(), rotatedEmissionVector.getZ());
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
    G4int npti = 200;
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
  G4int npti = 200;
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
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
    : polarFlag(false), fDirection(0., 0., 0.), fVerbose(0)
{
}

G4TwoPhotonTransition::~G4TwoPhotonTransition()
{
}

G4Fragment *
G4TwoPhotonTransition::SampleTransition(G4Fragment *nucleus,
                                        G4double newExcEnergy,
                                        G4double multipoleRatio,
                                        G4double angularRatio)
{
  G4Fragment *resultGamma1 = nullptr;

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

  // Do complete Lorentz computation
  G4LorentzVector lv = nucleus->GetMomentum();
  G4double mass = nucleus->GetGroundStateMass() + newExcEnergy;

  // select secondaries
  G4ParticleDefinition *gamma1 = G4Gamma::Gamma();
  G4ParticleDefinition *gamma2 = G4Gamma::Gamma();

  // find continuous energy distribution of photons
  if (!energySpectrumSampler)
  {
    if (fVerbose > 2)
    {
      G4cout << "## Initializing energy spectrum sampler ##" << G4endl;
    }

    SetUpEnergySpectrumSampler(totalTransEnergy, multipoleRatio);
  }

  if (energySpectrumSampler)
  {
    G4double eGamma1 = totalTransEnergy * energySpectrumSampler->shoot(G4Random::getTheEngine());
    G4double eGamma2 = totalTransEnergy - eGamma1; // keV

    if (fVerbose > 2)
    {
      G4cout << "G4TwoPhotonTransition::eGamma1 " << eGamma1 << " | "
             << "G4TwoPhotonTransition::eGamma2 " << eGamma2 << " | "
             << "G4TwoPhotonTransition::totalTransEnergy " << eGamma1 + eGamma2
             << G4endl;
    }
  }
  else
  {
    G4Exception("G4TwoPhotonTransition::SampleTransition()", "HAD_TWOPHOTON_001", FatalException, "No initialized energy spectrum sampler");
  }

  // finding angular distribution of photons
  if (!angularDistributionSampler)
  {
    if (fVerbose > 2)
    {
      G4cout << "## Initializing angular distribution sampler ##" << G4endl;
    }
    // angular ratio = alpha / chi
    // * there's probably a better way to do this but for now this works
    G4double alphaE1 = angularRatio;
    G4double chi = 1.0;

    SetUpAngularDistributionSampler(alphaE1, chi);
  }

  if (angularDistributionSampler)
  {
    G4double theta = CLHEP::pi * angularDistributionSampler->shoot(G4Random::getTheEngine());

    if (fVerbose > 2)
    {
      G4cout << "###### G4TwoPhotonTransition::angularRatio " << angularRatio << " | "
             << "G4TwoPhotonTransition::theta " << theta
             << G4endl;
    }
  }
  else
  {
    G4Exception("G4TwoPhotonTransition::SampleTransition()", "HAD_TWOPHOTON_002", FatalException, "No initialized angular distribution sampler");
  }

  energySpectrumSampler = NULL;
  angularDistributionSampler = NULL;

  return resultGamma1;
  /*
  if (polarFlag)
  {
    SampleDirection(nucleus, mpRatio, JP1, JP2, MP);
  }
  else
  {
    fDirection = G4RandomDirection();
  }

  G4double emass = part->GetPDGMass();

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
  result = new G4Fragment(res4mom, part);

  // G4cout << " DeltaE= " << e0 - lv.e() - res4mom.e() + emass
  //	 << "   Emass= " << emass << G4endl;
  if (fVerbose > 2)
  {
    G4cout << "G4TwoPhotonTransition::SampleTransition : " << *result << G4endl;
    G4cout << "       Left nucleus: " << *nucleus << G4endl;
  }
  return result;
  */
}

void G4TwoPhotonTransition::SampleDirection(G4Fragment *nuc, G4double ratio,
                                            G4int twoJ1, G4int twoJ2, G4int mp)
{
  G4double cosTheta, phi;
  G4NuclearPolarization *np = nuc->GetNuclearPolarization();
  if (fVerbose > 2)
  {
    G4cout << "G4TwoPhotonTransition::SampleDirection : 2J1= " << twoJ1
           << " 2J2= " << twoJ2 << " ratio= " << ratio
           << " mp= " << mp << G4endl;
    G4cout << "  Nucleus: " << *nuc << G4endl;
  }
  if (nullptr == np)
  {
    cosTheta = 2 * G4UniformRand() - 1.0;
    phi = CLHEP::twopi * G4UniformRand();
  }
  else
  {
    // PhotonEvaporation dataset:
    // The multipolarity number with 1,2,3,4,5,6,7 representing E0,E1,M1,E2,M2,E3,M3
    // monopole transition and 100*Nx+Ny representing multipolarity transition with
    // Ny and Ny taking the value 1,2,3,4,5,6,7 referring to E0,E1,M1,E2,M2,E3,M3,..
    // For example a M1+E2 transition would be written 304.
    // M1 is the primary transition (L) and E2 is the secondary (L')

    G4double mpRatio = ratio;

    G4int L0 = 0, Lp = 0;
    if (mp > 99)
    {
      L0 = mp / 200;
      Lp = (mp % 100) / 2;
    }
    else
    {
      L0 = mp / 2;
      Lp = 0;
      mpRatio = 0.;
    }
    fPolTrans.SampleGammaTransition(np, twoJ1, twoJ2, L0, Lp, mpRatio, cosTheta, phi);
  }

  G4double sinTheta = std::sqrt((1. - cosTheta) * (1. + cosTheta));
  fDirection.set(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
  if (fVerbose > 3)
  {
    G4cout << "G4TwoPhotonTransition::SampleDirection done: " << fDirection << G4endl;
    if (np)
    {
      G4cout << *np << G4endl;
    }
  }
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

  G4double norm = -8 / 3.;                                                        // normalizing factor
  G4double coeff = 4 * alphaE1 * chi / (std::pow(alphaE1, 2) + std::pow(chi, 2)); // interference term coefficient

  for (G4int ptn = 0; ptn <= npti; ptn++)
  {
    // Sample angular range
    cos_theta = 2.0 * (G4double(ptn) + 0.5) / G4double(npti) - 1.0;

    // Build numberical pdf
    // Normalized pdf for pure dipole transition
    // J. Kramp, D. Habs, R. Kroth, M. Music, J. Schirmer, D. Schwalm, and C. Broude, Nuclear Two-Photon Decay in 0+→0+ Transitions, Nuclear Physics, Section A 474, 412 (1987).
    w = 1 + coeff * cos_theta + std::pow(cos_theta, 2);
    f = norm * w;

    pdf[ptn] = f;
  }
  angularDistributionSampler = new G4RandGeneral(pdf, npti);
  delete[] pdf;
}
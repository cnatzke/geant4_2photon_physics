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

#include "G4TwoPhotonEvaporation.hh"

#include "G4NuclearPolarizationStore.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4LorentzVector.hh"
#include "G4FragmentVector.hh"
#include "G4TwoPhotonTransition.hh"
#include "G4Pow.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

G4float G4TwoPhotonEvaporation::GREnergy[] = {0.0f};
G4float G4TwoPhotonEvaporation::GRWidth[] = {0.0f};

#ifdef G4MULTITHREADED
G4Mutex G4TwoPhotonEvaporation::PhotonEvaporationMutex = G4MUTEX_INITIALIZER;
#endif

G4TwoPhotonEvaporation::G4TwoPhotonEvaporation(G4TwoPhotonTransition *p)
    : fLevelManager(nullptr), fTransition(p), fPolarization(nullptr),
      fVerbose(0), fPoints(0), fIndex(0), fMaxLifeTime(DBL_MAX), fRDM(false),
      fSampleTime(true), isInitialised(false)
{
    // G4cout << "### New G4TwoPhotonEvaporation() " << this << G4endl;
    fNuclearLevelData = G4NuclearLevelData::GetInstance();
    fNucPStore = G4NuclearPolarizationStore::GetInstance();
    Tolerance = 20 * CLHEP::eV;

    if (!fTransition)
    {
        fTransition = new G4TwoPhotonTransition();
    }

    theA = theZ = fCode = 0;
    fLevelEnergyMax = fStep = fExcEnergy = fProbability = 0.0;

    for (G4int i = 0; i < MAXDEPOINT; ++i)
    {
        fCummProbability[i] = 0.0;
    }
}

G4TwoPhotonEvaporation::~G4TwoPhotonEvaporation()
{
    delete fTransition;
}

void G4TwoPhotonEvaporation::Initialise()
{
    if (isInitialised)
    {
        return;
    }
    isInitialised = true;

    if (fVerbose > 0)
    {
        G4cout << "### G4TwoPhotonEvaporation is initialized " << this << G4endl;
    }
    G4DeexPrecoParameters *param = fNuclearLevelData->GetParameters();
    Tolerance = param->GetMinExcitation();
    fMaxLifeTime = param->GetMaxLifeTime();

    fTransition->SetVerbose(fVerbose);
}

G4FragmentVector *
G4TwoPhotonEvaporation::EmittedFragments(G4Fragment *nucleus)
{
    if (!isInitialised)
    {
        Initialise();
    }
    fSampleTime = (fRDM) ? false : true;

    if (fVerbose > 1)
    {
        G4cout << "G4TwoPhotonEvaporation::EmittedFragment: "
               << *nucleus << G4endl;
        if (fPolarization)
        {
            G4cout << "NucPolar: " << fPolarization << G4endl;
        }
        G4cout << " RDM: " << fRDM
               << " fPolarization: " << fPolarization << G4endl;
    }
    G4FragmentVector *gammas = GenerateGammas(nucleus);

    // remove G4NuclearPolarizaton when reach ground state
    if (fPolarization && 0 == fIndex)
    {
        if (fVerbose > 1)
        {
            G4cout << "G4TwoPhotonEvaporation::EmittedFragment: remove "
                   << fPolarization << G4endl;
        }
        fNucPStore->RemoveMe(fPolarization);
        fPolarization = nullptr;
        nucleus->SetNuclearPolarization(fPolarization);
    }

    if (fVerbose > 1)
    {
        G4cout << "G4TwoPhotonEvaporation::EmittedFragment: RDM= "
               << fRDM << " done:" << G4endl;
        if (gammas)
        {
            G4cout << gammas << G4endl;
        }
        G4cout << "   Residual: " << *nucleus << G4endl;
    }
    return gammas;
}

G4double
G4TwoPhotonEvaporation::GetFinalLevelEnergy(G4int Z, G4int A, G4double energy)
{
    G4double E = energy;
    InitialiseLevelManager(Z, A);
    if (fLevelManager)
    {
        E = fLevelManager->NearestLevelEnergy(energy, fIndex);
        if (E > fLevelEnergyMax + Tolerance)
        {
            E = energy;
        }
    }
    return E;
}

G4double G4TwoPhotonEvaporation::GetUpperLevelEnergy(G4int Z, G4int A)
{
    InitialiseLevelManager(Z, A);
    return fLevelEnergyMax;
}

G4FragmentVector *
G4TwoPhotonEvaporation::GenerateGammas(G4Fragment *nucleus)
{
    G4FragmentVector *products = new G4FragmentVector();
    std::vector<G4Fragment *> gammas;

    if (!isInitialised)
    {
        Initialise();
    }
    G4double eexc = nucleus->GetExcitationEnergy();
    if (eexc <= Tolerance)
    {
        return products;
    }

    InitialiseLevelManager(nucleus->GetZ_asInt(), nucleus->GetA_asInt());

    G4double time = nucleus->GetCreationTime();

    G4double efinal = 0.0;
    G4int JP1 = 0;
    G4int JP2 = 0;
    G4int multiP = 0;
    G4bool isGamma = true;
    G4bool isDiscrete = false;

    const G4NucLevel *level = nullptr;
    size_t ntrans = 0;

    // initial discrete state
    if (fLevelManager && eexc <= fLevelEnergyMax + Tolerance)
    {
        fIndex = fLevelManager->NearestLevelIndex(eexc, fIndex);
        if (0 < fIndex)
        {
            // for discrete transition
            level = fLevelManager->GetLevel(fIndex);
            if (level)
            {
                ntrans = level->NumberOfTransitions();
                JP1 = fLevelManager->SpinTwo(fIndex);
                if (ntrans > 0)
                {
                    isDiscrete = true;
                }
                else
                {
                    // if no transition available nothing is done for RDM
                    if (fRDM)
                    {
                        return products;
                    }
                    if (fLevelManager->FloatingLevel(fIndex) > 0)
                    {
                        --fIndex;
                        level = fLevelManager->GetLevel(fIndex);
                        ntrans = level->NumberOfTransitions();
                        JP1 = fLevelManager->SpinTwo(fIndex);
                        if (ntrans > 0)
                        {
                            isDiscrete = true;
                        }
                    }
                }
            }
        }
    }
    if (fVerbose > 1)
    {
        G4int prec = G4cout.precision(4);
        G4cout << "GenerateGammas: Z= " << nucleus->GetZ_asInt()
               << " A= " << nucleus->GetA_asInt()
               << " Exc= " << eexc << " Emax= "
               << fLevelEnergyMax << " idx= " << fIndex
               << " fCode= " << fCode
               << " Ntr= " << ntrans << " discrete: " << isDiscrete
               << " fProb= " << fProbability << G4endl;
        G4cout.precision(prec);
    }

    // No continuous decays
    if (!isDiscrete)
    {
        G4cerr << "Two-photon transition does not come from a discrete level" << G4endl;
    }
    else
    {
        if (fVerbose > 1)
        {
            G4cout << "Two-photon emission from level Index= " << fIndex
                   << " Elevel= " << fLevelManager->LevelEnergy(fIndex)
                   << "  RDM= " << fRDM << G4endl;
        }
        if (0 == fIndex || !level)
        {
            return products;
        }

        // stable fragment has life time -1
        // if called from radioactive decay the life time is not checked
        G4double ltime = fLevelManager->LifeTime(fIndex);
        if (ltime < 0.0 || (!fRDM && ltime > fMaxLifeTime))
        {
            return products;
        }

        // CRN chooses which gamma to emit, unneeded for 2-g decay
        size_t idx = 0;
        // it is discrete transition with possible gamma correlation
        isDiscrete = true;
        multiP = level->TransitionType(idx);
        fIndex = level->FinalExcitationIndex(idx);
        JP2 = fLevelManager->SpinTwo(fIndex);

        // final energy and time
        efinal = fLevelManager->LevelEnergy(fIndex);
        if (fSampleTime && ltime > 0.0)
        {
            time -= ltime * G4Log(G4UniformRand());
        }
        nucleus->SetFloatingLevelNumber(fLevelManager->FloatingLevel(fIndex));
    }
    // protection for floating levels
    if (std::abs(efinal - eexc) <= Tolerance)
    {
        return products;
    }

    G4double eTransTotal = std::abs(efinal - eexc);

    // Let transition class handle angular and energy distributions
    gammas = fTransition->SampleTransition(nucleus, efinal, fMultipoleMixing, fAngularRatio);

    if (gammas.size() != 0)
    {
        gammas.at(0)->SetCreationTime(time);
        gammas.at(1)->SetCreationTime(time);
        products->push_back(gammas.at(0));
        products->push_back(gammas.at(1));
    }

    // updated residual nucleus
    nucleus->SetCreationTime(time);
    nucleus->SetSpin(0.5 * JP2);
    if (fPolarization)
    {
        fPolarization->SetExcitationEnergy(efinal);
    }

    // ignore the floating levels with zero energy and create ground state
    if (efinal == 0.0 && fIndex > 0)
    {
        fIndex = 0;
        nucleus->SetFloatingLevelNumber(fLevelManager->FloatingLevel(fIndex));
    }

    if (fVerbose > 1)
    {
        G4cout << "Final level E= " << efinal << " time= " << time
               << " idxFinal= " << fIndex << " isDiscrete: " << isDiscrete
               << " isGamma: " << isGamma << " multiP= " << multiP
               << " JP1= " << JP1 << " JP2= " << JP2 << G4endl;
    }
    if (fVerbose > 1)
    {
        G4cout << G4endl;
        G4cout << "----------------> TWO PHOTON DECAY <---------------" << G4endl;
        // G4cout << "gamma energies: " << eGamma1 / CLHEP::keV << " | " << eGamma2 / CLHEP::keV << G4endl;
        G4cout << "total energy: " << eTransTotal / CLHEP::keV << G4endl;
        G4cout << "fVerbose : " << fVerbose << G4endl;
        G4cout << G4endl;
    }

    return products;
}

void G4TwoPhotonEvaporation::SetUpEnergySpectrumSampler(G4double transitionEnergy)
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
            e = transitionEnergy * G4double(ptn) / G4double(npti);

            // Build numberical pdf

            // Normalized pdf for pure dipole transition
            d = 140 * std::pow(e, 3) * pow((transitionEnergy - e), 3) / std::pow(transitionEnergy, 7);

            // Normalized pdf for pure quadrupole transition
            q = 2772 * std::pow(e, 5) * pow((transitionEnergy - e), 5) / std::pow(transitionEnergy, 11);

            // mixing ratio
            percentDipole = 1. / (1 + std::pow(fMultipoleMixing, 2));
            f = percentDipole * d + (1 - percentDipole) * q;

            pdf[ptn] = f;
        }
        energySpectrumSampler = new G4RandGeneral(pdf, npti);
        delete[] pdf;
    }
}

void G4TwoPhotonEvaporation::SetTwoPhotonTransition(G4TwoPhotonTransition *p)
{
    if (p != fTransition)
    {
        delete fTransition;
        fTransition = p;
    }
}
void G4TwoPhotonEvaporation::RDMForced(G4bool val)
{
    fRDM = val;
}

G4double
G4TwoPhotonEvaporation::GetEmissionProbability(G4Fragment *nucleus)
{
    if (!isInitialised)
    {
        Initialise();
    }
    fProbability = 0.0;
    fExcEnergy = nucleus->GetExcitationEnergy();
    G4int Z = nucleus->GetZ_asInt();
    G4int A = nucleus->GetA_asInt();
    fCode = 1000 * Z + A;
    if (fVerbose > 2)
    {
        G4cout << "G4PhotonEvaporation::GetEmissionProbability: Z="
               << Z << " A=" << A << " Eexc(MeV)= " << fExcEnergy << G4endl;
    }
    // ignore gamma de-excitation for exotic fragments
    // and for very low excitations
    if (0 >= Z || 1 >= A || Z == A || Tolerance >= fExcEnergy)
    {
        return fProbability;
    }

    // ignore gamma de-excitation for highly excited levels
    if (A >= MAXGRDATA)
    {
        A = MAXGRDATA - 1;
    }
    // G4cout<<" GREnergy= "<< GREnergy[A]<<" GRWidth= "<<GRWidth[A]<<G4endl;

    static const G4float GREfactor = 5.0f;
    if (fExcEnergy >= (G4double)(GREfactor * GRWidth[A] + GREnergy[A]))
    {
        return fProbability;
    }
    // probability computed assuming continium transitions
    // VI: continium transition are limited only to final states
    //     below Fermi energy (this approach needs further evaluation)
    G4double emax = std::max(0.0, nucleus->ComputeGroundStateMass(Z, A - 1) + CLHEP::neutron_mass_c2 - nucleus->GetGroundStateMass());

    // max energy level for continues transition
    emax = std::min(emax, fExcEnergy);
    const G4double eexcfac = 0.99;
    if (0.0 == emax || fExcEnergy * eexcfac <= emax)
    {
        emax = fExcEnergy * eexcfac;
    }

    fStep = emax;
    const G4double MaxDeltaEnergy = CLHEP::MeV;
    fPoints = std::min((G4int)(fStep / MaxDeltaEnergy) + 2, MAXDEPOINT);
    fStep /= ((G4double)(fPoints - 1));
    if (fVerbose > 2)
    {
        G4cout << "Emax= " << emax << " Npoints= " << fPoints
               << "  Eex= " << fExcEnergy << G4endl;
    }
    // integrate probabilities
    G4double eres = (G4double)GREnergy[A];
    G4double wres = (G4double)GRWidth[A];
    G4double eres2 = eres * eres;
    G4double wres2 = wres * wres;
    G4double levelDensity = fNuclearLevelData->GetLevelDensity(Z, A, fExcEnergy);
    G4double xsqr = std::sqrt(levelDensity * fExcEnergy);

    G4double egam = fExcEnergy;
    G4double gammaE2 = egam * egam;
    G4double gammaR2 = gammaE2 * wres2;
    G4double egdp2 = gammaE2 - eres2;

    G4double p0 = G4Exp(-2.0 * xsqr) * gammaR2 * gammaE2 / (egdp2 * egdp2 + gammaR2);
    G4double p1(0.0);

    for (G4int i = 1; i < fPoints; ++i)
    {
        egam -= fStep;
        gammaE2 = egam * egam;
        gammaR2 = gammaE2 * wres2;
        egdp2 = gammaE2 - eres2;
        p1 = G4Exp(2.0 * (std::sqrt(levelDensity * std::abs(fExcEnergy - egam)) - xsqr)) * gammaR2 * gammaE2 / (egdp2 * egdp2 + gammaR2);
        fProbability += (p1 + p0);
        fCummProbability[i] = fProbability;
        if (fVerbose > 3)
        {
            G4cout << "Egamma= " << egam << "  Eex= " << fExcEnergy
                   << "  p0= " << p0 << " p1= " << p1 << " sum= "
                   << fCummProbability[i] << G4endl;
        }
        p0 = p1;
    }

    static const G4double NormC = 1.25 * CLHEP::millibarn / (CLHEP::pi2 * CLHEP::hbarc * CLHEP::hbarc);
    fProbability *= fStep * NormC * A;
    if (fVerbose > 1)
    {
        G4cout << "prob= " << fProbability << G4endl;
    }
    return fProbability;
}
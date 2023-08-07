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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4TwoPhotonDecay.hh                                               //
//  Author: C. Natzke                                                         //
//  Based on file: G4TwoPhotonDecay.hh                                        //
//  Date:  2022-04-19                                                         //
//  Description: performs isomeric transition for excited states of           //
//               radioactive nuclei by emitting two gammas                    //
//               and returns daughter particles the rest frame                //
//               of the parent nucleus                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4TwoPhotonDecay.hh"
#include "G4IonTable.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4TwoPhotonEvaporation.hh"
#include "G4RadioactiveDecay.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShells.hh"
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4Fragment.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4TwoPhotonDecay::G4TwoPhotonDecay(const G4ParticleDefinition *theParentNucleus,
                                   const G4double &branch, const G4double &Qvalue,
                                   const G4double &excitationE, G4TwoPhotonEvaporation *aPhotoEvap, const G4String &dataFile)
    : G4NuclearDecay("Two-photon decay", TwoPhoton, excitationE, noFloat), transitionQ(Qvalue), applyARM(true), twoPhotonEvaporation(aPhotoEvap)
{

    SetParent(theParentNucleus); // Store name of parent nucleus, delete G4MT_parent
    SetBR(branch);

    parentZ = theParentNucleus->GetAtomicNumber();
    parentA = theParentNucleus->GetAtomicMass();

    SetNumberOfDaughters(1);
    G4IonTable *theIonTable =
        (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
    SetDaughter(0, theIonTable->GetIon(parentZ, parentA, excitationE, noFloat));

    ReadInTwoPhotonParameters(parentZ, parentA, dataFile);
}

G4TwoPhotonDecay::~G4TwoPhotonDecay()
{
}

G4DecayProducts *G4TwoPhotonDecay::DecayIt(G4double)
{
    // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)
    CheckAndFillParent();

    // Set up final state
    // parentParticle is set at rest here because boost with correct momentum
    // is done later
    G4LorentzVector atRest(G4MT_parent->GetPDGMass(),
                           G4ThreeVector(0., 0., 0.));
    G4DynamicParticle parentParticle(G4MT_parent, atRest);
    G4DecayProducts *products = new G4DecayProducts(parentParticle);

    // Let G4TwoPhotonEvaporation do the decay
    G4Fragment parentNucleus(parentA, parentZ, atRest);

    // twoPhotonEvaporation->SetVerboseLevel(2);
    twoPhotonEvaporation->SetRelativeBR(fBranchingRatio);
    twoPhotonEvaporation->SetMultipoleMixingRatio(fMultipoleMixing);
    twoPhotonEvaporation->SetAngularRatio(fAngularRatio);
    G4FragmentVector *decayProducts = twoPhotonEvaporation->EmittedFragments(&parentNucleus);

    // Modified nuclide is returned as dynDaughter
    G4IonTable *theIonTable =
        (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
    G4ParticleDefinition *daughterIon =
        theIonTable->GetIon(parentZ, parentA, parentNucleus.GetExcitationEnergy(),
                            G4Ions::FloatLevelBase(parentNucleus.GetFloatingLevelNumber()));
    G4DynamicParticle *dynDaughter = new G4DynamicParticle(daughterIon,
                                                           parentNucleus.GetMomentum());

    // Write gammas to products vector
    if (decayProducts)
    {
        if (G4DecayProducts->size() == 2)
        {
            G4DynamicParticle *gamma0Dyn =
                new G4DynamicParticle(decayProducts->at(0)->GetParticleDefinition(),
                                      decayProducts->at(0)->GetMomentum());
            gamma0Dyn->SetProperTime(decayProducts->at(0)->GetCreationTime());
            products->PushProducts(gamma0Dyn);

            G4DynamicParticle *gamma1Dyn =
                new G4DynamicParticle(decayProducts->at(1)->GetParticleDefinition(),
                                      decayProducts->at(1)->GetMomentum());
            gamma1Dyn->SetProperTime(decayProducts->at(1)->GetCreationTime());
            products->PushProducts(gamma1Dyn);

            delete decayProducts;
        }
        else
        {
            G4DynamicParticle *eOrGammaDyn =
                new G4DynamicParticle(decayProducts.at(0)->GetParticleDefinition(),
                                      decayProducts.at(0)->GetMomentum());
            eOrGammaDyn->SetProperTime(decayProducts.at(0)->GetCreationTime());
            products->PushProducts(eOrGammaDyn);
            delete decayProducts;

            // Now do atomic relaxation if e- is emitted
            if (applyARM)
            {
                G4int shellIndex = twoPhotonEvaporation->GetVacantShellNumber();
                if (shellIndex > -1)
                {
                    G4VAtomDeexcitation *atomDeex =
                        G4LossTableManager::Instance()->AtomDeexcitation();
                    if (atomDeex->IsFluoActive() && parentZ > 5 && parentZ < 100)
                    {
                        G4int nShells = G4AtomicShells::GetNumberOfShells(parentZ);
                        if (shellIndex >= nShells)
                            shellIndex = nShells;
                        G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIndex);
                        const G4AtomicShell *shell = atomDeex->GetAtomicShell(parentZ, as);
                        std::vector<G4DynamicParticle *> armProducts;

                        // VI, SI
                        // Allows fixing of Bugzilla 1727
                        G4double deexLimit = 0.1 * keV;
                        if (G4EmParameters::Instance()->DeexcitationIgnoreCut())
                            deexLimit = 0.;
                        //

                        atomDeex->GenerateParticles(&armProducts, shell, parentZ, deexLimit,
                                                    deexLimit);
                        G4double productEnergy = 0.;
                        for (G4int i = 0; i < G4int(armProducts.size()); i++)
                            productEnergy += armProducts[i]->GetKineticEnergy();

                        G4double deficit = shell->BindingEnergy() - productEnergy;
                        if (deficit > 0.0)
                        {
                            // Add a dummy electron to make up extra energy
                            G4double cosTh = 1. - 2. * G4UniformRand();
                            G4double sinTh = std::sqrt(1. - cosTh * cosTh);
                            G4double phi = twopi * G4UniformRand();

                            G4ThreeVector electronDirection(sinTh * std::sin(phi),
                                                            sinTh * std::cos(phi), cosTh);
                            G4DynamicParticle *extra =
                                new G4DynamicParticle(G4Electron::Electron(), electronDirection,
                                                      deficit);
                            armProducts.push_back(extra);
                        }

                        G4int nArm = armProducts.size();
                        if (nArm > 0)
                        {
                            G4ThreeVector bst = dynDaughter->Get4Momentum().boostVector();
                            for (G4int i = 0; i < nArm; ++i)
                            {
                                G4DynamicParticle *dp = armProducts[i];
                                G4LorentzVector lv = dp->Get4Momentum().boost(bst);
                                dp->Set4Momentum(lv);
                                products->PushProducts(dp);
                            }
                        }
                    }
                }
            } // if ARM on
        }     // end gamma/IC decay
    }

    products->PushProducts(dynDaughter);

    // Energy conservation check
    /*
    G4int newSize = products->entries();
    G4DynamicParticle* temp = 0;
    G4double KEsum = 0.0;
    for (G4int i = 0; i < newSize; i++) {
        temp = products->operator[](i);
        KEsum += temp->GetKineticEnergy();
    }
    G4double eCons = G4MT_parent->GetPDGMass() - dynDaughter->GetMass() - KEsum;
    G4cout << " Two-photon check: Ediff (keV) = " << eCons / CLHEP::keV << G4endl;
    */

    return products;
}

void G4TwoPhotonDecay::ReadInTwoPhotonParameters(G4int Z, G4int A, const G4String &filename)
{

    std::ifstream infile(filename, std::ios::in);

    // safety check
    if (!infile.is_open())
    {
        G4ExceptionDescription ed;
        ed << "User file for Z= " << Z << " A= " << A
           << " is not opened!";
        G4Exception("G4TwoPhotonDecay::ReadInTwoPhotonParameters(..)", "had014",
                    FatalException, ed, "");
        // return nullptr;
    }
    // file is open
    else
    {

        G4String line;
        // read in header line
        std::getline(infile, line);
        // read in data lines
        while (std::getline(infile, line))
        {
            std::istringstream sstr(line);

            sstr >> fLevelIndex >> fLevelEnergy >> branchingRatio >> fMultipoleMixing >> fAngularRatio;

            if (fVerbose > 1)
            {
                G4cout << "fLevelIndex " << fLevelIndex
                       << " | fLevelEnergy " << fLevelEnergy
                       << " | fMultipoleMixing " << fMultipoleMixing
                       << " | fAngularRatio " << fAngularRatio
                       << G4endl;
            }
        }

        if (fVerbose > 1)
        {
            G4cout << "G4TwoPhotonDecay: open file for Z= "
                   << Z << " A= " << A << G4endl;
        }
    }
    infile.close();
}

void G4TwoPhotonDecay::DumpNuclearInfo()
{
    G4cout << " G4TwoPhotonDecay for parent nucleus " << GetParentName() << G4endl;
    G4cout << " decays to " << GetDaughterName(0)
           << " + two gammas, with branching ratio " << GetBR()
           << "% and Q value " << transitionQ << G4endl;
}

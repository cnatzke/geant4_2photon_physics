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
                                   const G4double &excitationE, G4TwoPhotonEvaporation *aPhotoEvap)
    : G4NuclearDecay("Two-photon decay", TwoPhoton, excitationE, noFloat), transitionQ(Qvalue), twoPhotonEvaporation(aPhotoEvap)
{
    SetParent(theParentNucleus); // Store name of parent nucleus, delete G4MT_parent
    SetBR(branch);

    parentZ = theParentNucleus->GetAtomicNumber();
    parentA = theParentNucleus->GetAtomicMass();

    SetNumberOfDaughters(1);
    G4IonTable *theIonTable =
        (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
    SetDaughter(0, theIonTable->GetIon(parentZ, parentA, excitationE, noFloat));
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
    G4FragmentVector *emittedGammas = twoPhotonEvaporation->EmittedFragments(&parentNucleus);

    // Modified nuclide is returned as dynDaughter
    G4IonTable *theIonTable =
        (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
    G4ParticleDefinition *daughterIon =
        theIonTable->GetIon(parentZ, parentA, parentNucleus.GetExcitationEnergy(),
                            G4Ions::FloatLevelBase(parentNucleus.GetFloatingLevelNumber()));
    G4DynamicParticle *dynDaughter = new G4DynamicParticle(daughterIon,
                                                           parentNucleus.GetMomentum());

    // Write gammas to products vector
    if (emittedGammas)
    {
        for (size_t i = 0; i < emittedGammas->size(); i++)
        {
            G4DynamicParticle *emittedGammaDyn = new G4DynamicParticle(emittedGammas->at(i)->GetParticleDefinition(), emittedGammas->at(i)->GetMomentum());
            emittedGammaDyn->SetProperTime(emittedGammas->at(i)->GetCreationTime());
            products->PushProducts(emittedGammaDyn);
        }
        delete emittedGammas;
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

void G4TwoPhotonDecay::DumpNuclearInfo()
{
    G4cout << " G4TwoPhotonDecay for parent nucleus " << GetParentName() << G4endl;
    G4cout << " decays to " << GetDaughterName(0)
           << " + two gammas, with branching ratio " << GetBR()
           << "% and Q value " << transitionQ << G4endl;
}

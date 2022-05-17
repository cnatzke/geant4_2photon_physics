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

#ifndef G4TwoPhotonDecay_h
#define G4TwoPhotonDecay_h 1

#include "G4NuclearDecay.hh"
#include "G4Fragment.hh"

class G4TwoPhotonEvaporation;

class G4TwoPhotonDecay : public G4NuclearDecay
{
public:
G4TwoPhotonDecay(const G4ParticleDefinition* theParentNucleus,
          const G4double& theBR, const G4double& Qvalue,
          const G4double& excitation, G4TwoPhotonEvaporation* aPhotonEvap);

virtual ~G4TwoPhotonDecay();

virtual G4DecayProducts* DecayIt(G4double);

virtual void DumpNuclearInfo();

private:
const G4double transitionQ;
G4int parentZ;
G4int parentA;

G4TwoPhotonEvaporation* twoPhotonEvaporation;
};

#endif

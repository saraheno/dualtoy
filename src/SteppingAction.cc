#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "DetectorConstruction.hh"
#include "TString.h"
#include "TRandom3.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"



using namespace std;
using namespace CLHEP;

int to_int(string name)
{
  int Result; // int which will contain the result
  stringstream convert(name);
  string dummy;
  convert >> dummy;
  convert >> Result;
  return Result;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::SteppingAction(DetectorConstruction *detectorConstruction,
                               const G4int &scint, const G4int &cher) : fDetectorConstruction(detectorConstruction),
                                                                        propagateScintillation(scint),
                                                                        propagateCerenkov(cher)
{
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::~SteppingAction()
{
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void SteppingAction::UserSteppingAction(const G4Step *theStep)
{

  G4Track *theTrack = theStep->GetTrack();



  G4ParticleDefinition *particleType = theTrack->GetDefinition();


  G4StepPoint *thePrePoint = theStep->GetPreStepPoint();
  G4StepPoint *thePostPoint = theStep->GetPostStepPoint();
  const G4ThreeVector &thePrePosition = thePrePoint->GetPosition();
  G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
  G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
  G4String thePrePVName = "";
  if (thePrePV)
    thePrePVName = thePrePV->GetName();
  G4String thePostPVName = "";
  if (thePostPV)
    thePostPVName = thePostPV->GetName();

  //  G4VSolid* thePreS = thePrePV->GetLogicalVolume()->GetSolid();

  G4int nStep = theTrack->GetCurrentStepNumber();

  G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();

  //        cout << " step length = " << theStep->GetStepLength() << endl;
  //-------------

  // get position
  G4double global_x = thePrePosition.x() / mm;
  G4double global_y = thePrePosition.y() / mm;
  G4double global_z = thePrePosition.z() / mm;

  G4double energy = theStep->GetTotalEnergyDeposit();
  G4double energyIon = energy - theStep->GetNonIonizingEnergyDeposit();
  G4double energyElec = 0.;
  //total energy by particle types
  G4double energyPion_n = 0.;
  G4double energyPositron = 0.;
  G4double energyElectron = 0.;
  G4double energyPhoton = 0.;
  G4double energyPion_p = 0.;
  G4double energyKion = 0.;
  G4double energyNeutron = 0.;
  G4double energyProton = 0.;
  //ion energy by particle types
  G4double energyIonPion_n = 0.;
  G4double energyIonPositron = 0.;
  G4double energyIonElectron = 0.;
  G4double energyIonPhoton = 0.;
  G4double energyIonPion_p = 0.;
  G4double energyIonKion = 0.;
  G4double energyIonNeutron = 0.;
  G4double energyIonProton = 0.;

  if (TrPDGid == (-211))
  {
    energyPion_n = energy;
    energyIonPion_n = energyIon;
  }
  else if (TrPDGid == (-11))
  {
    energyPositron = energy;
    energyIonPositron = energyIon;
  }
  else if (TrPDGid == (11))
  {
    energyElectron = energy;
    energyIonElectron = energyIon;
  }
  else if (TrPDGid == (22))
  {
    energyPhoton = energy;
    energyIonPhoton = energyIon;
  }
  else if (TrPDGid == (211))
  {
    energyPion_p = energy;
    energyIonPion_p = energyIon;
  }
  else if (TrPDGid == (321))
  {
    energyKion = energy;
    energyIonKion = energyIon;
  }
  else if (TrPDGid == (2112))
  {
    energyNeutron = energy;
    energyIonNeutron = energyIon;
  }
  else if (TrPDGid == (2212))
  {
    energyProton = energy;
    energyIonProton = energyIon;
  }
  energyElec = energyIonPositron + energyIonElectron;

  //std::cout<<"TrPDGid energy energyIon enegyElec are "<<TrPDGid<<" "<<energy<<" "<<energyIon<<" "<<energyElec<<std::endl;

  CreateTree::Instance()->depositedEnergyTotal += energy / GeV;
  CreateTree::Instance()->depositedIonEnergyTotal += energyIon / GeV;
  CreateTree::Instance()->depositedElecEnergyTotal += energyElec / GeV;

  //    if(thePrePVName.contains("world")) {
  bool haha4 = ((theStep->GetPostStepPoint())->GetStepStatus()) == fWorldBoundary;
  if (haha4)
  {
    //std::cout<<"leaving "<<std::endl;
    CreateTree::Instance()->depositedEnergyEscapeWorld += (theStep->GetPostStepPoint())->GetKineticEnergy() / GeV;
  }
  //}

  // optical photon

  if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
  {

    G4String processName = theTrack->GetCreatorProcess()->GetProcessName();

    if (
        (nStep == 1) && (processName == "Cerenkov"))
    {


      TrackInformation *theTrackInfo = (TrackInformation *)(theTrack->GetUserInformation());
      //G4ThreeVector haha=theTrackInfo->GetParentMomentum();
      //G4double haha2=theTrackInfo->GetParentEnergy()/GeV;
      //G4double haha3=haha.mag()/GeV;
      //G4double betaa=0.;
      //if(haha2>0) betaa=haha3/haha2;

      G4int aapdgid = theTrackInfo->GetParentPDGid();


      //std::cout << " generated Cerenkov photon with parent " << theTrackInfo->GetParentName()<<" "<<aapdgid<<" with beta of "<<betaa<<" and energy "<<haha2<<std::endl;
      float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy() / eV);

      //kill very long wavelengths
      if (photWL > 1000 || photWL < 300)
        theTrack->SetTrackStatus(fKillTrackAndSecondaries);



      else if (thePrePVName.contains("ecalCrystalP_f_act2_cheren"))
      {
	
        CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_total += 1;
	//	std::cout<<"num cher "<<CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_total<<std::endl;
        if (aapdgid == (-211))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_particleID[0] += 1;
        if (aapdgid == (-11))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_particleID[1] += 1;
        if (aapdgid == (11))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_particleID[2] += 1;
        if (aapdgid == (22))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_particleID[3] += 1;
        if (aapdgid == (211))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_particleID[4] += 1;
        if (aapdgid == (321))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_particleID[5] += 1;
        if (aapdgid == (2112))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_particleID[6] += 1;
        if (aapdgid == (2212))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_particleID[7] += 1;

      }

      else if (thePrePVName.contains("ecalCrystalP_f_act1_scinti"))
      {
        CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_total += 1;
        if (aapdgid == (-211))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_particleID[0] += 1;
        if (aapdgid == (-11))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_particleID[1] += 1;
        if (aapdgid == (11))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_particleID[2] += 1;
        if (aapdgid == (22))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_particleID[3] += 1;
        if (aapdgid == (211))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_particleID[4] += 1;
        if (aapdgid == (321))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_particleID[5] += 1;
        if (aapdgid == (2112))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_particleID[6] += 1;
        if (aapdgid == (2212))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_f_particleID[7] += 1;

      }

      if (!propagateCerenkov)
        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }
  }

  else
  {



    //cal

    if (thePrePVName.contains("ecalCrystalP_f"))
    {

      if (thePrePVName.contains("ecalCrystalP_f_absorb"))
      {

        CreateTree::Instance()->depositedEnergyECAL_f[0] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_f[0] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_f[0] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[5] += energyKion / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[7] += energyIonProton / GeV;
      }

      if (thePrePVName.contains("ecalCrystalP_f_act1_scinti"))
      {
        CreateTree::Instance()->depositedEnergyECAL_f[1] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_f[1] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_f[1] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[5] += energyKion / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[7] += energyIonProton / GeV;

      }

      if (thePrePVName.contains("ecalCrystalP_f_act2_cheren"))
      {
        CreateTree::Instance()->depositedEnergyECAL_f[2] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_f[2] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_f[2] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[5] += energyKion / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[7] += energyIonProton / GeV;


      }
    }


    //    if (thePrePVName.contains("ecalCrystalL_f_absorb") || thePrePVName.contains("ecalCrystalL_r_absorb"))

    if (thePrePVName.contains("solenoid"))
    {
      CreateTree::Instance()->depositedEnergySolenoid += energy / GeV;
      CreateTree::Instance()->depositedIonEnergySolenoid += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergySolenoid += energyElec / GeV;
    }

    //G4cout << ">>> end non optical photon" << G4endl;
  } // non optical photon

  return;
}

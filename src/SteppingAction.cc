#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "DetectorConstruction.hh"
#include "TString.h"
#include "TRandom3.h"
//#include "TCint.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"

//long int CreateSeed();

using namespace std;
using namespace CLHEP;
/*
SteppingAction::SteppingAction (const string& configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  ConfigFile config (configFileName) ;

  config.readInto(core_material, "core_material"); 

  if (core_material == 0)
  {
	  config.readInto(toy_ly,	"toy_ly");
	  config.readInto(toy_decay,	"toy_decay");
	  config.readInto(toy_rise,	"toy_rise");
  }
  

}*/

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

  //  const G4ThreeVector& theTrackDirection = theTrack->GetMomentumDirection();
  //  const G4ThreeVector& theTrackVertexDirection = theTrack->GetVertexMomentumDirection();

  //  TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformation());

  G4ParticleDefinition *particleType = theTrack->GetDefinition();
  //G4int trackID = theTrack->GetTrackID();

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
      CreateTree::Instance()->h_phot_cer_parentID->Fill(aapdgid);

      //std::cout << " generated Cerenkov photon with parent " << theTrackInfo->GetParentName()<<" "<<aapdgid<<" with beta of "<<betaa<<" and energy "<<haha2<<std::endl;
      float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy() / eV);

      //kill very long wavelengths
      if (photWL > 1000 || photWL < 300)
        theTrack->SetTrackStatus(fKillTrackAndSecondaries);

      else if (thePrePVName.contains("corePV_front"))
      {
        CreateTree::Instance()->tot_phot_cer_Timing_f_total += 1;
        CreateTree::Instance()->h_phot_cer_lambda_Timing_f->Fill(photWL);
      }
      else if (thePrePVName.contains("corePV_rear"))
      {
        CreateTree::Instance()->tot_phot_cer_Timing_r_total += 1;
        CreateTree::Instance()->h_phot_cer_lambda_Timing_r->Fill(photWL);
      }

      else if (thePrePVName.contains("ecalCrystalP_f_fiber_cheren"))
      {
        CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_total += 1;
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
        CreateTree::Instance()->h_phot_cer_lambda_ECAL_f->Fill(photWL);
      }
      else if (thePrePVName.contains("ecalCrystalP_r_fiber_cheren"))
      {
        CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_total += 1;
        if (aapdgid == (-211))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_particleID[0] += 1;
        if (aapdgid == (-11))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_particleID[1] += 1;
        if (aapdgid == (11))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_particleID[2] += 1;
        if (aapdgid == (22))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_particleID[3] += 1;
        if (aapdgid == (211))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_particleID[4] += 1;
        if (aapdgid == (321))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_particleID[5] += 1;
        if (aapdgid == (2112))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_particleID[6] += 1;
        if (aapdgid == (2212))
          CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_particleID[7] += 1;
        CreateTree::Instance()->h_phot_cer_lambda_ECAL_r->Fill(photWL);
      }
      else if (thePrePVName.contains("ecalCrystalP_f_fiber_scinti"))
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
        CreateTree::Instance()->h_phot_cer_lambda_ECAL_f->Fill(photWL);
      }
      else if (thePrePVName.contains("ecalCrystalP_r_fiber_scinti"))
      {
        CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_total += 1;
        if (aapdgid == (-211))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_particleID[0] += 1;
        if (aapdgid == (-11))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_particleID[1] += 1;
        if (aapdgid == (11))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_particleID[2] += 1;
        if (aapdgid == (22))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_particleID[3] += 1;
        if (aapdgid == (211))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_particleID[4] += 1;
        if (aapdgid == (321))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_particleID[5] += 1;
        if (aapdgid == (2112))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_particleID[6] += 1;
        if (aapdgid == (2212))
          CreateTree::Instance()->tot_phot_cer_ECAL_scinti_r_particleID[7] += 1;
        CreateTree::Instance()->h_phot_cer_lambda_ECAL_r->Fill(photWL);
      }
      else if (thePrePVName.contains("hcalTile_layer"))
      {
        CreateTree::Instance()->tot_phot_cer_HCAL += 1;
        CreateTree::Instance()->h_phot_cer_lambda_HCAL->Fill(photWL);
      }

      if (!propagateCerenkov)
        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }
  }

  else
  {

    //count tracks before SCEPCAL at the tracker layers
    if (thePrePVName.contains("world") && thePostPVName.contains("trackerPV_Layer")) // interface between T1 and T2
    {
      for (int iLayer = 0; iLayer < 6; iLayer++)
      {
        if (thePostPVName == Form("trackerPV_Layer%d", iLayer)) // interface between T1 and T2
          CreateTree::Instance()->nTracksTRK[iLayer]++;         //counting tracks crossing the boundary
      }
    }

    //save primary particle position and momentum before timing layer T1 and before ECAL E1
    else if (thePrePVName.contains("world") && thePostPVName.contains("corePV_front")) // interface between world and T1
    {

      CreateTree::Instance()->nTracksT1++; //counting tracks crossing the boundary

      if (theTrack->GetParentID() == 0) // select only the primary particle
      {
        CreateTree::Instance()->primaryPosT1->at(0) = global_x;
        CreateTree::Instance()->primaryPosT1->at(1) = global_y;
        CreateTree::Instance()->primaryPosT1->at(2) = global_z;

        CreateTree::Instance()->primaryMomT1->at(0) = thePrePoint->GetMomentum().x() / GeV;
        CreateTree::Instance()->primaryMomT1->at(1) = thePrePoint->GetMomentum().y() / GeV;
        CreateTree::Instance()->primaryMomT1->at(2) = thePrePoint->GetMomentum().z() / GeV;
        CreateTree::Instance()->primaryMomT1->at(3) = thePrePoint->GetTotalEnergy() / GeV;
      }
    }

    else if (thePrePVName.contains("world") && thePostPVName.contains("corePV_rear")) // interface between T1 and T2
    {
      CreateTree::Instance()->nTracksT2++; //counting tracks crossing the boundary
    }

    else if ((thePrePVName.contains("world") || thePrePVName.contains("ecalGapP_f") || thePrePVName.contains("ecalDetP_f")) && thePostPVName.contains("ecalCrystalP_f") // interface between world and E1
    )
    {
      CreateTree::Instance()->nTracksE1++; //counting tracks crossing the boundary

      if (theTrack->GetParentID() == 0) // select only the primary particle
      {
        CreateTree::Instance()->primaryPosE1->at(0) = global_x;
        CreateTree::Instance()->primaryPosE1->at(1) = global_y;
        CreateTree::Instance()->primaryPosE1->at(2) = global_z;

        CreateTree::Instance()->primaryMomE1->at(0) = thePrePoint->GetMomentum().x() / GeV;
        CreateTree::Instance()->primaryMomE1->at(1) = thePrePoint->GetMomentum().y() / GeV;
        CreateTree::Instance()->primaryMomE1->at(2) = thePrePoint->GetMomentum().z() / GeV;
        CreateTree::Instance()->primaryMomE1->at(3) = thePrePoint->GetTotalEnergy() / GeV;
      }
    }

    else if ((thePrePVName.contains("ecalCrystalP_f") || thePrePVName.contains("world")) && thePostPVName.contains("ecalCrystalP_r")) // interface between E1 and E2
    {
      CreateTree::Instance()->nTracksE2++; //counting tracks crossing the boundary
    }

    //tracker
    if (thePrePVName.contains("trackerPV_Layer"))
    {
      for (int iLayer = 0; iLayer < 6; iLayer++)
      {
        if (thePrePVName == Form("trackerPV_Layer%d", iLayer))
          CreateTree::Instance()->Edep_Tracker_layer[iLayer] += energy / GeV;
      }
    }

    //timing
    if (thePrePVName.contains("corePV_front"))
    {
      CreateTree::Instance()->depositedEnergyTiming_f += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTiming_f += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTiming_f += energyElec / GeV;
      for (int iBar = 0; iBar < 18; iBar++)
      {
        if (thePrePVName == Form("corePV_front_%d", iBar))
          CreateTree::Instance()->Edep_Timing_f_ch[iBar] += energy / GeV;
      }
    }
    if (thePrePVName.contains("corePV_rear"))
    {
      CreateTree::Instance()->depositedEnergyTiming_r += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTiming_r += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTiming_r += energyElec / GeV;
      for (int iBar = 0; iBar < 18; iBar++)
      {
        if (thePrePVName == Form("corePV_rear_%d", iBar))
          CreateTree::Instance()->Edep_Timing_r_ch[iBar] += energy / GeV;
      }
    }

    //ecal
    int chn = 64;
    if (thePrePVName.contains("ecalCrystalP_f"))
    {
      for (int iCh = 0; iCh < chn; iCh++)
      {
        if (thePrePVName.contains(Form("_%d", iCh)))
        {
          CreateTree::Instance()->Edep_ECAL_f_ch[iCh] += energy / GeV;
          CreateTree::Instance()->IonEdep_ECAL_f_ch[iCh] += energyIon / GeV;
          CreateTree::Instance()->ElecEdep_ECAL_f_ch[iCh] += energyElec / GeV;
        }
      }

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

      if (thePrePVName.contains("ecalCrystalP_f_fiber_scinti"))
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

        for (int iCh = 0; iCh < chn; iCh++)
        {
          if (thePrePVName == Form("ecalCrystalP_f_fiber_scinti_%d", iCh))
          {
            CreateTree::Instance()->Edep_ECAL_f_scinti_ch[iCh] += energy / GeV;
            CreateTree::Instance()->IonEdep_ECAL_f_scinti_ch[iCh] += energyIon / GeV;
            CreateTree::Instance()->ElecEdep_ECAL_f_scinti_ch[iCh] += energyElec / GeV;
          }
        }
      }

      if (thePrePVName.contains("ecalCrystalP_f_fiber_cheren"))
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

        for (int iCh = 0; iCh < chn; iCh++)
        {
          if (thePrePVName == Form("ecalCrystalP_f_fiber_cherenp_%d", iCh) || thePrePVName == Form("ecalCrystalP_f_fiber_cherenc_%d", iCh))
          {
            CreateTree::Instance()->Edep_ECAL_f_cheren_ch[iCh] += energy / GeV;
            CreateTree::Instance()->IonEdep_ECAL_f_cheren_ch[iCh] += energyIon / GeV;
            CreateTree::Instance()->ElecEdep_ECAL_f_cheren_ch[iCh] += energyElec / GeV;
          }
        }
      }
    }

    if (thePrePVName.contains("ecalCrystalP_r"))
    {
      for (int iCh = 0; iCh < chn; iCh++)
      {
        if (thePrePVName.contains(Form("_%d", iCh)))
        {
          CreateTree::Instance()->Edep_ECAL_r_ch[iCh] += energy / GeV;
          CreateTree::Instance()->IonEdep_ECAL_r_ch[iCh] += energyIon / GeV;
          CreateTree::Instance()->ElecEdep_ECAL_r_ch[iCh] += energyElec / GeV;
        }
      }

      if (thePrePVName.contains("ecalCrystalP_r_absorb"))
      {
        CreateTree::Instance()->depositedEnergyECAL_r[0] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_r[0] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_r[0] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_r_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_r_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_r_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_r_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_r_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[5] += energyKion / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_r_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_r_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_r_particleID[7] += energyIonProton / GeV;
      }

      if (thePrePVName.contains("ecalCrystalP_r_fiber_scinti"))
      {
        CreateTree::Instance()->depositedEnergyECAL_r[1] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_r[1] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_r[1] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[5] += energyKion / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[7] += energyIonProton / GeV;
        for (int iCh = 0; iCh < chn; iCh++)
        {
          if (thePrePVName == Form("ecalCrystalP_r_fiber_scinti_%d", iCh))
          {
            CreateTree::Instance()->Edep_ECAL_r_scinti_ch[iCh] += energy / GeV;
            CreateTree::Instance()->IonEdep_ECAL_r_scinti_ch[iCh] += energyIon / GeV;
            CreateTree::Instance()->ElecEdep_ECAL_r_scinti_ch[iCh] += energyElec / GeV;
          }
        }
      }

      if (thePrePVName.contains("ecalCrystalP_r_fiber_cheren"))
      {
        CreateTree::Instance()->depositedEnergyECAL_r[2] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_r[2] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_r[2] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[5] += energyKion / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[7] += energyIonProton / GeV;
        for (int iCh = 0; iCh < chn; iCh++)
        {
          if (thePrePVName == Form("ecalCrystalP_r_fiber_cherenp_%d", iCh) || thePrePVName == Form("ecalCrystalP_r_fiber_cherenc_%d", iCh))
          {
            CreateTree::Instance()->Edep_ECAL_r_cheren_ch[iCh] += energy / GeV;
            CreateTree::Instance()->IonEdep_ECAL_r_cheren_ch[iCh] += energyIon / GeV;
            CreateTree::Instance()->ElecEdep_ECAL_r_cheren_ch[iCh] += energyElec / GeV;
          }
        }
      }
    }
    //    if (thePrePVName.contains("ecalCrystalL_f_absorb") || thePrePVName.contains("ecalCrystalL_r_absorb"))
    if (thePrePVName.contains("ecalCrystalP_r_fiber_scinti") || thePrePVName.contains("ecalCrystalP_f_fiber_scinti"))
    {

      CreateTree::Instance()->h_parentID_Edep->Fill(TrPDGid, energy / GeV);

      CreateTree::Instance()->h_time_z_st->Fill(global_z / mm, thePrePoint->GetGlobalTime() / ns - 20 / 3., energy / GeV);
      CreateTree::Instance()->h_time_z_mt->Fill(global_z / mm, thePrePoint->GetGlobalTime() / ns - 20 / 3., energy / GeV);
      CreateTree::Instance()->h_time_z_lt->Fill(global_z / mm, thePrePoint->GetGlobalTime() / ns - 20 / 3., energy / GeV);

      CreateTree::Instance()->h_localtime_z_st->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      CreateTree::Instance()->h_localtime_z_mt->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      CreateTree::Instance()->h_localtime_z_lt->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);

      if (TrPDGid == (-11))
      {
        CreateTree::Instance()->h_localtime_z_st_ep->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_mt_ep->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_lt_ep->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      }
      if (TrPDGid == (11))
      {
        CreateTree::Instance()->h_localtime_z_st_en->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_mt_en->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_lt_en->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      }
      if (TrPDGid == (22))
      {
        CreateTree::Instance()->h_localtime_z_st_gamma->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_mt_gamma->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_lt_gamma->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      }
      if (TrPDGid == (2112))
      {
        CreateTree::Instance()->h_localtime_z_st_n->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_mt_n->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_lt_n->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      }
      if (TrPDGid == (2212))
      {
        CreateTree::Instance()->h_localtime_z_st_p->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_mt_p->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_lt_p->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      }

      if (TrPDGid == (-211))
      {
        CreateTree::Instance()->h_localtime_z_st_pin->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_mt_pin->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_lt_pin->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      }
      if (TrPDGid == (211))
      {
        CreateTree::Instance()->h_localtime_z_st_pip->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_mt_pip->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
        CreateTree::Instance()->h_localtime_z_lt_pip->Fill(global_z / mm, (thePrePoint->GetGlobalTime() / ns) - global_z / mm / 300. - 19 / 3., energy / GeV);
      }
    }

    //hcal
    if (thePrePVName.contains("hcalTile_layer"))
    {
      CreateTree::Instance()->depositedEnergyHCALAct += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyHCALAct += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyHCALAct += energyElec / GeV;
    }
    if (thePrePVName.contains("hcalAbs"))
    {
      CreateTree::Instance()->depositedEnergyHCALPas += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyHCALPas += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyHCALPas += energyElec / GeV;
    }

    if (thePrePVName.contains("world"))
    {
      CreateTree::Instance()->depositedEnergyWorld += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyWorld += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyWorld += energyElec / GeV;
    }

    if (thePrePVName.contains("services"))
    {
      CreateTree::Instance()->depositedEnergyServices += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyServices += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyServices += energyElec / GeV;
    }

    if (thePrePVName.contains("TimingGap"))
    {
      CreateTree::Instance()->depositedEnergyTimingGap += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTimingGap += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTimingGap += energyElec / GeV;
    }

    if (thePrePVName.contains("ecalGap"))
    {
      CreateTree::Instance()->depositedEnergyEcalGap += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyEcalGap += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyEcalGap += energyElec / GeV;
    }

    if (thePrePVName.contains("ecalDet"))
    {
      CreateTree::Instance()->depositedEnergyEcalDet += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyEcalDet += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyEcalDet += energyElec / GeV;
    }

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

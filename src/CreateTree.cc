#include "CreateTree.hh"
#include <algorithm>

using namespace std;

CreateTree *CreateTree::fInstance = NULL;

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

CreateTree::CreateTree(TString name)
{
  if (fInstance)
  {
    return;
  }

  this->fInstance = this;
  this->fname = name;
  this->ftree = new TTree(name, name);

  this->GetTree()->Branch("Event", &this->Event, "Event/I");


  this->GetTree()->Branch("inputE1Thick", &this->inputE1Thick, "inputE1Thick/F");
  this->GetTree()->Branch("inputE2Thick", &this->inputE2Thick, "inputE2Thick/F");
  this->GetTree()->Branch("inputE1Width", &this->inputE1Width, "inputE1Width/F");


  inputInitialPosition = new vector<float>(3, 0.);
  inputMomentum = new vector<float>(4, 0.);
  primaryPosT1 = new vector<float>(3, 0.);
  primaryMomT1 = new vector<float>(4, 0.);
  primaryPosE1 = new vector<float>(3, 0.);
  primaryMomE1 = new vector<float>(4, 0.);

  this->GetTree()->Branch("inputInitialPosition", "vector<float>", &inputInitialPosition);
  this->GetTree()->Branch("inputMomentum", "vector<float>", &inputMomentum);
  this->GetTree()->Branch("primaryPosT1", "vector<float>", &primaryPosT1);
  this->GetTree()->Branch("primaryMomT1", "vector<float>", &primaryMomT1);
  this->GetTree()->Branch("primaryPosE1", "vector<float>", &primaryPosE1);
  this->GetTree()->Branch("primaryMomE1", "vector<float>", &primaryMomE1);




  //integrated per longitudinal layer
  this->GetTree()->Branch("depositedEnergyTotal", &this->depositedEnergyTotal, "depositedEnergyTotal/F");
  this->GetTree()->Branch("depositedEnergyEscapeWorld", &this->depositedEnergyEscapeWorld, "depositedEnergyEscapeWorld/F");
  this->GetTree()->Branch("depositedEnergyECAL_f", &this->depositedEnergyECAL_f, "depositedEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedEnergyWorld", &this->depositedEnergyWorld, "depositedEnergyWorld/F");
  this->GetTree()->Branch("depositedEnergySolenoid", &this->depositedEnergySolenoid, "depositedEnergySolenoid/F");

  this->GetTree()->Branch("depositedEnergyECAL_absorb_f_particleID", &this->depositedEnergyECAL_absorb_f_particleID, "depositedEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_scinti_f_particleID", &this->depositedEnergyECAL_scinti_f_particleID, "depositedEnergyECAL_scinti_f_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_cheren_f_particleID", &this->depositedEnergyECAL_cheren_f_particleID, "depositedEnergyECAL_cheren_f_particleID[8]/F");


  this->GetTree()->Branch("depositedIonEnergyTotal", &this->depositedIonEnergyTotal, "depositedIonEnergyTotal/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_f", &this->depositedIonEnergyECAL_f, "depositedIonEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedIonEnergyWorld", &this->depositedIonEnergyWorld, "depositedIonEnergyWorld/F");
  this->GetTree()->Branch("depositedIonEnergySolenoid", &this->depositedIonEnergySolenoid, "depositedIonEnergySolenoid/F");

  this->GetTree()->Branch("depositedIonEnergyECAL_absorb_f_particleID", &this->depositedIonEnergyECAL_absorb_f_particleID, "depositedIonEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_scinti_f_particleID", &this->depositedIonEnergyECAL_scinti_f_particleID, "depositedIonEnergyECAL_scinti_f_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_cheren_f_particleID", &this->depositedIonEnergyECAL_cheren_f_particleID, "depositedIonEnergyECAL_cheren_f_particleID[8]/F");

  this->GetTree()->Branch("depositedElecEnergyTotal", &this->depositedElecEnergyTotal, "depositedElecEnergyTotal/F");
  this->GetTree()->Branch("depositedElecEnergyECAL_f", &this->depositedElecEnergyECAL_f, "depositedElecEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedElecEnergyWorld", &this->depositedElecEnergyWorld, "depositedElecEnergyWorld/F");
  this->GetTree()->Branch("depositedElecEnergySolenoid", &this->depositedElecEnergySolenoid, "depositedElecEnergySolenoid/F");



  //Cerenkov photons
  //  this -> GetTree() -> Branch("tot_phot_sci_Timing",        &this->tot_phot_sci_Timing,               "tot_phot_sci_Timing/I");

  this->GetTree()->Branch("tot_phot_cer_ECAL_scinti_f_total", &this->tot_phot_cer_ECAL_scinti_f_total, "tot_phot_cer_ECAL_scinti_f_total/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_scinti_f_particleID", &this->tot_phot_cer_ECAL_scinti_f_particleID, "tot_phot_cer_ECAL_scinti_f_particleID[8]/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_cheren_f_total", &this->tot_phot_cer_ECAL_cheren_f_total, "tot_phot_cer_ECAL_cheren_f_total/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_cheren_f_particleID", &this->tot_phot_cer_ECAL_cheren_f_particleID, "tot_phot_cer_ECAL_cheren_f_particleID[8]/I");



  this->Clear();
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

CreateTree::~CreateTree()
{
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

int CreateTree::Fill()
{
  return this->GetTree()->Fill();
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

bool CreateTree::Write(TFile *outfile)
{
  outfile->cd();
  ftree->Write();


  return true;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void CreateTree::Clear()
{
  Event = 0;


  depositedEnergyEscapeWorld = 0.;

  depositedEnergyTotal = 0.;
  for (int i = 0; i < 3; i++){
    depositedEnergyECAL_f[i] = 0.;
  }
  depositedEnergyWorld = 0.;
  depositedEnergySolenoid = 0.;

  depositedIonEnergyTotal = 0.;
  for (int i = 0; i < 3; i++){
    depositedIonEnergyECAL_f[i] = 0.;
  }
  depositedIonEnergyWorld = 0.;
  depositedIonEnergySolenoid = 0.;

  depositedElecEnergyTotal = 0.;
  for (int i = 0; i < 3; i++){
    depositedElecEnergyECAL_f[i] = 0.;
  }
  depositedElecEnergyWorld = 0.;
  depositedElecEnergySolenoid = 0.;

  tot_phot_cer_ECAL_scinti_f_total = 0.;
  tot_phot_cer_ECAL_cheren_f_total = 0.;

  for (int iparticle = 0; iparticle < 8; iparticle++)
  {
    depositedEnergyECAL_absorb_f_particleID[iparticle] = 0.;
    depositedEnergyECAL_scinti_f_particleID[iparticle] = 0.;
    depositedEnergyECAL_cheren_f_particleID[iparticle] = 0.;

    depositedIonEnergyECAL_absorb_f_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_scinti_f_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_cheren_f_particleID[iparticle] = 0.;

    tot_phot_cer_ECAL_scinti_f_particleID[iparticle] = 0.;
    tot_phot_cer_ECAL_cheren_f_particleID[iparticle] = 0.;
  }


  for (int i = 0; i < 3; ++i)
  {
    inputInitialPosition->at(i) = 0.;
    primaryPosT1->at(i) = 0.;
    primaryPosE1->at(i) = 0.;
  }
  for (int i = 0; i < 4; ++i)
  {
    inputMomentum->at(i) = 0.;
    primaryMomT1->at(i) = 0.;
    primaryMomE1->at(i) = 0.;
  }

  for (int i = 0; i < 3; ++i)
  {
    inputInitialPosition->at(i) = 0.;
    primaryPosT1->at(i) = 0.;
    primaryPosE1->at(i) = 0.;
  }
}

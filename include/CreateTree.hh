#ifndef CreateTree_H
#define CreateTree_H 1

#include <iostream>
#include <vector>
#include "TString.h"
#include <map>

#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

class CreateTree
{
private:
  TTree *ftree;
  TString fname;

public:
  CreateTree(TString name);
  ~CreateTree();

  TTree *GetTree() const { return ftree; };
  TString GetName() const { return fname; };
  void AddEnergyDeposit(int index, float deposit);
  void AddScintillationPhoton(int index);
  void AddCerenkovPhoton(int index);
  int Fill();
  bool Write(TFile *);
  void Clear();

  static CreateTree *Instance() { return fInstance; };
  static CreateTree *fInstance;

  int Event;


  int inputE1Thick;
  int inputE2Thick;
  int inputE1Width;


  std::vector<float> *inputMomentum;        // Px Py Pz E
  std::vector<float> *inputInitialPosition; // x, y, z

  std::vector<float> *primaryMomT1; // Px Py Pz E
  std::vector<float> *primaryPosT1; // x, y, z

  std::vector<float> *primaryMomE1; // Px Py Pz E
  std::vector<float> *primaryPosE1; // x, y, z


  //integrated energy in each longitudinal layer
  float depositedEnergyEscapeWorld;

  float depositedEnergyTotal;
  float depositedEnergyECAL_f[3];
  float depositedEnergyECAL_r[3];
  float depositedEnergySolenoid;
  float depositedEnergyWorld;

  float depositedIonEnergyTotal;
  float depositedIonEnergyECAL_f[3];
  float depositedIonEnergyECAL_r[3];
  float depositedIonEnergySolenoid;
  float depositedIonEnergyWorld;

  float depositedElecEnergyTotal;
  float depositedElecEnergyECAL_f[3];
  float depositedElecEnergyECAL_r[3];
  float depositedElecEnergySolenoid;
  float depositedElecEnergyWorld;


  //store the energy deposition by components

  float depositedEnergyECAL_absorb_f_particleID[8];
  float depositedIonEnergyECAL_absorb_f_particleID[8];


  float depositedEnergyECAL_scinti_f_particleID[8];
  float depositedIonEnergyECAL_scinti_f_particleID[8];

  float depositedEnergyECAL_cheren_f_particleID[8];
  float depositedIonEnergyECAL_cheren_f_particleID[8];


  int tot_phot_cer_ECAL_scinti_f_total;
  int tot_phot_cer_ECAL_scinti_f_particleID[8]; 
  int tot_phot_cer_ECAL_cheren_f_total;
  int tot_phot_cer_ECAL_cheren_f_particleID[8]; 


};

#endif

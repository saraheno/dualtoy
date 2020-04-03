#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

const float sampling_fraction=100.;

void toyplot() {

  const char* inputfilename="electron_50GeV.root";
  const char* outputfilename="electron_50GeV_hist.root";

  TH1F *hTotalE = new TH1F("hTotalE","energy depisited /true",600,0.,1.1);
  TH1F *hWorldE = new TH1F("hWorldE","energy total world /true",600,0.,1.1);
  TH1F *hAbsE = new TH1F("hAbsE","energy absorber /true",600,0.,1.1);
  TH1F *hScintE = new TH1F("hScitE","energy scint /true",600,0.,0.1);
  TH1F *hQuartzE = new TH1F("hQuartzE","energy quartz /true",600,0.,0.1);
  TH1F *hScintQuartzE = new TH1F("hScintQuartzE","energy scint+quartz /true",600,0.,0.1);
  TH1F *hAbsOScintE = new TH1F("hAbsOScintE","energy absorber / energy scint",600,0.,200.);
  TH1F *hAbsOQuartzE = new TH1F("hAbsOQuartzE","energy absorber / energy scint",600,0.,200.);
  TH1F *hAbsOBothE = new TH1F("hAbsOBothE","energy absorber / energy scint",600,0.,200.);


  TFile *f = new TFile(inputfilename);
  TTree *t1 = (TTree*)f->Get("tree");

  vector<float> *inputMomentum = new vector<float>;
  float depositedEnergyECAL_f[3];  // this is total energy
  float depositedIonEnergyECAL_f[3];// this is total ionizing energy
  float depositedEnergyEscapeWorld;
  float depositedEnergyTotal,depositedEnergyWorld;
  int ncerk;

  t1->SetBranchAddress("inputMomentum",&inputMomentum);
  t1->SetBranchAddress("depositedEnergyEscapeWorld",&depositedEnergyEscapeWorld);
  t1->SetBranchAddress("depositedEnergyTotal",&depositedEnergyTotal);
  t1->SetBranchAddress("depositedEnergyECAL_f",&depositedEnergyECAL_f);
  t1->SetBranchAddress("depositedIonEnergyECAL_f",&depositedIonEnergyECAL_f);
  t1->SetBranchAddress("tot_phot_cer_ECAL_cheren_f_total",&ncerk);


  Int_t nentries = (Int_t)t1->GetEntries();
  for(Int_t i=0;i<nentries; i++) {
    t1->GetEntry(i);
    float trueE=9999999.;
    if((*inputMomentum)[3]>0) trueE=(*inputMomentum)[3];
    float Eabs=depositedEnergyECAL_f[0];
    float Escint=depositedIonEnergyECAL_f[1];
    float Equartz=depositedIonEnergyECAL_f[2];

    float EabsT=depositedIonEnergyECAL_f[0];
    float EscintT=depositedIonEnergyECAL_f[1];
    float EquartzT=depositedIonEnergyECAL_f[2];

    
    std::cout<<endl<<"event number "<<i<<std::endl;
    std::cout<<(*inputMomentum)[0]<<","<<(*inputMomentum)[1]<<","<<(*inputMomentum)[2]<<","<<(*inputMomentum)[3]<<std::endl;

    
    std::cout<<"total energy deposited is "<<depositedEnergyTotal<<std::endl;
    std::cout<<"world energy deposited is "<<depositedEnergyEscapeWorld<<std::endl;
    std::cout<<"Abs Scint Quartz "<<Eabs<<" "<<Escint<<" "<<Equartz<<std::endl;
    std::cout<<"AbsT ScintT QuartzT "<<EabsT<<" "<<EscintT<<" "<<EquartzT<<std::endl;
    std::cout<<" number of scint photos "<<ncerk<<std::endl;


    hTotalE->Fill(depositedEnergyTotal/trueE);
    hWorldE->Fill((depositedEnergyTotal+depositedEnergyEscapeWorld)/trueE);
    hAbsE->Fill(Eabs/trueE);
    hScintE->Fill(Escint/trueE);
    hQuartzE->Fill(Equartz/trueE);
    hScintQuartzE->Fill((Escint+Equartz)/trueE);
    float EbothScint=Escint+Equartz;
    if(Escint>0) hAbsOScintE->Fill(Eabs/Escint);
    if(Equartz>0) hAbsOQuartzE->Fill(Eabs/Equartz);
    if(EbothScint>0) hAbsOBothE->Fill(Eabs/EbothScint);
		    
  }

  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hTotalE->Write();
  hWorldE->Write();
  hAbsE->Write();
  hScintE->Write();
  hQuartzE->Write();
  hScintQuartzE->Write();
  hAbsOScintE->Write();
  hAbsOQuartzE->Write();
  hAbsOBothE->Write();
  out->Close();

}



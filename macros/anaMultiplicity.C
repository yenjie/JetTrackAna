#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>

#include 


void anaMultiplicity(char *infDataName= "../output-data.root",char *infMCName = "../output-pthat100.root")
{

   TFile *infData = new TFile(infDataName);
   TTree *tData = (TTree*) infData->Get("t");
   TFile *infMC = new TFile(infMCName);
   TTree *tMC = (TTree*) infMC->Get("t");

   

  


   TLegend *leg = new TLegend(0.2,0.2,0.7,0.6);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(pData,"PbPb 0-10%","");
   leg->AddEntry(pData,"PbPb Data","pl");
   leg->AddEntry(pMC,"Reco PYTHIA+HYDJET","pl");
   leg->AddEntry(pGen,"Gen PYTHIA+HYDJET","l");
   leg->Draw();
}

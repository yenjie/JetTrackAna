#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>

void plotMultVsMC(char *infDataName= "../output-data.root",char *infMCName = "../output-pthat100.root")
{

   TFile *infData = new TFile(infDataName);
   TTree *tData = (TTree*) infData->Get("t");
   TFile *infMC = new TFile(infMCName);
   TTree *tMC = (TTree*) infMC->Get("t");

   TCanvas *c = new TCanvas("c","",600,600);
   
   TProfile *pData = new TProfile("pData","",5,0,0.5);
   TProfile *pMC = new TProfile("pMC","",5,0,0.5);
   TProfile *pGen = new TProfile("pGen","",5,0,0.5);

  
   tData->Draw("Sum$(trkWt*(sqrt(acos(cos(subleadingJetPhi-trkPhi))**2+(subleadingJetEta-trkEta)**2)<0.8||sqrt(acos(cos(leadingJetPhi-trkPhi))**2+(leadingJetEta-trkEta)**2)<0.8)*(2*(acos(cos(leadingJetPhi-trkPhi))<3.1415926/2)-1)):Aj>>pData","hiBin<20&&leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>5*3.14159/6","prof");
   tMC->Draw("Sum$(trkWt*(sqrt(acos(cos(subleadingJetPhi-trkPhi))**2+(subleadingJetEta-trkEta)**2)<0.8||sqrt(acos(cos(leadingJetPhi-trkPhi))**2+(leadingJetEta-trkEta)**2)<0.8)*(2*(acos(cos(leadingJetPhi-trkPhi))<3.1415926/2)-1)):Aj>>pMC","hiBin<20&&leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>5*3.14159/6","prof");
   tMC->Draw("Sum$((sqrt(acos(cos(gensubleadingJetPhi-pPhi))**2+(gensubleadingJetEta-pEta)**2)<0.8||sqrt(acos(cos(genleadingJetPhi-pPhi))**2+(genleadingJetEta-pEta)**2)<0.8)*(2*(acos(cos(genleadingJetPhi-pPhi))<3.1415926/2)-1)):Aj>>pGen","hiBin<20&&genleadingJetPt>120&&gensubleadingJetPt>50&&acos(cos(genleadingJetPhi-gensubleadingJetPhi))>5*3.14159/6","prof");

   pData->SetXTitle("A_{J}");
   pData->SetYTitle("Multiplicity Difference");
   pData->SetAxisRange(-50,0,"Y");   
   pData->Draw();
   pMC->SetLineColor(2);
   pMC->SetMarkerColor(2);
   pGen->SetLineColor(2);
   pMC->Draw("same");
   pGen->Draw("hist same");

   TLegend *leg = new TLegend(0.2,0.2,0.7,0.6);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(pData,"PbPb 0-10%","");
   leg->AddEntry(pData,"PbPb Data","pl");
   leg->AddEntry(pMC,"Reco PYTHIA+HYDJET","pl");
   leg->AddEntry(pGen,"Gen PYTHIA+HYDJET","l");
   leg->Draw();
}

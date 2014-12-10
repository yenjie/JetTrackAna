#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include "../jettrkcorr/hiForest.h"

void plotSwap()
{
   TFile *inf = new TFile("../mptTree-merged.root");
   TTree *t = (TTree*)inf->Get("t");
   TFile *outf = new TFile("swapPlot.root","recreate");
   const int nCentBins=5;
   int centBins[nCentBins+1] = { 0,20,60,100,140,200};
   
   

   // initialize histograms   
   
   for (int i=0;i<nCentBins;i++) {
      TProfile *pLeading = new TProfile(Form("pLeading%d",i),"",50,0,1);
      TCanvas *c = new TCanvas(Form("c%i",i),Form("Leading Jet Matching Probability %d-%d",centBins[i],centBins[i+1]),600,600);
      pLeading->SetXTitle("Reco A_{J}");
      pLeading->SetYTitle("Leading Jet Matching Prob");
      TCut centCut = Form("%d<=hiBin&&hiBin<%d",centBins[i],centBins[i+1]);
      t->Draw(Form("leadingJetRefPt==genleadingJetPt:Aj>>pLeading%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3.","prof");
      // matching probability between leading reco jet and leading gen jet
      TF1 *fPLeading = new TF1(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]),"TMath::Erf(x*[0]+[1])");
      pLeading->Fit(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]));
      pLeading->Fit(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]));
      pLeading->Fit(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]));
      pLeading->Fit(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]));
      fPLeading->Write();
   } 
   // matching probability between subleading reco jet and subleading jet
   outf->Write();
}

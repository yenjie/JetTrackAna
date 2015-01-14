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
//   TFile *inf = new TFile("../mptTree-FFCor-merged.root");
//   TFile *inf = new TFile("../mptTree-80.root");
//   TFile *inf = new TFile("../mptTree-merged-jec1212.root");
//     TFile *inf = new TFile("../mptTree-80-R5.root");
     TFile *inf = new TFile("../mptTree-11080-R2.root");
   TTree *t = (TTree*)inf->Get("t");
   TFile *outf = new TFile("swapPlot.root","recreate");
   const int nCentBins=5;
   int centBins[nCentBins+1] = { 0,20,60,100,140,200};
   
   

   // initialize histograms   
   t->SetAlias("LD","leadingJetPt*0.00760+leadingJetNTrk4*(-0.00304)+subleadingJetPt*(-0.00883)+subleadingJetNTrk4*0.01045-0.0272");
   for (int i=0;i<nCentBins;i++) {
      TProfile *pLeading = new TProfile(Form("pLeading%d",i),"",50,0,1);
      TCanvas *c = new TCanvas(Form("c%i",i),Form("Leading Jet Matching Probability %d-%d",centBins[i],centBins[i+1]),600,600);
      pLeading->SetXTitle("Reco A_{J}");
      pLeading->SetYTitle("Leading Jet Matching Prob");
      TCut centCut = Form("%d<=hiBin&&hiBin<%d",centBins[i],centBins[i+1]);
      t->Draw(Form("(leadingJetRefPt==genleadingJetPt)*(LD>-0.2)+(LD<-0.2)*(subleadingJetRefPt==genleadingJetPt):Aj>>pLeading%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3.","prof");
      // matching probability between leading reco jet and leading gen jet
      TF1 *fPLeading = new TF1(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]),"TMath::Erf(x*[0])*0.5+0.5");
      pLeading->Fit(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]));
      pLeading->Fit(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]));
      pLeading->Fit(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]));
      pLeading->Fit(Form("fPLeading_%d_%d",centBins[i],centBins[i+1]));
      fPLeading->Write();
   } 
   // matching probability between subleading reco jet and subleading jet
   outf->Write();
}

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


void checkSwap()
{
//   TFile *inf = new TFile("../mptTree-FFCor-merged.root");
//   TFile *inf = new TFile("../mptTree-80.root");
//   TFile *inf = new TFile("../mptTree-merged-jec1212.root");
//     TFile *inf = new TFile("../mptTree-80-R3.root");
     TFile *inf = new TFile("../mptTreeFFCor-80-R3.root");
   TTree *t = (TTree*)inf->Get("t");
   TFile *outf = new TFile("swapPlot.root","recreate");
   const int nCentBins=5;
   int centBins[nCentBins+1] = { 0,20,60,100,140,200};
   
   

   // initialize histograms   
//   t->SetAlias("LD","leadingJetPt*0.00760+leadingJetNTrk4*(-0.00304)+subleadingJetPt*(-0.00883)+subleadingJetNTrk4*0.01045-0.0272");
   for (int i=0;i<nCentBins;i++) {
      TH1D *pLeading = new TH1D(Form("pLeading%d",i),"",20,0,1);
      TH1D *p123 = new TH1D(Form("p123_%d",i),"",20,0,1);
      TH1D *p132 = new TH1D(Form("p132_%d",i),"",20,0,1);
      TH1D *p213 = new TH1D(Form("p213_%d",i),"",20,0,1);
      TH1D *p231 = new TH1D(Form("p231_%d",i),"",20,0,1);
      TH1D *p312 = new TH1D(Form("p312_%d",i),"",20,0,1);
      TH1D *p321 = new TH1D(Form("p321_%d",i),"",20,0,1);
      TCanvas *c = new TCanvas(Form("c%i",i),Form("Leading Jet Matching Probability %d-%d",centBins[i],centBins[i+1]),600,600);
      pLeading->SetXTitle("Reco A_{J}");
      pLeading->SetYTitle("Entries");
      TCut centCut = Form("%d<=hiBin&&hiBin<%d",centBins[i],centBins[i+1]);
      TCut cut123 = "leadingJetRefPt==genleadingJetPt&&subleadingJetRefPt==gensubleadingJetPt&&((thirdleadingJetPt<50)||(thirdleadingJetRefPt==genthirdleadingJetPt))";
      TCut cut132 = "leadingJetRefPt==genleadingJetPt&&subleadingJetRefPt==genthirdleadingJetPt&&thirdleadingJetRefPt==gensubleadingJetPt";
      TCut cut213 = "leadingJetRefPt==gensubleadingJetPt&&subleadingJetRefPt==genleadingJetPt&&((thirdleadingJetPt<50)||(thirdleadingJetRefPt==genthirdleadingJetPt))";
      TCut cut231 = "leadingJetRefPt==gensubleadingJetPt&&subleadingJetRefPt==genthirdleadingJetPt&&thirdleadingJetRefPt==genleadingJetPt";
      TCut cut312 = "leadingJetRefPt==genthirdleadingJetPt&&subleadingJetRefPt==genleadingJetPt&&thirdleadingJetRefPt==gensubleadingJetPt";
      TCut cut321 = "leadingJetRefPt==genthirdleadingJetPt&&subleadingJetRefPt==gensubleadingJetPt&&thirdleadingJetRefPt==genleadingJetPt";
      t->Draw(Form("Aj>>pLeading%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3.");
      t->Draw(Form("Aj>>p123_%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3."&&cut123);
      t->Draw(Form("Aj>>p132_%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3."&&cut132);
      t->Draw(Form("Aj>>p213_%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3."&&cut213);
      t->Draw(Form("Aj>>p231_%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3."&&(cut231||cut312||cut321));
      t->Draw(Form("Aj>>p312_%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3."&&cut312);
      t->Draw(Form("Aj>>p321_%d",i),centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.14159/3."&&cut321);
      p123->SetLineColor(4);
      p213->SetLineColor(2);
      p132->SetLineColor(kGreen+2);
      p231->SetLineColor(6);
      pLeading->Draw();
      p123->Draw("same");
      p213->Draw("same");
      p132->Draw("same");
      p231->Draw("same");

      TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(pLeading,"All","l");
      leg->AddEntry(p123,"123","l");
      leg->AddEntry(p132,"132","l");
      leg->AddEntry(p213,"213","l");
      leg->AddEntry(p231,"231||312||321","l");
      leg->Draw();
   } 
   // matching probability between subleading reco jet and subleading jet
   outf->Write();
}

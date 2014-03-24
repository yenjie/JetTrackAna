#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>

#include "commonDefinition.h"
#include "commonTool.h"
#include "myTree.h"

#define PI 3.1415926535897932384626

void loop(treeData &data,TTree *tData)
{

   TH1D *hPtTmp = new TH1D("hPtTmp","",nPtBins,ptBins);
         
   TH1D *hMpt[10];
   
   for (int i=0; i<10; i++) {
      hMpt[i] = new TH1D(Form("hMpt%d",i),"",1000,-10000,10000);
      
   }
   
   
   int njet=0;
   for (int i=0;i<tData->GetEntries();i++) {
      if (i%1000==0) cout <<i<<" "<<tData->GetEntries()<<endl;
      tData->GetEntry(i);
      //cout <<data.hiBin<<endl;
            
      double AJ = (data.leadingJetPt-data.subleadingJetPt)/(data.leadingJetPt+data.subleadingJetPt);
      if (data.hiBin<20&&AJ>0.3&&data.leadingJetPt>120 && data.subleadingJetPt>50 && fabs(deltaPhi(data.leadingJetPhi,data.subleadingJetPhi))>5.*PI/6.) {
         for (int j=0;j<=nPtBins;j++){
	    data.missingPt[j]=0;
	 }
         for (int j=0;j<data.nTrk;j++) {
	    
	    // calculate the distance with respect to leading and subleading jet
	    double deta1 = data.leadingJetEta-data.trkEta[j];
 	    double dphi1 = fabs(deltaPhi(data.leadingJetPhi,data.trkPhi[j]));
	    double deta2 = data.subleadingJetEta-data.trkEta[j];
 	    double dphi2 = fabs(deltaPhi(data.subleadingJetPhi,data.trkPhi[j]));
//	    double dR   = sqrt(deta*deta+dphi*dphi);
	    
	    data.hMultVsEtaPhi->Fill(deta1,dphi1,data.trkWt[j]);
	    data.hMultVsEtaPhi->Fill(deta2,dphi2,-data.trkWt[j]);
	    data.hMultVsEtaPhi->Fill(deta1,-dphi1,data.trkWt[j]);
	    data.hMultVsEtaPhi->Fill(deta2,-dphi2,-data.trkWt[j]);
	    data.hMultVsEtaPhi->Fill(-deta1,dphi1,data.trkWt[j]);
	    data.hMultVsEtaPhi->Fill(-deta2,dphi2,-data.trkWt[j]);
	    data.hMultVsEtaPhi->Fill(-deta1,-dphi1,data.trkWt[j]);
	    data.hMultVsEtaPhi->Fill(-deta2,-dphi2,-data.trkWt[j]);
	    
	    Float_t avgPhi = getAvePhi(data.leadingJetPhi,data.subleadingJetPhi);
	    
	    int ptBin = hPtTmp->GetBin(data.trkPt[j]);
	    double mpt = data.trkPt[j] * data.trkWt[j] * cos(avgPhi-data.trkPhi[j]);
	    data.missingPt[ptBin]+= mpt;
	    data.missingPt[0]+=mpt;
	 }  
         njet++;
	 for (int j=0;j<10;j++){
	   hMpt[j]->Fill(data.missingPt[j]);
	 }
      }
   }
   data.hMultVsEtaPhi->Scale(1./njet);
}

void anaMultiplicity(char *infDataName= "../output-data.root",char *infMCName = "../output-pthat100.root")
{

   TFile *infData = new TFile(infDataName);
   TTree *tData = (TTree*) infData->Get("t");
   TFile *infMC = new TFile(infMCName);
   TTree *tMC = (TTree*) infMC->Get("t");

   treeData data;
   treeData mc;
   
   loadBranch(data,tData);
   loadBranch(mc,tMC);
   
   loop(data,tData);
   loop(mc,tMC);
   
   TCanvas *c = new TCanvas("c","Data",600,600);
   
   data.hMultVsEtaPhi->Draw("col");

   TCanvas *cMC = new TCanvas("cMC","MC",600,600);
   mc.hMultVsEtaPhi->Draw("col");
   
   TH2D *hDiff = (TH2D*)data.hMultVsEtaPhi->Clone("hDiff");
   hDiff->Add(mc.hMultVsEtaPhi,-1);
   TCanvas *cDiff = new TCanvas("cDiff","Data-MC",600,600);
   hDiff->Draw("col");
     
/*
   TLegend *leg = new TLegend(0.2,0.2,0.7,0.6);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(pData,"PbPb 0-10%","");
   leg->AddEntry(pData,"PbPb Data","pl");
   leg->AddEntry(pMC,"Reco PYTHIA+HYDJET","pl");
   leg->AddEntry(pGen,"Gen PYTHIA+HYDJET","l");
   */
//   leg->Draw();
}

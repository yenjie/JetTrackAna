#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TNtuple.h>

#include "commonDefinition.h"
#include "commonTool.h"
#include "myTree.h"

#define PI 3.1415926535897932384626

void loop(treeData &data,TTree *tData,char *tag="DATA")
{

   TH1D *hPtTmp = new TH1D(Form("%s_hPtTmp",tag),"",nPtBins,ptBins);
   data.nt = new TNtuple(Form("%s_nt",tag),"","pt:mpt:deta1:dphi1:deta2:dphi2:dR1:dR2:dphiAvg:hiBin:leadingJetPt:subleadingJetPt:leadingJetEta,subleadingJetEta");


   Float_t Aj=0;
   Float_t genAj=0;
   Float_t genDphi=0;
   Float_t dphi=0;
   Float_t avgPhi=0;
   Float_t genMultDiffPt[5];
   Float_t multDiffPt[5];
   
   TTree *t = new TTree(Form("%s_tree",tag),"");
   t->Branch("hiBin",&data.hiBin,"hiBin/I");
   t->Branch("leadingJetPt",&data.leadingJetPt,"leadingJetPt/F");
   t->Branch("leadingJetPhi",&data.leadingJetPhi,"leadingJetPhi/F");
   t->Branch("leadingJetEta",&data.leadingJetEta,"leadingJetEta/F");
   t->Branch("subleadingJetPt",&data.subleadingJetPt,"subleadingJetPt/F");
   t->Branch("subleadingJetPhi",&data.subleadingJetPhi,"subleadingJetPhi/F");
   t->Branch("subleadingJetEta",&data.subleadingJetEta,"subleadingJetEta/F");
   t->Branch("genleadingJetPt",&data.genleadingJetPt,"genleadingJetPt/F");
   t->Branch("genleadingJetPhi",&data.genleadingJetPhi,"genleadingJetPhi/F");
   t->Branch("genleadingJetEta",&data.genleadingJetEta,"genleadingJetEta/F");
   t->Branch("gensubleadingJetPt",&data.gensubleadingJetPt,"gensubleadingJetPt/F");
   t->Branch("gensubleadingJetPhi",&data.gensubleadingJetPhi,"gensubleadingJetPhi/F");
   t->Branch("gensubleadingJetEta",&data.gensubleadingJetEta,"gensubleadingJetEta/F");
   t->Branch("multDiff",&data.multDiff,"multDiff/F");
   t->Branch("genMultDiff",&data.genMultDiff,"genMultDiff/F");
   t->Branch("multDiffPt",multDiffPt,"multDiffPt[5]/F");
   t->Branch("genMultDiffPt",genMultDiffPt,"genMultDiffPt[5]/F");
   t->Branch("Aj",&Aj,"Aj/F");
   t->Branch("genAj",&genAj,"genAj/F");
   t->Branch("dphi",&dphi,"dphi/F");
   t->Branch("genDphi",&genDphi,"genDphi/F");

   t->Branch("coneMult",&data.coneMult,"coneMult/F");
   t->Branch("genconeMult",&data.genconeMult,"genconeMult/F");
   t->Branch("coneMultPt",coneMultPt,"coneMultPt[5]/F");
   t->Branch("genconeMultPt",genconeMultPt,"genconeMultPt[5]/F");


   int njet=0;
   for (int i=0;i<tData->GetEntries()/1.;i++) {
      if (i%1000==0) cout <<i<<" "<<tData->GetEntries()<<endl;
      tData->GetEntry(i);
      //cout <<data.hiBin<<endl;
            
      Aj    = (data.leadingJetPt-data.subleadingJetPt)/(data.leadingJetPt+data.subleadingJetPt);
      genAj = (data.genleadingJetPt-data.gensubleadingJetPt)/(data.genleadingJetPt+data.gensubleadingJetPt);
      dphi     = fabs(deltaPhi(data.leadingJetPhi,data.subleadingJetPhi)); 
      genDphi  = fabs(deltaPhi(data.genleadingJetPhi,data.gensubleadingJetPhi)); 

      
      // GEN study
            
      data.genMultDiff=0;
      for (int j=0;j<data.nP;j++) {
         avgPhi = getAvePhi(data.genleadingJetPhi,data.gensubleadingJetPhi);
         if (cos(avgPhi-data.pPhi[j])>0) {
            data.genMultDiff++;
         } else {  
            data.genMultDiff--;
         }
      }

      // RECO study
      data.multDiff=0;
      for (int j=0;j<data.nTrk;j++) {
         avgPhi = getAvePhi(data.leadingJetPhi,data.subleadingJetPhi);
         if (cos(avgPhi-data.trkPhi[j])>0) {
            data.multDiff+=data.trkWt[j];
         } else {  
            data.multDiff-=data.trkWt[j];
         }
      }  

      t->Fill();
   }
}

void anaMultiplicity(char *infDataName= "../output-data.root",char *infMCName = "merged.root", char *infPPName = "../output-ppData.root")
{

   TFile *infData = new TFile(infDataName);
   TTree *tData = (TTree*) infData->Get("t");
   TFile *infMC = new TFile(infMCName);
   TTree *tMC = (TTree*) infMC->Get("t");
   TFile *infPP = new TFile(infPPName);
   TTree *tPP = (TTree*) infPP->Get("t");

   TFile *outfile = new TFile("ana.root","recreate");

  
   treeData data;
   treeData mc;
   treeData pp;

   loadBranch(data,tData);
   loadBranch(mc,tMC);
   loadBranch(pp,tPP);

   loop(data,tData,"Data");
   loop(mc,tMC,"MC");
   loop(pp,tPP,"PP");

/*   
   TCanvas *c = new TCanvas("c","Data",600,600);
   
   data.hMultVsEtaPhi->Draw("col");

   TCanvas *cMC = new TCanvas("cMC","MC",600,600);
   mc.hMultVsEtaPhi->Draw("col");
   
   TH2D *hDiff = (TH2D*)data.hMultVsEtaPhi->Clone("hDiff");
   hDiff->Add(mc.hMultVsEtaPhi,-1);
   TCanvas *cDiff = new TCanvas("cDiff","Data-MC",600,600);
   hDiff->Draw("col");
     */
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

    outfile->Write();
}

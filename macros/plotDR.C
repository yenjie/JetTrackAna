#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TLegend.h>

void plotDR(char *infname = "output-pthat100.root")
{

   TFile *inf = new TFile(infname);
   TTree *t = (TTree*) inf->Get("t");
   
   TCanvas *c = new TCanvas("c","",600,600);
   
   TH1D *p = new TH1D("p","",20,0,5);
   TH1D *pGenReco = new TH1D("pGenReco","",20,0,5);
   TH1D *pRecoGen = new TH1D("pRecoGen","",20,0,5);
   TH1D *pGen = new TH1D("pGen","",20,0,5);

   t->Draw("sqrt(acos(cos(leadingJetPhi-trkPhi))**2+(leadingJetEta-trkEta)**2)>>p","trkWt*(trkPt<1&&leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>5*3.14159/6)","prof");
   t->Draw("sqrt(acos(cos(leadingJetPhi-pPhi))**2+(leadingJetEta-pEta)**2)>>pGenReco","pPt<1&&leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>5*3.14159/6","prof");
   t->Draw("sqrt(acos(cos(genleadingJetPhi-trkPhi))**2+(genleadingJetEta-trkEta)**2)>>pRecoGen","trkWt*(trkPt<1&&genleadingJetPt>120&&gensubleadingJetPt>50&&acos(cos(genleadingJetPhi-gensubleadingJetPhi)))>5*3.14159/6","prof");
   t->Draw("sqrt(acos(cos(genleadingJetPhi-pPhi))**2+(genleadingJetEta-pEta)**2)>>pGen","pPt<1&&genleadingJetPt>120&&gensubleadingJetPt>50&&acos(cos(genleadingJetPhi-gensubleadingJetPhi))>5*3.14159/6","prof");
   p->Sumw2();
   p->Scale(1/p->Integral(0,1000));
   pGenReco->Scale(1/pGenReco->Integral(0,1000));
   pRecoGen->Scale(1/pRecoGen->Integral(0,1000));
   pGen->Scale(1/pGen->Integral(0,1000));

   p->SetXTitle("#Delta R w/r leading jet");
   p->SetYTitle("Entries");
   p->SetAxisRange(0,1,"Y");   
   p->Draw("p");
   pGen->Draw("same");
   pGen->SetLineWidth(2);
   pRecoGen->SetLineColor(2);
   pRecoGen->Draw("hist same");
   pGenReco->SetLineColor(4);
   pGenReco->Draw("hist same");
   TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(p,"Reco T-Reco J","p");
   leg->AddEntry(pGenReco,"Gen T-Reco J","l");
   leg->AddEntry(pRecoGen,"Reco T-Gen J","l");
   leg->AddEntry(pGen,"Gen T-Gen J","l");

   leg->Draw();
}

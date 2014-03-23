#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TProfile.h>

void plotMpt(char *infname = "output-pthat100.root")
{

   TFile *inf = new TFile(infname);
   TTree *t = (TTree*) inf->Get("t");
   
   TCanvas *c = new TCanvas("c","",600,600);
   
   TProfile *p = new TProfile("p","",10,0,200);
   
   t->Draw("Sum$(trkPt*trkWt*:hiBin>>p","Aj>0.24&&leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>5*3.14159/6","prof");

   p->SetXTitle("Centrality Bin");
   p->SetYTitle("TrkPt Weighted Multiplicity Difference");
   p->SetAxisRange(-30,30,"Y");   
   p->Draw();

}

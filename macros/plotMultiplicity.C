#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TProfile.h>

void plotMultiplicity(char *infname = "output-pthat100.root")
{

   TFile *inf = new TFile(infname);
   TTree *t = (TTree*) inf->Get("t");
   
   TCanvas *c = new TCanvas("c","",600,600);
   
   TProfile *p = new TProfile("p","",10,0,0.5);
   TProfile *pGen = new TProfile("pGen","",10,0,0.5);

   t->Draw("Sum$(trkWt*(2*(acos(cos(leadingJetPhi-trkPhi))<3.1415926/2)-1)):Aj>>p","leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>5*3.14159/6","prof");
   t->Draw("Sum$((2*(acos(cos(genleadingJetPhi-pPhi))<3.1415926/2)-1)):Aj>>pGen","genleadingJetPt>120&&gensubleadingJetPt>50&&acos(cos(genleadingJetPhi-gensubleadingJetPhi))>5*3.14159/6","prof");

   p->SetXTitle("Centrality Bin");
   p->SetYTitle("Multiplicity Difference");
   p->SetAxisRange(-50,0,"Y");   
   p->Draw();
   pGen->Draw("hist same");

}

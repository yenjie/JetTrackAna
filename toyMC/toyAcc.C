#include <TH1D.h>
#include <TNtuple.h>
#include <TH2D.h>
#include <TFile.h>
#include <TRandom2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>

#define PI 3.141592653589793238462643

void toyAcc()
{

   TFile *outf = new TFile("toyMC.root","recreate");
   TNtuple *nt = new TNtuple("nt","","eta:phi");
   TRandom2 rnd;
   
   for (int i=0;i<10000000;i++){         
       nt->Fill(-10+20*rnd.Rndm(),-PI+2*PI*rnd.Rndm());
   }
   
   TProfile *p00 = new TProfile("p00","",100,0,10);
   TProfile *p05 = new TProfile("p05","",100,0,10);
   TProfile *p10 = new TProfile("p10","",100,0,10);
   TProfile *p15 = new TProfile("p15","",100,0,10);
   TProfile *p20 = new TProfile("p20","",100,0,10);

   TCanvas *c = new TCanvas("c","",600,600);

   nt->Draw("abs(eta)<2.4:sqrt((eta-0.0)**2+phi**2)>>p00","","prof");
   nt->Draw("abs(eta)<2.4:sqrt((eta-0.5)**2+phi**2)>>p05","","prof same");
   nt->Draw("abs(eta)<2.4:sqrt((eta-1.0)**2+phi**2)>>p10","","prof same");
   nt->Draw("abs(eta)<2.4:sqrt((eta-1.5)**2+phi**2)>>p15","","prof same");
   nt->Draw("abs(eta)<2.4:sqrt((eta-2.0)**2+phi**2)>>p20","","prof same");

   p00->SetXTitle("#Delta R");
   p00->SetYTitle("Acceptance");

   p00->SetMarkerColor(1);
   p05->SetMarkerColor(2);
   p10->SetMarkerColor(4);
   p15->SetMarkerColor(kGreen+2);
   p20->SetMarkerColor(kGray);
  
   TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(p00,"Acc for jet #eta=0.0","p");
   leg->AddEntry(p05,"Acc for jet #eta=0.5","p");
   leg->AddEntry(p10,"Acc for jet #eta=1.0","p");
   leg->AddEntry(p15,"Acc for jet #eta=1.5","p");
   leg->AddEntry(p20,"Acc for jet #eta=2.0","p");

   leg->Draw();

   outf->Write();
//   outf->Close();   
      

}

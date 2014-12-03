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

void plotMpt()
{
   // binning
   double AjBins[5] = { 0, 0.11, 0.22, 0.33, 0.499999};
   double trackPtBinL[10] = { 0.5,    0,   0.5, 1.0, 2.0, 4.0, 8.0};
   double trackPtBinH[10] = { 300, 0.5, 1.0, 2.0, 4.0, 8.0, 300};

   // open file
   TFile *inf = new TFile("mptTree-merged.root");
   TTree *t = (TTree*)inf->Get("t");
   TFile *outf = new TFile("mptPlot.root","recreate");

   // initialize histograms   
   TH1D *hProj = new TH1D("hProj","",10000,-500,500);
   TH1D *hMpt[10];
   for (int i=0;i<10;i++) hMpt[i] = new TH1D(Form("hMpt%d",i),"",4, AjBins);
   TH1D *hGenMptGenAj[10];
   for (int i=0;i<10;i++) hGenMptGenAj[i] = new TH1D(Form("hGenMptGenAj%d",i),"",4, AjBins);
   TH1D *hGenMptRecoAj[10];
   for (int i=0;i<10;i++) hGenMptRecoAj[i] = new TH1D(Form("hGenMptRecoAj%d",i),"",4, AjBins);
   TH1D *hGenMptWithRecoJetDirRecoAj[10];
   for (int i=0;i<10;i++) hGenMptWithRecoJetDirRecoAj[i] = new TH1D(Form("hGenMptWithRecoJetDirRecoAj%d",i),"",4, AjBins);

   // dijet selection cuts
   TCut centCut = "0<=hiBin&&hiBin<20";
//   TCut centCut = "hiBin>100";
//   TCut centCut = "hiBin<60";
   TCut recoCut = centCut&&"leadingJetPt>120&&subleadingJetPt>50&&acos(cos(leadingJetPhi-subleadingJetPhi))>2*3.1415927/3.";
   TCut genCut  = centCut&&"genleadingJetPt>120&&gensubleadingJetPt>50&&acos(cos(genleadingJetPhi-gensubleadingJetPhi))>2*3.1415927/3.";


   for (int iAj=0;iAj<4;iAj++) {
      TCut AjCut = Form("%.5f<=Aj&&Aj<%.5f",AjBins[iAj],AjBins[iAj+1]); 
      TCut genAjCut = Form("%.5f<=genAj&&genAj<%.5f",AjBins[iAj],AjBins[iAj+1]); 
      cout <<iAj<<" "<<AjCut.GetTitle()<<endl;
      for (int iPt=0; iPt<=6;iPt++) {

         t->Draw(Form("cormpt[%d]>>hProj",iPt),recoCut&&AjCut);
	 hMpt[iPt]->SetBinContent(iAj+1,hProj->GetMean());
	 hMpt[iPt]->SetBinError(iAj+1,hProj->GetMeanError());

         t->Draw(Form("genMpt[%d]>>hProj",iPt),genCut&&genAjCut);
	 hGenMptGenAj[iPt]->SetBinContent(iAj+1,hProj->GetMean());
	 hGenMptGenAj[iPt]->SetBinError(iAj+1,hProj->GetMeanError());

         t->Draw(Form("genMpt[%d]>>hProj",iPt),recoCut&&AjCut);
	 hGenMptRecoAj[iPt]->SetBinContent(iAj+1,hProj->GetMean());
	 hGenMptRecoAj[iPt]->SetBinError(iAj+1,hProj->GetMeanError());

         t->Draw(Form("genPMpt[%d]>>hProj",iPt),recoCut&&AjCut);
	 hGenMptWithRecoJetDirRecoAj[iPt]->SetBinContent(iAj+1,hProj->GetMean());
	 hGenMptWithRecoJetDirRecoAj[iPt]->SetBinError(iAj+1,hProj->GetMeanError());

      }
   }
  
  
  
   for (int i=0;i<=6;i++)
   {
      TCanvas *c = new TCanvas(Form("c%d",i),Form("%.5f<pT<%.5f",trackPtBinL[i],trackPtBinH[i]),600,600);
      TH2D *h = new TH2D(Form("h%d",i),"",100,0,0.5,100,-10,10);
      TLine *l = new TLine(0,0,0.5,0);
      l->SetLineStyle(2);
      h->SetXTitle("A_{J}");
      h->SetYTitle("Mpt");
      h->Draw();
      l->Draw();
      hMpt[i]->Draw("same");
      hMpt[i]->SetMarkerStyle(25);
      hGenMptWithRecoJetDirRecoAj[i]->SetMarkerStyle(24);
      hGenMptGenAj[i]->SetLineColor(2);
      hGenMptRecoAj[i]->SetLineColor(4);
      hGenMptWithRecoJetDirRecoAj[i]->SetLineColor(4);
      hGenMptGenAj[i]->SetMarkerColor(2);
      hGenMptWithRecoJetDirRecoAj[i]->SetMarkerColor(4);
      hGenMptRecoAj[i]->SetMarkerColor(4);
      hGenMptGenAj[i]->Draw("same");
      hGenMptRecoAj[i]->Draw("same");
      hGenMptWithRecoJetDirRecoAj[i]->Draw("same");
   
      TLegend *leg =new TLegend(0.3,0.7,0.9,0.9);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(hGenMptGenAj[i], "PYTHIA+HYDJET 0-10%","t");
      leg->AddEntry(hGenMptGenAj[i],Form("%.1f<pT<%.1f",trackPtBinL[i],trackPtBinH[i]),"t");
      leg->AddEntry(hGenMptGenAj[i], "Gen MpT vs Gen Aj, Gen Jet Sel","pl");
      leg->AddEntry(hGenMptRecoAj[i], "Gen MpT vs Reco Aj, Reco Jet Sel","pl");
      leg->AddEntry(hGenMptWithRecoJetDirRecoAj[i], "Gen MpT (recojet dir) vs Reco Aj, Reco Jet Sel","pl");
      leg->AddEntry(hMpt[i], "Reco MpT vs Reco Aj, Reco Jet Sel","pl");
      leg->Draw();
      c->Write();  
   }
   
   
}

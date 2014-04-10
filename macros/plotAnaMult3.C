#include "commonTool.h"

#include <TCanvas.h>
#include <TColor.h>


void plotAnaMult3(char *infname="ana.root")
{    
   TFile *inf = new TFile(infname);
   
   TTree *tData = (TTree*) inf->Get("Data_tree");
   TTree *tMC   = (TTree*) inf->Get("MC_tree");
   TTree *tPP   = (TTree*) inf->Get("PP_tree");

   TCanvas *c1 = new TCanvas("c","",1200,700);
//   c->Divide(4,1);
   TCut recoCut = "leadingJetPt>120&&subleadingJetPt>50&&dphi>2*3.14159265358979/3.&&abs(leadingJetEta)<1.6&&abs(subleadingJetEta)<1.6";
   TCut genCut = "genleadingJetPt>120&&gensubleadingJetPt>50&&genDphi>2*3.14159265358979/3.&&abs(genleadingJetEta)<1.6&&abs(gensubleadingJetEta)<1.6";

   
   const int nPtBin=5;
   Float_t PtBins[nPtBin+1] = {2.1,2.6,2.8,2.9,3.0,3.15};
   
   Int_t centBin[5] = {200,100,60,20,0};
   makeMultiPanelCanvas(c1,4,2,0.0,0.0,0.2,0.2,0.02);
   
   for (int i=0;i<4;i++) {
      c1->cd(i+1);
      TH1D * empty=new TH1D("empty","",100,2.1,3.15);
      empty->Fill(0.5,1000); 
      empty->SetMaximum(39.99); 
      empty->SetMinimum(0.001); 
      empty->SetNdivisions(105,"X");
      empty->GetXaxis()->SetTitleSize(28);
      empty->GetXaxis()->SetTitleFont(43); 
      empty->GetXaxis()->SetTitleOffset(1.8);
      empty->GetXaxis()->SetLabelSize(22);
      empty->GetXaxis()->SetLabelFont(43);
      empty->GetYaxis()->SetTitleSize(28);
      empty->GetYaxis()->SetTitleFont(43); 
      empty->GetYaxis()->SetTitleOffset(1.8);
      empty->GetYaxis()->SetLabelSize(22);
      empty->GetYaxis()->SetLabelFont(43);
      empty->GetXaxis()->CenterTitle();
      empty->GetYaxis()->CenterTitle();

      empty->SetXTitle("#Delta#phi_{1,2}");
      empty->SetYTitle("Multiplicity Difference");
   
      TProfile *pData = new TProfile(Form("pData%d",i),"",nPtBin,PtBins);
      TProfile *pMC = new TProfile(Form("pMC%d",i),"",nPtBin,PtBins);
      TProfile *pPP = new TProfile(Form("pPP%d",i),"",nPtBin,PtBins);
      TProfile *pGen = new TProfile(Form("pGen%d",i),"",nPtBin,PtBins);
      TCut centCut = Form("hiBin>=%d&&hiBin<%d",centBin[i+1],centBin[i]);
      tData->Draw(Form("-multDiff:dphi>>pData%d",i),recoCut&&centCut);
      tPP->Draw(Form("-multDiff:dphi>>pPP%d",i),recoCut);
      tMC->Draw(Form("-multDiff:dphi>>pMC%d",i),recoCut&&centCut);  
      tMC->Draw(Form("-genMultDiff:genDphi>>pGen%d",i),genCut&&centCut);  

      pMC->SetLineColor(2);
      pMC->SetMarkerColor(2);
      pPP->SetLineColor(4);
      pPP->SetMarkerColor(4);
      
//      pData->SetAxisRange(0,50,"Y");
//      pData->SetAxisRange(0,0.49,"X");
      empty->Draw();

      double diff=0;
      if (i==0) diff=0.1;
      drawText(Form("%d-%d %%",(int)(0.5*centBin[i+1]),(int)(0.5*centBin[i])),0.22+diff,0.65);
      if (i==0) drawText("PbPb #sqrt{s_{NN}}=2.76 TeV 150/#mub",0.22+diff,0.85);
      if (i==0) drawText("pp      #sqrt{s_{NN}}=2.76 TeV 5.3/pb",0.22+diff,0.75);
      if (i==3) drawText("CMS Preliminary",0.3+diff,0.85);
      
      Float_t sys[4]={1,1,2.5,3};
      
      
      for (int j=1;j<=pData->GetNbinsX();j++)
      {
         TBox *b = new TBox(pData->GetBinLowEdge(j),pData->GetBinContent(j)-sys[i],pData->GetBinLowEdge(j+1),pData->GetBinContent(j)+sys[i]);
	 //b->SetFillColor(kGray);
	 b->SetFillStyle(0);
	 b->SetLineColor(1);
	 b->Draw();
         TBox *b2 = new TBox(pPP->GetBinLowEdge(j),pPP->GetBinContent(j)-1,pPP->GetBinLowEdge(j+1),pPP->GetBinContent(j)+1);
	 //b2->SetFillColor(65);
	 b2->SetFillStyle(0);
	 b2->SetLineColor(4);
	 b2->Draw();
      }
      pData->Draw("same");
      pMC->Draw("same");
      pPP->Draw("same");
      pGen->SetMarkerColor(4);
      pGen->SetMarkerStyle(24);
      pData->Draw("same");
      
//      pGen->Draw(" same");

      c1->cd(5+i);
      TH1D *empty2 = (TH1D*)empty->Clone("empty2");
      empty2->SetYTitle("PbPb - pp");
      empty2->SetMinimum(-5);
      empty2->SetMaximum(+29.99);
      empty2->Draw();
//      TProfile *pDiff = (TProfile*)pData->Clone("pDiff");
      TH1D *pDiff = new TH1D("pDiff","",nPtBin,PtBins);
      
      for (int j=1;j<=pData->GetNbinsX();j++)
      {
	 pDiff->SetBinContent(j,pData->GetBinContent(j)-pPP->GetBinContent(j));
  	 pDiff->SetBinError(j,sqrt(pData->GetBinError(j)*pData->GetBinError(j)+pPP->GetBinError(j)*pPP->GetBinError(j)));
         TBox *b = new TBox(pDiff->GetBinLowEdge(j),pDiff->GetBinContent(j)-sys[i],pDiff->GetBinLowEdge(j+1),pDiff->GetBinContent(j)+sys[i]);
         TBox *b2 = new TBox(pDiff->GetBinLowEdge(j),pDiff->GetBinContent(j)-sys[i],pDiff->GetBinLowEdge(j+1),pDiff->GetBinContent(j)+sys[i]);
         b->SetFillColor(TColor::GetColor("#FFEE00"));
	 b2->SetLineColor(1);
	 b2->SetFillStyle(0);
	 b->Draw();
	 b2->Draw();
      }
      
      pDiff->Draw("p same");


      TLegend *leg = new TLegend(0.3,0.6,0.9,0.9);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(pData,"PbPb","pl");
      leg->AddEntry(pPP,"pp","pl");
      leg->AddEntry(pMC,"PYTHIA+HYDJET","pl");
      if (i==0) leg->Draw();

      if (i==1) drawText("p_{T,1} > 120 GeV/c",0.22+diff,0.85);
      if (i==1) drawText("p_{T,2} >  50 GeV/c",0.22+diff,0.75);
      if (i==1) drawText("#Delta#phi_{1,2} > 5#pi/6",0.22+diff,0.65);

   }
/*
   TCanvas *c1 = new TCanvas("c1","",(ncent+1)*300,700);
   makeMultiPanelCanvas(c1,ncent+1,2,0.0,0.0,0.2,0.2,0.02);
   TH1D * empty=new TH1D("empty",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV)",axistitle[index_var].Data()),nalpha/2+1,frac);
   TH1D * empty2=new TH1D("empty2",Form(";%s;<#slash{p}_{T}^{#parallel}> (GeV)",axistitle[index_var].Data()),nalpha/2+1,frac);
   empty->Fill(0.5,1000); 
   empty2->Fill(0.5,1000);   
   if(doIntegrate){
      if(index_var==0){
         empty->SetMaximum(30); 
         empty->SetMinimum(-70); 
      }else{
         empty->SetMaximum(35); 
         empty->SetMinimum(-45); 
      }
   }else{
      empty->SetMaximum(15); 
      empty->SetMinimum(-10); 
      empty2->SetMaximum(15); 
      empty2->SetMinimum(-10); 
   }
   empty->GetXaxis()->SetTitleSize(28);
   empty->GetXaxis()->SetTitleFont(43); 
   empty->GetXaxis()->SetTitleOffset(2.2);
   empty->GetXaxis()->SetLabelSize(22);
   empty->GetXaxis()->SetLabelFont(43);
   empty->GetYaxis()->SetTitleSize(28);
   empty->GetYaxis()->SetTitleFont(43); 
   empty->GetYaxis()->SetTitleOffset(2.2);
   empty->GetYaxis()->SetLabelSize(22);
   empty->GetYaxis()->SetLabelFont(43);
   empty2->GetXaxis()->SetTitleSize(28);
   empty2->GetXaxis()->SetTitleFont(43); 
   empty2->GetXaxis()->SetTitleOffset(2.2);
   empty2->GetXaxis()->SetLabelSize(22);
   empty2->GetXaxis()->SetLabelFont(43);
   empty2->GetYaxis()->SetTitleSize(28);
   empty2->GetYaxis()->SetTitleFont(43); 
   empty2->GetYaxis()->SetTitleOffset(2.2);
   empty2->GetYaxis()->SetLabelSize(22);
   empty2->GetYaxis()->SetLabelFont(43);
   
   
   c1->cd(ncent+2); 
*/

   c1->SaveAs("results/MultiplicityDifference-DeltaPhi.C");
   c1->SaveAs("results/MultiplicityDifference-DeltaPhi.gif");
   c1->SaveAs("results/MultiplicityDifference-DeltaPhi.eps");
   c1->SaveAs("results/MultiplicityDifference-DeltaPhi.pdf");


}

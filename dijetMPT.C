#include "jettrkcorr/hiForest.h"
#include "utilities.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>
#include "../fragmentation_JEC/fragmentation_JEC.h"


#define PI 3.141592653589793238462643

void dijetMPT(double tag=0, char *infName = "/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet21_STARTHI53_LV1/merged3/HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root")
{
   // Define the input file and HiForest
   
   HiForest *c = new HiForest(infName);
 //  c->hasPFTree=0;
   c->hasPhotonTree=0;
   c->hasTowerTree=0;
   c->hasHbheTree=0;
   c->hasEbTree=0;
   c->hasGenpTree=0;
   c->hasGenParticleTree=0;   
   c->hasAk5CaloJetTree=0;
   c->hasAkPu2CaloJetTree=0;
//   c->hasAkPu3CaloJetTree=0;
   c->hasAkPu4CaloJetTree=0;
   c->hasAkPu5CaloJetTree=0;
   c->hasAkPu2JetTree=0;
   c->hasAkPu3JetTree=0;
   c->hasAkPu4JetTree=0;
   c->hasAkPu5JetTree=0;
   c->hasAkVs2PFJetTree=0;
   c->hasAkVs3PFJetTree=0;
   c->hasAkVs4PFJetTree=0;
   c->hasAkVs5PFJetTree=0;
//   c->doTrackCorrections=1;
//   c->InitTree();
   fragmentation_JEC *FF_JEC=new fragmentation_JEC(3, 1, 0, 1,1); //3rd variable is only for when do_PbPb is false
   FF_JEC->set_correction();

   
   // Output file
   TFile *output = new TFile(Form("mptTree-%.0f.root",tag),"recreate");
   
   // Output

   TTree * t = new TTree("t","gammajet");
   
   JetData data(t,1);

   HistoData histos_MergedGeneralCalo("MergedGeneral");
   HistoData histos2_MergedGeneral("MergedGeneral2");   // phi dependent corr
   
   TH1D *hWeight = new TH1D("hWeight","",1000,0,100);
   TH1D *hWeight2 = new TH1D("hWeight2","",1000,0,100);
   TH1D *hPt = new TH1D("hPt","",150,0,150);
   TH1D *hGenPt = new TH1D("hGenPt","",150,0,150);
   TH1D *hNoWPt = new TH1D("hNoWPt","",150,0,150);
   
   double trackPtBinL[10] = { 0.5,    0,   0.5, 1.0, 2.0, 4.0, 8.0};
   double trackPtBinH[10] = { 300, 0.5, 1.0, 2.0, 4.0, 8.0, 300};
   
   // Main loop
   for (int i=0;i<c->GetEntries();i++) {
      c->GetEntry(i);
      data.hiBin = c->evt.hiBin;
      if (i % 1000 == 0) cout <<i<<" / "<<c->GetEntries()<<endl;

      data.leadingJetPt = -1;      						// reconstructed leading jet
      data.subleadingJetPt = -1;						// reconstructed subleading jet
      data.leadingJetIt = -1;							// reconstructed leading jet index
      data.subleadingJetIt = -1;						// reconstructed subleading jet index
      data.genleadingJetPt = -1;						// gen level leading jet
      data.gensubleadingJetPt = -1;						// gen level subleading jet
      
      // Event selection
      if (fabs(c->evt.vz)>15) continue;
//      if (!c->selectEvent()) continue;

      // Select leading and subleading jet
      for (int j=0;j<c->akVs3Calo.nref;j++) {
         if (fabs(c->akVs3Calo.jteta[j])>2) continue;

 
         int npf=0;
         for(int ipf=0;ipf<c->pf.nPFpart;ipf++){
            if(0){
//               if(FF_JEC->passes_PF_selection(c->pf.pfPt[ipf], c->pf.pfEta[ipf], c->pf.pfPhi[ipf], c->pf.pfId[ipf], ->jteta[j], fjet->jtphi[j])) npf++;
            }else{
               if(FF_JEC->passes_PF_selection(c->pf.pfVsPt[ipf], c->pf.pfEta[ipf], c->pf.pfPhi[ipf], c->pf.pfId[ipf], c->akVs3Calo.jteta[j], c->akVs3Calo.jtphi[j])) npf++;
            }
         }

        //Then get the corrected pt for each jet before dijet selection by doing

        double corrected_pt=FF_JEC->get_corrected_pt(c->akVs3Calo.jtpt[j], npf, c->evt.hiBin);
   
        double residual_corrected_pt=FF_JEC->get_residual_corrected_pt(corrected_pt,c->evt.hiBin);
        double ptUnCor=c->akVs3Calo.jtpt[j];
	c->akVs3Calo.jtpt[j]=residual_corrected_pt;
	
         if (c->akVs3Calo.jtpt[j]>data.leadingJetPt) {
	    data.leadingJetPt = c->akVs3Calo.jtpt[j];
	    data.leadingJetPtUnCor = ptUnCor;
	    data.leadingJetEta = c->akVs3Calo.jteta[j];
	    data.leadingJetPhi = c->akVs3Calo.jtphi[j];
	    data.leadingJetRefPt = c->akVs3Calo.refpt[j];
	    data.leadingJetTrkMax = c->akVs3Calo.trackMax[j];
	    data.leadingJetIt = j;
	 }   
	 if (c->akVs3Calo.jtpt[j]>data.subleadingJetPt && c->akVs3Calo.jtpt[j] < data.leadingJetPt) {
	    data.subleadingJetPt = c->akVs3Calo.jtpt[j];
	    data.subleadingJetPtUnCor = ptUnCor;
	    data.subleadingJetEta = c->akVs3Calo.jteta[j];
	    data.subleadingJetPhi = c->akVs3Calo.jtphi[j];
	    data.subleadingJetRefPt = c->akVs3Calo.refpt[j];
	    data.subleadingJetTrkMax = c->akVs3Calo.trackMax[j];
	    data.subleadingJetIt = j;
         }
//	 if (c->akVs3Calo.jtpt[j]<data.subleadingJetPt) break;	 
      } 

      // Select generator level leading and subleading jet
      for (int j=0;j<c->akVs3Calo.ngen;j++) {
         if (fabs(c->akVs3Calo.geneta[j])>2) continue;
         if (c->akVs3Calo.genpt[j]>data.genleadingJetPt) {
	    data.genleadingJetPt = c->akVs3Calo.genpt[j];
	    data.genleadingJetEta = c->akVs3Calo.geneta[j];
	    data.genleadingJetPhi = c->akVs3Calo.genphi[j];
	 }   
	 if (c->akVs3Calo.genpt[j]>data.gensubleadingJetPt && c->akVs3Calo.genpt[j] < data.genleadingJetPt) {
	    data.gensubleadingJetPt = c->akVs3Calo.genpt[j];
	    data.gensubleadingJetEta = c->akVs3Calo.geneta[j];
	    data.gensubleadingJetPhi = c->akVs3Calo.genphi[j];
         }
//	 if (c->akVs3Calo.genpt[j]<data.gensubleadingJetPt) break;	 
      } 
      
      //if (data.subleadingJetPt<50||data.subleadingJetPt>80) continue;
      double dijetAvgPhi = getAvePhi(data.leadingJetPhi,data.subleadingJetPhi);
      double gendijetAvgPhi = getAvePhi(data.genleadingJetPhi,data.gensubleadingJetPhi);
      data.dijetPhi = dijetAvgPhi;
      data.genDijetPhi = gendijetAvgPhi;
      
      // MPT calculation
      for (int k=0;k<10;k++)
      {
         data.mpt[k] = 0;
         data.cormpt[k] = 0;
         data.cormpt2[k] = 0;
         data.genMpt[k] = 0;
         data.genPMpt[k] = 0;
      }

/*
      for (int j=0;j<c->track.nTrk;j++) {
         if (fabs(c->track.trkEta[j])>2.4) continue;
	 if (fabs(c->track.trkPt[j])<0.5) continue;
	 double mptTrk = -c->track.trkPt[j]*cos(c->track.trkPhi[j]-data.leadingJetPhi);
	 for (k=0;k<10;k++)
   	    if (c->track.trkPt[j]>trackPtBinL[j]&&c->track.trkPt[j]<=trackPtBinH[j]) data.mpt[j]+=mptTrk;
	 }
      }
*/

      data.leadingJetNTrk1=0;
      data.leadingJetNTrk2=0;
      data.leadingJetNTrk4=0;
      data.subleadingJetNTrk1=0;
      data.subleadingJetNTrk2=0;
      data.subleadingJetNTrk4=0;
      
      for (int j=0;j<c->track.nTrk;j++) {
         if (fabs(c->track.trkEta[j])>2.4) continue;
	 if ((c->track.trkPt[j])<0.5) continue;
         if (!(c->track.highPurity[j] &&
              fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3.0 &&
              fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3.0 && 
              (c->track.trkPtError[j]/c->track.trkPt[j])<0.1)) continue;
	 double dphi1 = acos(cos(c->track.trkPhi[j]-data.leadingJetPhi));
         double deta1 = fabs(c->track.trkEta[j]-data.leadingJetEta);
         double dphi2 = acos(cos(c->track.trkPhi[j]-data.subleadingJetPhi));
         double deta2 = fabs(c->track.trkEta[j]-data.subleadingJetEta);
	 
	 double dr1 = sqrt(dphi1*dphi1+deta1*deta1);
	 double dr2 = sqrt(dphi2*dphi2+deta2*deta2);
	 
	 if (dr1<0.3&&c->track.trkPt[j]>1) data.leadingJetNTrk1++;
	 if (dr2<0.3&&c->track.trkPt[j]>1) data.subleadingJetNTrk1++;
	 if (dr1<0.3&&c->track.trkPt[j]>2) data.leadingJetNTrk2++;
	 if (dr2<0.3&&c->track.trkPt[j]>2) data.subleadingJetNTrk2++;
	 if (dr1<0.3&&c->track.trkPt[j]>4) data.leadingJetNTrk4++;
	 if (dr2<0.3&&c->track.trkPt[j]>4) data.subleadingJetNTrk4++;
	 
	 double trkWt=c->getTrackCorrection(j);
         double trkWt2=0;         

	 //cout <<trkWt<<endl;
	 double mptTrk = -c->track.trkPt[j]*cos(c->track.trkPhi[j]-dijetAvgPhi)*trkWt;
	 for (int k=0;k<10;k++) {
   	    if (c->track.trkPt[j]>trackPtBinL[k]&&c->track.trkPt[j]<=trackPtBinH[k]) data.cormpt[k]+=mptTrk;
	 }

	 
	 double ptWeight = c->track.trkPt[j];
	 histos_MergedGeneralCalo.hRecoPt->Fill(c->track.trkPt[j],ptWeight);
	 histos_MergedGeneralCalo.hCorrectedPt->Fill(c->track.trkPt[j],ptWeight*trkWt);
	 histos_MergedGeneralCalo.hRecoEta->Fill(c->track.trkEta[j],ptWeight);
	 histos_MergedGeneralCalo.hCorrectedEta->Fill(c->track.trkEta[j],ptWeight*trkWt);
	 histos_MergedGeneralCalo.hRecoPhi->Fill(c->track.trkPhi[j],ptWeight);
	 histos_MergedGeneralCalo.hCorrectedPhi->Fill(c->track.trkPhi[j],ptWeight*trkWt);
	 histos_MergedGeneralCalo.hRecoDR->Fill(dr1,ptWeight);
	 histos_MergedGeneralCalo.hCorrectedDR->Fill(dr1,ptWeight*trkWt);

         hWeight->Fill(trkWt);
	 hPt->Fill(c->track.trkPt[j],trkWt);
      }

      for (int j=0;j<c->track.nParticle;j++) {
         if (fabs(c->track.pEta[j])>2.4) continue;
	 //if (fabs(c->track.pPt[j])<0.5) continue;

	       
	 double mptPTrk = -c->track.pPt[j]*cos(c->track.pPhi[j]-dijetAvgPhi);
	 
	 double dphi1 = acos(cos(c->track.pPhi[j]-data.leadingJetPhi));
         double deta1 = fabs(c->track.pEta[j]-data.leadingJetEta);
	 
	 double dr1 = sqrt(dphi1*dphi1+deta1*deta1);
	 double ptWeight = c->track.pPt[j];
	 
	 histos_MergedGeneralCalo.hGenPt->Fill(c->track.pPt[j],ptWeight);
	 histos_MergedGeneralCalo.hGenEta->Fill(c->track.pEta[j],ptWeight);
	 histos_MergedGeneralCalo.hGenPhi->Fill(c->track.pPhi[j],ptWeight);
	 histos_MergedGeneralCalo.hGenDR->Fill(dr1,ptWeight);
	 histos2_MergedGeneral.hGenPt->Fill(c->track.pPt[j],ptWeight);
	 histos2_MergedGeneral.hGenEta->Fill(c->track.pEta[j],ptWeight);
	 histos2_MergedGeneral.hGenPhi->Fill(c->track.pPhi[j],ptWeight);
	 histos2_MergedGeneral.hGenDR->Fill(dr1,ptWeight);


	 double mptTrk = -c->track.pPt[j]*cos(c->track.pPhi[j]-gendijetAvgPhi);

	 for (int k=0;k<10;k++) {
   	    if (c->track.pPt[j]>trackPtBinL[k]&&c->track.pPt[j]<=trackPtBinH[k]){
	       data.genMpt[k]+=mptTrk;
	       data.genPMpt[k]+=mptPTrk;
            }
	 }


	 hGenPt->Fill(c->track.pPt[j]);

      }

      //cout <<data.mpt<<endl;
      t->Fill();
   }
  // t->Write();
   histos_MergedGeneralCalo.calcEff();
   histos2_MergedGeneral.calcEff();
   output->Write();
   output->Close();
}

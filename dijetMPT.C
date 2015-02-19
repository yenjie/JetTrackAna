#include "jettrkcorr/hiForest.h"
#include "utilities.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>
#include "../fragmentation_JEC/fragmentation_JEC.h"


#define PI 3.141592653589793238462643

void dijetMPT(int tag=0, char *infName =
"/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet21_STARTHI53_LV1/merged3/HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root", bool
doPbPb=0, bool doPPTracking=0, int conesize=3, int jetType=0)
{
   // Define the input file and HiForest
   collisionType cMode = cPbPb;
   if (!doPbPb) cMode= cPP;
   HiForest *c = new HiForest(infName,"",cMode);
 //  c->hasPFTree=0;
   c->hasPhotonTree=0;
   c->hasTowerTree=0;
   c->hasHbheTree=0;
   c->hasEbTree=0;
   c->hasGenpTree=0;
   c->hasGenParticleTree=0;   
 //  c->hasAk5CaloJetTree=0;
   c->hasAkPu2CaloJetTree=0;
//   c->hasAkPu3CaloJetTree=0;
   c->hasAkPu4CaloJetTree=0;
   c->hasAkPu5CaloJetTree=0;
   c->hasAkPu2JetTree=0;
   //c->hasAkPu3JetTree=0;
   c->hasAkPu4JetTree=0;
   c->hasAkPu5JetTree=0;
   c->hasAkVs2PFJetTree=0;
   //c->hasAkVs3PFJetTree=0;
   c->hasAkVs4PFJetTree=0;
   c->hasAkVs5PFJetTree=0;
//   c->doTrackCorrections=1;
   c->InitTree();
   fragmentation_JEC *FF_JEC=new fragmentation_JEC(3, doPbPb, doPPTracking, 1,1,2,jetType); //3rd variable is only for when do_PbPb is false
//   FF_JEC->set_effcorrection();
   FF_JEC->set_correction();

   
   // Output file
   TFile *output = new TFile(Form("mptTreeFFCor-%d-R%d-type%d.root",tag,conesize,jetType),"recreate");
   
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
   
   Jets jets;
   
   // Main loop
   for (int i=0;i<c->GetEntries();i++) {
      c->GetEntry(i);
      data.hiBin = c->evt.hiBin;
      if (doPbPb==0) data.hiBin=0;
      if (i % 1000 == 0) cout <<i<<" / "<<c->GetEntries()<<endl;

      data.leadingJetPt = -1;      						// reconstructed leading jet with ff cor
      data.leadingJetPtUnCor = -1;      					// reconstructed leading jet
      data.leadingJetIt = -1;							// reconstructed leading jet index
      data.subleadingJetPt = -1;						// reconstructed subleading jet with ff cor
      data.subleadingJetPtUnCor = -1;						// reconstructed subleading jet
      data.subleadingJetIt = -1;						// reconstructed subleading jet index
      data.thirdleadingJetPt = -1;						// reconstructed thirdleading jet with ff cor
      data.thirdleadingJetPtUnCor = -1;						// reconstructed thirdleading jet
      data.thirdleadingJetIt = -1;						// reconstructed thirdleading jet index
      data.genleadingJetPt = -1;						// gen level leading jet
      data.gensubleadingJetPt = -1;						// gen level subleading jet
      data.genthirdleadingJetPt = -1;						// gen level thirdleading jet

      if (jetType==0) {
         // calojet
         if (conesize==2) jets = c->akVs2Calo;    
         if (conesize==3) jets = c->akVs3Calo;    
         if (conesize==4) jets = c->akVs4Calo;    
         if (conesize==5) jets = c->akVs5Calo;    
      } else if (jetType==1) {
         if (conesize==2) jets = c->akVs2PF;    
         if (conesize==3) jets = c->akVs3PF;    
         if (conesize==4) jets = c->akVs4PF;    
         if (conesize==5) jets = c->akVs5PF;    
      } else {
         if (conesize==2) jets = c->akPu2Calo;    
         if (conesize==3) jets = c->akPu3Calo;    
         if (conesize==4) jets = c->akPu4Calo;    
         if (conesize==5) jets = c->akPu5Calo;    
      }
//      if (doPbPb==1) jets = c->ak3Calo; else jets = c->akVs3Calo;
      
      // Event selection
      if (fabs(c->evt.vz)>15) continue;
//      if (!c->selectEvent()) continue;

      // Select leading and subleading jet
      for (int j=0;j<jets.nref;j++) {
         if (fabs(jets.jteta[j])>2) continue;

         double eff = 1;
         if (doPbPb) {
//	    eff = FF_JEC->get_efficiency(jets.jtpt[j],jets.jteta[j], c->evt.hiBin);
	 }
         int npf=0;
         for(int ipf=0;ipf<c->pf.nPFpart;ipf++){
            if(doPbPb==0){
               if(FF_JEC->passes_PF_selection(c->pf.pfPt[ipf], c->pf.pfEta[ipf], c->pf.pfPhi[ipf], c->pf.pfId[ipf], jets.jteta[j], jets.jtphi[j])) npf++;
            }else{
               if(FF_JEC->passes_PF_selection(c->pf.pfVsPt[ipf], c->pf.pfEta[ipf], c->pf.pfPhi[ipf], c->pf.pfId[ipf], jets.jteta[j], jets.jtphi[j])) npf++;
            }
         }
         npf/=eff;
	 
        //Then get the corrected pt for each jet before dijet selection by doing

        double corrected_pt=0;
	if (doPbPb) corrected_pt=FF_JEC->get_corrected_pt(jets.jtpt[j], npf, c->evt.hiBin);
             else  corrected_pt=FF_JEC->get_corrected_pt(jets.jtpt[j], npf);
        double residual_corrected_pt=FF_JEC->get_residual_corrected_pt(corrected_pt,c->evt.hiBin);
        double ptUnCor=jets.jtpt[j];
	jets.jtpt[j]=residual_corrected_pt;
	
         if (jets.jtpt[j]>data.leadingJetPt) {
	    data.leadingJetPt = jets.jtpt[j];
	    data.leadingJetPtUnCor = ptUnCor;
	    data.leadingJetEta = jets.jteta[j];
	    data.leadingJetPhi = jets.jtphi[j];
	    data.leadingJetRefPt = jets.refpt[j];
	    data.leadingJetTrkMax = jets.trackMax[j];
	    data.leadingJetIt = j;
	    data.leadingJetNPF = npf;
	 }   
	 if (jets.jtpt[j]>data.subleadingJetPt && jets.jtpt[j] < data.leadingJetPt) {
	    data.subleadingJetPt = jets.jtpt[j];
	    data.subleadingJetPtUnCor = ptUnCor;
	    data.subleadingJetEta = jets.jteta[j];
	    data.subleadingJetPhi = jets.jtphi[j];
	    data.subleadingJetRefPt = jets.refpt[j];
	    data.subleadingJetTrkMax = jets.trackMax[j];
	    data.subleadingJetIt = j;
	    data.subleadingJetNPF = npf;
         }
	 if (jets.jtpt[j]>data.thirdleadingJetPt && jets.jtpt[j] < data.subleadingJetPt) {
	    data.thirdleadingJetPt = jets.jtpt[j];
	    data.thirdleadingJetPtUnCor = ptUnCor;
	    data.thirdleadingJetEta = jets.jteta[j];
	    data.thirdleadingJetPhi = jets.jtphi[j];
	    data.thirdleadingJetRefPt = jets.refpt[j];
	    data.thirdleadingJetTrkMax = jets.trackMax[j];
	    data.thirdleadingJetIt = j;
	    data.thirdleadingJetNPF = npf;
         }
//	 if (jets.jtpt[j]<data.subleadingJetPt) break;	 
      } 
//      cout <<data.leadingJetPt<<endl;
      // Select generator level leading and subleading jet
      for (int j=0;j<jets.ngen;j++) {
         if (fabs(jets.geneta[j])>2) continue;
         if (jets.genpt[j]>data.genleadingJetPt) {
	    data.genleadingJetPt = jets.genpt[j];
	    data.genleadingJetEta = jets.geneta[j];
	    data.genleadingJetPhi = jets.genphi[j];
	 }   
	 if (jets.genpt[j]>data.gensubleadingJetPt && jets.genpt[j] < data.genleadingJetPt) {
	    data.gensubleadingJetPt = jets.genpt[j];
	    data.gensubleadingJetEta = jets.geneta[j];
	    data.gensubleadingJetPhi = jets.genphi[j];
         }
	 if (jets.genpt[j]>data.genthirdleadingJetPt && jets.genpt[j] < data.gensubleadingJetPt) {
	    data.genthirdleadingJetPt = jets.genpt[j];
	    data.genthirdleadingJetEta = jets.geneta[j];
	    data.genthirdleadingJetPhi = jets.genphi[j];
         }
//	 if (jets.genpt[j]<data.gensubleadingJetPt) break;	 
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
      data.leadingJetNTrk8=0;
      data.subleadingJetNTrk1=0;
      data.subleadingJetNTrk2=0;
      data.subleadingJetNTrk4=0;
      data.subleadingJetNTrk8=0;
      data.thirdleadingJetNTrk1=0;
      data.thirdleadingJetNTrk2=0;
      data.thirdleadingJetNTrk4=0;
      data.thirdleadingJetNTrk8=0;
      
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
         double dphi3 = acos(cos(c->track.trkPhi[j]-data.thirdleadingJetPhi));
         double deta3 = fabs(c->track.trkEta[j]-data.thirdleadingJetEta);
	 
	 double dr1 = sqrt(dphi1*dphi1+deta1*deta1);
	 double dr2 = sqrt(dphi2*dphi2+deta2*deta2);
	 double dr3 = sqrt(dphi3*dphi3+deta3*deta3);
	 
	 double trkWt=c->getTrackCorrection(j);
         double trkWt2=0;         

	 if (dr1<0.3&&c->track.trkPt[j]>1) data.leadingJetNTrk1++;
	 if (dr2<0.3&&c->track.trkPt[j]>1) data.subleadingJetNTrk1++;
	 if (dr3<0.3&&c->track.trkPt[j]>1) data.thirdleadingJetNTrk1++;
	 if (dr1<0.3&&c->track.trkPt[j]>2) data.leadingJetNTrk2++;
	 if (dr2<0.3&&c->track.trkPt[j]>2) data.subleadingJetNTrk2++;
	 if (dr3<0.3&&c->track.trkPt[j]>2) data.thirdleadingJetNTrk2++;
	 if (dr1<0.3&&c->track.trkPt[j]>4) data.leadingJetNTrk4++;
	 if (dr2<0.3&&c->track.trkPt[j]>4) data.subleadingJetNTrk4++;
	 if (dr3<0.3&&c->track.trkPt[j]>4) data.thirdleadingJetNTrk4++;
	 if (dr1<0.3&&c->track.trkPt[j]>8) data.leadingJetNTrk8++;
	 if (dr2<0.3&&c->track.trkPt[j]>8) data.subleadingJetNTrk8++;
	 if (dr3<0.3&&c->track.trkPt[j]>8) data.thirdleadingJetNTrk8++;

	 if (dr1<0.3&&c->track.trkPt[j]>1) data.leadingJetSumPt1+=c->track.trkPt[j]*trkWt;
	 if (dr2<0.3&&c->track.trkPt[j]>1) data.subleadingJetSumPt1+=c->track.trkPt[j]*trkWt;
	 if (dr3<0.3&&c->track.trkPt[j]>1) data.thirdleadingJetSumPt1+=c->track.trkPt[j]*trkWt;
	 if (dr1<0.3&&c->track.trkPt[j]>2) data.leadingJetSumPt2+=c->track.trkPt[j]*trkWt;
	 if (dr2<0.3&&c->track.trkPt[j]>2) data.subleadingJetSumPt2+=c->track.trkPt[j]*trkWt;
	 if (dr3<0.3&&c->track.trkPt[j]>2) data.thirdleadingJetSumPt2+=c->track.trkPt[j]*trkWt;
	 if (dr1<0.3&&c->track.trkPt[j]>4) data.leadingJetSumPt4+=c->track.trkPt[j]*trkWt;
	 if (dr2<0.3&&c->track.trkPt[j]>4) data.subleadingJetSumPt4+=c->track.trkPt[j]*trkWt;
	 if (dr3<0.3&&c->track.trkPt[j]>4) data.thirdleadingJetSumPt4+=c->track.trkPt[j]*trkWt;
	 if (dr1<0.3&&c->track.trkPt[j]>8) data.leadingJetSumPt8+=c->track.trkPt[j]*trkWt;
	 if (dr2<0.3&&c->track.trkPt[j]>8) data.subleadingJetSumPt8+=c->track.trkPt[j]*trkWt;
	 if (dr3<0.3&&c->track.trkPt[j]>8) data.thirdleadingJetSumPt8+=c->track.trkPt[j]*trkWt;
	 

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

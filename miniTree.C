#include "jettrkcorr/hiForest.h"
#include "utilities.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>

#define PI 3.141592653589793238462643

void miniTree(double tag=0, char *infName = "/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet21_STARTHI53_LV1/merged3/HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root")
{
   // Define the input file and HiForest
   HiForest *c = new HiForest(infName);
   c->hasHltTree=0;
   c->hasPFTree=0;
   c->hasPhotonTree=0;
   c->hasTowerTree=0;
   c->hasHbheTree=0;
   c->hasEbTree=0;
   c->hasGenpTree=0;
   c->hasGenParticleTree=0;   
   c->hasAk5CaloJetTree=0;
   c->hasAkPu2CaloJetTree=0;
   c->hasAkPu3CaloJetTree=0;
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
   
   // Output file
   TFile *output = new TFile(Form("output-%.0f.root",tag),"recreate");
   
   // Output
   TTree * t = new TTree("t","gammajet");
   
   JetData data(t,1);

//   HistoData histos_MergedGeneralCalo("MergedGeneral");
//   HistoData histos2_MergedGeneral("MergedGeneral2");   // phi dependent corr
   
   TH1D *hWeight = new TH1D("hWeight","",1000,0,100);
   TH1D *hWeight2 = new TH1D("hWeight2","",1000,0,100);
   TH1D *hPt = new TH1D("hPt","",150,0,150);
   TH1D *hGenPt = new TH1D("hGenPt","",150,0,150);
   TH1D *hNoWPt = new TH1D("hNoWPt","",150,0,150);
   TH1D *hRmin = new TH1D("hRmin","",100,0,10);
   // Main loop
   for (int i=0;i<c->GetEntries();i++) {
      c->GetEntry(i);
      data.hiBin = c->evt.hiBin;
      if (i % 1000 == 0) cout <<i<<" / "<<c->GetEntries()<<endl;
      data.leadingJetPt = -1;
      data.subleadingJetPt = -1;
      data.leadingJetIt = -1;
      data.subleadingJetIt = -1;
      data.genleadingJetPt = -1;
      data.gensubleadingJetPt = -1;
//      if (data.hiBin>=4) continue;
      
      // Event selection
//      if (fabs(c->evt.vz)>15) continue;
      
      
      // Select leading and subleading jet
      for (int j=0;j<c->akVs3Calo.nref;j++) {
         if (fabs(c->akVs3Calo.jteta[j])>2) continue;
         if (c->akVs3Calo.jtpt[j]>data.leadingJetPt) {
	    data.leadingJetPt = c->akVs3Calo.jtpt[j];
	    data.leadingJetEta = c->akVs3Calo.jteta[j];
	    data.leadingJetPhi = c->akVs3Calo.jtphi[j];
	    data.leadingJetIt = j;
	 }   
	 if (c->akVs3Calo.jtpt[j]>data.subleadingJetPt && c->akVs3Calo.jtpt[j] < data.leadingJetPt) {
	    data.subleadingJetPt = c->akVs3Calo.jtpt[j];
	    data.subleadingJetEta = c->akVs3Calo.jteta[j];
	    data.subleadingJetPhi = c->akVs3Calo.jtphi[j];
	    data.subleadingJetIt = j;
         }
	 if (c->akVs3Calo.jtpt[j]<data.subleadingJetPt) break;	 
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
	 if (c->akVs3Calo.genpt[j]<data.gensubleadingJetPt) break;	 
      } 
      
      if (!(data.subleadingJetPt>50&&data.leadingJetPt>120||data.gensubleadingJetPt>50&&data.genleadingJetPt>120)) continue;
      
      // MPT calculation
      data.mpt = 0;
      data.cormpt = 0;
      data.cormpt2 = 0;
      data.genMpt = 0;
      data.genPMpt = 0;
      /*
      for (int j=0;j<c->track.nTrk;j++) {
         if (fabs(c->track.trkEta[j])>2.4) continue;
	 if (fabs(c->track.trkPt[j])<0.5) continue;
	 double mptTrk = -c->track.trkPt[j]*cos(c->track.trkPhi[j]-data.leadingJetPhi);
	 data.mpt+=mptTrk;
      }
      */
      data.nTrk=0;

      for (int j=0;j<c->track.nTrk;j++) {
         if (fabs(c->track.trkEta[j])>2.4) continue;
	 if ((c->track.trkPt[j])<0.5) continue;
         if (!(c->track.highPurity[j] &&
              fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3.0 &&
              fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3.0 && 
              (c->track.trkPtError[j]/c->track.trkPt[j])<0.1)) continue;
	 data.trkWt[data.nTrk]=c->getTrackCorrection(j);
	 data.trkPt[data.nTrk]=c->track.trkPt[j];
         data.trkEta[data.nTrk]=c->track.trkEta[j];
	 data.trkPhi[data.nTrk]=c->track.trkPhi[j];
	 //data.trkRmin[data.nTrk]=c->getTrkRMin(c->track.trkPhi[j],c->track.trkEta[j],c->akVs3Calo);
	 //hPt->Fill(data.trkPt[data.nTrk],data.trkWt[data.nTrk]);
         //hRmin->Fill(data.trkRmin[data.nTrk]);
         data.nTrk++;
      }


      data.nP=0;

      for (int j=0;j<c->track.nParticle;j++) {
         if (fabs(c->track.pEta[j])>2.4) continue;
	 if (fabs(c->track.pPt[j])<0.5) continue;
//	 if (fabs(c->track.pPt[j])>1) continue;
/*
	 double mptPTrk = -c->track.pPt[j]*cos(c->track.pPhi[j]-(data.leadingJetPhi+(PI-data.subleadingJetPhi))/2);
	 data.genPMpt+=mptPTrk;

	 double dphi1 = acos(cos(c->track.pPhi[j]-data.leadingJetPhi));
         double deta1 = fabs(c->track.pEta[j]-data.leadingJetEta);
	 
	 double dr1 = sqrt(dphi1*dphi1+deta1*deta1);
	 double ptWeight = c->track.pPt[j];
	 
	 double mptTrk = -c->track.pPt[j]*cos(c->track.pPhi[j]-data.genleadingJetPhi);
	 data.genMpt+=mptTrk;
	 hGenPt->Fill(c->track.pPt[j]);
	*/ 
	 data.pPt[data.nP]=c->track.pPt[j];
         data.pEta[data.nP]=c->track.pEta[j];
	 data.pPhi[data.nP]=c->track.pPhi[j];
         data.nP++;

      }

      //cout <<data.mpt<<endl;
      t->Fill();
   }
  // t->Write();
//   histos_MergedGeneralCalo.calcEff();
//   histos2_MergedGeneral.calcEff();
   output->Write();
   output->Close();
}

#include "jettrkcorr/hiForest.h"
#include "partonFlavorTree.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>

#define PI 3.141592653589793238462643

void partonFlavorTree(double tag=0, char *infName = "/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet21_STARTHI53_LV1/merged3/HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root",collisionType cType = cPbPb)
{
   // Define the input file and HiForest
   HiForest *c = new HiForest(infName,"",cType);
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
      
      // Event selection
//      if (fabs(c->evt.vz)>15) continue;
//      if (!c->selectEvent()) continue;
      // Select leading and subleading jet
      
      data.nJet=0;
      for (int j=0;j<c->akVs3Calo.nref;j++) {
         if (fabs(c->akVs3Calo.jteta[j])>2) continue;
	 if (c->akVs3Calo.jtpt[j]<100) continue;
	 data.jetPt[data.nJet]=c->akVs3Calo.jtpt[j];
	 data.jetEta[data.nJet]=c->akVs3Calo.jteta[j];
	 data.jetPhi[data.nJet]=c->akVs3Calo.jtphi[j];
	 data.jetFlavor[data.nJet]=c->akVs3Calo.refparton_flavor[j];
	 data.refPt[data.nJet]=c->akVs3Calo.refpt[j];
//	 cout <<data.nJet<<endl;
	 
	 data.jetTrkMult1[data.nJet]=0;
	 data.jetTrkMult2[data.nJet]=0;
	 data.jetTrkMult3[data.nJet]=0;
	 data.jetTrkMult4[data.nJet]=0;
	 data.jetTrkMult5[data.nJet]=0;
	 data.jetTrkMult6[data.nJet]=0;
	 data.jetTrkMult7[data.nJet]=0;
	 data.jetTrkMult8[data.nJet]=0;
	 data.jetTrkSum1[data.nJet]=0;
	 data.jetTrkSum2[data.nJet]=0;
	 data.jetTrkSum3[data.nJet]=0;
	 data.jetTrkSum4[data.nJet]=0;
	 data.jetTrkSum5[data.nJet]=0;
	 data.jetTrkSum6[data.nJet]=0;
	 data.jetTrkSum7[data.nJet]=0;
	 data.jetTrkSum8[data.nJet]=0;
	 data.jetTrkW1[data.nJet]=0;
	 data.jetTrkW2[data.nJet]=0;
	 data.jetTrkW3[data.nJet]=0;
	 data.jetTrkW4[data.nJet]=0;
	 data.jetTrkW5[data.nJet]=0;
	 data.jetTrkW6[data.nJet]=0;
	 data.jetTrkW7[data.nJet]=0;
	 data.jetTrkW8[data.nJet]=0;
	 
	 
         for (int k=0;k<c->track.nTrk;k++) {
            if (fabs(c->track.trkEta[k])>2.4) continue;
   	    if ((c->track.trkPt[k])<4) continue;
            if (!(c->track.highPurity[k] &&
              fabs(c->track.trkDxy1[k]/c->track.trkDxyError1[k])<3.0 &&
              fabs(c->track.trkDz1[k]/c->track.trkDzError1[k])<3.0 && 
              (c->track.trkPtError[k]/c->track.trkPt[k])<0.1)) continue;

            
            double dR = c->getDR( c->track.trkEta[k], c->track.trkPhi[k], data.jetEta[data.nJet], data.jetPhi[data.nJet]);
            double dR2 = c->getDR( c->track.trkEta[k], c->track.trkPhi[k], -data.jetEta[data.nJet], data.jetPhi[data.nJet]);
	    if (dR<0.1) data.jetTrkMult1[data.nJet]++;
	    if (dR<0.2) data.jetTrkMult2[data.nJet]++;
	    if (dR<0.3) data.jetTrkMult3[data.nJet]++;
	    if (dR<0.4) data.jetTrkMult4[data.nJet]++;
	    if (dR<0.5) data.jetTrkMult5[data.nJet]++;
	    if (dR<0.6) data.jetTrkMult6[data.nJet]++;
	    if (dR<0.7) data.jetTrkMult7[data.nJet]++;
	    if (dR<0.8) data.jetTrkMult8[data.nJet]++;

	    if (dR2<0.1) data.jetTrkMult1[data.nJet]--;
	    if (dR2<0.2) data.jetTrkMult2[data.nJet]--;
	    if (dR2<0.3) data.jetTrkMult3[data.nJet]--;
	    if (dR2<0.4) data.jetTrkMult4[data.nJet]--;
	    if (dR2<0.5) data.jetTrkMult5[data.nJet]--;
	    if (dR2<0.6) data.jetTrkMult6[data.nJet]--;
	    if (dR2<0.7) data.jetTrkMult7[data.nJet]--;
	    if (dR2<0.8) data.jetTrkMult8[data.nJet]--;

	    if (dR<0.1) data.jetTrkW1[data.nJet]+=c->track.trkPt[k]*dR;
	    if (dR<0.2) data.jetTrkW2[data.nJet]+=c->track.trkPt[k]*dR;
	    if (dR<0.3) data.jetTrkW3[data.nJet]+=c->track.trkPt[k]*dR;
	    if (dR<0.4) data.jetTrkW4[data.nJet]+=c->track.trkPt[k]*dR;
	    if (dR<0.5) data.jetTrkW5[data.nJet]+=c->track.trkPt[k]*dR;
	    if (dR<0.6) data.jetTrkW6[data.nJet]+=c->track.trkPt[k]*dR;
	    if (dR<0.7) data.jetTrkW7[data.nJet]+=c->track.trkPt[k]*dR;
	    if (dR<0.8) data.jetTrkW8[data.nJet]+=c->track.trkPt[k]*dR;
/*
	    if (dR2<0.1) data.jetTrkW1[data.nJet]-=c->track.trkPt[k]*dR;
	    if (dR2<0.2) data.jetTrkW2[data.nJet]-=c->track.trkPt[k]*dR;
	    if (dR2<0.3) data.jetTrkW3[data.nJet]-=c->track.trkPt[k]*dR;
	    if (dR2<0.4) data.jetTrkW4[data.nJet]-=c->track.trkPt[k]*dR;
	    if (dR2<0.5) data.jetTrkW5[data.nJet]-=c->track.trkPt[k]*dR;
	    if (dR2<0.6) data.jetTrkW6[data.nJet]-=c->track.trkPt[k]*dR;
	    if (dR2<0.7) data.jetTrkW7[data.nJet]-=c->track.trkPt[k]*dR;
	    if (dR2<0.8) data.jetTrkW8[data.nJet]-=c->track.trkPt[k]*dR;
*/
	    if (dR<0.1) data.jetTrkSum1[data.nJet]+=c->track.trkPt[k];
	    if (dR<0.2) data.jetTrkSum2[data.nJet]+=c->track.trkPt[k];
	    if (dR<0.3) data.jetTrkSum3[data.nJet]+=c->track.trkPt[k];
	    if (dR<0.4) data.jetTrkSum4[data.nJet]+=c->track.trkPt[k];
	    if (dR<0.5) data.jetTrkSum5[data.nJet]+=c->track.trkPt[k];
	    if (dR<0.6) data.jetTrkSum6[data.nJet]+=c->track.trkPt[k];
	    if (dR<0.7) data.jetTrkSum7[data.nJet]+=c->track.trkPt[k];
	    if (dR<0.8) data.jetTrkSum8[data.nJet]+=c->track.trkPt[k];
/*

	    if (dR2<0.1) data.jetTrkSum1[data.nJet]-=c->track.trkPt[k];
	    if (dR2<0.2) data.jetTrkSum2[data.nJet]-=c->track.trkPt[k];
	    if (dR2<0.3) data.jetTrkSum3[data.nJet]-=c->track.trkPt[k];
	    if (dR2<0.4) data.jetTrkSum4[data.nJet]-=c->track.trkPt[k];
	    if (dR2<0.5) data.jetTrkSum5[data.nJet]-=c->track.trkPt[k];
	    if (dR2<0.6) data.jetTrkSum6[data.nJet]-=c->track.trkPt[k];
	    if (dR2<0.7) data.jetTrkSum7[data.nJet]-=c->track.trkPt[k];
	    if (dR2<0.8) data.jetTrkSum8[data.nJet]-=c->track.trkPt[k];
*/
	 }       
	 
	 if (data.jetTrkSum1[data.nJet]!=0) data.jetTrkW1[data.nJet]/=data.jetTrkSum1[data.nJet];
	 if (data.jetTrkSum2[data.nJet]!=0) data.jetTrkW2[data.nJet]/=data.jetTrkSum2[data.nJet];
	 if (data.jetTrkSum3[data.nJet]!=0) data.jetTrkW3[data.nJet]/=data.jetTrkSum3[data.nJet];
	 if (data.jetTrkSum4[data.nJet]!=0) data.jetTrkW4[data.nJet]/=data.jetTrkSum4[data.nJet];
	 if (data.jetTrkSum5[data.nJet]!=0) data.jetTrkW5[data.nJet]/=data.jetTrkSum5[data.nJet];
	 if (data.jetTrkSum6[data.nJet]!=0) data.jetTrkW6[data.nJet]/=data.jetTrkSum6[data.nJet];
	 if (data.jetTrkSum7[data.nJet]!=0) data.jetTrkW7[data.nJet]/=data.jetTrkSum7[data.nJet];
	 if (data.jetTrkSum8[data.nJet]!=0) data.jetTrkW8[data.nJet]/=data.jetTrkSum8[data.nJet];
	 
	 data.nJet++;
      } 

/*            
      // Select generator level leading and subleading jet
      for (int j=0;j<c->akVs3Calo.ngen;j++) {
         if (fabs(c->akVs3Calo.geneta[j])>2) continue;
      } 
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

	 data.pPt[data.nP]=c->track.pPt[j];
         data.pEta[data.nP]=c->track.pEta[j];
	 data.pPhi[data.nP]=c->track.pPhi[j];
         data.nP++;

      }
*/
      //cout <<data.mpt<<endl;
      t->Fill();
   }

   output->Write();
   output->Close();
}

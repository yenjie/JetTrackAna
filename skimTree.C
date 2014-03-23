#include "jettrkcorr/hiForest.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>

// Example of forest skim

void skimTree(char *infname = "../JetSample/hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0.root")
{
   // Define the input file and HiForest
   HiForest *c = new HiForest(infname);
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
   c->SetOutputFile("skim_jet.root");

   int filtered=0;
   // Main loop
   for (int i=0;i<c->GetEntries();i++)
   {
      c->GetEntry(i);
      if (i%1000==0) cout <<filtered<<" "<<i<<" / "<<c->GetEntries()<<endl;
      //if (c->evt.hiBin>=20) continue;
      int flag=0;
      int flag2=0;
      for (int j=0;j<c->akVs3Calo.nref;j++) {
        if (fabs(c->akVs3Calo.jteta[j])>2) continue;
	if (c->akVs3Calo.jtpt[j]>120) flag=1;
	if (c->akVs3Calo.jtpt[j]>50) flag2++;
	if (flag>=1&&flag2>=2) break;
      }
        if (flag>=1&&flag2>=2) {
	   c->FillOutput(); // Write output forest
	   filtered++;
	}  
   }

   delete c;
}

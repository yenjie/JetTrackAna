#include "../jettrkcorr/hiForest.h"
#include <TGraphAsymmErrors.h>

void plotTriggerEff(char *infname = "/mnt/hadoop/cms/store/user/luck/pp_minbiasSkim_forest_53x_2013-08-15-0155/pp_minbiasSkim_forest_53x_2013-08-15-0155.root")
{
    HiForest *c = new HiForest(infname);

    TH1D *h = new TH1D("h","",10,0,100);
    c->Draw("ak3Calo.jtpt>>h","(abs(ak3Calo.jteta)<2)");

    TH1D *h2 = new TH1D("h2","",10,0,100);
    c->Draw("ak3Calo.jtpt>>h2","HLT_PAJet80_NoJetID_v1&&(abs(ak3Calo.jteta)<2)","same");

    TGraphAsymmErrors *g = new TGraphAsymmErrors;

    g->Divide(h2,h);

    g->Draw("ap");
}

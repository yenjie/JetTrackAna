#include <TH1D.h>
#include <TH2D.h>

class treeData{

public:
   Int_t           hiBin;
   Float_t         leadingJetPt;
   Float_t         leadingJetPhi;
   Float_t         leadingJetEta;
   Float_t         subleadingJetPt;
   Float_t         subleadingJetPhi;
   Float_t         subleadingJetEta;
   Float_t         genleadingJetPt;
   Float_t         genleadingJetPhi;
   Float_t         genleadingJetEta;
   Float_t         gensubleadingJetPt;
   Float_t         gensubleadingJetPhi;
   Float_t         gensubleadingJetEta;
   Int_t           leadingJetIt;
   Int_t           subleadingJetIt;
   Float_t         mpt;
   Float_t         cormpt;
   Float_t         cormpt2;
   Float_t         genMpt;
   Float_t         genPMpt;
   Int_t           nTrk;
   Float_t         trkPt[10000];
   Float_t         trkRmin[10000];
   Float_t         trkPhi[10000];
   Float_t         trkEta[10000];
   Float_t         trkWt[10000];
   Int_t           nP;
   Float_t         pPt[10000];
   Float_t         pPhi[10000];
   Float_t         pEta[10000];

   Float_t         missingPt[10];
   Float_t         multDiff;
   Float_t         genMultDiff;
   
   TH1D *hMultVsAJ = new TH1D("hMultVsAJ","",10,0,0.5);
   TH1D *hMultVsAJGen = new TH1D("hMultVsAJGen","",10,0,0.5);

   TH2D *hMultVsEtaPhi = new TH2D("hMultVsEtaPhi",";#Delta#eta;#Delta#phi",10,-1,1,10,-1,1);
   TH1D *hMpt[10];
   
   
   TNtuple *nt;
   
   TTree *t;
   
};

void loadBranch(treeData &data, TTree *t)
{
   t->SetBranchAddress("hiBin",&data.hiBin);
   t->SetBranchAddress("leadingJetPt",&data.leadingJetPt);
   t->SetBranchAddress("leadingJetPhi",&data.leadingJetPhi);
   t->SetBranchAddress("leadingJetEta",&data.leadingJetEta);
   t->SetBranchAddress("subleadingJetPt",&data.subleadingJetPt);
   t->SetBranchAddress("subleadingJetPhi",&data.subleadingJetPhi);
   t->SetBranchAddress("subleadingJetEta",&data.subleadingJetEta);
   t->SetBranchAddress("genleadingJetPt",&data.genleadingJetPt);
   t->SetBranchAddress("genleadingJetPhi",&data.genleadingJetPhi);
   t->SetBranchAddress("genleadingJetEta",&data.genleadingJetEta);
   t->SetBranchAddress("gensubleadingJetPt",&data.gensubleadingJetPt);
   t->SetBranchAddress("gensubleadingJetPhi",&data.gensubleadingJetPhi);
   t->SetBranchAddress("gensubleadingJetEta",&data.gensubleadingJetEta);
   t->SetBranchAddress("leadingJetIt",&data.leadingJetIt);
   t->SetBranchAddress("subleadingJetIt",&data.subleadingJetIt);
   t->SetBranchAddress("mpt",&data.mpt);
   t->SetBranchAddress("cormpt",&data.cormpt);
   t->SetBranchAddress("cormpt2",&data.cormpt2);
   t->SetBranchAddress("genMpt",&data.genMpt);
   t->SetBranchAddress("genPMpt",&data.genPMpt);
   t->SetBranchAddress("nTrk",&data.nTrk);
   t->SetBranchAddress("trkPt",data.trkPt);
   t->SetBranchAddress("trkRmin",data.trkRmin);
   t->SetBranchAddress("trkPhi",data.trkPhi);
   t->SetBranchAddress("trkEta",data.trkEta);
   t->SetBranchAddress("trkWt",data.trkWt);
   t->SetBranchAddress("nP",&data.nP);
   t->SetBranchAddress("pPt",&data.pPt);
   t->SetBranchAddress("pPhi",&data.pPhi);
   t->SetBranchAddress("pEta",&data.pEta);
}

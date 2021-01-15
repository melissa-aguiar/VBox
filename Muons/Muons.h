//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 11 19:00:20 2020 by ROOT version 6.18/04
// from TChain Data/
//////////////////////////////////////////////////////////

#ifndef Muons_h
#define Muons_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Riostream.h"
#include "TRobustEstimator.h"
#include "TPrincipal.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TMatrixD.h"
#include "TVector.h"
#include "TLegend.h"
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <bitset>
#include <stdio.h>
#include <string>
#include "TGraph.h"
#include "TMultiGraph.h"

// Header file for the classes stored in the TTree if any.

class Muons {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         SmpMuon[2][64][4][7];
   Float_t         TileTrue[2][64];
   Float_t         SmpNoise[2][64][4][7];
   Float_t         TileNoise[2][64];
   Float_t         TMDBTrue[2][64];
   Float_t         TMDBNoise[2][64];

   // List of branches
   TBranch        *b_SmpMuon;   //!
   TBranch        *b_TileTrue;   //!
   TBranch        *b_SmpNoise;   //!
   TBranch        *b_TileNoise;   //!
   TBranch        *b_TMDBTrue;   //!
   TBranch        *b_TMDBNoise;   //!

   Muons(TTree *tree=0);
   virtual ~Muons();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void SalveMuons();

   Float_t mf[2][64][4];
   Float_t peso[2][64][4][7];
   Double_t caltxt[5];
   Float_t cal[2][64][4][2];
   Float_t mf_muon[2][64][4];
   Float_t mf_ruido[2][64][4];
   Float_t mev_muon[2][64][4];
   Float_t mev_ruido[2][64][4];
   Float_t est_s_muon[2][64][2];
   Float_t est_s_ruido[2][64][2];
   Float_t coefs[2][64][4];
   Float_t est_d_muon[2][64][2];
   Float_t est_d_ruido[2][64][2];
   TH1F** h0;
   TH1F** h1;
   TH1F** h2;
   TH1F** h3;
   char buf0[100],buf1[100], buf2[100], buf3[100], buf4[100], buf5[100];
   TCanvas** c;
   char leg0[100], leg1[100], leg2[100], leg3[100];
   int media0, variancia0, media1, variancia1, media2, variancia2, media3, variancia3;
   char filename[100];
   Float_t pd_s[2][64][2][1000];
   Float_t fa_s[2][64][2][1000];
   TCanvas** c2;
   //TGraph *roc;
   TMultiGraph *mroc;
   Float_t ruido[2][64];
   Float_t r_den[2][64];
   Float_t muon[2][64];
   Float_t m_den[2][64];
   Float_t ROC_pds[2][64][1000];
   Float_t ROC_fas[2][64][1000];
   Float_t ROC_pdd[2][64][1000];
   Float_t ROC_fad[2][64][1000];
   Long64_t nentries;
};

#endif

#ifdef Muons_cxx
Muons::Muons(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("Data",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("Data","");
      chain->Add("/home/mel/VBox/Muons/DataTileMuon_EB.root/Data");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

Muons::~Muons()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Muons::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Muons::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Muons::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("SmpMuon", SmpMuon, &b_SmpMuon);
   fChain->SetBranchAddress("TileTrue", TileTrue, &b_TileTrue);
   fChain->SetBranchAddress("SmpNoise", SmpNoise, &b_SmpNoise);
   fChain->SetBranchAddress("TileNoise", TileNoise, &b_TileNoise);
   fChain->SetBranchAddress("TMDBTrue", TMDBTrue, &b_TMDBTrue);
   fChain->SetBranchAddress("TMDBNoise", TMDBNoise, &b_TMDBNoise);
   Notify();
}

Bool_t Muons::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Muons::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Muons::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Muons_cxx

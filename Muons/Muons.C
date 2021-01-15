#define Muons_cxx
#include "Muons.h"
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

#include "Denoising_Muons.h"

void Muons::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Muons.C
//      root> Muons t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   gROOT->SetBatch(1);

   if (fChain == 0) return;

   nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      
   }
   Muons::SalveMuons();
   
   gROOT->SetBatch(0);
}

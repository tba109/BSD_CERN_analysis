//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan  5 19:09:12 2013 by ROOT version 5.32/01
// from TChain t/
//////////////////////////////////////////////////////////

#ifndef Efrac_by_electron_energy_h
#define Efrac_by_electron_energy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Efrac_by_electron_energy {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           feventno;
   Int_t           fspillno;
   Int_t           fspillid;
   Int_t           fcurtime;
   Int_t           fCHA[12];
   Int_t           fCHB[12];
   Int_t           trig;
   Double_t        CALSum;
   Double_t        CALProfile[20];
   Double_t        cal[2][10][50];
   Double_t        sbt[4][64];
   Float_t         BSDLate;
   Float_t         BSDLatePE;
   Float_t         BSDEarly;
   Int_t           IsWithCal;
   Int_t           Xribbon;
   Int_t           Yribbon;
   Int_t           ptype;
   Float_t         penergy;
   Int_t           windowDelay;
   Int_t           windowWidth;
   Char_t          BSDfile[47];
   Char_t          CALfile[21];

   // List of branches
   TBranch        *b_feventno;   //!
   TBranch        *b_fspillno;   //!
   TBranch        *b_fspillid;   //!
   TBranch        *b_fcurtime;   //!
   TBranch        *b_fCHA;   //!
   TBranch        *b_fCHB;   //!
   TBranch        *b_fevent;   //!
   TBranch        *b_CALSum;   //!
   TBranch        *b_CALProfile;   //!
   TBranch        *b_cal;   //!
   TBranch        *b_sbt;   //!
   TBranch        *b_BSDLate;   //!
   TBranch        *b_BSDLatePE;   //!
   TBranch        *b_BSDEarly;   //!
   TBranch        *b_IsWithCal;   //!
   TBranch        *b_Xribbon;   //!
   TBranch        *b_Yribbon;   //!
   TBranch        *b_ptype;   //!
   TBranch        *b_penergy;   //!
   TBranch        *b_windowDelay;   //!
   TBranch        *b_windowWidth;   //!
   TBranch        *b_BSDfile;   //!
   TBranch        *b_CALfile;   //!

   Efrac_by_electron_energy(TTree *tree=0);
   virtual ~Efrac_by_electron_energy();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop1(Double_t myEnergy);
   virtual void     Loop2();
   virtual void     Loop3(Double_t myEnergyE, Double_t myEnergyP);
   virtual void     Loop4(Double_t myEnergyE, Double_t myEnergyP, Int_t CalLayer);
   virtual void     Loop5(Double_t myEnergyE, Double_t myEnergyP, Int_t CalLayer);
   virtual void     Loop6(Double_t myEnergyE, Double_t myEnergyP, Int_t CalLayer);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Efrac_by_electron_energy_cxx
Efrac_by_electron_energy::Efrac_by_electron_energy(TTree *tree) : fChain(0) 
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
      f->GetObject("t",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("t","");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-18-24-46.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-18-33-03.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-18-43-53.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-18-51-40.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-18-59-26.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-19-07-23.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-19-20-32.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-19-28-15.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-19-36-26.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-19-55-23.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-20-03-16.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-20-13-18.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-20-28-18.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-20-51-09.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-20-36-11.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-21-02-29.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-13-01-07.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-13-09-10.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-13-17-02.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-13-32-47.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-13-37-35.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-15-09-26.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-15-31-00.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-15-35-09.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-15-44-33.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-15-53-21.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-15-57-47.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-17-31-19.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-17-38-59.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-17-45-45.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-17-51-53.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-17-59-01.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-21-12-45.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-21-18-37.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-21-25-05.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-21-31-44.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-30-2012-21-38-17.root/t");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

Efrac_by_electron_energy::~Efrac_by_electron_energy()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Efrac_by_electron_energy::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Efrac_by_electron_energy::LoadTree(Long64_t entry)
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

void Efrac_by_electron_energy::Init(TTree *tree)
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

   fChain->SetBranchAddress("feventno", &feventno, &b_feventno);
   fChain->SetBranchAddress("fspillno", &fspillno, &b_fspillno);
   fChain->SetBranchAddress("fspillid", &fspillid, &b_fspillid);
   fChain->SetBranchAddress("fcurtime", &fcurtime, &b_fcurtime);
   fChain->SetBranchAddress("fCHA", fCHA, &b_fCHA);
   fChain->SetBranchAddress("fCHB", fCHB, &b_fCHB);
   fChain->SetBranchAddress("trig", &trig, &b_fevent);
   fChain->SetBranchAddress("CALSum", &CALSum, &b_CALSum);
   fChain->SetBranchAddress("CALProfile", CALProfile, &b_CALProfile);
   fChain->SetBranchAddress("cal", cal, &b_cal);
   fChain->SetBranchAddress("sbt", sbt, &b_sbt);
   fChain->SetBranchAddress("BSDLate", &BSDLate, &b_BSDLate);
   fChain->SetBranchAddress("BSDLatePE", &BSDLatePE, &b_BSDLatePE);
   fChain->SetBranchAddress("BSDEarly", &BSDEarly, &b_BSDEarly);
   fChain->SetBranchAddress("IsWithCal", &IsWithCal, &b_IsWithCal);
   fChain->SetBranchAddress("Xribbon", &Xribbon, &b_Xribbon);
   fChain->SetBranchAddress("Yribbon", &Yribbon, &b_Yribbon);
   fChain->SetBranchAddress("ptype", &ptype, &b_ptype);
   fChain->SetBranchAddress("penergy", &penergy, &b_penergy);
   fChain->SetBranchAddress("windowDelay", &windowDelay, &b_windowDelay);
   fChain->SetBranchAddress("windowWidth", &windowWidth, &b_windowWidth);
   fChain->SetBranchAddress("BSDfile", BSDfile, &b_BSDfile);
   fChain->SetBranchAddress("CALfile", CALfile, &b_CALfile);
   Notify();
}

Bool_t Efrac_by_electron_energy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Efrac_by_electron_energy::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Efrac_by_electron_energy::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Efrac_by_electron_energy_cxx

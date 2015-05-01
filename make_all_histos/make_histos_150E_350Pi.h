//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 18 13:03:39 2012 by ROOT version 5.32/01
// from TChain t/
//////////////////////////////////////////////////////////

#ifndef make_histos_150E_350Pi_h
#define make_histos_150E_350Pi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class make_histos_150E_350Pi {
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

   make_histos_150E_350Pi(TTree *tree=0);
   virtual ~make_histos_150E_350Pi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef make_histos_150E_350Pi_cxx
make_histos_150E_350Pi::make_histos_150E_350Pi(TTree *tree) : fChain(0) 
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
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-28-2012-00-00-18.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-28-2012-00-08-38.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-28-2012-00-13-56.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-28-2012-00-18-29.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-28-2012-00-24-21.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-28-2012-00-28-39.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-28-2012-00-33-04.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-13-24-15.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-13-24-15.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-10-16-35.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-10-24-38.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-10-30-54.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-13-46-34.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-13-37-37.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-19-36-38.root/t");
      chain->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/BSDCALCombined_11-27-2012-19-49-55.root/t");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

make_histos_150E_350Pi::~make_histos_150E_350Pi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t make_histos_150E_350Pi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t make_histos_150E_350Pi::LoadTree(Long64_t entry)
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

void make_histos_150E_350Pi::Init(TTree *tree)
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

Bool_t make_histos_150E_350Pi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void make_histos_150E_350Pi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t make_histos_150E_350Pi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef make_histos_150E_350Pi_cxx

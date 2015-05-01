#define make_histos_150E_350Pi_cxx
#include "make_histos_150E_350Pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>

void make_histos_150E_350Pi::Loop()
{
    
  if (fChain == 0) return;
  
  TH1F * hPi_CALSum = new TH1F("hPi_CALSum","",200,-1000,30000);
  TH1F * hE_CALSum = new TH1F("hE_CALSum","",500,-1000,30000);
  TH1F * hPi_myBSDLatePE = new TH1F("hPi_myBSDLatePE","",200,-100,8000);
  TH1F * hE_myBSDLatePE = new TH1F("hE_myBSDLatePE","",200,-100,8000);
  TH1F * hE_Efrac = new TH1F("hE_Efrac","",500,-0.2,3);
  TH1F * hPi_Efrac = new TH1F("hPi_Efrac","",500,-0.2,3);

  Int_t NE150 = 0;
  Int_t NPi350 = 0;    
    
  Double_t gainA[12] = {4.23,4.23,3.92,3.09,5.81,5.03,3.17,5.30,4.31,3.90,4.04,3.68};
  Double_t myfCHA[12];
  Double_t myBSDLatePE = 0;

  fChain->SetBranchStatus("*",1);
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
      if( !IsWithCal ) continue;

      myBSDLatePE = 0;
      for( int i = 0; i < 12; i++ )
	myfCHA[i] = -1;
            
      for( int i = 1; i <= 10; i++ )
	{
	  if( fCHA[i] <= 0 )
	    {
	      myfCHA[i] = 2048;
	    }
	  else
	    {
	      myfCHA[i] = fCHA[i];
	    }
	}
      
      // Simulate CAL trigger by requiring 6 consecutive layers with > 40 MeV
      Int_t Ntrig_layers = 0;
      for( int i = 0; i < 20; i++ )
	{
	  if( CALProfile[i]/9 > 40 )
	    Ntrig_layers++;
	  else
	    Ntrig_layers = 0;
	  if( Ntrig_layers >= 6 )
	    break;
	}
      if( Ntrig_layers < 6 ) continue;

      
      for( int i = 1; i <= 10; i++ )
	myBSDLatePE += myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19);
      
      // Check if we've got enough of that particular flavor, i.e. all electrons 
      // need to be weighted equally relative to pions
      if( ptype == 2 && penergy == 350 && NPi350 >= 9000) continue; 
      if( ptype == 1 && penergy == 150 && NE150 >= 9000 ) continue;
     
      // Fill Efrac for electrons
      if( ptype == 1 && penergy == 150 )
	{
	  hE_Efrac->Fill(myBSDLatePE/CALSum);
	  hE_CALSum->Fill(CALSum);
	  hE_myBSDLatePE->Fill(myBSDLatePE);
	}
      
      // Fill Efrac for pions
      if( ptype == 2 && penergy == 350 )
	{
	  hPi_Efrac->Fill(myBSDLatePE/CALSum);
	  hPi_CALSum->Fill(CALSum);
	  hPi_myBSDLatePE->Fill(myBSDLatePE);
	}
	  
      // Increment counts, since everything has to be equally weighted
      if( ptype == 2 && penergy == 350 ) NPi350++;
      if( ptype == 1 && penergy == 150 ) NE150++;
      
    }
  
  TCanvas * cv1 = new TCanvas("cv1");
  hE_CALSum->Draw();

  TCanvas * cv2 = new TCanvas("cv2");
  hPi_CALSum->Draw();

  TCanvas * cv3 = new TCanvas("cv3");
  hE_myBSDLatePE->Draw();

  TCanvas * cv4 = new TCanvas("cv4");
  hPi_myBSDLatePE->Draw();

  TCanvas * cv5 = new TCanvas("cv5");
  hE_Efrac->Draw();

  TCanvas * cv6 = new TCanvas("cv6");
  hPi_Efrac->Draw();

  TFile f("histos_150E_350Pi.root","new");
  hE_CALSum->Write();
  hPi_CALSum->Write();
  hE_myBSDLatePE->Write();
  hPi_myBSDLatePE->Write();
  hE_Efrac->Write();
  hPi_Efrac->Write();
  f.Close();

}

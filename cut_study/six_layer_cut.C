#define six_layer_cut_cxx
#include "six_layer_cut.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLine.h>
#include <TF1.h>

using namespace std;

void six_layer_cut::Loop1()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  fChain->SetBranchStatus("*",1);
  
  Long64_t nbytes = 0, nb = 0;
  
  Double_t gainA[12] = {4.23,4.23,3.92,3.09,5.81,5.03,3.17,5.30,4.31,3.90,4.04,3.68};
  Double_t myfCHA[12];
  Double_t myBSDLatePE = 0;
  for( int i = 0; i < 12; i++ )
    myfCHA[i] = -1;
  
  TH1F * hE_no6_CAL = new TH1F("hE_no6_CAL","",150,-1000,20000);
  TH1F * hE_no6_myBSDLatePE = new TH1F("hE_no6_myBSDLatePE","",150,-1000,8000);
  TH1F * hE_with6_CAL = new TH1F("hE_with6_CAL","",150,-1000,20000);
  TH1F * hE_with6_myBSDLatePE = new TH1F("hE_with6_myBSDLatePE","",150,-1000,8000);
  TH1F * hPi_no6_CAL = new TH1F("hPi_no6_CAL","",150,-1000,20000);
  TH1F * hPi_no6_myBSDLatePE = new TH1F("hPi_no6_myBSDLatePE","",150,-1000,8000);
  TH1F * hPi_with6_CAL = new TH1F("hPi_with6_CAL","",150,-1000,20000);
  TH1F * hPi_with6_myBSDLatePE = new TH1F("hPi_with6_myBSDLatePE","",150,-1000,8000);
  
  Int_t NE50_no6 = 0;
  Int_t NE75_no6 = 0;
  Int_t NE100_no6 = 0;
  Int_t NE125_no6 = 0;
  Int_t NE150_no6 = 0;
  Int_t NE175_no6 = 0;
  
  Int_t NPi250_no6 = 0;
  Int_t NPi300_no6 = 0;
  Int_t NPi350_no6 = 0;    
  
  Int_t NE50_with6 = 0;
  Int_t NE75_with6 = 0;
  Int_t NE100_with6 = 0;
  Int_t NE125_with6 = 0;
  Int_t NE150_with6 = 0;
  Int_t NE175_with6 = 0;

  Int_t NPi250_with6 = 0;
  Int_t NPi300_with6 = 0;
  Int_t NPi350_with6 = 0;    
  
  Double_t Ncut_Total = 0;
  Double_t Ncut_BSDLatePE = 0;
  Double_t Ncut_IsWithCal = 0;
  Double_t Ncut_CALSum = 0;
  Double_t Ncut_trig_layers = 0;
  Double_t Ncut_neg = 0;
  Double_t Nsurvive = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
      Ncut_Total++;
      
      if( !IsWithCal ) 
	{
	  Ncut_IsWithCal++;
	  continue;
	}
      
      myBSDLatePE = 0;
      for( int i = 1; i <= 10; i++ )
	myfCHA[i] = -1;
      
      Bool_t neg = 0;
      for( int i = 1; i <= 10; i++ )
	{
	  if( fCHA[i] <= 0 )
	    {
	      neg = kTRUE;
	      myfCHA[i] = 2048;
	    }
	  else
	    myfCHA[i] = fCHA[i];
	}
      if( neg == kTRUE ) 
	{
	  Ncut_neg++;
	  // continue;
	}
      
      if( myBSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  // continue;
	}
      
      if( CALSum < 0 ) 
	{
	  Ncut_CALSum++;
	  continue;
	}
      
      myBSDLatePE = 0;
      // Calculate myBSDLatePE
      for( int i = 1; i <= 10; i++ )
	myBSDLatePE += myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19); 
      
      ///////////////////////////////////////////////////////////////////////////////
      // Now we start the weighting game!
      // For Loop1, we want to see how the 6 layer cut down selects from a set of
      // otherwise viable events for pions and electrons
      // Since the 6 layer cuts much more strongly on lower energies
      //////////////////////////////////////////////////////////////////////////////
      
      // Want to weight all energies equally before the 6 layer cut
      if( (ptype == 1) && (penergy == 50) && (NE50_no6 >= 2585) ) continue;
      if( (ptype == 1) && (penergy == 75) && (NE75_no6 >= 2585) ) continue;
      if( (ptype == 1) && (penergy == 100) && (NE100_no6 >= 2585) ) continue;
      if( (ptype == 1) && (penergy == 125) && (NE125_no6 >= 2585) ) continue;
      if( (ptype == 1) && (penergy == 150) && (NE150_no6 >= 2585) ) continue;
      if( (ptype == 1) && (penergy == 175) && (NE175_no6 >= 2585) ) continue;
      
      if( (ptype == 2) && (penergy == 250) && (NPi250_no6 >= 18000) ) continue; 
      if( (ptype == 2) && (penergy == 300) && (NPi300_no6 >= 18000) ) continue; 
      if( (ptype == 2) && (penergy == 350) && (NPi350_no6 >= 18000) ) continue; 
      
      // Fill no six layer histo for electrons
      if( ptype == 1 )
	{
	  hE_no6_CAL->Fill(CALSum);
	  hE_no6_myBSDLatePE->Fill(myBSDLatePE);
	  if( penergy == 175 ) 
	    {
	      NE175_no6++;
	    }
	  if( penergy == 150 ) 
	    {
	      NE150_no6++;
	    }
	  if( penergy == 125 ) 
	    {
	      NE125_no6++;
	    }
	  if( penergy == 100 ) 
	    {
	      NE100_no6++;
	    }
	  if( penergy == 75 ) 
	    {
	      NE75_no6++;
	    }
	  if( penergy == 50 )
	    {
	      NE50_no6++;
	    }
	}
      
      // Fill no six layer histo for pions
      if( ptype == 2 )
	{
	  hPi_no6_CAL->Fill(CALSum);
	  hPi_no6_myBSDLatePE->Fill(myBSDLatePE);
	  if( penergy == 350 ) 
	    {
	      NPi350_no6++;
	    }
	  if( penergy == 300 ) 
	    {
	      NPi300_no6++;
	    }
	  if( penergy == 250 ) 
	    {
	      NPi250_no6++;
	    }
	}
      
      /////////////////////////////////////////////////////////////////////////////
      // Now apply cut to downselected events
      // Simulate CAL trigger by requiring 6 consecutive layers with > 40 MeV
      /////////////////////////////////////////////////////////////////////////////
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
      if( Ntrig_layers < 6 )
	{
	  Ncut_trig_layers++;
	  continue;
	}
      else
	{
	  // Fill no six layer histo for electrons
	  if( ptype == 1 )
	    {
	      hE_with6_CAL->Fill(CALSum);
	      hE_with6_myBSDLatePE->Fill(myBSDLatePE);
	      
	      if( penergy == 175 ) 
		{
		  NE175_with6++;
		}
	      if( penergy == 150 ) 
		{
		  NE150_with6++;
		}
	      if( penergy == 125 ) 
		{
		  NE125_with6++;
		}
	      if( penergy == 100 ) 
		{
		  NE100_with6++;
		}
	      if( penergy == 75 ) 
		{
		  NE75_with6++;
		}
	      if( penergy == 50 )
		{
		  NE50_with6++;
		}
	    }
	  
	  // Fill no six layer histo for pions
	  if( ptype == 2 )
	    {
	      hPi_with6_CAL->Fill(CALSum);
	      hPi_with6_myBSDLatePE->Fill(myBSDLatePE);
	      
	      if( penergy == 350 ) 
		{
		  NPi350_with6++;
		}
	      if( penergy == 300 ) 
		{
		  NPi300_with6++;
		}
	      if( penergy == 250 ) 
		{ 
		  NPi250_with6++;
		}
	    }
	}
      Nsurvive++;
    }
  
  // Report on cuts
  Int_t Nremain = 0;

  cout << "===========================================================" << endl;
  cout << "Ncut_Total = " << Ncut_Total << endl;
  Nremain = Ncut_Total;
  
  cout << "Ncut_IsWithCal = " << Ncut_IsWithCal << ", " << 100*Ncut_IsWithCal/Nremain << endl;
  Nremain -= Ncut_IsWithCal;
  
  cout << "Ncut_neg = " << Ncut_neg << ", " << 100*Ncut_neg/Nremain << endl;
  // Nremain -= Ncut_neg;
  
  cout << "Ncut_BSDLatePE = " << Ncut_BSDLatePE << ", " << 100*Ncut_BSDLatePE/Nremain << endl;
  Nremain -= Ncut_BSDLatePE;
  
  cout << "Ncut_CALSum = " << Ncut_CALSum << ", " << 100*Ncut_CALSum/Nremain << endl;
  Nremain -= Ncut_CALSum;
  
  cout << "Ncut_trig_layers = " << Ncut_trig_layers << ", " << 100*Ncut_trig_layers/Ncut_Total << endl;
  Nremain -= Ncut_trig_layers;
  
  cout << "Nsurvive = " << Nsurvive << ", " << 100*Nsurvive/Ncut_Total << endl;
  cout << "Nremain = " << Nremain << ", " << 100*Nremain/Ncut_Total << endl;
  
  cout << "============================================================" << endl;
  cout << "NE50_no6 = " << NE50_no6 << endl;
  cout << "NE75_no6 = " << NE75_no6 << endl;
  cout << "NE75_with6 = " << NE75_with6 << endl;
  cout << "NE100_no6 = " << NE100_no6 << endl;
  cout << "NE100_with6 = " << NE100_with6 << endl;
  cout << "NE125_no6 = " << NE125_no6 << endl;
  cout << "NE125_with6 = " << NE125_with6 << endl;
  cout << "NE150_no6 = " << NE150_no6 << endl;
  cout << "NE150_with6 = " << NE150_with6 << endl;
  cout << "NE175_no6 = " << NE175_no6 << endl;
  cout << "NE175_with6 = " << NE175_with6 << endl;
  
  cout << "=============================================================" << endl;
  cout << "NPi250_no6 = " << NPi250_no6 << endl;
  cout << "NPi250_with6 = " << NPi250_with6 << endl;
  cout << "NPi300_no6 = " << NPi300_no6 << endl;
  cout << "NPi300_with6 = " << NPi300_with6 << endl;
  cout << "NPi350_no6 = " << NPi350_no6 << endl;
  cout << "NPi350_with6 = " << NPi350_with6 << endl;
  
  
  // Draw
  TCanvas * cv1 = new TCanvas("cv1");
  cv1->Draw();
  
  hE_no6_CAL->SetLineColor(kBlue);
  hE_with6_CAL->SetLineColor(kBlack);
  hE_no6_CAL->GetXaxis()->SetTitle("CALSum for Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_CAL->GetXaxis()->CenterTitle();
  
  hE_no6_CAL->Draw();
  hE_with6_CAL->Draw("sames");
  
  TCanvas * cv2 = new TCanvas("cv2");
  cv2->Draw();
  
  hPi_no6_CAL->SetLineColor(kBlue);
  hPi_with6_CAL->SetLineColor(kBlack);
  hPi_no6_CAL->GetXaxis()->SetTitle("CALSum for Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_CAL->GetXaxis()->CenterTitle();
  
  hPi_no6_CAL->Draw();
  hPi_with6_CAL->Draw("sames");
  
  TCanvas * cv3 = new TCanvas("cv3");
  cv3->Draw();

  hE_no6_myBSDLatePE->SetLineColor(kBlue);
  hE_with6_myBSDLatePE->SetLineColor(kBlack);
  hE_no6_myBSDLatePE->GetXaxis()->SetTitle("myBSDLatePE for Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE->GetXaxis()->CenterTitle();
  
  hE_no6_myBSDLatePE->Draw();
  hE_with6_myBSDLatePE->Draw("sames");
  
  TCanvas * cv4 = new TCanvas("cv4");
  cv4->Draw();
  
  hPi_no6_myBSDLatePE->SetLineColor(kBlue);
  hPi_with6_myBSDLatePE->SetLineColor(kBlack);
  hPi_no6_myBSDLatePE->GetXaxis()->SetTitle("myBSDLatePE for Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_myBSDLatePE->GetXaxis()->CenterTitle();
  
  hPi_no6_myBSDLatePE->Draw();
  hPi_with6_myBSDLatePE->Draw("sames");
  
}

// Make plots for each energy so we can see how it effects each energy
void six_layer_cut::Loop2()
{
if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  fChain->SetBranchStatus("*",1);
  
  Long64_t nbytes = 0, nb = 0;
  
  Double_t gainA[12] = {4.23,4.23,3.92,3.09,5.81,5.03,3.17,5.30,4.31,3.90,4.04,3.68};
  Double_t myfCHA[12];
  Double_t myBSDLatePE = 0;
  for( int i = 0; i < 12; i++ )
    myfCHA[i] = -1;
  
  
  // 50 GeV electrons
  TH1F * hE_no6_CAL_50 = new TH1F("hE_no6_CAL_50","",150,-1000,20000);
  TH1F * hE_with6_CAL_50 = new TH1F("hE_with6_CAL_50","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_50 = new TH1F("hE_no6_myBSDLatePE_50","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_50 = new TH1F("hE_with6_myBSDLatePE_50","",150,-1000,8000);
  
  // 75 GeV electrons
  TH1F * hE_no6_CAL_75 = new TH1F("hE_no6_CAL_75","",150,-1000,20000);
  TH1F * hE_with6_CAL_75 = new TH1F("hE_with6_CAL_75","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_75 = new TH1F("hE_no6_myBSDLatePE_75","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_75 = new TH1F("hE_with6_myBSDLatePE_75","",150,-1000,8000);
  
  // 100 GeV electrons
  TH1F * hE_no6_CAL_100 = new TH1F("hE_no6_CAL_100","",150,-1000,20000);
  TH1F * hE_with6_CAL_100 = new TH1F("hE_with6_CAL_100","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_100 = new TH1F("hE_no6_myBSDLatePE_100","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_100 = new TH1F("hE_with6_myBSDLatePE_100","",150,-1000,8000);
  
  // 125 GeV electrons
  TH1F * hE_no6_CAL_125 = new TH1F("hE_no6_CAL_125","",150,-1000,20000);
  TH1F * hE_with6_CAL_125 = new TH1F("hE_with6_CAL_125","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_125 = new TH1F("hE_no6_myBSDLatePE_125","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_125 = new TH1F("hE_with6_myBSDLatePE_125","",150,-1000,8000);
  
  // 150 GeV electrons
  TH1F * hE_no6_CAL_150 = new TH1F("hE_no6_CAL_150","",150,-1000,20000);
  TH1F * hE_with6_CAL_150 = new TH1F("hE_with6_CAL_150","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_150 = new TH1F("hE_no6_myBSDLatePE_150","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_150 = new TH1F("hE_with6_myBSDLatePE_150","",150,-1000,8000);
  
  // 175 GeV electrons
  TH1F * hE_no6_CAL_175 = new TH1F("hE_no6_CAL_175","",150,-1000,20000);
  TH1F * hE_with6_CAL_175 = new TH1F("hE_with6_CAL_175","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_175 = new TH1F("hE_no6_myBSDLatePE_175","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_175 = new TH1F("hE_with6_myBSDLatePE_175","",150,-1000,8000);

  // 250 GeV pions
  TH1F * hPi_no6_CAL_250 = new TH1F("hPi_no6_CAL_250","",150,-1000,20000);
  TH1F * hPi_with6_CAL_250 = new TH1F("hPi_with6_CAL_250","",150,-1000,20000);

  TH1F * hPi_no6_myBSDLatePE_250 = new TH1F("hPi_no6_myBSDLatePE_250","",150,-1000,8000);
  TH1F * hPi_with6_myBSDLatePE_250 = new TH1F("hPi_with6_myBSDLatePE_250","",150,-1000,8000);

  // 300 GeV pions
  TH1F * hPi_no6_CAL_300 = new TH1F("hPi_no6_CAL_300","",150,-1000,20000);
  TH1F * hPi_with6_CAL_300 = new TH1F("hPi_with6_CAL_300","",150,-1000,20000);

  TH1F * hPi_no6_myBSDLatePE_300 = new TH1F("hPi_no6_myBSDLatePE_300","",150,-1000,8000);
  TH1F * hPi_with6_myBSDLatePE_300 = new TH1F("hPi_with6_myBSDLatePE_300","",150,-1000,8000);

  // 350 GeV pions
  TH1F * hPi_no6_CAL_350 = new TH1F("hPi_no6_CAL_350","",150,-1000,20000);
  TH1F * hPi_with6_CAL_350 = new TH1F("hPi_with6_CAL_350","",150,-1000,20000);

  TH1F * hPi_no6_myBSDLatePE_350 = new TH1F("hPi_no6_myBSDLatePE_350","",150,-1000,8000);
  TH1F * hPi_with6_myBSDLatePE_350 = new TH1F("hPi_with6_myBSDLatePE_350","",150,-1000,8000);

  
  Int_t NE50_no6 = 0;
  Int_t NE75_no6 = 0;
  Int_t NE100_no6 = 0;
  Int_t NE125_no6 = 0;
  Int_t NE150_no6 = 0;
  Int_t NE175_no6 = 0;
  
  Int_t NPi250_no6 = 0;
  Int_t NPi300_no6 = 0;
  Int_t NPi350_no6 = 0;    
  
  Int_t NE50_with6 = 0;
  Int_t NE75_with6 = 0;
  Int_t NE100_with6 = 0;
  Int_t NE125_with6 = 0;
  Int_t NE150_with6 = 0;
  Int_t NE175_with6 = 0;

  Int_t NPi250_with6 = 0;
  Int_t NPi300_with6 = 0;
  Int_t NPi350_with6 = 0;    
  
  Double_t Ncut_Total = 0;
  Double_t Ncut_BSDLatePE = 0;
  Double_t Ncut_IsWithCal = 0;
  Double_t Ncut_CALSum = 0;
  Double_t Ncut_trig_layers = 0;
  Double_t Ncut_neg = 0;
  Double_t Nsurvive = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
      Ncut_Total++;
      
      if( !IsWithCal ) 
	{
	  Ncut_IsWithCal++;
	  continue;
	}
      
      myBSDLatePE = 0;
      for( int i = 1; i <= 10; i++ )
	myfCHA[i] = -1;
      
      Bool_t neg = 0;
      for( int i = 1; i <= 10; i++ )
	{
	  if( fCHA[i] <= 0 )
	    {
	      neg = kTRUE;
	      myfCHA[i] = 2048;
	    }
	  else
	    myfCHA[i] = fCHA[i];
	}
      if( neg == kTRUE ) 
	{
	  Ncut_neg++;
	  // continue;
	}
      
      if( myBSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  // continue;
	}
      
      if( CALSum < 0 ) 
	{
	  Ncut_CALSum++;
	  continue;
	}
      
      myBSDLatePE = 0;
      // Calculate myBSDLatePE
      for( int i = 1; i <= 10; i++ )
	myBSDLatePE += myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19); 
      
      ///////////////////////////////////////////////////////////////////////////////
      // Now we start the weighting game!
      // For Loop1, we want to see how the 6 layer cut down selects from a set of
      // otherwise viable events for pions and electrons
      // Since the 6 layer cuts much more strongly on lower energies
      //////////////////////////////////////////////////////////////////////////////
      
      // Want to weight all energies equally before the 6 layer cut
      if( (ptype == 1) && (penergy == 50) && (NE50_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 75) && (NE75_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 100) && (NE100_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 125) && (NE125_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 150) && (NE150_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 175) && (NE175_no6 > 7000) ) continue;
      
      if( (ptype == 2) && (penergy == 250) && (NPi250_no6 > 18000) ) continue; 
      if( (ptype == 2) && (penergy == 300) && (NPi300_no6 > 18000) ) continue; 
      if( (ptype == 2) && (penergy == 350) && (NPi350_no6 > 18000) ) continue; 
      
      // Fill no six layer histo for electrons
      if( ptype == 1 )
	{
	  if( penergy == 175 ) 
	    {
	      hE_no6_CAL_175->Fill(CALSum);
	      hE_no6_myBSDLatePE_175->Fill(myBSDLatePE);
	      NE175_no6++;
	    }
	  if( penergy == 150 ) 
	    {
	      hE_no6_CAL_150->Fill(CALSum);
	      hE_no6_myBSDLatePE_150->Fill(myBSDLatePE);
	      NE150_no6++;
	    }
	  if( penergy == 125 ) 
	    {
	      hE_no6_CAL_125->Fill(CALSum);
	      hE_no6_myBSDLatePE_125->Fill(myBSDLatePE);
	      NE125_no6++;
	    }
	  if( penergy == 100 ) 
	    {
	      hE_no6_CAL_100->Fill(CALSum);
	      hE_no6_myBSDLatePE_100->Fill(myBSDLatePE);
	      NE100_no6++;
	    }
	  if( penergy == 75 ) 
	    {
	      hE_no6_CAL_75->Fill(CALSum);
	      hE_no6_myBSDLatePE_75->Fill(myBSDLatePE);
	      NE75_no6++;
	    }
	  if( penergy == 50 ) 
	    {
	      hE_no6_CAL_50->Fill(CALSum);
	      hE_no6_myBSDLatePE_50->Fill(myBSDLatePE);
	      NE50_no6++;
	    }
	}
      
      // Fill no six layer histo for pions
      if( ptype == 2 )
	{
	  if( penergy == 350 ) 
	    {
	      hPi_no6_CAL_350->Fill(CALSum);
	      hPi_no6_myBSDLatePE_350->Fill(myBSDLatePE);
	      NPi350_no6++;
	    }
	  if( penergy == 300 ) 
	    {
	      hPi_no6_CAL_300->Fill(CALSum);
	      hPi_no6_myBSDLatePE_300->Fill(myBSDLatePE);
	      NPi300_no6++;
	    }
	  if( penergy == 250 ) 
	    {
	      hPi_no6_CAL_250->Fill(CALSum);
	      hPi_no6_myBSDLatePE_250->Fill(myBSDLatePE);
	      NPi250_no6++;
	    }
	}
      
      /////////////////////////////////////////////////////////////////////////////
      // Now apply cut to downselected events
      // Simulate CAL trigger by requiring 6 consecutive layers with > 40 MeV
      /////////////////////////////////////////////////////////////////////////////
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
      if( Ntrig_layers < 6 )
	{
	  Ncut_trig_layers++;
	  continue;
	}
      else
	{
	  // Fill no six layer histo for electrons
	  if( ptype == 1 )
	    {
	      if( penergy == 175 ) 
		{
		  hE_with6_CAL_175->Fill(CALSum);
		  hE_with6_myBSDLatePE_175->Fill(myBSDLatePE);
		  NE175_with6++;
		}
	      if( penergy == 150 ) 
		{
		  hE_with6_CAL_150->Fill(CALSum);
		  hE_with6_myBSDLatePE_150->Fill(myBSDLatePE);
		  NE150_with6++;
		}
	      if( penergy == 125 ) 
		{ hE_with6_CAL_125->Fill(CALSum);
		  hE_with6_myBSDLatePE_125->Fill(myBSDLatePE);
		  NE125_with6++;
		}
	      if( penergy == 100 ) 
		{ hE_with6_CAL_100->Fill(CALSum);
		  hE_with6_myBSDLatePE_100->Fill(myBSDLatePE);
		  NE100_with6++;
		}
	      if( penergy == 75 ) 
		{ hE_with6_CAL_75->Fill(CALSum);
		  hE_with6_myBSDLatePE_75->Fill(myBSDLatePE);
		  NE75_with6++;
		}
	      if( penergy == 50 ) 
		{ hE_with6_CAL_50->Fill(CALSum);
		  hE_with6_myBSDLatePE_50->Fill(myBSDLatePE);
		  NE75_with6++;
		}
	    }
	  
	  // Fill no six layer histo for pions
	  if( ptype == 2 )
	    {
	      if( penergy == 350 ) 
		{
		  hPi_with6_CAL_350->Fill(CALSum);
		  hPi_with6_myBSDLatePE_350->Fill(myBSDLatePE);
		  NPi350_with6++;
		}
	      if( penergy == 300 ) 
		{
		  hPi_with6_CAL_300->Fill(CALSum);
		  hPi_with6_myBSDLatePE_300->Fill(myBSDLatePE);
	     	  NPi300_with6++;
		}
	      if( penergy == 250 ) 
		{ 
		  hPi_with6_CAL_250->Fill(CALSum);
		  hPi_with6_myBSDLatePE_250->Fill(myBSDLatePE);
	     	  NPi250_with6++;
		}
	    }
	}
      Nsurvive++;
    }
  
  // Report on cuts
  Int_t Nremain = 0;

  cout << "===========================================================" << endl;
  cout << "Ncut_Total = " << Ncut_Total << endl;
  Nremain = Ncut_Total;
  
  cout << "Ncut_IsWithCal = " << Ncut_IsWithCal << ", " << 100*Ncut_IsWithCal/Nremain << endl;
  Nremain -= Ncut_IsWithCal;
  
  cout << "Ncut_neg = " << Ncut_neg << ", " << 100*Ncut_neg/Nremain << endl;
  // Nremain -= Ncut_neg;
  
  cout << "Ncut_BSDLatePE = " << Ncut_BSDLatePE << ", " << 100*Ncut_BSDLatePE/Nremain << endl;
  Nremain -= Ncut_BSDLatePE;
  
  cout << "Ncut_CALSum = " << Ncut_CALSum << ", " << 100*Ncut_CALSum/Nremain << endl;
  Nremain -= Ncut_CALSum;
  
  cout << "Ncut_trig_layers = " << Ncut_trig_layers << ", " << 100*Ncut_trig_layers/Ncut_Total << endl;
  Nremain -= Ncut_trig_layers;
  
  cout << "Nsurvive = " << Nsurvive << ", " << 100*Nsurvive/Ncut_Total << endl;
  cout << "Nremain = " << Nremain << ", " << 100*Nremain/Ncut_Total << endl;
  
  cout << "============================================================" << endl;
  cout << "NE75_no6 = " << NE75_no6 << endl;
  cout << "NE75_with6 = " << NE75_with6 << endl;
  cout << "NE100_no6 = " << NE100_no6 << endl;
  cout << "NE100_with6 = " << NE100_with6 << endl;
  cout << "NE125_no6 = " << NE125_no6 << endl;
  cout << "NE125_with6 = " << NE125_with6 << endl;
  cout << "NE150_no6 = " << NE150_no6 << endl;
  cout << "NE150_with6 = " << NE150_with6 << endl;
  cout << "NE175_no6 = " << NE175_no6 << endl;
  cout << "NE175_with6 = " << NE175_with6 << endl;
  
  cout << "=============================================================" << endl;
  cout << "NPi250_no6 = " << NPi250_no6 << endl;
  cout << "NPi250_with6 = " << NPi250_with6 << endl;
  cout << "NPi300_no6 = " << NPi300_no6 << endl;
  cout << "NPi300_with6 = " << NPi300_with6 << endl;
  cout << "NPi350_no6 = " << NPi350_no6 << endl;
  cout << "NPi350_with6 = " << NPi350_with6 << endl;
  
  
  // Draw
  TCanvas * cv1 = new TCanvas("cv1");
  cv1->Draw();
  cv1->Divide(3,2);
  
  // 50 GeV electrons
  cv1->cd(1)->SetLogy();
  hE_no6_CAL_50->SetLineColor(kBlue);
  hE_with6_CAL_50->SetLineColor(kBlack);
  hE_no6_CAL_50->GetXaxis()->SetTitle("CALSum for 50 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_50->GetXaxis()->CenterTitle();
  hE_no6_CAL_50->Draw();
  hE_with6_CAL_50->Draw("sames");

  // 75 GeV electrons
  cv1->cd(2)->SetLogy();
  hE_no6_CAL_75->SetLineColor(kBlue);
  hE_with6_CAL_75->SetLineColor(kBlack);
  hE_no6_CAL_75->GetXaxis()->SetTitle("CALSum for 75 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_75->GetXaxis()->CenterTitle();
  hE_no6_CAL_75->Draw();
  hE_with6_CAL_75->Draw("sames");

  // 100 GeV electrons
  cv1->cd(3)->SetLogy();
  hE_no6_CAL_100->SetLineColor(kBlue);
  hE_with6_CAL_100->SetLineColor(kBlack);
  hE_no6_CAL_100->GetXaxis()->SetTitle("CALSum for 100 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_100->GetXaxis()->CenterTitle();
  hE_no6_CAL_100->Draw();
  hE_with6_CAL_100->Draw("sames");

  // 125 GeV electrons
  cv1->cd(4)->SetLogy();
  hE_no6_CAL_125->SetLineColor(kBlue);
  hE_with6_CAL_125->SetLineColor(kBlack);
  hE_no6_CAL_125->GetXaxis()->SetTitle("CALSum for 125 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_125->GetXaxis()->CenterTitle();
  hE_no6_CAL_125->Draw();
  hE_with6_CAL_125->Draw("sames");

  // 150 GeV electrons
  cv1->cd(5)->SetLogy();
  hE_no6_CAL_150->SetLineColor(kBlue);
  hE_with6_CAL_150->SetLineColor(kBlack);
  hE_no6_CAL_150->GetXaxis()->SetTitle("CALSum for 150 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_150->GetXaxis()->CenterTitle();
  hE_no6_CAL_150->Draw();
  hE_with6_CAL_150->Draw("sames");

  // 175 GeV electrons
  cv1->cd(6)->SetLogy();
  hE_no6_CAL_175->SetLineColor(kBlue);
  hE_with6_CAL_175->SetLineColor(kBlack);
  hE_no6_CAL_175->GetXaxis()->SetTitle("CALSum for 175 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_175->GetXaxis()->CenterTitle();
  hE_no6_CAL_175->Draw();
  hE_with6_CAL_175->Draw("sames");
  
  TCanvas * cv2 = new TCanvas("cv2");
  cv2->Draw();
  cv2->Divide(3,1);

  // 250 GeV pions
  cv2->cd(1)->SetLogy();
  hPi_no6_CAL_250->SetLineColor(kBlue);
  hPi_with6_CAL_250->SetLineColor(kBlack);
  hPi_no6_CAL_250->GetXaxis()->SetTitle("CALSum for 250 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_CAL_250->GetXaxis()->CenterTitle();
  hPi_no6_CAL_250->Draw();
  hPi_with6_CAL_250->Draw("sames");
  
  // 300 GeV pions
  cv2->cd(2)->SetLogy();
  hPi_no6_CAL_300->SetLineColor(kBlue);
  hPi_with6_CAL_300->SetLineColor(kBlack);
  hPi_no6_CAL_300->GetXaxis()->SetTitle("CALSum for 300 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_CAL_300->GetXaxis()->CenterTitle();
  hPi_no6_CAL_300->Draw();
  hPi_with6_CAL_300->Draw("sames");
  
  // 350 GeV pions
  cv2->cd(3)->SetLogy();
  hPi_no6_CAL_350->SetLineColor(kBlue);
  hPi_with6_CAL_350->SetLineColor(kBlack);
  hPi_no6_CAL_350->GetXaxis()->SetTitle("CALSum for 350 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_CAL_350->GetXaxis()->CenterTitle();
  hPi_no6_CAL_350->Draw();
  hPi_with6_CAL_350->Draw("sames");

  TCanvas * cv3 = new TCanvas("cv3");
  cv3->Draw();
  cv3->Divide(3,2);
  
  // 50 GeV electrons
  cv3->cd(1)->SetLogy();
  hE_no6_myBSDLatePE_50->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_50->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_50->GetXaxis()->SetTitle("myBSDLatePE for 50 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_50->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_50->Draw();
  hE_with6_myBSDLatePE_50->Draw("sames");

  // 75 GeV electrons
  cv3->cd(2)->SetLogy();
  hE_no6_myBSDLatePE_75->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_75->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_75->GetXaxis()->SetTitle("myBSDLatePE for 75 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_75->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_75->Draw();
  hE_with6_myBSDLatePE_75->Draw("sames");

  // 100 GeV electrons
  cv3->cd(3)->SetLogy();
  hE_no6_myBSDLatePE_100->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_100->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_100->GetXaxis()->SetTitle("myBSDLatePE for 100 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_100->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_100->Draw();
  hE_with6_myBSDLatePE_100->Draw("sames");
  
  // 125 GeV electrons
  cv3->cd(4)->SetLogy();
  hE_no6_myBSDLatePE_125->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_125->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_125->GetXaxis()->SetTitle("myBSDLatePE for 125 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_125->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_125->Draw();
  hE_with6_myBSDLatePE_125->Draw("sames");

  // 150 GeV electrons
  cv3->cd(5)->SetLogy();
  hE_no6_myBSDLatePE_150->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_150->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_150->GetXaxis()->SetTitle("myBSDLatePE for 150 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_150->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_150->Draw();
  hE_with6_myBSDLatePE_150->Draw("sames");

  // 175 GeV electrons
  cv3->cd(6)->SetLogy();
  hE_no6_myBSDLatePE_175->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_175->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_175->GetXaxis()->SetTitle("myBSDLatePE for 175 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_175->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_175->Draw();
  hE_with6_myBSDLatePE_175->Draw("sames");
  
  TCanvas * cv4 = new TCanvas("cv4");
  cv4->Draw();
  cv4->Divide(3,1);
  
  // 250 GeV pions
  cv4->cd(1)->SetLogy();
  hPi_no6_myBSDLatePE_250->SetLineColor(kBlue);
  hPi_with6_myBSDLatePE_250->SetLineColor(kBlack);
  hPi_no6_myBSDLatePE_250->GetXaxis()->SetTitle("myBSDLatePE for 250 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_myBSDLatePE_250->GetXaxis()->CenterTitle();
  hPi_no6_myBSDLatePE_250->Draw();
  hPi_with6_myBSDLatePE_250->Draw("sames");
  
  // 300 GeV pions
  cv4->cd(2)->SetLogy();
  hPi_no6_myBSDLatePE_300->SetLineColor(kBlue);
  hPi_with6_myBSDLatePE_300->SetLineColor(kBlack);
  hPi_no6_myBSDLatePE_300->GetXaxis()->SetTitle("myBSDLatePE for 300 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_myBSDLatePE_300->GetXaxis()->CenterTitle();
  hPi_no6_myBSDLatePE_300->Draw();
  hPi_with6_myBSDLatePE_300->Draw("sames");
  
  // 350 GeV pions
  cv4->cd(3)->SetLogy();
  hPi_no6_myBSDLatePE_350->SetLineColor(kBlue);
  hPi_with6_myBSDLatePE_350->SetLineColor(kBlack);
  hPi_no6_myBSDLatePE_350->GetXaxis()->SetTitle("myBSDLatePE for 350 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_myBSDLatePE_350->GetXaxis()->CenterTitle();
  hPi_no6_myBSDLatePE_350->Draw();
  hPi_with6_myBSDLatePE_350->Draw("sames");
    
}


// Make plots for each energy so we can see how it effects each energy (same as Loop2)
// Now fit cut histo to a gaussian and see how many points fall in the 2 sigma and
// 3 sigma range
void six_layer_cut::Loop3()
{
if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  fChain->SetBranchStatus("*",1);
  
  Long64_t nbytes = 0, nb = 0;
  
  Double_t gainA[12] = {4.23,4.23,3.92,3.09,5.81,5.03,3.17,5.30,4.31,3.90,4.04,3.68};
  Double_t myfCHA[12];
  Double_t myBSDLatePE = 0;
  for( int i = 0; i < 12; i++ )
    myfCHA[i] = -1;
  
  
  // 50 GeV electrons
  TH1F * hE_no6_CAL_50 = new TH1F("hE_no6_CAL_50","",150,-1000,20000);
  TH1F * hE_with6_CAL_50 = new TH1F("hE_with6_CAL_50","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_50 = new TH1F("hE_no6_myBSDLatePE_50","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_50 = new TH1F("hE_with6_myBSDLatePE_50","",150,-1000,8000);
  
  // 75 GeV electrons
  TH1F * hE_no6_CAL_75 = new TH1F("hE_no6_CAL_75","",150,-1000,20000);
  TH1F * hE_with6_CAL_75 = new TH1F("hE_with6_CAL_75","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_75 = new TH1F("hE_no6_myBSDLatePE_75","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_75 = new TH1F("hE_with6_myBSDLatePE_75","",150,-1000,8000);
  
  // 100 GeV electrons
  TH1F * hE_no6_CAL_100 = new TH1F("hE_no6_CAL_100","",150,-1000,20000);
  TH1F * hE_with6_CAL_100 = new TH1F("hE_with6_CAL_100","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_100 = new TH1F("hE_no6_myBSDLatePE_100","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_100 = new TH1F("hE_with6_myBSDLatePE_100","",150,-1000,8000);
  
  // 125 GeV electrons
  TH1F * hE_no6_CAL_125 = new TH1F("hE_no6_CAL_125","",150,-1000,20000);
  TH1F * hE_with6_CAL_125 = new TH1F("hE_with6_CAL_125","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_125 = new TH1F("hE_no6_myBSDLatePE_125","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_125 = new TH1F("hE_with6_myBSDLatePE_125","",150,-1000,8000);
  
  // 150 GeV electrons
  TH1F * hE_no6_CAL_150 = new TH1F("hE_no6_CAL_150","",150,-1000,20000);
  TH1F * hE_with6_CAL_150 = new TH1F("hE_with6_CAL_150","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_150 = new TH1F("hE_no6_myBSDLatePE_150","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_150 = new TH1F("hE_with6_myBSDLatePE_150","",150,-1000,8000);
  
  // 175 GeV electrons
  TH1F * hE_no6_CAL_175 = new TH1F("hE_no6_CAL_175","",150,-1000,20000);
  TH1F * hE_with6_CAL_175 = new TH1F("hE_with6_CAL_175","",150,-1000,20000);
  
  TH1F * hE_no6_myBSDLatePE_175 = new TH1F("hE_no6_myBSDLatePE_175","",150,-1000,8000);
  TH1F * hE_with6_myBSDLatePE_175 = new TH1F("hE_with6_myBSDLatePE_175","",150,-1000,8000);

  // 250 GeV pions
  TH1F * hPi_no6_CAL_250 = new TH1F("hPi_no6_CAL_250","",150,-1000,20000);
  TH1F * hPi_with6_CAL_250 = new TH1F("hPi_with6_CAL_250","",150,-1000,20000);

  TH1F * hPi_no6_myBSDLatePE_250 = new TH1F("hPi_no6_myBSDLatePE_250","",150,-1000,8000);
  TH1F * hPi_with6_myBSDLatePE_250 = new TH1F("hPi_with6_myBSDLatePE_250","",150,-1000,8000);

  // 300 GeV pions
  TH1F * hPi_no6_CAL_300 = new TH1F("hPi_no6_CAL_300","",150,-1000,20000);
  TH1F * hPi_with6_CAL_300 = new TH1F("hPi_with6_CAL_300","",150,-1000,20000);

  TH1F * hPi_no6_myBSDLatePE_300 = new TH1F("hPi_no6_myBSDLatePE_300","",150,-1000,8000);
  TH1F * hPi_with6_myBSDLatePE_300 = new TH1F("hPi_with6_myBSDLatePE_300","",150,-1000,8000);

  // 350 GeV pions
  TH1F * hPi_no6_CAL_350 = new TH1F("hPi_no6_CAL_350","",150,-1000,20000);
  TH1F * hPi_with6_CAL_350 = new TH1F("hPi_with6_CAL_350","",150,-1000,20000);

  TH1F * hPi_no6_myBSDLatePE_350 = new TH1F("hPi_no6_myBSDLatePE_350","",150,-1000,8000);
  TH1F * hPi_with6_myBSDLatePE_350 = new TH1F("hPi_with6_myBSDLatePE_350","",150,-1000,8000);

  
  Int_t NE50_no6 = 0;
  Int_t NE75_no6 = 0;
  Int_t NE100_no6 = 0;
  Int_t NE125_no6 = 0;
  Int_t NE150_no6 = 0;
  Int_t NE175_no6 = 0;
  
  Int_t NPi250_no6 = 0;
  Int_t NPi300_no6 = 0;
  Int_t NPi350_no6 = 0;    
  
  Int_t NE50_with6 = 0;
  Int_t NE75_with6 = 0;
  Int_t NE100_with6 = 0;
  Int_t NE125_with6 = 0;
  Int_t NE150_with6 = 0;
  Int_t NE175_with6 = 0;

  Int_t NPi250_with6 = 0;
  Int_t NPi300_with6 = 0;
  Int_t NPi350_with6 = 0;    
  
  Double_t Ncut_Total = 0;
  Double_t Ncut_BSDLatePE = 0;
  Double_t Ncut_IsWithCal = 0;
  Double_t Ncut_CALSum = 0;
  Double_t Ncut_trig_layers = 0;
  Double_t Ncut_neg = 0;
  Double_t Nsurvive = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
      Ncut_Total++;
      
      if( !IsWithCal ) 
	{
	  Ncut_IsWithCal++;
	  continue;
	}
      
      myBSDLatePE = 0;
      for( int i = 1; i <= 10; i++ )
	myfCHA[i] = -1;
      
      Bool_t neg = 0;
      for( int i = 1; i <= 10; i++ )
	{
	  if( fCHA[i] <= 0 )
	    {
	      neg = kTRUE;
	      myfCHA[i] = 2048;
	    }
	  else
	    myfCHA[i] = fCHA[i];
	}
      if( neg == kTRUE ) 
	{
	  Ncut_neg++;
	  // continue;
	}
      
      if( myBSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  // continue;
	}
      
      if( CALSum < 0 ) 
	{
	  Ncut_CALSum++;
	  // continue;
	}
      
      myBSDLatePE = 0;
      // Calculate myBSDLatePE
      for( int i = 1; i <= 10; i++ )
	myBSDLatePE += myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19); 
      
      ///////////////////////////////////////////////////////////////////////////////
      // Now we start the weighting game!
      // For Loop1, we want to see how the 6 layer cut down selects from a set of
      // otherwise viable events for pions and electrons
      // Since the 6 layer cuts much more strongly on lower energies
      //////////////////////////////////////////////////////////////////////////////
      
      // Want to weight all energies equally before the 6 layer cut
      // TBA Thu Jan  3 20:36:00 EST 2013
      // Don't actually need to do this for the analysis I'm doing now
      /*
      if( (ptype == 1) && (penergy == 50) && (NE50_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 75) && (NE75_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 100) && (NE100_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 125) && (NE125_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 150) && (NE150_no6 > 7000) ) continue;
      if( (ptype == 1) && (penergy == 175) && (NE175_no6 > 7000) ) continue;
      
      if( (ptype == 2) && (penergy == 250) && (NPi250_no6 > 18000) ) continue; 
      if( (ptype == 2) && (penergy == 300) && (NPi300_no6 > 18000) ) continue; 
      if( (ptype == 2) && (penergy == 350) && (NPi350_no6 > 18000) ) continue; 
      */
 
      // Fill no six layer histo for electrons
      if( ptype == 1 )
	{
	  if( penergy == 175 ) 
	    {
	      hE_no6_CAL_175->Fill(CALSum);
	      hE_no6_myBSDLatePE_175->Fill(myBSDLatePE);
	      NE175_no6++;
	    }
	  if( penergy == 150 ) 
	    {
	      hE_no6_CAL_150->Fill(CALSum);
	      hE_no6_myBSDLatePE_150->Fill(myBSDLatePE);
	      NE150_no6++;
	    }
	  if( penergy == 125 ) 
	    {
	      hE_no6_CAL_125->Fill(CALSum);
	      hE_no6_myBSDLatePE_125->Fill(myBSDLatePE);
	      NE125_no6++;
	    }
	  if( penergy == 100 ) 
	    {
	      hE_no6_CAL_100->Fill(CALSum);
	      hE_no6_myBSDLatePE_100->Fill(myBSDLatePE);
	      NE100_no6++;
	    }
	  if( penergy == 75 ) 
	    {
	      hE_no6_CAL_75->Fill(CALSum);
	      hE_no6_myBSDLatePE_75->Fill(myBSDLatePE);
	      NE75_no6++;
	    }
	  if( penergy == 50 ) 
	    {
	      hE_no6_CAL_50->Fill(CALSum);
	      hE_no6_myBSDLatePE_50->Fill(myBSDLatePE);
	      NE50_no6++;
	    }
	}
      
      // Fill no six layer histo for pions
      if( ptype == 2 )
	{
	  if( penergy == 350 ) 
	    {
	      hPi_no6_CAL_350->Fill(CALSum);
	      hPi_no6_myBSDLatePE_350->Fill(myBSDLatePE);
	      NPi350_no6++;
	    }
	  if( penergy == 300 ) 
	    {
	      hPi_no6_CAL_300->Fill(CALSum);
	      hPi_no6_myBSDLatePE_300->Fill(myBSDLatePE);
	      NPi300_no6++;
	    }
	  if( penergy == 250 ) 
	    {
	      hPi_no6_CAL_250->Fill(CALSum);
	      hPi_no6_myBSDLatePE_250->Fill(myBSDLatePE);
	      NPi250_no6++;
	    }
	}
      
      /////////////////////////////////////////////////////////////////////////////
      // Now apply cut to downselected events
      // Simulate CAL trigger by requiring 6 consecutive layers with > 40 MeV
      /////////////////////////////////////////////////////////////////////////////
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
      if( Ntrig_layers < 6 )
	{
	  Ncut_trig_layers++;
	  continue;
	}
      else
	{
	  // Fill no six layer histo for electrons
	  if( ptype == 1 )
	    {
	      if( penergy == 175 ) 
		{
		  hE_with6_CAL_175->Fill(CALSum);
		  hE_with6_myBSDLatePE_175->Fill(myBSDLatePE);
		  NE175_with6++;
		}
	      if( penergy == 150 ) 
		{
		  hE_with6_CAL_150->Fill(CALSum);
		  hE_with6_myBSDLatePE_150->Fill(myBSDLatePE);
		  NE150_with6++;
		}
	      if( penergy == 125 ) 
		{ hE_with6_CAL_125->Fill(CALSum);
		  hE_with6_myBSDLatePE_125->Fill(myBSDLatePE);
		  NE125_with6++;
		}
	      if( penergy == 100 ) 
		{ hE_with6_CAL_100->Fill(CALSum);
		  hE_with6_myBSDLatePE_100->Fill(myBSDLatePE);
		  NE100_with6++;
		}
	      if( penergy == 75 ) 
		{ hE_with6_CAL_75->Fill(CALSum);
		  hE_with6_myBSDLatePE_75->Fill(myBSDLatePE);
		  NE75_with6++;
		}
	      if( penergy == 50 ) 
		{ hE_with6_CAL_50->Fill(CALSum);
		  hE_with6_myBSDLatePE_50->Fill(myBSDLatePE);
		  NE75_with6++;
		}
	    }
	  
	  // Fill no six layer histo for pions
	  if( ptype == 2 )
	    {
	      if( penergy == 350 ) 
		{
		  hPi_with6_CAL_350->Fill(CALSum);
		  hPi_with6_myBSDLatePE_350->Fill(myBSDLatePE);
		  NPi350_with6++;
		}
	      if( penergy == 300 ) 
		{
		  hPi_with6_CAL_300->Fill(CALSum);
		  hPi_with6_myBSDLatePE_300->Fill(myBSDLatePE);
	     	  NPi300_with6++;
		}
	      if( penergy == 250 ) 
		{ 
		  hPi_with6_CAL_250->Fill(CALSum);
		  hPi_with6_myBSDLatePE_250->Fill(myBSDLatePE);
	     	  NPi250_with6++;
		}
	    }
	}
      Nsurvive++;
    }
  
  // Report on cuts
  Int_t Nremain = 0;

  cout << "===========================================================" << endl;
  cout << "Ncut_Total = " << Ncut_Total << endl;
  Nremain = Ncut_Total;
  
  cout << "Ncut_IsWithCal = " << Ncut_IsWithCal << ", " << 100*Ncut_IsWithCal/Nremain << endl;
  Nremain -= Ncut_IsWithCal;
  
  cout << "Ncut_neg = " << Ncut_neg << ", " << 100*Ncut_neg/Nremain << endl;
  // Nremain -= Ncut_neg;
  
  cout << "Ncut_BSDLatePE = " << Ncut_BSDLatePE << ", " << 100*Ncut_BSDLatePE/Nremain << endl;
  Nremain -= Ncut_BSDLatePE;
  
  cout << "Ncut_CALSum = " << Ncut_CALSum << ", " << 100*Ncut_CALSum/Nremain << endl;
  Nremain -= Ncut_CALSum;
  
  cout << "Ncut_trig_layers = " << Ncut_trig_layers << ", " << 100*Ncut_trig_layers/Ncut_Total << endl;
  Nremain -= Ncut_trig_layers;
  
  cout << "Nsurvive = " << Nsurvive << ", " << 100*Nsurvive/Ncut_Total << endl;
  cout << "Nremain = " << Nremain << ", " << 100*Nremain/Ncut_Total << endl;
  
  cout << "============================================================" << endl;
  cout << "NE75_no6 = " << NE75_no6 << endl;
  cout << "NE75_with6 = " << NE75_with6 << endl;
  cout << "NE100_no6 = " << NE100_no6 << endl;
  cout << "NE100_with6 = " << NE100_with6 << endl;
  cout << "NE125_no6 = " << NE125_no6 << endl;
  cout << "NE125_with6 = " << NE125_with6 << endl;
  cout << "NE150_no6 = " << NE150_no6 << endl;
  cout << "NE150_with6 = " << NE150_with6 << endl;
  cout << "NE175_no6 = " << NE175_no6 << endl;
  cout << "NE175_with6 = " << NE175_with6 << endl;
  
  cout << "=============================================================" << endl;
  cout << "NPi250_no6 = " << NPi250_no6 << endl;
  cout << "NPi250_with6 = " << NPi250_with6 << endl;
  cout << "NPi300_no6 = " << NPi300_no6 << endl;
  cout << "NPi300_with6 = " << NPi300_with6 << endl;
  cout << "NPi350_no6 = " << NPi350_no6 << endl;
  cout << "NPi350_with6 = " << NPi350_with6 << endl;
  
  
  // Draw
  TCanvas * cv1 = new TCanvas("cv1");
  cv1->Draw();
  cv1->Divide(3,2);
  
  Double_t means[5];
  Double_t sigmas[5];

  TF1 * f1 = new TF1("f1","gaus"); // dummy pointer for the first pass fit
  
  // 50 GeV electrons
  cv1->cd(1)->SetLogy();
  hE_no6_CAL_50->SetLineColor(kBlue);
  hE_with6_CAL_50->SetLineColor(kBlack);
  hE_no6_CAL_50->GetXaxis()->SetTitle("CALSum for 50 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_50->GetXaxis()->CenterTitle();
  hE_no6_CAL_50->Draw();
  hE_with6_CAL_50->Draw("sames");
  
  ////////////////////////////////////////////////////////////////////////////////////
  // 75 GeV electrons
  cv1->cd(2)->SetLogy();
  hE_no6_CAL_75->SetLineColor(kBlue);
  hE_with6_CAL_75->SetLineColor(kBlack);
  hE_no6_CAL_75->GetXaxis()->SetTitle("CALSum for 75 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_75->GetXaxis()->CenterTitle();
  hE_no6_CAL_75->Draw();
  hE_with6_CAL_75->Draw("sames");
  hE_with6_CAL_75->Fit("f1","","same");
  means[0] = f1->GetParameter(1);
  sigmas[0] = f1->GetParameter(2);
  f1->SetRange(means[0]-2*sigmas[0],means[0]+2*sigmas[0]);
  hE_with6_CAL_75->Fit("f1","R","same");
  means[0] = f1->GetParameter(1);
  sigmas[0] = f1->GetParameter(2);
  
  ////////////////////////////////////////////////////////////////////////////////////
  // 100 GeV electrons
  cv1->cd(3)->SetLogy();
  hE_no6_CAL_100->SetLineColor(kBlue);
  hE_with6_CAL_100->SetLineColor(kBlack);
  hE_no6_CAL_100->GetXaxis()->SetTitle("CALSum for 100 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_100->GetXaxis()->CenterTitle();
  hE_no6_CAL_100->Draw();
  hE_with6_CAL_100->Draw("sames");
  hE_with6_CAL_100->Fit("f1","","same");
  means[1] = f1->GetParameter(1);
  sigmas[1] = f1->GetParameter(2);
  f1->SetRange(means[1]-2*sigmas[1],means[1]+2*sigmas[1]);
  hE_with6_CAL_100->Fit("f1","R","same");
  means[1] = f1->GetParameter(1);
  sigmas[1] = f1->GetParameter(2);
  
  ////////////////////////////////////////////////////////////////////////////////////
  // 125 GeV electrons
  cv1->cd(4)->SetLogy();
  hE_no6_CAL_125->SetLineColor(kBlue);
  hE_with6_CAL_125->SetLineColor(kBlack);
  hE_no6_CAL_125->GetXaxis()->SetTitle("CALSum for 125 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_125->GetXaxis()->CenterTitle();
  hE_no6_CAL_125->Draw();
  hE_with6_CAL_125->Draw("sames");
  hE_with6_CAL_125->Fit("f1","","same");
  means[2] = f1->GetParameter(1);
  sigmas[2] = f1->GetParameter(2);
  f1->SetRange(means[2]-2*sigmas[2],means[2]+2*sigmas[2]);
  hE_with6_CAL_125->Fit("f1","R","same");
  means[2] = f1->GetParameter(1);
  sigmas[2] = f1->GetParameter(2);

  ////////////////////////////////////////////////////////////////////////////////////
  // 150 GeV electrons
  cv1->cd(5)->SetLogy();
  hE_no6_CAL_150->SetLineColor(kBlue);
  hE_with6_CAL_150->SetLineColor(kBlack);
  hE_no6_CAL_150->GetXaxis()->SetTitle("CALSum for 150 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_150->GetXaxis()->CenterTitle();
  hE_no6_CAL_150->Draw();
  hE_with6_CAL_150->Draw("sames");
  hE_with6_CAL_150->Fit("f1","","same");
  means[3] = f1->GetParameter(1);
  sigmas[3] = f1->GetParameter(2);
  f1->SetRange(means[3]-2*sigmas[3],means[3]+2*sigmas[3]);
  hE_with6_CAL_150->Fit("f1","R","same");
  means[3] = f1->GetParameter(1);
  sigmas[3] = f1->GetParameter(2);

  ////////////////////////////////////////////////////////////////////////////////////
  // 175 GeV electrons
  cv1->cd(6)->SetLogy();
  hE_no6_CAL_175->SetLineColor(kBlue);
  hE_with6_CAL_175->SetLineColor(kBlack);
  hE_no6_CAL_175->GetXaxis()->SetTitle("CALSum for 175 GeV Electrons with (black) and without (blue) six layer cut");
  hE_no6_CAL_175->GetXaxis()->CenterTitle();
  hE_no6_CAL_175->Draw();
  hE_with6_CAL_175->Draw("sames");
  hE_with6_CAL_175->Fit("f1","","same");
  means[4] = f1->GetParameter(1);
  sigmas[4] = f1->GetParameter(2);
  f1->SetRange(means[4]-2*sigmas[4],means[4]+2*sigmas[4]);
  hE_with6_CAL_175->Fit("f1","R","same");
  means[4] = f1->GetParameter(1);
  sigmas[4] = f1->GetParameter(2);

  cout << "means[5] = {" 
       << means[0] << ", "
       << means[1] << ", "
       << means[2] << ", "
       << means[3] << ", "
       << means[4] << "};"
       << endl;
  
  cout << "sigmas[5] = {" 
       << sigmas[0] << ", "
       << sigmas[1] << ", "
       << sigmas[2] << ", "
       << sigmas[3] << ", "
       << sigmas[4] << "};"
       << endl;
  
  return;

  TCanvas * cv2 = new TCanvas("cv2");
  cv2->Draw();
  cv2->Divide(3,1);

  // 250 GeV pions
  cv2->cd(1)->SetLogy();
  hPi_no6_CAL_250->SetLineColor(kBlue);
  hPi_with6_CAL_250->SetLineColor(kBlack);
  hPi_no6_CAL_250->GetXaxis()->SetTitle("CALSum for 250 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_CAL_250->GetXaxis()->CenterTitle();
  hPi_no6_CAL_250->Draw();
  hPi_with6_CAL_250->Draw("sames");
    
  // 300 GeV pions
  cv2->cd(2)->SetLogy();
  hPi_no6_CAL_300->SetLineColor(kBlue);
  hPi_with6_CAL_300->SetLineColor(kBlack);
  hPi_no6_CAL_300->GetXaxis()->SetTitle("CALSum for 300 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_CAL_300->GetXaxis()->CenterTitle();
  hPi_no6_CAL_300->Draw();
  hPi_with6_CAL_300->Draw("sames");
    
  // 350 GeV pions
  cv2->cd(3)->SetLogy();
  hPi_no6_CAL_350->SetLineColor(kBlue);
  hPi_with6_CAL_350->SetLineColor(kBlack);
  hPi_no6_CAL_350->GetXaxis()->SetTitle("CALSum for 350 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_CAL_350->GetXaxis()->CenterTitle();
  hPi_no6_CAL_350->Draw();
  hPi_with6_CAL_350->Draw("sames");
  
  TCanvas * cv3 = new TCanvas("cv3");
  cv3->Draw();
  cv3->Divide(3,2);
  
  // 50 GeV electrons
  cv3->cd(1)->SetLogy();
  hE_no6_myBSDLatePE_50->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_50->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_50->GetXaxis()->SetTitle("myBSDLatePE for 50 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_50->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_50->Draw();
  hE_with6_myBSDLatePE_50->Draw("sames");

  // 75 GeV electrons
  cv3->cd(2)->SetLogy();
  hE_no6_myBSDLatePE_75->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_75->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_75->GetXaxis()->SetTitle("myBSDLatePE for 75 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_75->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_75->Draw();
  hE_with6_myBSDLatePE_75->Draw("sames");

  // 100 GeV electrons
  cv3->cd(3)->SetLogy();
  hE_no6_myBSDLatePE_100->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_100->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_100->GetXaxis()->SetTitle("myBSDLatePE for 100 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_100->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_100->Draw();
  hE_with6_myBSDLatePE_100->Draw("sames");
  
  // 125 GeV electrons
  cv3->cd(4)->SetLogy();
  hE_no6_myBSDLatePE_125->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_125->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_125->GetXaxis()->SetTitle("myBSDLatePE for 125 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_125->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_125->Draw();
  hE_with6_myBSDLatePE_125->Draw("sames");

  // 150 GeV electrons
  cv3->cd(5)->SetLogy();
  hE_no6_myBSDLatePE_150->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_150->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_150->GetXaxis()->SetTitle("myBSDLatePE for 150 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_150->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_150->Draw();
  hE_with6_myBSDLatePE_150->Draw("sames");

  // 175 GeV electrons
  cv3->cd(6)->SetLogy();
  hE_no6_myBSDLatePE_175->SetLineColor(kBlue);
  hE_with6_myBSDLatePE_175->SetLineColor(kBlack);
  hE_no6_myBSDLatePE_175->GetXaxis()->SetTitle("myBSDLatePE for 175 GeV Electrons with (black) and without (blue) six layer cut (detector units)");
  hE_no6_myBSDLatePE_175->GetXaxis()->CenterTitle();
  hE_no6_myBSDLatePE_175->Draw();
  hE_with6_myBSDLatePE_175->Draw("sames");
  
  TCanvas * cv4 = new TCanvas("cv4");
  cv4->Draw();
  cv4->Divide(3,1);
  
  // 250 GeV pions
  cv4->cd(1)->SetLogy();
  hPi_no6_myBSDLatePE_250->SetLineColor(kBlue);
  hPi_with6_myBSDLatePE_250->SetLineColor(kBlack);
  hPi_no6_myBSDLatePE_250->GetXaxis()->SetTitle("myBSDLatePE for 250 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_myBSDLatePE_250->GetXaxis()->CenterTitle();
  hPi_no6_myBSDLatePE_250->Draw();
  hPi_with6_myBSDLatePE_250->Draw("sames");
  
  // 300 GeV pions
  cv4->cd(2)->SetLogy();
  hPi_no6_myBSDLatePE_300->SetLineColor(kBlue);
  hPi_with6_myBSDLatePE_300->SetLineColor(kBlack);
  hPi_no6_myBSDLatePE_300->GetXaxis()->SetTitle("myBSDLatePE for 300 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_myBSDLatePE_300->GetXaxis()->CenterTitle();
  hPi_no6_myBSDLatePE_300->Draw();
  hPi_with6_myBSDLatePE_300->Draw("sames");
  
  // 350 GeV pions
  cv4->cd(3)->SetLogy();
  hPi_no6_myBSDLatePE_350->SetLineColor(kBlue);
  hPi_with6_myBSDLatePE_350->SetLineColor(kBlack);
  hPi_no6_myBSDLatePE_350->GetXaxis()->SetTitle("myBSDLatePE for 350 GeV Pions with (black) and without (blue) six layer cut (detector units)");
  hPi_no6_myBSDLatePE_350->GetXaxis()->CenterTitle();
  hPi_no6_myBSDLatePE_350->Draw();
  hPi_with6_myBSDLatePE_350->Draw("sames");
    
}


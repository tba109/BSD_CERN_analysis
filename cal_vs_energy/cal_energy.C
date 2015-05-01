#define cal_energy_cxx
#include "cal_energy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void cal_energy::Loop1()
{
  if (fChain == 0) return;

  TH1F * hPi350 = new TH1F("hPi350","",200,-1000,30000);
  TH1F * hPi300 = new TH1F("hPi300","",200,-1000,30000);
  TH1F * hPi250 = new TH1F("hPi250","",200,-1000,30000);

  TH1F * hE175 = new TH1F("hE175","",200,-1000,30000);
  TH1F * hE150 = new TH1F("hE150","",200,-1000,30000);
  TH1F * hE125 = new TH1F("hE125","",200,-1000,30000);
  TH1F * hE100 = new TH1F("hE100","",200,-1000,30000);
  TH1F * hE75 = new TH1F("hE75","",200,-1000,30000);
  TH1F * hE50 = new TH1F("hE50","",200,-1000,30000);
  
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

      // 350 GeV pions
      if( ptype == 2 && penergy == 350) hPi350->Fill(CALSum);
      
      // 300 GeV pions
      if( ptype == 2 && penergy == 300) hPi300->Fill(CALSum);
      
      // 250 GeV pions
      if( ptype == 2 && penergy == 250) hPi250->Fill(CALSum);
      
      // 175 GeV electrons
      if( ptype == 1 && penergy == 175) hE175->Fill(CALSum);

      // 150 GeV electrons
      if( ptype == 1 && penergy == 150) hE150->Fill(CALSum);
      
      // 125 GeV electrons
      if( ptype == 1 && penergy == 125) hE125->Fill(CALSum);

      // 100 GeV electrons
      if( ptype == 1 && penergy == 100) hE100->Fill(CALSum);
      
      // 75 GeV electrons
      if( ptype == 1 && penergy == 75) hE75->Fill(CALSum);
      
      // 50 GeV electrons
      if( ptype == 1 && penergy == 50) hE50->Fill(CALSum);
    }
  
  TCanvas * cvPi = new TCanvas("cvPi");
  hPi350->SetLineColor(kRed);
  hPi300->SetLineColor(kGreen);
  hPi250->SetLineColor(kBlue);

  if( hPi350->GetEntries() > 0 )
    {
      if( hPi300->GetEntries() > 0 )
	hPi300->Scale(hPi350->GetEntries()/hPi300->GetEntries());
      
      if( hPi250->GetEntries() > 0 )
	hPi250->Scale(hPi350->GetEntries()/hPi250->GetEntries());
    }
  
  TH2F * haxesPi = new TH2F("haxesPi","",100,-1000,30000,100,0,1000);
  haxesPi->SetStats(kFALSE);
  haxesPi->Draw();
  
  hPi350->Draw("sames");
  hPi300->Draw("sames");
  hPi250->Draw("sames");

  TCanvas * cvE = new TCanvas("cvE");
  hE175->SetLineColor(kRed);
  hE150->SetLineColor(kOrange);
  hE125->SetLineColor(kYellow);
  hE100->SetLineColor(kGreen);
  hE75->SetLineColor(kBlue);
  hE50->SetLineColor(kViolet);

  if( hE175->GetEntries() > 0 )
    {
      if( hE150->GetEntries() > 0 )
	hE150->Scale(hE175->GetEntries()/hE150->GetEntries());
      
      if( hE125->GetEntries() > 0 )
	hE125->Scale(hE175->GetEntries()/hE125->GetEntries());
      
      if( hE100->GetEntries() > 0 )
	hE100->Scale(hE175->GetEntries()/hE100->GetEntries());
      
      if( hE75->GetEntries() > 0 )
	hE75->Scale(hE175->GetEntries()/hE75->GetEntries());
      
      if( hE50->GetEntries() > 0 )
	hE50->Scale(hE175->GetEntries()/hE50->GetEntries());
    }

  TH2F * haxesE = new TH2F("haxesE","",100,-1000,30000,100,0,800);
  haxesE->SetStats(kFALSE);
  haxesE->Draw();

  hE175->Draw("sames");
  hE150->Draw("sames");
  hE125->Draw("sames");
  hE100->Draw("sames");
  hE75->Draw("sames");
  // hE50->Draw("sames");
  
}

void cal_energy::Loop2()
{
  if (fChain == 0) return;

  TH1F * hPi350 = new TH1F("hPi350","",200,-1000,30000);
  TH1F * hPi300 = new TH1F("hPi300","",200,-1000,30000);
  TH1F * hPi250 = new TH1F("hPi250","",200,-1000,30000);

  TH1F * hE175 = new TH1F("hE175","",200,-1000,30000);
  TH1F * hE150 = new TH1F("hE150","",200,-1000,30000);
  TH1F * hE125 = new TH1F("hE125","",200,-1000,30000);
  TH1F * hE100 = new TH1F("hE100","",200,-1000,30000);
  TH1F * hE75 = new TH1F("hE75","",200,-1000,30000);
  TH1F * hE50 = new TH1F("hE50","",200,-1000,30000);
  
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

      // 350 GeV pions
      if( ptype == 2 && penergy == 350) hPi350->Fill(CALSum);
      
      // 300 GeV pions
      if( ptype == 2 && penergy == 300) hPi300->Fill(CALSum);
      
      // 250 GeV pions
      if( ptype == 2 && penergy == 250) hPi250->Fill(CALSum);
      
      // 175 GeV electrons
      if( ptype == 1 && penergy == 175) hE175->Fill(CALSum);

      // 150 GeV electrons
      if( ptype == 1 && penergy == 150) hE150->Fill(CALSum);
      
      // 125 GeV electrons
      if( ptype == 1 && penergy == 125) hE125->Fill(CALSum);

      // 100 GeV electrons
      if( ptype == 1 && penergy == 100) hE100->Fill(CALSum);
      
      // 75 GeV electrons
      if( ptype == 1 && penergy == 75) hE75->Fill(CALSum);
      
      // 50 GeV electrons
      if( ptype == 1 && penergy == 50) hE50->Fill(CALSum);
    }
  
  TCanvas * cvPi = new TCanvas("cvPi");
  hPi350->SetLineColor(kRed);
  hPi300->SetLineColor(kGreen);
  hPi250->SetLineColor(kBlue);

  if( hPi350->GetEntries() > 0 )
    hPi350->Scale(1./hPi350->GetEntries());
  
  if( hPi300->GetEntries() > 0 )
    hPi300->Scale(1./hPi300->GetEntries());
  
  if( hPi250->GetEntries() > 0 )
    hPi250->Scale(1./hPi250->GetEntries());
  
  
  TH2F * haxesPi = new TH2F("haxesPi","",100,-1000,30000,100,0,0.03);
  haxesPi->SetStats(kFALSE);
  haxesPi->Draw();
  
  hPi350->Draw("sames");
  hPi300->Draw("sames");
  hPi250->Draw("sames");

  TCanvas * cvE = new TCanvas("cvE");
  hE175->SetLineColor(kRed);
  hE150->SetLineColor(kOrange);
  hE125->SetLineColor(kYellow);
  hE100->SetLineColor(kGreen);
  hE75->SetLineColor(kBlue);
  hE50->SetLineColor(kViolet);

  if( hE175->GetEntries() > 0 )
    hE175->Scale(1./hE175->GetEntries());
  
  if( hE150->GetEntries() > 0 )
    hE150->Scale(1./hE150->GetEntries());
  
  if( hE125->GetEntries() > 0 )
    hE125->Scale(1./hE125->GetEntries());
  
  if( hE100->GetEntries() > 0 )
    hE100->Scale(1./hE100->GetEntries());
  
  if( hE75->GetEntries() > 0 )
    hE75->Scale(1./hE75->GetEntries());
  
  if( hE50->GetEntries() > 0 )
    hE50->Scale(1./hE50->GetEntries());
  
  TH2F * haxesE = new TH2F("haxesE","",100,-1000,30000,100,0,0.16);
  haxesE->SetStats(kFALSE);
  haxesE->Draw();

  hE175->Draw("sames");
  hE150->Draw("sames");
  hE125->Draw("sames");
  hE100->Draw("sames");
  hE75->Draw("sames");
  // hE50->Draw("sames");
  
}

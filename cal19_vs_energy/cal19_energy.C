#define cal19_energy_cxx
#include "cal19_energy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void cal19_energy::Loop1()
{
  if (fChain == 0) return;

  TH1F * hPi350 = new TH1F("hPi350","",2000,-1000,30000);
  TH1F * hPi300 = new TH1F("hPi300","",2000,-1000,30000);
  TH1F * hPi250 = new TH1F("hPi250","",2000,-1000,30000);

  TH1F * hE175 = new TH1F("hE175","",2000,-1000,30000);
  TH1F * hE150 = new TH1F("hE150","",2000,-1000,30000);
  TH1F * hE125 = new TH1F("hE125","",2000,-1000,30000);
  TH1F * hE100 = new TH1F("hE100","",2000,-1000,30000);
  TH1F * hE75 = new TH1F("hE75","",2000,-1000,30000);
  TH1F * hE50 = new TH1F("hE50","",2000,-1000,30000);
  
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
      if( ptype == 2 && penergy == 350) hPi350->Fill(CALProfile[19]);
      
      // 300 GeV pions
      if( ptype == 2 && penergy == 300) hPi300->Fill(CALProfile[19]);
      
      // 250 GeV pions
      if( ptype == 2 && penergy == 250) hPi250->Fill(CALProfile[19]);
      
      // 175 GeV electrons
      if( ptype == 1 && penergy == 175) hE175->Fill(CALProfile[19]);

      // 150 GeV electrons
      if( ptype == 1 && penergy == 150) hE150->Fill(CALProfile[19]);
      
      // 125 GeV electrons
      if( ptype == 1 && penergy == 125) hE125->Fill(CALProfile[19]);

      // 100 GeV electrons
      if( ptype == 1 && penergy == 100) hE100->Fill(CALProfile[19]);
      
      // 75 GeV electrons
      if( ptype == 1 && penergy == 75) hE75->Fill(CALProfile[19]);
      
      // 50 GeV electrons
      if( ptype == 1 && penergy == 50) hE50->Fill(CALProfile[19]);
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

// This one looks over the entire last layer of the calorimeter
void cal19_energy::Loop2()
{
  if (fChain == 0) return;

  TH1F * hPi350 = new TH1F("hPi350","",2000,-1000,30000);
  TH1F * hPi300 = new TH1F("hPi300","",2000,-1000,30000);
  TH1F * hPi250 = new TH1F("hPi250","",2000,-1000,30000);

  TH1F * hE175 = new TH1F("hE175","",2000,-1000,30000);
  TH1F * hE150 = new TH1F("hE150","",2000,-1000,30000);
  TH1F * hE125 = new TH1F("hE125","",2000,-1000,30000);
  TH1F * hE100 = new TH1F("hE100","",2000,-1000,30000);
  TH1F * hE75 = new TH1F("hE75","",2000,-1000,30000);
  TH1F * hE50 = new TH1F("hE50","",2000,-1000,30000);
  
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


      Double_t sum = 0;
      for( int k = 0; k < 50; k++ )
	{
	  if( cal[1][9][k] > 0 )
	    sum += cal[1][9][k];
	}

      printf("sum = %f\n",sum);

      // 350 GeV pions
      if( ptype == 2 && penergy == 350) hPi350->Fill(sum);
      
      // 300 GeV pions
      if( ptype == 2 && penergy == 300) hPi300->Fill(sum);
      
      // 250 GeV pions
      if( ptype == 2 && penergy == 250) hPi250->Fill(sum);
      
      // 175 GeV electrons
      if( ptype == 1 && penergy == 175) hE175->Fill(sum);

      // 150 GeV electrons
      if( ptype == 1 && penergy == 150) hE150->Fill(sum);
      
      // 125 GeV electrons
      if( ptype == 1 && penergy == 125) hE125->Fill(sum);

      // 100 GeV electrons
      if( ptype == 1 && penergy == 100) hE100->Fill(sum);
      
      // 75 GeV electrons
      if( ptype == 1 && penergy == 75) hE75->Fill(sum);
      
      // 50 GeV electrons
      if( ptype == 1 && penergy == 50) hE50->Fill(sum);
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

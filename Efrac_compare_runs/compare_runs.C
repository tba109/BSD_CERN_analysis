#define compare_runs_cxx
#include "compare_runs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>

using namespace std;

void compare_runs::Loop()
{
  if (fChain == 0) return;
  
  TH1F * hPi = new TH1F("hPi","",500,-0.2,2);
  TH1F * hE = new TH1F("hE","",500,-0.2,2);

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
       
      // Fill Efrac for electrons
      if( ptype == 1 )
	hE->Fill(myBSDLatePE/CALSum);
      
      // Fill Efrac for pions
      if( ptype == 2 )
	hPi->Fill(myBSDLatePE/CALSum);
    }
  
  hE->SetLineColor(kRed);
  hE->SetTitle("Efrac for 350 GeV Pions (Blue) and 150 GeV Electrons (Red)");
  hE->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE->GetXaxis()->CenterTitle();
  hE->Draw();
  hPi->Scale(hE->GetEntries()/hPi->GetEntries());
  hPi->Draw("sames");

  Double_t reg_power[1000];
  Double_t e_accept[1000];
  Double_t pi_accept[1000];
  Double_t endbin[1000];
  Int_t Npts = 0;
  Int_t Nskip = 0;

  Int_t start = hE->FindBin(0);
  Int_t nbins = 40;
  // ToDo: Add error bars
  cout << "Start = " << start << endl;

  Double_t x_err[1000];
  Double_t y_err[1000];

  Double_t num = 0;
  Double_t den = 0;
  Double_t num_err = 0;
  Double_t den_err = 0;
  Double_t Ne = 0;
  Double_t Ne_total = 0;
  Double_t Npi = 0;
  Double_t Npi_total = 0;
  for( int i = 0; i < nbins; i++ )
    {
      cout << "Efrac = " << hE->GetBinLowEdge(start+i) << endl;
      if( hPi->Integral(start,start+i) > 50 && hE->Integral(start,start+i) > 50 )
	{
	  Ne = hE->Integral(start,start+i);
	  Ne_total = hE->Integral();
	  Npi = hPi->Integral(start,start+i);
	  Npi_total = hPi->Integral();
	  
	  reg_power[i-Nskip] = (Ne/Ne_total) / (Npi/Npi_total);
	  e_accept[i-Nskip] = Ne/Ne_total;
	  pi_accept[i-Nskip] = Npi/Npi_total;
	  endbin[i-Nskip] = start+i;
	  
	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);
	  num_err = (num) * sqrt(1/Ne + 1/Ne_total);
	  den_err = (den) * sqrt(1/Npi + 1/Npi_total);
	  y_err[i-Nskip] = (num/den) * sqrt( (num_err*num_err)/(num*num) + (den_err*den_err)/(den*den) );
	  x_err[i-Nskip] = num_err;
	
	  cout << i-Nskip 
	       << " " << e_accept[i-Nskip]
	       << " " << 1/pi_accept[i-Nskip]
	       << " " << reg_power[i-Nskip] 
	       << " " << hPi->Integral(start,start+i)
	       << " " << hE->Integral(start,start+i)
	       << endl;
	  Npts++;
	}
      else
	Nskip++;
    }
  
  
  // For quick and easy copy/paste graphing
  cout << "Double_t rej[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << reg_power[i] << ", " << endl;
  cout << endl;
  
  // For quick and easy copy/paste graphing
  cout << "Double_t e_accept[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << e_accept[i] << ", " << endl;
  cout << endl;
  
  cout << "Double_t x_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << x_err[i] << ", " << endl;
  cout << endl;

  cout << "Double_t y_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << y_err[i] << ", " << endl;
  cout << endl;

  TCanvas * cv2 = new TCanvas("cv2");
  TGraphErrors * gr = new TGraphErrors(Npts,e_accept,reg_power,x_err,y_err);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Electron acceptance");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("Rejection Power");
  gr->GetYaxis()->CenterTitle();
  gr->SetMarkerStyle(20);
  gr->Draw("APL");
}

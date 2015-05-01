#define Efrac_by_electron_energy_cxx
#include "Efrac_by_electron_energy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

using namespace std;

// Basically the same, but only choose a certain energy
void Efrac_by_electron_energy::Loop1(Double_t myEnergy)
{
  // The means and sigmas by energy index 0 for 75, etc, up to 4 for 175
  Int_t index = 0;
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
 
  Double_t myMean = 0;
  Double_t mySigma = 0;

  if( myEnergy == 75 )
    {
      myMean = means[0];
      mySigma = sigmas[0];
    }
  else if( myEnergy == 100 )
    {
      myMean = means[1];
      mySigma = sigmas[1];
    }
   else if( myEnergy == 125 )
    {
      myMean = means[2];
      mySigma = sigmas[2];
    }
  else if( myEnergy == 150 )
    {
      myMean = means[3];
      mySigma = sigmas[3];
    }
  else if( myEnergy == 175 )
    {
      myMean = means[4];
      mySigma = sigmas[4];
    }
  else
    {
      cout << "Wrong energy choice! You have " << myEnergy  << endl;
      return;
    }

  if (fChain == 0) return;
  
  TH1F * hPi = new TH1F("hPi","",500,-0.2,3);
  TH1F * hE = new TH1F("hE","",500,-0.2,3);
  
  Int_t NE75 = 0;
  Int_t NE100 = 0;
  Int_t NE125 = 0;
  Int_t NE150 = 0;
  Int_t NE175 = 0;

  Int_t NPi250 = 0;
  Int_t NPi300 = 0;
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
      if( ptype == 2 && penergy == 250 && NPi250 > 6617) continue; 
      if( ptype == 2 && penergy == 300 && NPi300 > 6617) continue; 
      if( ptype == 2 && penergy == 350 && NPi350 > 6617) continue; 
      
      // Don't need these, since we aren't mixing electron energies
      // if( ptype == 1 && penergy == 75 && NE75 > 2442 ) continue;
      // if( ptype == 1 && penergy == 100 && NE100 > 2442 ) continue;
      // if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      // if( ptype == 1 && penergy == 150 && NE150 > 2442 ) continue;
      // if( ptype == 1 && penergy == 175 && NE175 > 2442 ) continue;
      
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == myEnergy )
	{
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    hE->Fill(myBSDLatePE/CALSum);
	}
      
      // Fill Efrac for pions
      if( ptype == 2 && (penergy >= 250 || penergy <= 350) )
	{
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    hPi->Fill(myBSDLatePE/CALSum);
	}
      
      // Increment counts, since everything has to be equally weighted
      if( ptype == 2 && penergy == 250 ) NPi250++;
      if( ptype == 2 && penergy == 300 ) NPi300++;
      if( ptype == 2 && penergy == 350 ) NPi350++;
      
      if( ptype == 1 && penergy == 75 ) NE75++;
      if( ptype == 1 && penergy == 100 ) NE100++;
      if( ptype == 1 && penergy == 125 ) NE125++;
      if( ptype == 1 && penergy == 150 ) NE150++;
      if( ptype == 1 && penergy == 175 ) NE175++;
    }
  
  hE->SetLineColor(kRed);
  hE->SetTitle("Efrac for 250-350 GeV Pions (Blue) and 75-175 GeV Electrons (Red) and CALSum Between 6000 and 9000");
  hE->GetXaxis()->SetTitle("Efrac (arbitrary units)");
  hE->GetXaxis()->CenterTitle();
  hE->Scale(1./hE->GetEntries());
  hE->Draw();
  hPi->Scale(1./hPi->GetEntries());
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
  Double_t Ne_entries = hE->GetEntries();
  Double_t Npi_entries = hPi->GetEntries();
  for( int i = 0; i < nbins; i++ )
    {
      cout << "Efrac = " << hE->GetBinLowEdge(start+i) << endl;
      if( hPi->Integral(start,start+i) > 10./Ne_entries && 
	  hE->Integral(start,start+i) > 10./Npi_entries )
	{
	  num_err = 0;
	  den_err = 0;
	  Ne = 0;
	  Npi = 0;
	  num = 0;
	  den = 0;
	  Npi_total = 0;
	  Ne_total = 0;
	  
	  Ne = hE->Integral(start,start+i);
	  Ne_total = hE->Integral();
	  Npi = hPi->Integral(start,start+i);
	  Npi_total = hPi->Integral();
	  
	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);

	  reg_power[i-Nskip] = num/den;
	  e_accept[i-Nskip] = num;
	  pi_accept[i-Nskip] = den;
	  endbin[i-Nskip] = start+i;
	  
	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);
	  num_err = (num) * sqrt(1/(Ne*Ne_entries) + 1/(Ne_total*Ne_entries));
	  den_err = (den) * sqrt(1/(Npi*Npi_entries) + 1/(Npi_total*Npi_entries));
	  y_err[i-Nskip] = (num/den) * sqrt( (num_err*num_err)/(num*num) + (den_err*den_err)/(den*den) );
	  x_err[i-Nskip] = num_err;
	
	  cout << i-Nskip 
	       << " " << e_accept[i-Nskip]
	       << " " << pi_accept[i-Nskip]
	       << " " << reg_power[i-Nskip] 
	       << " " << hPi->Integral(start,start+i)
	       << " " << hE->Integral(start,start+i)
	       << " " << num_err
	       << " " << den_err
	       << endl;
	  Npts++;
	}
      else
	Nskip++;
    }
  
  /*  
  // For quick and easy copy/paste graphing
  cout << "rej[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << reg_power[i] << ", " << endl;
  cout << endl;
  
  // For quick and easy copy/paste graphing
  cout << "e_accept[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << e_accept[i] << ", " << endl;
  cout << endl;
  
  cout << "x_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << x_err[i] << endl;
  cout << endl;

  cout << "y_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << y_err[i] << endl;
  cout << endl;
  */
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

// Basically the same, but only choose a certain energy
void Efrac_by_electron_energy::Loop2()
{
  // The means and sigmas by energy index 0 for 75, etc, up to 4 for 175
  Int_t index = 0;
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
 
  Double_t myMean = 0;
  Double_t mySigma = 0;

  if (fChain == 0) return;
  
  TH1F * hPi = new TH1F("hPi","",500,-0.2,3);
  TH1F * hE75 = new TH1F("hE75","",500,-0.2,3);
  TH1F * hE100 = new TH1F("hE100","",500,-0.2,3);
  TH1F * hE125 = new TH1F("hE125","",500,-0.2,3);
  TH1F * hE150 = new TH1F("hE150","",500,-0.2,3);
  TH1F * hE175 = new TH1F("hE175","",500,-0.2,3);

  Int_t NE75 = 0;
  Int_t NE100 = 0;
  Int_t NE125 = 0;
  Int_t NE150 = 0;
  Int_t NE175 = 0;

  Int_t NPi250 = 0;
  Int_t NPi300 = 0;
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
      if( ptype == 2 && penergy == 250 && NPi250 > 6617) continue; 
      if( ptype == 2 && penergy == 300 && NPi300 > 6617) continue; 
      if( ptype == 2 && penergy == 350 && NPi350 > 6617) continue; 
      
      // Don't need these, since we aren't mixing electron energies
      // if( ptype == 1 && penergy == 75 && NE75 > 2442 ) continue;
      // if( ptype == 1 && penergy == 100 && NE100 > 2442 ) continue;
      // if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      // if( ptype == 1 && penergy == 150 && NE150 > 2442 ) continue;
      // if( ptype == 1 && penergy == 175 && NE175 > 2442 ) continue;
      
      if( ptype == 1 && penergy == 75 )
	{
	  myMean = means[0];
	  mySigma = sigmas[0];
	}
      else if( ptype == 1 && penergy == 100 )
	{
	  myMean = means[1];
	  mySigma = sigmas[1];
	}
      else if( ptype == 1 && penergy == 125 )
	{
	  myMean = means[2];
	  mySigma = sigmas[2];
	}
      else if( ptype == 1 && penergy == 150 )
	{
	  myMean = means[3];
	  mySigma = sigmas[3];
	}
      else if( ptype == 1 && penergy == 175 )
	{
	  myMean = means[4];
	  mySigma = sigmas[4];
	}
      else if( ptype == 1 )
	{
      continue;
	}
      
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == 75 )
	{
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    hE75->Fill(myBSDLatePE/CALSum);
	}
      
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == 100 )
	{
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    hE100->Fill(myBSDLatePE/CALSum);
	}

      // See if we're looking at the right energy
      if( ptype == 1 && penergy == 125 )
	{
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    hE125->Fill(myBSDLatePE/CALSum);
	}
      
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == 150 )
	{
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    hE150->Fill(myBSDLatePE/CALSum);
	}
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == 175 )
	{
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    hE175->Fill(myBSDLatePE/CALSum);
	}

      // Fill Efrac for pions
      if( ptype == 2 && (penergy >= 250 || penergy <= 350) )
	{
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    hPi->Fill(myBSDLatePE/CALSum);
	}
      
      // Increment counts, since everything has to be equally weighted
      if( ptype == 2 && penergy == 250 ) NPi250++;
      if( ptype == 2 && penergy == 300 ) NPi300++;
      if( ptype == 2 && penergy == 350 ) NPi350++;
      
      if( ptype == 1 && penergy == 75 ) NE75++;
      if( ptype == 1 && penergy == 100 ) NE100++;
      if( ptype == 1 && penergy == 125 ) NE125++;
      if( ptype == 1 && penergy == 150 ) NE150++;
      if( ptype == 1 && penergy == 175 ) NE175++;
    }

  ///////////////////////////////////////////////////////////////////////////////////////
  // 75 GeV electrons
  ///////////////////////////////////////////////////////////////////////////////////////
  
  Double_t reg_power75[1000];
  Double_t e_accept75[1000];
  Double_t pi_accept75[1000];
  Double_t endbin75[1000];
  Int_t Npts75 = 0;
  Int_t Nskip75 = 0;

  Int_t start75 = hE75->FindBin(0);
  Int_t nbins75 = 40;
  
  Double_t x_err75[1000];
  Double_t y_err75[1000];

  Double_t num75 = 0;
  Double_t den75 = 0;
  Double_t num_err75 = 0;
  Double_t den_err75 = 0;
  Double_t Ne75 = 0;
  Double_t Ne_total75 = 0;
  Double_t Npi75 = 0;
  Double_t Npi_total75 = 0;
  Double_t Ne_entries75 = hE75->GetEntries();
  Double_t Npi_entries75 = hPi->GetEntries();
  for( int i = 0; i < nbins75; i++ )
    {
      cout << "Efrac = " << hE75->GetBinLowEdge(start75+i) << endl;
      if( hPi->Integral(start75,start75+i) > 100./Ne_entries75 && 
	  hE75->Integral(start75,start75+i) > 100./Npi_entries75 )
	{
	  num_err75 = 0;
	  den_err75 = 0;
	  Ne75 = 0;
	  Npi75 = 0;
	  num75 = 0;
	  den75 = 0;
	  Npi_total75 = 0;
	  Ne_total75 = 0;
	  
	  Ne75 = hE75->Integral(start75,start75+i);
	  Ne_total75 = hE75->Integral();
	  Npi75 = hPi->Integral(start75,start75+i);
	  Npi_total75 = hPi->Integral();
	  
	  num75 = (Ne75/Ne_total75);
	  den75 = (Npi75/Npi_total75);

	  reg_power75[i-Nskip75] = num75/den75;
	  e_accept75[i-Nskip75] = num75;
	  pi_accept75[i-Nskip75] = den75;
	  endbin75[i-Nskip75] = start75+i;
	  
	  num75 = (Ne75/Ne_total75);
	  den75 = (Npi75/Npi_total75);
	  num_err75 = (num75) * sqrt(1/(Ne75*Ne_entries75) + 1/(Ne_total75*Ne_entries75));
	  den_err75 = (den75) * sqrt(1/(Npi75*Npi_entries75) + 1/(Npi_total75*Npi_entries75));
	  y_err75[i-Nskip75] = (num75/den75) * sqrt( (num_err75*num_err75)/(num75*num75) + (den_err75*den_err75)/(den75*den75) );
	  x_err75[i-Nskip75] = num_err75;
	
	  cout << i-Nskip75 
	       << " " << e_accept75[i-Nskip75]
	       << " " << pi_accept75[i-Nskip75]
	       << " " << reg_power75[i-Nskip75] 
	       << " " << hPi->Integral(start75,start75+i)
	       << " " << hE75->Integral(start75,start75+i)
	       << " " << num_err75
	       << " " << den_err75
	       << endl;
	  Npts75++;
	}
      else
	{
	  Nskip75++;
	  cout << "No go" << " "
	       << hPi->Integral(start75,start75+i) << " " 
	       << hE75->Integral(start75,start75+i) << " "
	       << endl;
	}
    }
  
  TCanvas * cv75 = new TCanvas("cv75");
  TGraphErrors * gr75 = new TGraphErrors(Npts75,e_accept75,reg_power75,x_err75,y_err75);
  gr75->SetTitle("");
  gr75->GetXaxis()->SetTitle("Electron acceptance");
  gr75->GetXaxis()->CenterTitle();
  gr75->GetYaxis()->SetTitle("Rejection Power");
  gr75->GetYaxis()->CenterTitle();
  gr75->SetMarkerStyle(20);
  gr75->Draw("APL");


  ///////////////////////////////////////////////////////////////////////////////////////
  // 100 GeV electrons
  ///////////////////////////////////////////////////////////////////////////////////////

  Double_t reg_power100[1000];
  Double_t e_accept100[1000];
  Double_t pi_accept100[1000];
  Double_t endbin100[1000];
  Int_t Npts100 = 0;
  Int_t Nskip100 = 0;

  Int_t start100 = hE100->FindBin(0);
  Int_t nbins100 = 40;
  
  Double_t x_err100[1000];
  Double_t y_err100[1000];

  Double_t num100 = 0;
  Double_t den100 = 0;
  Double_t num_err100 = 0;
  Double_t den_err100 = 0;
  Double_t Ne100 = 0;
  Double_t Ne_total100 = 0;
  Double_t Npi100 = 0;
  Double_t Npi_total100 = 0;
  Double_t Ne_entries100 = hE100->GetEntries();
  Double_t Npi_entries100 = hPi->GetEntries();
  for( int i = 0; i < nbins100; i++ )
    {
      cout << "Efrac = " << hE100->GetBinLowEdge(start100+i) << endl;
      if( hPi->Integral(start100,start100+i) > 100./Ne_entries100 && 
	  hE100->Integral(start100,start100+i) > 100./Npi_entries100 )
	{
	  num_err100 = 0;
	  den_err100 = 0;
	  Ne100 = 0;
	  Npi100 = 0;
	  num100 = 0;
	  den100 = 0;
	  Npi_total100 = 0;
	  Ne_total100 = 0;
	  
	  Ne100 = hE100->Integral(start100,start100+i);
	  Ne_total100 = hE100->Integral();
	  Npi100 = hPi->Integral(start100,start100+i);
	  Npi_total100 = hPi->Integral();
	  
	  num100 = (Ne100/Ne_total100);
	  den100 = (Npi100/Npi_total100);

	  reg_power100[i-Nskip100] = num100/den100;
	  e_accept100[i-Nskip100] = num100;
	  pi_accept100[i-Nskip100] = den100;
	  endbin100[i-Nskip100] = start100+i;
	  
	  num100 = (Ne100/Ne_total100);
	  den100 = (Npi100/Npi_total100);
	  num_err100 = (num100) * sqrt(1/(Ne100*Ne_entries100) + 1/(Ne_total100*Ne_entries100));
	  den_err100 = (den100) * sqrt(1/(Npi100*Npi_entries100) + 1/(Npi_total100*Npi_entries100));
	  y_err100[i-Nskip100] = (num100/den100) * sqrt( (num_err100*num_err100)/(num100*num100) + (den_err100*den_err100)/(den100*den100) );
	  x_err100[i-Nskip100] = num_err100;
	
	  cout << i-Nskip100 
	       << " " << e_accept100[i-Nskip100]
	       << " " << pi_accept100[i-Nskip100]
	       << " " << reg_power100[i-Nskip100] 
	       << " " << hPi->Integral(start100,start100+i)
	       << " " << hE100->Integral(start100,start100+i)
	       << " " << num_err100
	       << " " << den_err100
	       << endl;
	  Npts100++;
	}
      else
	Nskip100++;
    }
  
  TCanvas * cv100 = new TCanvas("cv100");
  TGraphErrors * gr100 = new TGraphErrors(Npts100,e_accept100,reg_power100,x_err100,y_err100);
  gr100->SetTitle("");
  gr100->GetXaxis()->SetTitle("Electron acceptance");
  gr100->GetXaxis()->CenterTitle();
  gr100->GetYaxis()->SetTitle("Rejection Power");
  gr100->GetYaxis()->CenterTitle();
  gr100->SetMarkerStyle(20);
  gr100->Draw("APL");


  ///////////////////////////////////////////////////////////////////////////////////////
  // 125 GeV electrons
  ///////////////////////////////////////////////////////////////////////////////////////
  
  Double_t reg_power125[1000];
  Double_t e_accept125[1000];
  Double_t pi_accept125[1000];
  Double_t endbin125[1000];
  Int_t Npts125 = 0;
  Int_t Nskip125 = 0;

  Int_t start125 = hE125->FindBin(0);
  Int_t nbins125 = 40;
  
  Double_t x_err125[1000];
  Double_t y_err125[1000];

  Double_t num125 = 0;
  Double_t den125 = 0;
  Double_t num_err125 = 0;
  Double_t den_err125 = 0;
  Double_t Ne125 = 0;
  Double_t Ne_total125 = 0;
  Double_t Npi125 = 0;
  Double_t Npi_total125 = 0;
  Double_t Ne_entries125 = hE125->GetEntries();
  Double_t Npi_entries125 = hPi->GetEntries();
  for( int i = 0; i < nbins125; i++ )
    {
      cout << "Efrac = " << hE125->GetBinLowEdge(start125+i) << endl;
      if( hPi->Integral(start125,start125+i) > 100./Ne_entries125 && 
	  hE125->Integral(start125,start125+i) > 100./Npi_entries125 )
	{
	  num_err125 = 0;
	  den_err125 = 0;
	  Ne125 = 0;
	  Npi125 = 0;
	  num125 = 0;
	  den125 = 0;
	  Npi_total125 = 0;
	  Ne_total125 = 0;
	  
	  Ne125 = hE125->Integral(start125,start125+i);
	  Ne_total125 = hE125->Integral();
	  Npi125 = hPi->Integral(start125,start125+i);
	  Npi_total125 = hPi->Integral();
	  
	  num125 = (Ne125/Ne_total125);
	  den125 = (Npi125/Npi_total125);

	  reg_power125[i-Nskip125] = num125/den125;
	  e_accept125[i-Nskip125] = num125;
	  pi_accept125[i-Nskip125] = den125;
	  endbin125[i-Nskip125] = start125+i;
	  
	  num125 = (Ne125/Ne_total125);
	  den125 = (Npi125/Npi_total125);
	  num_err125 = (num125) * sqrt(1/(Ne125*Ne_entries125) + 1/(Ne_total125*Ne_entries125));
	  den_err125 = (den125) * sqrt(1/(Npi125*Npi_entries125) + 1/(Npi_total125*Npi_entries125));
	  y_err125[i-Nskip125] = (num125/den125) * sqrt( (num_err125*num_err125)/(num125*num125) + (den_err125*den_err125)/(den125*den125) );
	  x_err125[i-Nskip125] = num_err125;
	
	  cout << i-Nskip125 
	       << " " << e_accept125[i-Nskip125]
	       << " " << pi_accept125[i-Nskip125]
	       << " " << reg_power125[i-Nskip125] 
	       << " " << hPi->Integral(start125,start125+i)
	       << " " << hE125->Integral(start125,start125+i)
	       << " " << num_err125
	       << " " << den_err125
	       << endl;
	  Npts125++;
	}
      else
	Nskip125++;
    }
  
  TCanvas * cv125 = new TCanvas("cv125");
  TGraphErrors * gr125 = new TGraphErrors(Npts125,e_accept125,reg_power125,x_err125,y_err125);
  gr125->SetTitle("");
  gr125->GetXaxis()->SetTitle("Electron acceptance");
  gr125->GetXaxis()->CenterTitle();
  gr125->GetYaxis()->SetTitle("Rejection Power");
  gr125->GetYaxis()->CenterTitle();
  gr125->SetMarkerStyle(20);
  gr125->Draw("APL");

  ///////////////////////////////////////////////////////////////////////////////////////
  // 150 GeV electrons
  ///////////////////////////////////////////////////////////////////////////////////////
  
  Double_t reg_power150[1000];
  Double_t e_accept150[1000];
  Double_t pi_accept150[1000];
  Double_t endbin150[1000];
  Int_t Npts150 = 0;
  Int_t Nskip150 = 0;

  Int_t start150 = hE150->FindBin(0);
  Int_t nbins150 = 40;
  
  Double_t x_err150[1000];
  Double_t y_err150[1000];

  Double_t num150 = 0;
  Double_t den150 = 0;
  Double_t num_err150 = 0;
  Double_t den_err150 = 0;
  Double_t Ne150 = 0;
  Double_t Ne_total150 = 0;
  Double_t Npi150 = 0;
  Double_t Npi_total150 = 0;
  Double_t Ne_entries150 = hE150->GetEntries();
  Double_t Npi_entries150 = hPi->GetEntries();
  for( int i = 0; i < nbins150; i++ )
    {
      cout << "Efrac = " << hE150->GetBinLowEdge(start150+i) << endl;
      if( hPi->Integral(start150,start150+i) > 100./Ne_entries150 && 
	  hE150->Integral(start150,start150+i) > 100./Npi_entries150 )
	{
	  num_err150 = 0;
	  den_err150 = 0;
	  Ne150 = 0;
	  Npi150 = 0;
	  num150 = 0;
	  den150 = 0;
	  Npi_total150 = 0;
	  Ne_total150 = 0;
	  
	  Ne150 = hE150->Integral(start150,start150+i);
	  Ne_total150 = hE150->Integral();
	  Npi150 = hPi->Integral(start150,start150+i);
	  Npi_total150 = hPi->Integral();
	  
	  num150 = (Ne150/Ne_total150);
	  den150 = (Npi150/Npi_total150);

	  reg_power150[i-Nskip150] = num150/den150;
	  e_accept150[i-Nskip150] = num150;
	  pi_accept150[i-Nskip150] = den150;
	  endbin150[i-Nskip150] = start150+i;
	  
	  num150 = (Ne150/Ne_total150);
	  den150 = (Npi150/Npi_total150);
	  num_err150 = (num150) * sqrt(1/(Ne150*Ne_entries150) + 1/(Ne_total150*Ne_entries150));
	  den_err150 = (den150) * sqrt(1/(Npi150*Npi_entries150) + 1/(Npi_total150*Npi_entries150));
	  y_err150[i-Nskip150] = (num150/den150) * sqrt( (num_err150*num_err150)/(num150*num150) + (den_err150*den_err150)/(den150*den150) );
	  x_err150[i-Nskip150] = num_err150;
	
	  cout << i-Nskip150 
	       << " " << e_accept150[i-Nskip150]
	       << " " << pi_accept150[i-Nskip150]
	       << " " << reg_power150[i-Nskip150] 
	       << " " << hPi->Integral(start150,start150+i)
	       << " " << hE150->Integral(start150,start150+i)
	       << " " << num_err150
	       << " " << den_err150
	       << endl;
	  Npts150++;
	}
      else
	Nskip150++;
    }
  
  TCanvas * cv150 = new TCanvas("cv150");
  TGraphErrors * gr150 = new TGraphErrors(Npts150,e_accept150,reg_power150,x_err150,y_err150);
  gr150->SetTitle("");
  gr150->GetXaxis()->SetTitle("Electron acceptance");
  gr150->GetXaxis()->CenterTitle();
  gr150->GetYaxis()->SetTitle("Rejection Power");
  gr150->GetYaxis()->CenterTitle();
  gr150->SetMarkerStyle(20);
  gr150->Draw("APL");

  ///////////////////////////////////////////////////////////////////////////////////////
  // 175 GeV electrons
  ///////////////////////////////////////////////////////////////////////////////////////
  
  Double_t reg_power175[1000];
  Double_t e_accept175[1000];
  Double_t pi_accept175[1000];
  Double_t endbin175[1000];
  Int_t Npts175 = 0;
  Int_t Nskip175 = 0;

  Int_t start175 = hE175->FindBin(0);
  Int_t nbins175 = 40;
  
  Double_t x_err175[1000];
  Double_t y_err175[1000];

  Double_t num175 = 0;
  Double_t den175 = 0;
  Double_t num_err175 = 0;
  Double_t den_err175 = 0;
  Double_t Ne175 = 0;
  Double_t Ne_total175 = 0;
  Double_t Npi175 = 0;
  Double_t Npi_total175 = 0;
  Double_t Ne_entries175 = hE175->GetEntries();
  Double_t Npi_entries175 = hPi->GetEntries();
  for( int i = 0; i < nbins175; i++ )
    {
      cout << "Efrac = " << hE175->GetBinLowEdge(start175+i) << endl;
      if( hPi->Integral(start175,start175+i) > 100./Ne_entries175 && 
	  hE175->Integral(start175,start175+i) > 100./Npi_entries175 )
	{
	  num_err175 = 0;
	  den_err175 = 0;
	  Ne175 = 0;
	  Npi175 = 0;
	  num175 = 0;
	  den175 = 0;
	  Npi_total175 = 0;
	  Ne_total175 = 0;
	  
	  Ne175 = hE175->Integral(start175,start175+i);
	  Ne_total175 = hE175->Integral();
	  Npi175 = hPi->Integral(start175,start175+i);
	  Npi_total175 = hPi->Integral();
	  
	  num175 = (Ne175/Ne_total175);
	  den175 = (Npi175/Npi_total175);

	  reg_power175[i-Nskip175] = num175/den175;
	  e_accept175[i-Nskip175] = num175;
	  pi_accept175[i-Nskip175] = den175;
	  endbin175[i-Nskip175] = start175+i;
	  
	  num175 = (Ne175/Ne_total175);
	  den175 = (Npi175/Npi_total175);
	  num_err175 = (num175) * sqrt(1/(Ne175*Ne_entries175) + 1/(Ne_total175*Ne_entries175));
	  den_err175 = (den175) * sqrt(1/(Npi175*Npi_entries175) + 1/(Npi_total175*Npi_entries175));
	  y_err175[i-Nskip175] = (num175/den175) * sqrt( (num_err175*num_err175)/(num175*num175) + (den_err175*den_err175)/(den175*den175) );
	  x_err175[i-Nskip175] = num_err175;
	
	  cout << i-Nskip175 
	       << " " << e_accept175[i-Nskip175]
	       << " " << pi_accept175[i-Nskip175]
	       << " " << reg_power175[i-Nskip175] 
	       << " " << hPi->Integral(start175,start175+i)
	       << " " << hE175->Integral(start175,start175+i)
	       << " " << num_err175
	       << " " << den_err175
	       << endl;
	  Npts175++;
	}
      else
	Nskip175++;
    }
  
  TCanvas * cv175 = new TCanvas("cv175");
  TGraphErrors * gr175 = new TGraphErrors(Npts175,e_accept175,reg_power175,x_err175,y_err175);
  gr175->SetTitle("");
  gr175->GetXaxis()->SetTitle("Electron acceptance");
  gr175->GetXaxis()->CenterTitle();
  gr175->GetYaxis()->SetTitle("Rejection Power");
  gr175->GetYaxis()->CenterTitle();
  gr175->SetMarkerStyle(20);
  gr175->Draw("APL");

}

// Basically the same, but only choose a certain energy
void Efrac_by_electron_energy::Loop3(Double_t myEnergyE, Double_t myEnergyP)
{
  // The means and sigmas by energy index 0 for 75, etc, up to 4 for 175
  Int_t index = 0;
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
 
  Double_t myMean = 0;
  Double_t mySigma = 0;

  if( myEnergyE == 75 )
    {
      myMean = means[0];
      mySigma = sigmas[0];
    }
  else if( myEnergyE == 100 )
    {
      myMean = means[1];
      mySigma = sigmas[1];
    }
   else if( myEnergyE == 125 )
    {
      myMean = means[2];
      mySigma = sigmas[2];
    }
  else if( myEnergyE == 150 )
    {
      myMean = means[3];
      mySigma = sigmas[3];
    }
  else if( myEnergyE == 175 )
    {
      myMean = means[4];
      mySigma = sigmas[4];
    }
  else
    {
      cout << "Wrong energy choice! You have " << myEnergyE  << endl;
      return;
    }

  if (fChain == 0) return;
  
  TH1F * hPi = new TH1F("hPi","",500,-0.2,3);
  TH1F * hE = new TH1F("hE","",500,-0.2,3);

  TH1F * hPiCal = new TH1F("hPiCal","",200,-1000,30000);
  TH1F * hECal = new TH1F("hECal","",200,-1000,30000);

  TH1F * hPiBSD = new TH1F("hPiBSD","",200,-100,8000);
  TH1F * hEBSD = new TH1F("hEBSD","",200,-100,8000);
    
  Int_t NE75 = 0;
  Int_t NE100 = 0;
  Int_t NE125 = 0;
  Int_t NE150 = 0;
  Int_t NE175 = 0;

  Int_t NPi250 = 0;
  Int_t NPi300 = 0;
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
      // if( ptype == 2 && penergy == 250 && NPi250 > 6617) continue; 
      // if( ptype == 2 && penergy == 300 && NPi300 > 6617) continue; 
      // if( ptype == 2 && penergy == 350 && NPi350 > 6617) continue; 
      
      // Don't need these, since we aren't mixing electron energies
      // if( ptype == 1 && penergy == 75 && NE75 > 2442 ) continue;
      // if( ptype == 1 && penergy == 100 && NE100 > 2442 ) continue;
      // if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      // if( ptype == 1 && penergy == 150 && NE150 > 2442 ) continue;
      // if( ptype == 1 && penergy == 175 && NE175 > 2442 ) continue;
      
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == myEnergyE )
	{
	  hECal->Fill(CALSum);
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    {
	      hE->Fill(myBSDLatePE/CALSum);
	      hEBSD->Fill(myBSDLatePE);
	    }
	}
      
      // Fill Efrac for pions
      if( ptype == 2 && penergy == myEnergyP )
	{
	  hPiCal->Fill(CALSum);
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    {
	      hPi->Fill(myBSDLatePE/CALSum);
	      hPiBSD->Fill(myBSDLatePE);
	    }
	}


      // Increment counts, since everything has to be equally weighted
      if( ptype == 2 && penergy == 250 ) NPi250++;
      if( ptype == 2 && penergy == 300 ) NPi300++;
      if( ptype == 2 && penergy == 350 ) NPi350++;
      
      if( ptype == 1 && penergy == 75 ) NE75++;
      if( ptype == 1 && penergy == 100 ) NE100++;
      if( ptype == 1 && penergy == 125 ) NE125++;
      if( ptype == 1 && penergy == 150 ) NE150++;
      if( ptype == 1 && penergy == 175 ) NE175++;
    }

  TCanvas * cv1 = new TCanvas("cv1");
  
  gStyle->SetOptStat(0);

  hPiCal->SetLineColor(kBlue);
  hECal->SetLineColor(kRed);
  hPiCal->Scale(1./hPiCal->GetEntries());
  hECal->Scale(1./hECal->GetEntries());
  hECal->GetXaxis()->SetTitle("Calorimeter Signal (arbitrary units)");
  hECal->GetXaxis()->CenterTitle();
  hECal->Draw();
  hPiCal->Draw("sames");
  
  TF1 * f1 = new TF1("f1","gaus",myMean - 2*mySigma, myMean + 2*mySigma);
  hECal->Fit("f1","R","same");

  cv1->Update();

  TLine * l1 = new TLine(myMean - 2*mySigma, 0, myMean - 2*mySigma, gPad->GetUymax());
  TLine * l2 = new TLine(myMean + 2*mySigma, 0, myMean + 2*mySigma, gPad->GetUymax());

  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  
  l1->Draw("same");
  l2->Draw("same");
  
  TCanvas * cv2 = new TCanvas("cv2");
  hEBSD->SetLineColor(kRed);
  hPiBSD->SetLineColor(kBlue);
  hEBSD->Scale(1./hEBSD->GetEntries());
  hPiBSD->Scale(1./hPiBSD->GetEntries());
  hEBSD->GetXaxis()->SetTitle("BSD Signal (arbitrary units)");
  hEBSD->GetXaxis()->CenterTitle();
  hEBSD->Draw("");
  hPiBSD->Draw("sames");
  


  TCanvas * cv3 = new TCanvas("cv3");
  
  hE->SetLineColor(kRed);
  hE->SetTitle("");
  hE->GetXaxis()->SetTitle("E_{frac} (\"detector\" units)");
  hE->GetXaxis()->CenterTitle();
  hE->Scale(1./hE->GetEntries());
  hE->Draw();
  hPi->Scale(1./hPi->GetEntries());
  hPi->Draw("sames");

  Double_t reg_power[1000];
  Double_t e_accept[1000];
  Double_t pi_accept[1000];
  Double_t endbin[1000];
  Int_t Npts = 0;
  Int_t Nskip = 0;
  Double_t efrac[1000];

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
  Double_t Ne_entries = hE->GetEntries();
  Double_t Npi_entries = hPi->GetEntries();
  Double_t Efrac_temp = 0;

  for( int i = 0; i < nbins; i++ )
    {
      cout << "Efrac = " << hE->GetBinLowEdge(start+i) << endl;
      Efrac_temp = hE->GetBinLowEdge(start+i);
      
      if( hPi->Integral(start,start+i) > 10./Ne_entries && 
	  hE->Integral(start,start+i) > 10./Npi_entries )
	{
	  num_err = 0;
	  den_err = 0;
	  Ne = 0;
	  Npi = 0;
	  num = 0;
	  den = 0;
	  Npi_total = 0;
	  Ne_total = 0;
	  
	  Ne = hE->Integral(start,start+i);
	  Ne_total = hE->Integral();
	  Npi = hPi->Integral(start,start+i);
	  Npi_total = hPi->Integral();
	  
	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);

	  reg_power[i-Nskip] = num/den;
	  e_accept[i-Nskip] = num;
	  pi_accept[i-Nskip] = den;
	  endbin[i-Nskip] = start+i;
	  efrac[i-Nskip] = Efrac_temp;

	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);
	  num_err = (num) * sqrt(1/(Ne*Ne_entries) + 1/(Ne_total*Ne_entries));
	  den_err = (den) * sqrt(1/(Npi*Npi_entries) + 1/(Npi_total*Npi_entries));
	  y_err[i-Nskip] = (num/den) * sqrt( (num_err*num_err)/(num*num) + (den_err*den_err)/(den*den) );
	  x_err[i-Nskip] = num_err;
	
	  cout << i-Nskip 
	       << " " << e_accept[i-Nskip]
	       << " " << pi_accept[i-Nskip]
	       << " " << reg_power[i-Nskip] 
	       << " " << hPi->Integral(start,start+i)
	       << " " << hE->Integral(start,start+i)
	       << " " << num_err
	       << " " << den_err
	       << endl;
	  Npts++;
	}
      else
	Nskip++;
    }
  
  /*  
  // For quick and easy copy/paste graphing
  cout << "rej[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << reg_power[i] << ", " << endl;
  cout << endl;
  
  // For quick and easy copy/paste graphing
  cout << "e_accept[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << e_accept[i] << ", " << endl;
  cout << endl;
  
  cout << "x_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << x_err[i] << endl;
  cout << endl;

  cout << "y_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << y_err[i] << endl;
  cout << endl;
  */
  TCanvas * cv4 = new TCanvas("cv4");
  TGraphErrors * gr = new TGraphErrors(Npts,e_accept,reg_power,x_err,y_err);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Electron Acceptance");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("Rejection Power");
  gr->GetYaxis()->CenterTitle();
  gr->SetMarkerStyle(20);
  gr->Draw("APL");

  TCanvas * cv5 = new TCanvas("cv5");
  TGraphErrors * gr5 = new TGraphErrors(10,efrac,e_accept,0,x_err);
  gr5->SetTitle("");
  gr5->SetMarkerStyle(20);
  gr5->GetXaxis()->SetTitle("E_{frac} (arbitrary units)");
  gr5->GetXaxis()->CenterTitle();
  gr5->GetYaxis()->SetTitle("Electron Acceptance");
  gr5->GetYaxis()->CenterTitle();
  gr5->Draw("APL");
}


// myBSDLatePE and CALProfile[19] seem to be somewhat correlated...
// That's a bit disturbing from the perspective of neutron captures
// I wonder how much the BSD actually adds compared to the last
// layer of the CAL...
void Efrac_by_electron_energy::Loop4(Double_t myEnergyE, Double_t myEnergyP, Int_t CalLayer)
{
  // The means and sigmas by energy index 0 for 75, etc, up to 4 for 175
  Int_t index = 0;
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
 
  Double_t myMean = 0;
  Double_t mySigma = 0;

  TH2F * hCAL19vsBSDE = new TH2F("hCAL19vsBSDE","",200,-100,1000,200,-100,8000);
  TH2F * hCAL19vsBSDPi = new TH2F("hCAL19vsBSDPi","",200,-100,1000,200,-100,8000);
  TH2F * hCAL19vsBSDE2 = new TH2F("hCAL19vsBSDE2","",100,-100,1000,100,-100,8000);
  // TH2F * hCAL19vsBSDPi2 = new TH2F("hCAL19vsBSDPi2","",100,-100,1000,100,-100,8000);
  TH2F * hCAL19vsBSDPi2 = new TH2F("hCAL19vsBSDPi2","",50,-5,80,50,-100,8000);

  if( myEnergyE == 75 )
    {
      myMean = means[0];
      mySigma = sigmas[0];
    }
  else if( myEnergyE == 100 )
    {
      myMean = means[1];
      mySigma = sigmas[1];
    }
   else if( myEnergyE == 125 )
    {
      myMean = means[2];
      mySigma = sigmas[2];
    }
  else if( myEnergyE == 150 )
    {
      myMean = means[3];
      mySigma = sigmas[3];
    }
  else if( myEnergyE == 175 )
    {
      myMean = means[4];
      mySigma = sigmas[4];
    }
  else
    {
      cout << "Wrong energy choice! You have " << myEnergyE  << endl;
      return;
    }

  if (fChain == 0) return;
  
  TH1F * hPi = new TH1F("hPi","",500,-0.2,3);
  TH1F * hE = new TH1F("hE","",500,-0.2,3);

  TH1F * hPiCal = new TH1F("hPiCal","",200,-1000,30000);
  TH1F * hECal = new TH1F("hECal","",200,-1000,30000);

  TH1F * hPiBSD = new TH1F("hPiBSD","",200,-100,8000);
  TH1F * hEBSD = new TH1F("hEBSD","",200,-100,8000);

  TH1F * hPiCAL_19 = new TH1F("hPiCAL_19","",200,-100,1000);
  TH1F * hECAL_19 = new TH1F("hECAL_19","",200,-100,1000);
    
  Int_t NE75 = 0;
  Int_t NE100 = 0;
  Int_t NE125 = 0;
  Int_t NE150 = 0;
  Int_t NE175 = 0;

  Int_t NPi250 = 0;
  Int_t NPi300 = 0;
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
      // if( ptype == 2 && penergy == 250 && NPi250 > 6617) continue; 
      // if( ptype == 2 && penergy == 300 && NPi300 > 6617) continue; 
      // if( ptype == 2 && penergy == 350 && NPi350 > 6617) continue; 
      
      // Don't need these, since we aren't mixing electron energies
      // if( ptype == 1 && penergy == 75 && NE75 > 2442 ) continue;
      // if( ptype == 1 && penergy == 100 && NE100 > 2442 ) continue;
      // if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      // if( ptype == 1 && penergy == 150 && NE150 > 2442 ) continue;
      // if( ptype == 1 && penergy == 175 && NE175 > 2442 ) continue;
      
      // if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      // if( ptype == 2 && penergy == 350 && NPi350 > 2442) continue; 
      
      // if( ptype == 2 && penergy == 350 && NPi350 > 4134) continue; 
      
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == myEnergyE )
	{
	  hECal->Fill(CALSum);
	  hCAL19vsBSDE2->Fill(CALProfile[CalLayer],myBSDLatePE);
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    {
	      hE->Fill(myBSDLatePE/CALSum);
	      hEBSD->Fill(myBSDLatePE);
	      hECAL_19->Fill(CALProfile[CalLayer]);
	      hCAL19vsBSDE->Fill(CALProfile[CalLayer],myBSDLatePE);
	    }
	}
      
      // Fill Efrac for pions
      if( ptype == 2 && penergy == myEnergyP )
	{
	  hPiCal->Fill(CALSum);
	  // hCAL19vsBSDPi2->Fill(CALProfile[CalLayer],myBSDLatePE);
	  // plot in MeV
	  hCAL19vsBSDPi2->Fill(1.34*CALProfile[CalLayer]/9,myBSDLatePE);
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    {
	      hPi->Fill(myBSDLatePE/CALSum);
	      hPiBSD->Fill(myBSDLatePE);
	      hPiCAL_19->Fill(CALProfile[CalLayer]);
	      hCAL19vsBSDPi->Fill(CALProfile[CalLayer],myBSDLatePE);
	    }
	}


      // Increment counts, since everything has to be equally weighted
      if( ptype == 2 && penergy == 250 ) NPi250++;
      if( ptype == 2 && penergy == 300 ) NPi300++;
      if( ptype == 2 && penergy == 350 ) NPi350++;
      
      if( ptype == 1 && penergy == 75 ) NE75++;
      if( ptype == 1 && penergy == 100 ) NE100++;
      if( ptype == 1 && penergy == 125 ) NE125++;
      if( ptype == 1 && penergy == 150 ) NE150++;
      if( ptype == 1 && penergy == 175 ) NE175++;
    }

  TCanvas * cv1 = new TCanvas("cv1");
  
  // gStyle->SetOptStat(0);

  hPiCal->SetLineColor(kBlue);
  hECal->SetLineColor(kRed);
  hPiCal->Scale(1./hPiCal->GetEntries());
  hECal->Scale(1./hECal->GetEntries());
  hECal->GetXaxis()->SetTitle("Calorimeter Signal (arbitrary units)");
  hECal->GetXaxis()->CenterTitle();
  hECal->Draw();
  hPiCal->Draw("sames");
  
  TF1 * f1 = new TF1("f1","gaus",myMean - 2*mySigma, myMean + 2*mySigma);
  hECal->Fit("f1","R","same");

  cv1->Update();

  TLine * l1 = new TLine(myMean - 2*mySigma, 0, myMean - 2*mySigma, gPad->GetUymax());
  TLine * l2 = new TLine(myMean + 2*mySigma, 0, myMean + 2*mySigma, gPad->GetUymax());

  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  
  l1->Draw("same");
  l2->Draw("same");
  
  TCanvas * cv2 = new TCanvas("cv2");
  hEBSD->SetLineColor(kRed);
  hPiBSD->SetLineColor(kBlue);
  hEBSD->Scale(1./hEBSD->GetEntries());
  hPiBSD->Scale(1./hPiBSD->GetEntries());
  hEBSD->GetXaxis()->SetTitle("BSD Signal (arbitrary units)");
  hEBSD->GetXaxis()->CenterTitle();
  hEBSD->Draw("");
  hPiBSD->Draw("sames");
  
  TCanvas * cv3 = new TCanvas("cv3");
  
  hE->SetLineColor(kRed);
  hE->SetTitle("");
  hE->GetXaxis()->SetTitle("E_{frac} (\"detector\" units)");
  hE->GetXaxis()->CenterTitle();
  hE->Scale(1./hE->GetEntries());
  hE->Draw();
  hPi->Scale(1./hPi->GetEntries());
  hPi->Draw("sames");

  Double_t reg_power[1000];
  Double_t e_accept[1000];
  Double_t pi_accept[1000];
  Double_t endbin[1000];
  Int_t Npts = 0;
  Int_t Nskip = 0;
  Double_t efrac[1000];

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
  Double_t Ne_entries = hE->GetEntries();
  Double_t Npi_entries = hPi->GetEntries();
  Double_t Efrac_temp = 0;

  for( int i = 0; i < nbins; i++ )
    {
      cout << "Efrac = " << hE->GetBinLowEdge(start+i) << endl;
      Efrac_temp = hE->GetBinLowEdge(start+i);
      
      if( hPi->Integral(start,start+i) > 10./Ne_entries && 
	  hE->Integral(start,start+i) > 10./Npi_entries )
	{
	  num_err = 0;
	  den_err = 0;
	  Ne = 0;
	  Npi = 0;
	  num = 0;
	  den = 0;
	  Npi_total = 0;
	  Ne_total = 0;
	  
	  Ne = hE->Integral(start,start+i);
	  Ne_total = hE->Integral();
	  Npi = hPi->Integral(start,start+i);
	  Npi_total = hPi->Integral();
	  
	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);

	  reg_power[i-Nskip] = num/den;
	  e_accept[i-Nskip] = num;
	  pi_accept[i-Nskip] = den;
	  endbin[i-Nskip] = start+i;
	  efrac[i-Nskip] = Efrac_temp;

	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);
	  num_err = (num) * sqrt(1/(Ne*Ne_entries) + 1/(Ne_total*Ne_entries));
	  den_err = (den) * sqrt(1/(Npi*Npi_entries) + 1/(Npi_total*Npi_entries));
	  y_err[i-Nskip] = (num/den) * sqrt( (num_err*num_err)/(num*num) + (den_err*den_err)/(den*den) );
	  x_err[i-Nskip] = num_err;
	
	  cout << i-Nskip 
	       << " " << e_accept[i-Nskip]
	       << " " << pi_accept[i-Nskip]
	       << " " << reg_power[i-Nskip] 
	       << " " << hPi->Integral(start,start+i)
	       << " " << hE->Integral(start,start+i)
	       << " " << num_err
	       << " " << den_err
	       << endl;
	  Npts++;
	}
      else
	Nskip++;
    }
  
  /*  
  // For quick and easy copy/paste graphing
  cout << "rej[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << reg_power[i] << ", " << endl;
  cout << endl;
  
  // For quick and easy copy/paste graphing
  cout << "e_accept[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << e_accept[i] << ", " << endl;
  cout << endl;
  
  cout << "x_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << x_err[i] << endl;
  cout << endl;

  cout << "y_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << y_err[i] << endl;
  cout << endl;
  */
  TCanvas * cv4 = new TCanvas("cv4");
  TGraphErrors * gr = new TGraphErrors(Npts,e_accept,reg_power,x_err,y_err);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Electron Acceptance");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("Rejection Power");
  gr->GetYaxis()->CenterTitle();
  gr->SetMarkerStyle(20);
  gr->Draw("APL");

  TCanvas * cv5 = new TCanvas("cv5");
  TGraphErrors * gr5 = new TGraphErrors(10,efrac,e_accept,0,x_err);
  gr5->SetTitle("");
  gr5->SetMarkerStyle(20);
  gr5->GetXaxis()->SetTitle("E_{frac} (arbitrary units)");
  gr5->GetXaxis()->CenterTitle();
  gr5->GetYaxis()->SetTitle("Electron Acceptance");
  gr5->GetYaxis()->CenterTitle();
  gr5->Draw("APL");

  // This is looking at how much rejection power the bottom layer of the calorimeter
  // has, since it seems to be correlated to the BSD signal...
  TCanvas * cv6 = new TCanvas("cv6");
  hECAL_19->SetLineColor(kRed);
  hPiCAL_19->SetLineColor(kBlue);
  hECAL_19->Scale(1./hECAL_19->GetEntries());
  hPiCAL_19->Scale(1./hPiCAL_19->GetEntries());

  hECAL_19->Draw();
  hPiCAL_19->Draw("sames");
  
  TCanvas * cv7 = new TCanvas("cv7");
  hCAL19vsBSDE->SetMarkerColor(kRed);
  hCAL19vsBSDPi->SetMarkerColor(kBlue);
  // hCAL19vsBSDE->SetMarkerStyle(20);
  // hCAL19vsBSDPi->SetMarkerStyle(20);
  hCAL19vsBSDE->Draw();
  hCAL19vsBSDPi->Draw("sames");

  TCanvas * cv8 = new TCanvas("cv8");
  hCAL19vsBSDE2->Draw("COLZ");

  gStyle->SetOptStat(1);
  TCanvas * cv9 = new TCanvas("cv9");
  hCAL19vsBSDPi2->Draw("COLZ");
}

// Use this one to try and do rej power vs electron acceptance on BSD signal
void Efrac_by_electron_energy::Loop5(Double_t myEnergyE, Double_t myEnergyP, Int_t CalLayer)
{
  // The means and sigmas by energy index 0 for 75, etc, up to 4 for 175
  Int_t index = 0;
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
 
  Double_t myMean = 0;
  Double_t mySigma = 0;

  TH2F * hCAL19vsBSDE = new TH2F("hCAL19vsBSDE","",200,-100,1000,200,-100,8000);
  TH2F * hCAL19vsBSDPi = new TH2F("hCAL19vsBSDPi","",200,-100,1000,200,-100,8000);
  TH2F * hCAL19vsBSDE2 = new TH2F("hCAL19vsBSDE2","",100,-100,1000,100,-100,8000);
  // TH2F * hCAL19vsBSDPi2 = new TH2F("hCAL19vsBSDPi2","",100,-100,1000,100,-100,8000);
  TH2F * hCAL19vsBSDPi2 = new TH2F("hCAL19vsBSDPi2","",50,-5,80,50,-100,8000);

  if( myEnergyE == 75 )
    {
      myMean = means[0];
      mySigma = sigmas[0];
    }
  else if( myEnergyE == 100 )
    {
      myMean = means[1];
      mySigma = sigmas[1];
    }
   else if( myEnergyE == 125 )
    {
      myMean = means[2];
      mySigma = sigmas[2];
    }
  else if( myEnergyE == 150 )
    {
      myMean = means[3];
      mySigma = sigmas[3];
    }
  else if( myEnergyE == 175 )
    {
      myMean = means[4];
      mySigma = sigmas[4];
    }
  else
    {
      cout << "Wrong energy choice! You have " << myEnergyE  << endl;
      return;
    }

  if (fChain == 0) return;
  
  TH1F * hPi = new TH1F("hPi","",500,-0.2,3);
  TH1F * hE = new TH1F("hE","",500,-0.2,3);

  TH1F * hPiCal = new TH1F("hPiCal","",200,-1000,30000);
  TH1F * hECal = new TH1F("hECal","",200,-1000,30000);

  TH1F * hPiBSD = new TH1F("hPiBSD","",200,-100,8000);
  TH1F * hEBSD = new TH1F("hEBSD","",200,-100,8000);

  TH1F * hPiCAL_19 = new TH1F("hPiCAL_19","",200,-100,1000);
  TH1F * hECAL_19 = new TH1F("hECAL_19","",200,-100,1000);
    
  Int_t NE75 = 0;
  Int_t NE100 = 0;
  Int_t NE125 = 0;
  Int_t NE150 = 0;
  Int_t NE175 = 0;

  Int_t NPi250 = 0;
  Int_t NPi300 = 0;
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
      // if( ptype == 2 && penergy == 250 && NPi250 > 6617) continue; 
      // if( ptype == 2 && penergy == 300 && NPi300 > 6617) continue; 
      // if( ptype == 2 && penergy == 350 && NPi350 > 6617) continue; 
      
      // Don't need these, since we aren't mixing electron energies
      // if( ptype == 1 && penergy == 75 && NE75 > 2442 ) continue;
      // if( ptype == 1 && penergy == 100 && NE100 > 2442 ) continue;
      // if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      // if( ptype == 1 && penergy == 150 && NE150 > 2442 ) continue;
      // if( ptype == 1 && penergy == 175 && NE175 > 2442 ) continue;
      
      if( ptype == 2 && penergy == 350 && NPi350 >= 6000) continue; 
      if( ptype == 1 && penergy == 125 && NE125 >= 6000 ) continue;
      
      
      // if( ptype == 2 && penergy == 350 && NPi350 > 4134) continue; 
      
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == myEnergyE )
	{
	  hECal->Fill(CALSum);
	  hCAL19vsBSDE2->Fill(CALProfile[CalLayer],myBSDLatePE);
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    {
	      hE->Fill(myBSDLatePE/CALSum);
	      hEBSD->Fill(myBSDLatePE);
	      hECAL_19->Fill(CALProfile[CalLayer]);
	      hCAL19vsBSDE->Fill(CALProfile[CalLayer],myBSDLatePE);
	      
	      if( ptype == 1 && penergy == 125 ) NE125++;
	    }
	}
      
      // Fill Efrac for pions
      if( ptype == 2 && penergy == myEnergyP )
	{
	  hPiCal->Fill(CALSum);
	  // hCAL19vsBSDPi2->Fill(CALProfile[CalLayer],myBSDLatePE);
	  // plot in MeV
	  hCAL19vsBSDPi2->Fill(1.34*CALProfile[CalLayer]/9,myBSDLatePE);
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    {
	      hPi->Fill(myBSDLatePE/CALSum);
	      hPiBSD->Fill(myBSDLatePE);
	      hPiCAL_19->Fill(CALProfile[CalLayer]);
	      hCAL19vsBSDPi->Fill(CALProfile[CalLayer],myBSDLatePE);
	    
	      if( ptype == 2 && penergy == 350 ) NPi350++;
	      	    
	    }
	}
    }

  TCanvas * cv1 = new TCanvas("cv1");
  
  // gStyle->SetOptStat(0);

  hPiCal->SetLineColor(kBlue);
  hECal->SetLineColor(kRed);
  hPiCal->Scale(1./hPiCal->GetEntries());
  hECal->Scale(1./hECal->GetEntries());
  hECal->GetXaxis()->SetTitle("Calorimeter Signal (arbitrary units)");
  hECal->GetXaxis()->CenterTitle();
  hECal->Draw();
  hPiCal->Draw("sames");
  
  TF1 * f1 = new TF1("f1","gaus",myMean - 2*mySigma, myMean + 2*mySigma);
  hECal->Fit("f1","R","same");

  cv1->Update();

  TLine * l1 = new TLine(myMean - 2*mySigma, 0, myMean - 2*mySigma, gPad->GetUymax());
  TLine * l2 = new TLine(myMean + 2*mySigma, 0, myMean + 2*mySigma, gPad->GetUymax());

  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  
  l1->Draw("same");
  l2->Draw("same");
  
  TCanvas * cv2 = new TCanvas("cv2");
  hEBSD->SetLineColor(kRed);
  hPiBSD->SetLineColor(kBlue);
  hEBSD->Scale(1./hEBSD->GetEntries());
  hPiBSD->Scale(1./hPiBSD->GetEntries());
  hEBSD->GetXaxis()->SetTitle("BSD Signal (arbitrary units)");
  hEBSD->GetXaxis()->CenterTitle();
  hEBSD->Draw("");
  hPiBSD->Draw("sames");
  
  TCanvas * cv3 = new TCanvas("cv3");
  
  hE->SetLineColor(kRed);
  hE->SetTitle("");
  hE->GetXaxis()->SetTitle("E_{frac} (\"detector\" units)");
  hE->GetXaxis()->CenterTitle();
  hE->Scale(1./hE->GetEntries());
  hE->Draw();
  hPi->Scale(1./hPi->GetEntries());
  hPi->Draw("sames");

  Double_t reg_power[1000];
  Double_t e_accept[1000];
  Double_t pi_accept[1000];
  Double_t endbin[1000];
  Int_t Npts = 0;
  Int_t Nskip = 0;
  Double_t efrac[1000];

  Int_t start = hEBSD->FindBin(0);
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
  Double_t Ne_entries = hEBSD->GetEntries();
  Double_t Npi_entries = hPiBSD->GetEntries();
  Double_t Efrac_temp = 0;

  // This analysis is actually in terms of BSD signal, not Efrac

  for( int i = 0; i < nbins; i++ )
    {
      cout << "Efrac = " << hEBSD->GetBinLowEdge(start+i) << endl;
      Efrac_temp = hEBSD->GetBinLowEdge(start+i);
      
      if( hPiBSD->Integral(start,start+i) > 10./Ne_entries && 
	  hEBSD->Integral(start,start+i) > 10./Npi_entries )
	{
	  num_err = 0;
	  den_err = 0;
	  Ne = 0;
	  Npi = 0;
	  num = 0;
	  den = 0;
	  Npi_total = 0;
	  Ne_total = 0;
	  
	  Ne = hEBSD->Integral(start,start+i);
	  Ne_total = hEBSD->Integral();
	  Npi = hPiBSD->Integral(start,start+i);
	  Npi_total = hPiBSD->Integral();
	  
	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);

	  reg_power[i-Nskip] = num/den;
	  e_accept[i-Nskip] = num;
	  pi_accept[i-Nskip] = den;
	  endbin[i-Nskip] = start+i;
	  efrac[i-Nskip] = Efrac_temp;

	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);
	  num_err = (num) * sqrt(1/(Ne*Ne_entries) + 1/(Ne_total*Ne_entries));
	  den_err = (den) * sqrt(1/(Npi*Npi_entries) + 1/(Npi_total*Npi_entries));
	  y_err[i-Nskip] = (num/den) * sqrt( (num_err*num_err)/(num*num) + (den_err*den_err)/(den*den) );
	  x_err[i-Nskip] = num_err;
	
	  cout << i-Nskip 
	       << " " << e_accept[i-Nskip]
	       << " " << pi_accept[i-Nskip]
	       << " " << reg_power[i-Nskip] 
	       << " " << hPiBSD->Integral(start,start+i)
	       << " " << hEBSD->Integral(start,start+i)
	       << " " << num_err
	       << " " << den_err
	       << endl;
	  Npts++;
	}
      else
	Nskip++;
    }
  
  /*  
  // For quick and easy copy/paste graphing
  cout << "rej[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << reg_power[i] << ", " << endl;
  cout << endl;
  
  // For quick and easy copy/paste graphing
  cout << "e_accept[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << e_accept[i] << ", " << endl;
  cout << endl;
  
  cout << "x_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << x_err[i] << endl;
  cout << endl;

  cout << "y_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << y_err[i] << endl;
  cout << endl;
  */
  TCanvas * cv4 = new TCanvas("cv4");
  TGraphErrors * gr = new TGraphErrors(Npts,e_accept,reg_power,x_err,y_err);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Electron Acceptance");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("Rejection Power");
  gr->GetYaxis()->CenterTitle();
  gr->SetMarkerStyle(20);
  gr->Draw("APL");

  TCanvas * cv5 = new TCanvas("cv5");
  TGraphErrors * gr5 = new TGraphErrors(10,efrac,e_accept,0,x_err);
  gr5->SetTitle("");
  gr5->SetMarkerStyle(20);
  gr5->GetXaxis()->SetTitle("E_{frac} (arbitrary units)");
  gr5->GetXaxis()->CenterTitle();
  gr5->GetYaxis()->SetTitle("Electron Acceptance");
  gr5->GetYaxis()->CenterTitle();
  gr5->Draw("APL");

  // This is looking at how much rejection power the bottom layer of the calorimeter
  // has, since it seems to be correlated to the BSD signal...
  TCanvas * cv6 = new TCanvas("cv6");
  hECAL_19->SetLineColor(kRed);
  hPiCAL_19->SetLineColor(kBlue);
  hECAL_19->Scale(1./hECAL_19->GetEntries());
  hPiCAL_19->Scale(1./hPiCAL_19->GetEntries());

  hECAL_19->Draw();
  hPiCAL_19->Draw("sames");
  
  TCanvas * cv7 = new TCanvas("cv7");
  hCAL19vsBSDE->SetMarkerColor(kRed);
  hCAL19vsBSDPi->SetMarkerColor(kBlue);
  // hCAL19vsBSDE->SetMarkerStyle(20);
  // hCAL19vsBSDPi->SetMarkerStyle(20);
  hCAL19vsBSDE->GetXaxis()->SetTitle("CERN Last Calorimeter Layer");
  hCAL19vsBSDE->GetXaxis()->CenterTitle("");
  hCAL19vsBSDE->GetYaxis()->SetTitle("CERN BSD Late Photoelectron Signal Summed Over Ten PMTs");
  hCAL19vsBSDE->GetYaxis()->SetTitleOffset(1.4);
  hCAL19vsBSDE->GetYaxis()->CenterTitle("");
  
  gStyle->SetOptStat(0);
  TH2F * haxes7 = new TH2F("hCAL19vsBSDE2","",100,-100,1000,100,-100,8000);
  haxes7->GetXaxis()->SetTitle("CERN Last Calorimeter Layer");
  haxes7->GetXaxis()->CenterTitle("");
  haxes7->GetYaxis()->SetTitle("CERN BSD Late Photoelectron Signal Summed Over Ten PMTs");
  haxes7->GetYaxis()->SetTitleOffset(1.4);
  haxes7->GetYaxis()->CenterTitle("");
  
  haxes7->Draw();

  gStyle->SetOptStat(1);
  hCAL19vsBSDE->Draw("sames");
  hCAL19vsBSDPi->Draw("sames");

  TCanvas * cv8 = new TCanvas("cv8");
  hCAL19vsBSDE2->Draw("COLZ");

  gStyle->SetOptStat(1);
  TCanvas * cv9 = new TCanvas("cv9");
  hCAL19vsBSDPi2->Draw("COLZ");
}


// myBSDLatePE and cal layer, but this time total cal layer, not 5 ribbon sum
void Efrac_by_electron_energy::Loop6(Double_t myEnergyE, Double_t myEnergyP, Int_t CalLayer)
{
  // The means and sigmas by energy index 0 for 75, etc, up to 4 for 175
  Int_t index = 0;
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
 
  Double_t myMean = 0;
  Double_t mySigma = 0;

  TH2F * hCAL19vsBSDE = new TH2F("hCAL19vsBSDE","",200,-100,1000,200,-100,8000);
  TH2F * hCAL19vsBSDPi = new TH2F("hCAL19vsBSDPi","",200,-100,1000,200,-100,8000);
  TH2F * hCAL19vsBSDE2 = new TH2F("hCAL19vsBSDE2","",100,-100,1000,100,-100,8000);
  TH2F * hCAL19vsBSDPi2 = new TH2F("hCAL19vsBSDPi2","",100,-100,1000,100,-100,8000);
  // TH2F * hCAL19vsBSDPi2 = new TH2F("hCAL19vsBSDPi2","",50,-5,80,50,-100,8000);

  if( myEnergyE == 75 )
    {
      myMean = means[0];
      mySigma = sigmas[0];
    }
  else if( myEnergyE == 100 )
    {
      myMean = means[1];
      mySigma = sigmas[1];
    }
   else if( myEnergyE == 125 )
    {
      myMean = means[2];
      mySigma = sigmas[2];
    }
  else if( myEnergyE == 150 )
    {
      myMean = means[3];
      mySigma = sigmas[3];
    }
  else if( myEnergyE == 175 )
    {
      myMean = means[4];
      mySigma = sigmas[4];
    }
  else
    {
      cout << "Wrong energy choice! You have " << myEnergyE  << endl;
      return;
    }

  if (fChain == 0) return;
  
  TH1F * hPi = new TH1F("hPi","",500,-0.2,3);
  TH1F * hE = new TH1F("hE","",500,-0.2,3);

  TH1F * hPiCal = new TH1F("hPiCal","",200,-1000,30000);
  TH1F * hECal = new TH1F("hECal","",200,-1000,30000);

  TH1F * hPiBSD = new TH1F("hPiBSD","",200,-100,8000);
  TH1F * hEBSD = new TH1F("hEBSD","",200,-100,8000);

  TH1F * hPiCAL_19 = new TH1F("hPiCAL_19","",200,-100,1000);
  TH1F * hECAL_19 = new TH1F("hECAL_19","",200,-100,1000);
    
  Int_t NE75 = 0;
  Int_t NE100 = 0;
  Int_t NE125 = 0;
  Int_t NE150 = 0;
  Int_t NE175 = 0;

  Int_t NPi250 = 0;
  Int_t NPi300 = 0;
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
      // if( ptype == 2 && penergy == 250 && NPi250 > 6617) continue; 
      // if( ptype == 2 && penergy == 300 && NPi300 > 6617) continue; 
      // if( ptype == 2 && penergy == 350 && NPi350 > 6617) continue; 
      
      // Don't need these, since we aren't mixing electron energies
      // if( ptype == 1 && penergy == 75 && NE75 > 2442 ) continue;
      // if( ptype == 1 && penergy == 100 && NE100 > 2442 ) continue;
      // if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      // if( ptype == 1 && penergy == 150 && NE150 > 2442 ) continue;
      // if( ptype == 1 && penergy == 175 && NE175 > 2442 ) continue;

      if( ptype == 1 && penergy == 125 && NE125 >= 6000 ) continue;
      if( ptype == 2 && penergy == 350 && NPi350 >= 6000) continue; 
      
      Double_t sum = 0;
      Int_t iXY = 0;
      if( CalLayer % 2 == 1 )
	iXY = 1;
      Int_t layers_to_10 = (Int_t)CalLayer/2;
      for( int i = 0; i < 50; i++ )
	if( cal[iXY][layers_to_10][i] > 0)
	  sum += cal[iXY][layers_to_10][i];
      
      // if( ptype == 2 && penergy == 350 && NPi350 > 4134) continue; 
      
      // See if we're looking at the right energy
      if( ptype == 1 && penergy == myEnergyE )
	{
	  hECal->Fill(CALSum);
	  // hCAL19vsBSDE2->Fill(CALProfile[CalLayer],myBSDLatePE);
	  hCAL19vsBSDE2->Fill(sum,myBSDLatePE);
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    {
	      hE->Fill(myBSDLatePE/CALSum);
	      hEBSD->Fill(myBSDLatePE);
	      // hECAL_19->Fill(CALProfile[CalLayer]);
	      hECAL_19->Fill(sum);
	      // hCAL19vsBSDE->Fill(CALProfile[CalLayer],myBSDLatePE);
	      hCAL19vsBSDE->Fill(sum,myBSDLatePE);
	      NE125++;
	    }
	}
      
      // Fill Efrac for pions
      if( ptype == 2 && penergy == myEnergyP )
	{
	  hPiCal->Fill(CALSum);
	  // hCAL19vsBSDPi2->Fill(CALProfile[CalLayer],myBSDLatePE);
	  // plot in MeV
	  // hCAL19vsBSDPi2->Fill(1.34*CALProfile[CalLayer]/9,myBSDLatePE);
	  hCAL19vsBSDPi2->Fill(sum,myBSDLatePE);
	  if( (CALSum > (myMean - 2*mySigma)) && (CALSum < (myMean + 2*mySigma)) )
	    {
	      hPi->Fill(myBSDLatePE/CALSum);
	      hPiBSD->Fill(myBSDLatePE);
	      // hPiCAL_19->Fill(CALProfile[CalLayer]);
	      hPiCAL_19->Fill(sum);
	      hCAL19vsBSDPi->Fill(sum,myBSDLatePE);
	      NPi350++;
	    }
	}


      // Increment counts, since everything has to be equally weighted
      // if( ptype == 2 && penergy == 250 ) NPi250++;
      // if( ptype == 2 && penergy == 300 ) NPi300++;
      // if( ptype == 2 && penergy == 350 ) NPi350++;
      
      // if( ptype == 1 && penergy == 75 ) NE75++;
      // if( ptype == 1 && penergy == 100 ) NE100++;
      // if( ptype == 1 && penergy == 125 ) NE125++;
      // if( ptype == 1 && penergy == 150 ) NE150++;
      // if( ptype == 1 && penergy == 175 ) NE175++;
    }

  TCanvas * cv1 = new TCanvas("cv1");
  
  // gStyle->SetOptStat(0);

  hPiCal->SetLineColor(kBlue);
  hECal->SetLineColor(kRed);
  hPiCal->Scale(1./hPiCal->GetEntries());
  hECal->Scale(1./hECal->GetEntries());
  hECal->GetXaxis()->SetTitle("Calorimeter Signal (arbitrary units)");
  hECal->GetXaxis()->CenterTitle();
  hECal->Draw();
  hPiCal->Draw("sames");
  
  TF1 * f1 = new TF1("f1","gaus",myMean - 2*mySigma, myMean + 2*mySigma);
  hECal->Fit("f1","R","same");

  cv1->Update();

  TLine * l1 = new TLine(myMean - 2*mySigma, 0, myMean - 2*mySigma, gPad->GetUymax());
  TLine * l2 = new TLine(myMean + 2*mySigma, 0, myMean + 2*mySigma, gPad->GetUymax());

  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  
  l1->Draw("same");
  l2->Draw("same");
  
  TCanvas * cv2 = new TCanvas("cv2");
  hEBSD->SetLineColor(kRed);
  hPiBSD->SetLineColor(kBlue);
  hEBSD->Scale(1./hEBSD->GetEntries());
  hPiBSD->Scale(1./hPiBSD->GetEntries());
  hEBSD->GetXaxis()->SetTitle("BSD Signal (arbitrary units)");
  hEBSD->GetXaxis()->CenterTitle();
  hEBSD->Draw("");
  hPiBSD->Draw("sames");
  
  TCanvas * cv3 = new TCanvas("cv3");
  
  hE->SetLineColor(kRed);
  hE->SetTitle("");
  hE->GetXaxis()->SetTitle("E_{frac} (\"detector\" units)");
  hE->GetXaxis()->CenterTitle();
  hE->Scale(1./hE->GetEntries());
  hE->Draw();
  hPi->Scale(1./hPi->GetEntries());
  hPi->Draw("sames");

  Double_t reg_power[1000];
  Double_t e_accept[1000];
  Double_t pi_accept[1000];
  Double_t endbin[1000];
  Int_t Npts = 0;
  Int_t Nskip = 0;
  Double_t efrac[1000];

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
  Double_t Ne_entries = hE->GetEntries();
  Double_t Npi_entries = hPi->GetEntries();
  Double_t Efrac_temp = 0;

  for( int i = 0; i < nbins; i++ )
    {
      cout << "Efrac = " << hE->GetBinLowEdge(start+i) << endl;
      Efrac_temp = hE->GetBinLowEdge(start+i);
      
      if( hPi->Integral(start,start+i) > 10./Ne_entries && 
	  hE->Integral(start,start+i) > 10./Npi_entries )
	{
	  num_err = 0;
	  den_err = 0;
	  Ne = 0;
	  Npi = 0;
	  num = 0;
	  den = 0;
	  Npi_total = 0;
	  Ne_total = 0;
	  
	  Ne = hE->Integral(start,start+i);
	  Ne_total = hE->Integral();
	  Npi = hPi->Integral(start,start+i);
	  Npi_total = hPi->Integral();
	  
	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);

	  reg_power[i-Nskip] = num/den;
	  e_accept[i-Nskip] = num;
	  pi_accept[i-Nskip] = den;
	  endbin[i-Nskip] = start+i;
	  efrac[i-Nskip] = Efrac_temp;

	  num = (Ne/Ne_total);
	  den = (Npi/Npi_total);
	  num_err = (num) * sqrt(1/(Ne*Ne_entries) + 1/(Ne_total*Ne_entries));
	  den_err = (den) * sqrt(1/(Npi*Npi_entries) + 1/(Npi_total*Npi_entries));
	  y_err[i-Nskip] = (num/den) * sqrt( (num_err*num_err)/(num*num) + (den_err*den_err)/(den*den) );
	  x_err[i-Nskip] = num_err;
	
	  cout << i-Nskip 
	       << " " << e_accept[i-Nskip]
	       << " " << pi_accept[i-Nskip]
	       << " " << reg_power[i-Nskip] 
	       << " " << hPi->Integral(start,start+i)
	       << " " << hE->Integral(start,start+i)
	       << " " << num_err
	       << " " << den_err
	       << endl;
	  Npts++;
	}
      else
	Nskip++;
    }
  
  /*  
  // For quick and easy copy/paste graphing
  cout << "rej[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << reg_power[i] << ", " << endl;
  cout << endl;
  
  // For quick and easy copy/paste graphing
  cout << "e_accept[" << Npts << "] = {";
  for(int i = 0; i < Npts; i++ )
    cout << e_accept[i] << ", " << endl;
  cout << endl;
  
  cout << "x_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << x_err[i] << endl;
  cout << endl;

  cout << "y_err[" << Npts << "] = {";
  for( int i = 0; i < Npts; i++ )
    cout << y_err[i] << endl;
  cout << endl;
  */
  TCanvas * cv4 = new TCanvas("cv4");
  TGraphErrors * gr = new TGraphErrors(Npts,e_accept,reg_power,x_err,y_err);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Electron Acceptance");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("Rejection Power");
  gr->GetYaxis()->CenterTitle();
  gr->SetMarkerStyle(20);
  gr->Draw("APL");

  TCanvas * cv5 = new TCanvas("cv5");
  TGraphErrors * gr5 = new TGraphErrors(10,efrac,e_accept,0,x_err);
  gr5->SetTitle("");
  gr5->SetMarkerStyle(20);
  gr5->GetXaxis()->SetTitle("E_{frac} (arbitrary units)");
  gr5->GetXaxis()->CenterTitle();
  gr5->GetYaxis()->SetTitle("Electron Acceptance");
  gr5->GetYaxis()->CenterTitle();
  gr5->Draw("APL");

  // This is looking at how much rejection power the bottom layer of the calorimeter
  // has, since it seems to be correlated to the BSD signal...
  TCanvas * cv6 = new TCanvas("cv6");
  hECAL_19->SetLineColor(kRed);
  hPiCAL_19->SetLineColor(kBlue);
  hECAL_19->Scale(1./hECAL_19->GetEntries());
  hPiCAL_19->Scale(1./hPiCAL_19->GetEntries());

  hECAL_19->Draw();
  hPiCAL_19->Draw("sames");
  
  TCanvas * cv7 = new TCanvas("cv7");
  hCAL19vsBSDE->SetMarkerColor(kRed);
  hCAL19vsBSDPi->SetMarkerColor(kBlue);
  // hCAL19vsBSDE->SetMarkerStyle(20);
  // hCAL19vsBSDPi->SetMarkerStyle(20);
  hCAL19vsBSDE->Draw();
  hCAL19vsBSDPi->Draw("sames");

  TCanvas * cv8 = new TCanvas("cv8");
  hCAL19vsBSDE2->Draw("COLZ");

  gStyle->SetOptStat(1);
  TCanvas * cv9 = new TCanvas("cv9");
  hCAL19vsBSDPi2->Draw("COLZ");
}

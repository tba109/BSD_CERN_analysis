#define Efrac_electron_sigma_cut_cxx
#include "Efrac_electron_sigma_cut.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>

using namespace std;

void Efrac_electron_sigma_cut::Loop1()
{
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
  
  // Let's draw Efrac for this range
  if (fChain == 0) return;
  

  Int_t index = 0;
  TH1F * hE_2sigma = new TH1F("hE_2sigma","",500,-0.2,3);
  TH1F * hE_3sigma = new TH1F("hE_3sigma","",500,-0.2,3);

  Int_t NE75 = 0;
  Int_t NE100 = 0;
  Int_t NE125 = 0;
  Int_t NE150 = 0;
  Int_t NE175 = 0;
    
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
      if( ptype == 1 && penergy == 75 && NE75 > 2442 ) continue;
      if( ptype == 1 && penergy == 100 && NE100 > 2442 ) continue;
      if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      if( ptype == 1 && penergy == 150 && NE150 > 2442 ) continue;
      if( ptype == 1 && penergy == 175 && NE175 > 2442 ) continue;

      // Need to find the index for mean and sigma
      if( ptype == 1 && penergy == 75 )  index = 0;
      if( ptype == 1 && penergy == 100 ) index = 1;
      if( ptype == 1 && penergy == 125 ) index = 2;
      if( ptype == 1 && penergy == 150 ) index = 3;
      if( ptype == 1 && penergy == 175 ) index = 4;


      // Fill Efrac for 2 sigma electrons
      if( ptype == 1 && (penergy >=75 || penergy <= 175) )
	if( (CALSum < (means[index] - 2*sigmas[index])) || (CALSum > (means[index] + 2*sigmas[index])) )
	   hE_2sigma->Fill(myBSDLatePE/CALSum);
      
      // Fill Efrac for 3 sigma electrons
      if( ptype == 1 && (penergy >=75 || penergy <= 175) )
	if( (CALSum < (means[index] - 3*sigmas[index])) || (CALSum > (means[index] + 3*sigmas[index])) )
	   hE_3sigma->Fill(myBSDLatePE/CALSum);
      
      // Increment counts, since everything has to be equally weighted
      if( ptype == 1 && penergy == 75 ) NE75++;
      if( ptype == 1 && penergy == 100 ) NE100++;
      if( ptype == 1 && penergy == 125 ) NE125++;
      if( ptype == 1 && penergy == 150 ) NE150++;
      if( ptype == 1 && penergy == 175 ) NE175++;
    }
  
  TCanvas * cv2sigma = new TCanvas();
  hE_2sigma->SetLineColor(kRed);
  hE_2sigma->SetTitle("Efrac for 75-175 GeV Electrons outside 2#sigma (Red)");
  hE_2sigma->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_2sigma->GetXaxis()->CenterTitle();
  hE_2sigma->Scale(1./hE_2sigma->GetEntries());
  hE_2sigma->Draw();

  TCanvas * cv3sigma = new TCanvas();
  hE_3sigma->SetLineColor(kRed);
  hE_3sigma->SetTitle("Efrac for 75-175 GeV Electrons outside 3#sigma (Red)");
  hE_3sigma->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_3sigma->GetXaxis()->CenterTitle();
  hE_3sigma->Scale(1./hE_3sigma->GetEntries());
  hE_3sigma->Draw();


}

// This one does on a per energy run basis
void Efrac_electron_sigma_cut::Loop2()
{
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
  
  // Let's draw Efrac for this range
  if (fChain == 0) return;

  Int_t index = 0;
  TH1F * hE_2sigma75 = new TH1F("hE_2sigma75","",500,-0.2,3);
  TH1F * hE_2sigma100 = new TH1F("hE_2sigma100","",500,-0.2,3);
  TH1F * hE_2sigma125 = new TH1F("hE_2sigma125","",500,-0.2,3);
  TH1F * hE_2sigma150 = new TH1F("hE_2sigma150","",500,-0.2,3);
  TH1F * hE_2sigma175 = new TH1F("hE_2sigma175","",500,-0.2,3);
  TH1F * hE_3sigma75 = new TH1F("hE_3sigma75","",500,-0.2,3);
  TH1F * hE_3sigma100 = new TH1F("hE_3sigma100","",500,-0.2,3);
  TH1F * hE_3sigma125 = new TH1F("hE_3sigma125","",500,-0.2,3);
  TH1F * hE_3sigma150 = new TH1F("hE_3sigma150","",500,-0.2,3);
  TH1F * hE_3sigma175 = new TH1F("hE_3sigma175","",500,-0.2,3);
 
  Int_t NE75_2sigma = 0;
  Int_t NE100_2sigma = 0;
  Int_t NE125_2sigma = 0;
  Int_t NE150_2sigma = 0;
  Int_t NE175_2sigma = 0;
   
  Int_t NE75_3sigma = 0;
  Int_t NE100_3sigma = 0;
  Int_t NE125_3sigma = 0;
  Int_t NE150_3sigma = 0;
  Int_t NE175_3sigma = 0;
    
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
            
      // Need to find the index for mean and sigma
      if( ptype == 1 && penergy == 75 )  index = 0;
      if( ptype == 1 && penergy == 100 ) index = 1;
      if( ptype == 1 && penergy == 125 ) index = 2;
      if( ptype == 1 && penergy == 150 ) index = 3;
      if( ptype == 1 && penergy == 175 ) index = 4;


      // Fill Efrac for 2 sigma electrons
      if( ptype == 1 && (penergy >=75 || penergy <= 175) )
	if( (CALSum < (means[index] - 2*sigmas[index])) || (CALSum > (means[index] + 2*sigmas[index])) )
	  {
	    // Increment counts, since everything has to be equally weighted
	    if( ptype == 1 && penergy == 75 ) 
	      {
		hE_2sigma75->Fill(myBSDLatePE/CALSum);
		NE75_2sigma++;
	      }
	    if( ptype == 1 && penergy == 100 ) 
	      {
		hE_2sigma100->Fill(myBSDLatePE/CALSum);
		NE100_2sigma++;
	      }
	    if( ptype == 1 && penergy == 125 )
	      {
		hE_2sigma125->Fill(myBSDLatePE/CALSum);
		NE125_2sigma++;
	      }
	    if( ptype == 1 && penergy == 150 ) 
	      {
		hE_2sigma150->Fill(myBSDLatePE/CALSum);
		NE150_2sigma++;
	      }
	    if( ptype == 1 && penergy == 175 ) 
	      {
		hE_2sigma175->Fill(myBSDLatePE/CALSum);
		NE175_2sigma++;
	      }
	  }
      
      // Fill Efrac for 3 sigma electrons
      if( ptype == 1 && (penergy >=75 || penergy <= 175) )
	if( (CALSum < (means[index] - 3*sigmas[index])) || (CALSum > (means[index] + 3*sigmas[index])) )
	  {
	    // Increment counts, since everything has to be equally weighted
	    if( ptype == 1 && penergy == 75 ) 
	      {
		 hE_3sigma75->Fill(myBSDLatePE/CALSum);
		 NE75_3sigma++;
	      }
	    if( ptype == 1 && penergy == 100 ) 
	      {
		hE_3sigma100->Fill(myBSDLatePE/CALSum);
		NE100_3sigma++;
	      }
	    if( ptype == 1 && penergy == 125 ) 
	      {
		hE_3sigma125->Fill(myBSDLatePE/CALSum);
		NE125_3sigma++;
	      }
	    if( ptype == 1 && penergy == 150 ) 
	      {
		hE_3sigma150->Fill(myBSDLatePE/CALSum);
		NE150_3sigma++;
	      }
	    if( ptype == 1 && penergy == 175 ) 
	      {
		hE_3sigma175->Fill(myBSDLatePE/CALSum);
		NE175_3sigma++;
	      }
	  }
    }
  
  // Do 2 sigma plots
  TCanvas * cv2sigma = new TCanvas();
  cv2sigma->Divide(3,2);
  
  // 75 GeV electrons
  cv2sigma->cd(2);
  hE_2sigma75->SetLineColor(kRed);
  hE_2sigma75->SetTitle("Efrac for 75 GeV Electrons outside 2#sigma (Red)");
  hE_2sigma75->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_2sigma75->GetXaxis()->CenterTitle();
  hE_2sigma75->Scale(1./hE_2sigma75->GetEntries());
  hE_2sigma75->Draw();

  // 100 GeV electrons
  cv2sigma->cd(3);
  hE_2sigma100->SetLineColor(kRed);
  hE_2sigma100->SetTitle("Efrac for 100 GeV Electrons outside 2#sigma (Red)");
  hE_2sigma100->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_2sigma100->GetXaxis()->CenterTitle();
  hE_2sigma100->Scale(1./hE_2sigma100->GetEntries());
  hE_2sigma100->Draw();

  // 125 GeV electrons
  cv2sigma->cd(4);
  hE_2sigma125->SetLineColor(kRed);
  hE_2sigma125->SetTitle("Efrac for 125 GeV Electrons outside 2#sigma (Red)");
  hE_2sigma125->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_2sigma125->GetXaxis()->CenterTitle();
  hE_2sigma125->Scale(1./hE_2sigma125->GetEntries());
  hE_2sigma125->Draw();

  // 150 GeV electrons
  cv2sigma->cd(5);
  hE_2sigma150->SetLineColor(kRed);
  hE_2sigma150->SetTitle("Efrac for 150 GeV Electrons outside 2#sigma (Red)");
  hE_2sigma150->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_2sigma150->GetXaxis()->CenterTitle();
  hE_2sigma150->Scale(1./hE_2sigma150->GetEntries());
  hE_2sigma150->Draw();
  
  // 175 GeV electrons
  cv2sigma->cd(6);
  hE_2sigma175->SetLineColor(kRed);
  hE_2sigma175->SetTitle("Efrac for 175 GeV Electrons outside 2#sigma (Red)");
  hE_2sigma175->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_2sigma175->GetXaxis()->CenterTitle();
  hE_2sigma175->Scale(1./hE_2sigma175->GetEntries());
  hE_2sigma175->Draw();

  // Do 3 sigma plots
  TCanvas * cv3sigma = new TCanvas();
  cv3sigma->Divide(3,2);

  // 75 GeV electrons
  cv3sigma->cd(2);
  hE_3sigma75->SetLineColor(kRed);
  hE_3sigma75->SetTitle("Efrac for 75 GeV Electrons outside 3#sigma (Red)");
  hE_3sigma75->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_3sigma75->GetXaxis()->CenterTitle();
  hE_3sigma75->Scale(1./hE_3sigma75->GetEntries());
  hE_3sigma75->Draw();

  // 100 GeV electrons
  cv3sigma->cd(3);
  hE_3sigma100->SetLineColor(kRed);
  hE_3sigma100->SetTitle("Efrac for 100 GeV Electrons outside 3#sigma (Red)");
  hE_3sigma100->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_3sigma100->GetXaxis()->CenterTitle();
  hE_3sigma100->Scale(1./hE_3sigma100->GetEntries());
  hE_3sigma100->Draw();

  // 125 GeV electrons
  cv3sigma->cd(4);
  hE_3sigma125->SetLineColor(kRed);
  hE_3sigma125->SetTitle("Efrac for 125 GeV Electrons outside 3#sigma (Red)");
  hE_3sigma125->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_3sigma125->GetXaxis()->CenterTitle();
  hE_3sigma125->Scale(1./hE_3sigma125->GetEntries());
  hE_3sigma125->Draw();

  // 150 GeV electrons
  cv3sigma->cd(5);
  hE_3sigma150->SetLineColor(kRed);
  hE_3sigma150->SetTitle("Efrac for 150 GeV Electrons outside 3#sigma (Red)");
  hE_3sigma150->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_3sigma150->GetXaxis()->CenterTitle();
  hE_3sigma150->Scale(1./hE_3sigma150->GetEntries());
  hE_3sigma150->Draw();

  // 175 GeV electrons
  cv3sigma->cd(6);
  hE_3sigma175->SetLineColor(kRed);
  hE_3sigma175->SetTitle("Efrac for 175 GeV Electrons outside 3#sigma (Red)");
  hE_3sigma175->GetXaxis()->SetTitle("Efrac (\"detector\" units)");
  hE_3sigma175->GetXaxis()->CenterTitle();
  hE_3sigma175->Scale(1./hE_3sigma175->GetEntries());
  hE_3sigma175->Draw();

}

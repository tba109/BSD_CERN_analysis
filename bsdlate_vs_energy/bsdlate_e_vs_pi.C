#define bsdlate_e_vs_pi_cxx
#include "bsdlate_e_vs_pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void bsdlate_e_vs_pi::Loop()
{
  if (fChain == 0) return;
  
  TH1F * hPi = new TH1F("hPi","",200,-100,8000);
  TH1F * hE = new TH1F("hE","",200,-100,8000);

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
      
      if( ptype == 1 && penergy == 75 && NE75 > 2442 ) continue;
      if( ptype == 1 && penergy == 100 && NE100 > 2442 ) continue;
      if( ptype == 1 && penergy == 125 && NE125 > 2442 ) continue;
      if( ptype == 1 && penergy == 150 && NE150 > 2442 ) continue;
      if( ptype == 1 && penergy == 175 && NE175 > 2442 ) continue;
      
      // 250-350 GeV pions
      if( ptype == 2 && (penergy <=  350 || penergy >= 250) ) hPi->Fill(myBSDLatePE);
      
      // 75-175 GeV electrons
      if( ptype == 1 && (penergy <= 175 || penergy >= 75) ) hE->Fill(myBSDLatePE);

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
  hPi->SetLineColor(kBlue);
  hE->SetLineColor(kRed);

  if( hPi->GetEntries() > 0 )
    hPi->Scale(1./hPi->GetEntries());

  if( hE->GetEntries() > 0 )
    hE->Scale(1./hE->GetEntries());
  
  TH2F * haxes = new TH2F("haxes","",100,-100,8000,100,0,0.07);
  haxes->SetStats(kFALSE);
  haxes->GetXaxis()->SetTitle("myBSDLatePE for 75-175 GeV Electrons (Red) and 250-350 GeV Pions (Blue)");
  haxes->GetXaxis()->CenterTitle();
  haxes->Draw();
  
  hPi->Draw("sames");
  hE->Draw("sames");
}

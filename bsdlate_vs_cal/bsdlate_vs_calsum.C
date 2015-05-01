#define bsdlate_vs_calsum_cxx
#include "bsdlate_vs_calsum.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>

void bsdlate_vs_calsum::Loop()
{
    
  if (fChain == 0) return;
  
  Double_t aE[2][100000];
  Double_t aPi[2][100000];
  Int_t NE = 0;
  Int_t NPi = 0;

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
      if( ptype == 2 && NPi > 10000 ) continue;
      if( ptype == 2 && penergy == 250 && NPi250 > 3333) continue; 
      if( ptype == 2 && penergy == 300 && NPi300 > 3333) continue; 
      if( ptype == 2 && penergy == 350 && NPi350 > 3333) continue; 
      
      if( ptype == 1 && NE > 10000 ) continue;
      if( ptype == 1 && penergy == 75 && NE75 > 2000 ) continue;
      if( ptype == 1 && penergy == 100 && NE100 > 2000 ) continue;
      if( ptype == 1 && penergy == 125 && NE125 > 2000 ) continue;
      if( ptype == 1 && penergy == 150 && NE150 > 2000 ) continue;
      if( ptype == 1 && penergy == 175 && NE175 > 2000 ) continue;

      // Fill Efrac for electrons
      if( ptype == 1 && (penergy >=75 || penergy <= 175) )
	{
	  aE[0][NE] = CALSum;
	  aE[1][NE] = myBSDLatePE;
	}
      
      // Fill Efrac for pions
      if( ptype == 2 && (penergy >= 250 || penergy <= 350) )
	{
	  aPi[0][NPi] = CALSum;
	  aPi[1][NPi] = myBSDLatePE;
	}
      
      // Increment counts, since everything has to be equally weighted
      if( ptype == 2 ) NPi++;
      if( ptype == 2 && penergy == 250 ) NPi250++;
      if( ptype == 2 && penergy == 300 ) NPi300++;
      if( ptype == 2 && penergy == 350 ) NPi350++;
      
      if( ptype == 1 ) NE++;
      if( ptype == 1 && penergy == 75 ) NE75++;
      if( ptype == 1 && penergy == 100 ) NE100++;
      if( ptype == 1 && penergy == 125 ) NE125++;
      if( ptype == 1 && penergy == 150 ) NE150++;
      if( ptype == 1 && penergy == 175 ) NE175++;
    }
  
  TGraph * grE = new TGraph(NE,aE[0],aE[1]);
  TGraph * grPi = new TGraph(NPi,aPi[0],aPi[1]);
  grE->SetMarkerColor(kRed);
  grE->GetXaxis()->SetTitle("CALSum (\"detector\" units)");
  grE->GetXaxis()->CenterTitle();
  grE->GetYaxis()->SetTitle("myBSDLatePE (PE)");
  grE->GetYaxis()->CenterTitle();
  grE->SetTitle("");
  grE->Draw("AP");
  grPi->SetMarkerColor(kBlue);
  grPi->Draw("Psame");
}

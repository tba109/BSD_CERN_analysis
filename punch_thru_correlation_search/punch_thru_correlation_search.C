#define punch_thru_correlation_search_cxx
#include "punch_thru_correlation_search.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <iostream>

using namespace std;

void punch_thru_correlation_search::Loop(Int_t myType, Int_t myEnergy)
{
  TH2F * hPi2D[20];
  
  char name[100];
  for(int i = 0; i < 20; i++ )
    {
      sprintf(name,"hPi2D_%d",i);
      hPi2D[i] = new TH2F(name,"",100,-100,1000,100,-100,7000);
    }    

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
      
      Double_t sum0 = 0;
      Double_t sum1 = 0;
      if( ptype == myType && penergy == myEnergy )
	{
	  for( int i = 0; i < 10; i++ )
	    {
	      sum0 = 0;
	      sum1 = 0;
	      for( int j = 0; j < 50; j++ )
		{
		  if( cal[0][i][j] > 0 )
		    sum0 += cal[0][i][j];
		  if( cal[1][i][j] > 0 )
		    sum1 += cal[1][i][j];
		}
	      hPi2D[2*i]->Fill(sum0,myBSDLatePE);
	      hPi2D[2*i+1]->Fill(sum1,myBSDLatePE);
	      // for( int i = 0; i < 20; i++ )
	      // hPi2D[i]->Fill(CALProfile[i],myBSDLatePE);
	    }
	}
    }
  
  char go;
  for( int j = 0; j < 20; j++ )
    {
      hPi2D[j]->Draw("COLZ");
      gPad->Update();
      // cin >> go;
      usleep(250000);
    }
}

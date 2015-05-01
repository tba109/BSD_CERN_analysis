#define make_all_histos_cxx
#include "make_all_histos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>

void make_all_histos::Loop()
{
  if (fChain == 0) return;


  char name[100];
  Int_t pi_evals[3] = {250,300,350};
  Int_t e_evals[5] = {75,100,125,150,175};

  TH1F * hPi_CAL19[3];
  TH1F * hPi_CALSum[3];
  TH1F * hPi_myBSDLatePE[3];
  TH1F * hPi_Efrac[3];
  
  TH1F * hE_CAL19[5];
  TH1F * hE_CALSum[5];
  TH1F * hE_myBSDLatePE[5];
  TH1F * hE_Efrac[5];

  for( int i = 0; i < 3; i++ )
    {
      sprintf(name,"hPi_CAL19_%d",pi_evals[i]);
      hPi_CAL19[i] = new TH1F(name,"",200,-1000,30000);

      sprintf(name,"hPi_CALSum_%d",pi_evals[i]);
      hPi_CALSum[i] = new TH1F(name,"",200,-1000,30000);

      sprintf(name,"hPi_myBSDLatePE_%d",pi_evals[i]);
      hPi_myBSDLatePE[i] = new TH1F(name,"",200,-100,8000);
      
      sprintf(name,"hPi_Efrac_%d",pi_evals[i]);
      hPi_Efrac[i] = new TH1F(name,"",500,-0.2,3);
    }
  
  for( int i = 0; i < 5; i++ )
    {
      sprintf(name,"hE_CAL19_%d",e_evals[i]);
      hE_CAL19[i] = new TH1F(name,"",200,-1000,30000);

      sprintf(name,"hE_CALSum_%d",e_evals[i]);
      hE_CALSum[i] = new TH1F(name,"",200,-1000,30000);

      sprintf(name,"hE_myBSDLatePE_%d",e_evals[i]);
      hE_myBSDLatePE[i] = new TH1F(name,"",200,-100,8000);
      
      sprintf(name,"hE_Efrac_%d",e_evals[i]);
      hE_Efrac[i] = new TH1F(name,"",500,-0.2,3);
    }

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

      // The BSD Late Time PEs
      for( int i = 1; i <= 10; i++ )
	myBSDLatePE += myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19);
      
      Double_t sum = 0;
      for( int k = 0; k < 50; k++ )
	{
	  if( cal[1][9][k] > 0 )
	    sum += cal[1][9][k];
	}


      // Fill Efrac for electrons
      for( int i = 0; i < 5; i++ )
	{
	  if( ptype == 1 && penergy == e_evals[i] )
	    {
	      hE_Efrac[i]->Fill(myBSDLatePE/CALSum);
	      hE_CALSum[i]->Fill(CALSum);
	      hE_myBSDLatePE[i]->Fill(myBSDLatePE);
	      hE_CAL19[i]->Fill(sum);
	    }
	}

      for( int i = 0; i < 3; i++ )
	{
	  // Fill Efrac for pions
	  if( ptype == 2 && penergy == pi_evals[i] )
	    {
	      hPi_Efrac[i]->Fill(myBSDLatePE/CALSum);
	      hPi_CALSum[i]->Fill(CALSum);
	      hPi_myBSDLatePE[i]->Fill(myBSDLatePE);
	      hPi_CAL19[i]->Fill(sum);
	    }
	}
	 
    }

  printf("here\n");
  
  TFile f("histos_E_Pi.root","recreate");
  TCanvas * cv1 = new TCanvas("cv1");
  
  for( int i = 0; i < 5; i++ )
    {
      hE_CALSum[i]->Draw();
      hE_CALSum[i]->Write();
      gPad->Update();
      sleep(1);

      hE_CAL19[i]->Draw();
      hE_CAL19[i]->Write();
      gPad->Update();
      sleep(1);
      
      hE_myBSDLatePE[i]->Draw();
      hE_myBSDLatePE[i]->Write();
      gPad->Update();
      sleep(1);
      
      hE_Efrac[i]->Draw();
      hE_Efrac[i]->Write();
      gPad->Update();
      sleep(1);
    }
  
  for( int i = 0; i < 3; i++ )
    {
      hPi_CALSum[i]->Draw();
      hPi_CALSum[i]->Write();
      gPad->Update();
      sleep(1);

      hPi_CAL19[i]->Draw();
      hPi_CAL19[i]->Write();
      gPad->Update();
      sleep(1);
      
      hPi_myBSDLatePE[i]->Draw();
      hPi_myBSDLatePE[i]->Write();
      gPad->Update();
      sleep(1);
      
      hPi_Efrac[i]->Draw();
      hPi_Efrac[i]->Write();
      gPad->Update();
      sleep(1);
    }

  printf("here\n");
  f.Close();
  
}

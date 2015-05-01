#define bsdlate_pe_energies_cxx
#include "bsdlate_pe_energies.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

// These have all the standard cuts
void bsdlate_pe_energies::Loop1()
{
  if (fChain == 0) return;
  
  TH1F * hPi350[12];
  TH1F * hPi300[12];
  TH1F * hPi250[12];
  
  TH1F * hE175[12];
  TH1F * hE150[12];
  TH1F * hE125[12];
  TH1F * hE100[12];
  TH1F * hE75[12];
    
  char name[1000];
  
  for( int i = 0; i < 12; i++ )
    {
      sprintf(name,"hPi350_%d",i);
      hPi350[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hPi300_%d",i);
      hPi300[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hPi250_%d",i);
      hPi250[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE75_%d",i);
      hE75[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE100_%d",i);
      hE100[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE125_%d",i);
      hE125[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE150_%d",i);
      hE150[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE175_%d",i);
      hE175[i] = new TH1F(name,"",200,-100,1000);
    }

  Double_t gainA[12] = {4.23,4.23,3.92,3.09,5.81,5.03,3.17,5.30,4.31,3.90,4.04,3.68};
  Double_t myfCHA[12];
  Double_t myBSDLatePE[12];
  
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
      
      for( int i = 0; i < 12; i++ )
	{
	  myfCHA[i] = -1;
	  myBSDLatePE[i] = 0;
	}

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
	{
	  myBSDLatePE[i] = myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19);
	  
	  if( ptype == 2 && penergy == 350 ) hPi350[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 2 && penergy == 300 ) hPi300[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 2 && penergy == 250 ) hPi250[i]->Fill(myBSDLatePE[i]);
	  
	  if( ptype == 1 && penergy == 75 ) hE75[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 1 && penergy == 100 ) hE100[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 1 && penergy == 125 ) hE125[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 1 && penergy == 150 ) hE150[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 1 && penergy == 175 ) hE175[i]->Fill(myBSDLatePE[i]);
	}
    }
  
  TCanvas * cvPi350 = new TCanvas("cvPi350");
  cvPi350->Divide(5,2);
  
  TCanvas * cvPi300 = new TCanvas("cvPi300");
  cvPi300->Divide(5,2);
  
  TCanvas * cvPi250 = new TCanvas("cvPi250");
  cvPi250->Divide(5,2);
  
  TCanvas * cvE175 = new TCanvas("cvE175");
  cvE175->Divide(5,2);
  
  TCanvas * cvE150 = new TCanvas("cvE150");
  cvE150->Divide(5,2);
  
  TCanvas * cvE125 = new TCanvas("cvE125");
  cvE125->Divide(5,2);
  
  TCanvas * cvE100 = new TCanvas("cvE100");
  cvE100->Divide(5,2);
  
  TCanvas * cvE75 = new TCanvas("cvE75");
  cvE75->Divide(5,2);
  
  for( int i = 1; i <= 10; i++ )
    {
      cvPi350->cd(i);
      hPi350[i]->Draw();
      
      cvPi300->cd(i);
      hPi300[i]->Draw();
      
      cvPi250->cd(i);
      hPi250[i]->Draw();
      
      cvE75->cd(i);
      hE75[i]->Draw();
      
      cvE100->cd(i);
      hE100[i]->Draw();
      
      cvE125->cd(i);
      hE125[i]->Draw();
      
      cvE150->cd(i);
      hE150[i]->Draw();
      
      cvE175->cd(i);
      hE175[i]->Draw();
    }

  Double_t mean = 0;

  cout << "=======================================" << endl;
  cout << "75 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE75[i]->GetMean() << endl;
      mean += hE75[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "100 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE100[i]->GetMean() << endl;
      mean += hE100[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "125 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE125[i]->GetMean() << endl;
      mean += hE125[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "150 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE150[i]->GetMean() << endl;
      mean += hE150[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "175 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE175[i]->GetMean() << endl;
      mean += hE175[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "250 GeV Pions:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hPi250[i]->GetMean() << endl;
      mean += hPi250[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "300 GeV Pions:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hPi300[i]->GetMean() << endl;
      mean += hPi300[i]->GetMean();
    }
  cout << mean/10 << endl;;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "350 GeV Pions:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hPi350[i]->GetMean() << endl;
      mean += hPi350[i]->GetMean();
    }
  cout << mean/10 << endl;;
  mean = 0;

  cout << "=======================================" << endl;
}


// These have all the standard cuts
void bsdlate_pe_energies::Loop2()
{
  if (fChain == 0) return;
  
  TH1F * hPi350[12];
  TH1F * hPi300[12];
  TH1F * hPi250[12];
  
  TH1F * hE175[12];
  TH1F * hE150[12];
  TH1F * hE125[12];
  TH1F * hE100[12];
  TH1F * hE75[12];
    
  char name[1000];
  
  for( int i = 0; i < 12; i++ )
    {
      sprintf(name,"hPi350_%d",i);
      hPi350[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hPi300_%d",i);
      hPi300[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hPi250_%d",i);
      hPi250[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE75_%d",i);
      hE75[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE100_%d",i);
      hE100[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE125_%d",i);
      hE125[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE150_%d",i);
      hE150[i] = new TH1F(name,"",200,-100,1000);
   
      sprintf(name,"hE175_%d",i);
      hE175[i] = new TH1F(name,"",200,-100,1000);
    }

  Double_t gainA[12] = {4.23,4.23,3.92,3.09,5.81,5.03,3.17,5.30,4.31,3.90,4.04,3.68};
  Double_t myfCHA[12];
  Double_t myBSDLatePE[12];
  
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
            
      for( int i = 0; i < 12; i++ )
	{
	  myfCHA[i] = -1;
	  myBSDLatePE[i] = 0;
	}

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
      
      for( int i = 1; i <= 10; i++ )
	{
	  myBSDLatePE[i] = myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19);
	  
	  if( ptype == 2 && penergy == 350 ) hPi350[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 2 && penergy == 300 ) hPi300[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 2 && penergy == 250 ) hPi250[i]->Fill(myBSDLatePE[i]);
	  
	  if( ptype == 1 && penergy == 75 ) hE75[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 1 && penergy == 100 ) hE100[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 1 && penergy == 125 ) hE125[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 1 && penergy == 150 ) hE150[i]->Fill(myBSDLatePE[i]);
	  if( ptype == 1 && penergy == 175 ) hE175[i]->Fill(myBSDLatePE[i]);
	}
    }
  
  TCanvas * cvPi350 = new TCanvas("cvPi350");
  cvPi350->Divide(5,2);
  
  TCanvas * cvPi300 = new TCanvas("cvPi300");
  cvPi300->Divide(5,2);
  
  TCanvas * cvPi250 = new TCanvas("cvPi250");
  cvPi250->Divide(5,2);
  
  TCanvas * cvE175 = new TCanvas("cvE175");
  cvE175->Divide(5,2);
  
  TCanvas * cvE150 = new TCanvas("cvE150");
  cvE150->Divide(5,2);
  
  TCanvas * cvE125 = new TCanvas("cvE125");
  cvE125->Divide(5,2);
  
  TCanvas * cvE100 = new TCanvas("cvE100");
  cvE100->Divide(5,2);
  
  TCanvas * cvE75 = new TCanvas("cvE75");
  cvE75->Divide(5,2);
  
  for( int i = 1; i <= 10; i++ )
    {
      cvPi350->cd(i);
      hPi350[i]->Draw();
      
      cvPi300->cd(i);
      hPi300[i]->Draw();
      
      cvPi250->cd(i);
      hPi250[i]->Draw();
      
      cvE75->cd(i);
      hE75[i]->Draw();
      
      cvE100->cd(i);
      hE100[i]->Draw();
      
      cvE125->cd(i);
      hE125[i]->Draw();
      
      cvE150->cd(i);
      hE150[i]->Draw();
      
      cvE175->cd(i);
      hE175[i]->Draw();
    }

  Double_t mean = 0;

  cout << "=======================================" << endl;
  cout << "75 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE75[i]->GetMean() << endl;
      mean += hE75[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "100 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE100[i]->GetMean() << endl;
      mean += hE100[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "125 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE125[i]->GetMean() << endl;
      mean += hE125[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "150 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE150[i]->GetMean() << endl;
      mean += hE150[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "175 GeV Electrons:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hE175[i]->GetMean() << endl;
      mean += hE175[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "250 GeV Pions:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hPi250[i]->GetMean() << endl;
      mean += hPi250[i]->GetMean();
    }
  cout << mean/10 << endl;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "300 GeV Pions:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hPi300[i]->GetMean() << endl;
      mean += hPi300[i]->GetMean();
    }
  cout << mean/10 << endl;;
  mean = 0;

  cout << "=======================================" << endl;
  cout << "350 GeV Pions:" << endl;
  for( int i = 1; i <= 10; i++ )
    {
      cout << hPi350[i]->GetMean() << endl;
      mean += hPi350[i]->GetMean();
    }
  cout << mean/10 << endl;;
  mean = 0;

  cout << "=======================================" << endl;
}

#define Efrac_vs_calsum_sigma_cut_cxx
#include "Efrac_vs_calsum_sigma_cut.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>

using namespace std;

void Efrac_vs_calsum_sigma_cut::Loop()
{
  Int_t index = 0;
  Double_t means[5] = {5134.9, 6757.48, 8413.58, 10199.2, 11840.9};
  Double_t sigmas[5] = {390.96, 555.482, 685.362, 777.684, 875.024};
  
  if (fChain == 0) return;
    
  Double_t aE[5][2][10000];
  Double_t aE_2sigma[5][2][10000];
  Double_t aE_3sigma[5][2][10000];
    
  Int_t NE[5] = {0,0,0,0,0};
  Int_t NE_2sigma[5] = {0,0,0,0,0};
  Int_t NE_3sigma[5] = {0,0,0,0,0};

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

      // Need to find the index for mean and sigma
      if( ptype == 1 && penergy == 75 )  index = 0;
      else if( ptype == 1 && penergy == 100 ) index = 1;
      else if( ptype == 1 && penergy == 125 ) index = 2;
      else if( ptype == 1 && penergy == 150 ) index = 3;
      else if( ptype == 1 && penergy == 175 ) index = 4;
      else continue;

      for( int i = 1; i <= 10; i++ )
	myBSDLatePE += myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19);

      // Fill variables for electrons
      if( ptype == 1 && (penergy >=75 || penergy <= 175) )
	{
	  // all CALSum values
	  aE[index][0][NE[index]] = CALSum;
	  aE[index][1][NE[index]] = myBSDLatePE/CALSum;
	  NE[index]++;

	  // 2 sigma CALSum values
	  if( (CALSum < (means[index] - 2*sigmas[index])) ||
	      (CALSum > (means[index] + 2*sigmas[index])) )
	    {
	      aE_2sigma[index][0][NE_2sigma[index]] = CALSum;
	      aE_2sigma[index][1][NE_2sigma[index]] = myBSDLatePE/CALSum;
	      NE_2sigma[index]++;
	    }
	  
	  // 3 sigma CALSum values
	  if( (CALSum < (means[index] - 3*sigmas[index])) ||
	      (CALSum > (means[index] + 3*sigmas[index])) )
	    {
	      aE_3sigma[index][0][NE_3sigma[index]] = CALSum;
	      aE_3sigma[index][1][NE_3sigma[index]] = myBSDLatePE/CALSum;
	      NE_3sigma[index]++;
	    }
	}
    }
 
  cout << endl;
  
  // Make 75 GeV Plots
  TCanvas * cv75 = new TCanvas("cv75");
 
  cout << "====================================================" << endl;
  cout << "NE[0] = " << NE[0] << endl;
  cout << "NE_2sigma[0] = " << NE_2sigma[0] << endl;
  cout << "NE_3sigma[0] = " << NE_3sigma[0] << endl;
 
  TGraph * grE75 = new TGraph(NE[0],aE[0][0],aE[0][1]);
  grE75->SetMarkerColor(kRed);
  grE75->SetMarkerStyle(6);
  grE75->GetXaxis()->SetTitle("75 GeV CALSum (\"detector\" units)");
  grE75->GetXaxis()->CenterTitle();
  grE75->GetYaxis()->SetTitle("75 GeV Efrac (\"detector\" units)");
  grE75->GetYaxis()->CenterTitle();
  grE75->GetYaxis()->SetTitleOffset(1.2);
  grE75->SetTitle("");
  grE75->Draw("AP");
  
  TGraph * grE75_2sigma = new TGraph(NE_2sigma[0],aE_2sigma[0][0],aE_2sigma[0][1]);
  grE75_2sigma->SetMarkerStyle(6);
  grE75_2sigma->SetMarkerColor(kGreen);
  grE75_2sigma->Draw("Psame");

  TGraph * grE75_3sigma = new TGraph(NE_3sigma[0],aE_3sigma[0][0],aE_3sigma[0][1]);
  grE75_3sigma->SetMarkerStyle(6);
  grE75_3sigma->SetMarkerColor(kViolet);
  grE75_3sigma->Draw("Psame");
  
  // Make 100 GeV Plots
  TCanvas * cv100 = new TCanvas("cv100");
 
  cout << "====================================================" << endl;
  cout << "NE[1] = " << NE[1] << endl;
  cout << "NE_2sigma[1] = " << NE_2sigma[1] << endl;
  cout << "NE_3sigma[1] = " << NE_3sigma[1] << endl;
 
  TGraph * grE100 = new TGraph(NE[1],aE[1][0],aE[1][1]);
  grE100->SetMarkerColor(kRed);
  grE100->SetMarkerStyle(6);
  grE100->GetXaxis()->SetTitle("100 GeV CALSum (\"detector\" units)");
  grE100->GetXaxis()->CenterTitle();
  grE100->GetYaxis()->SetTitle("100 GeV Efrac (\"detector\" units)");
  grE100->GetYaxis()->CenterTitle();
  grE100->GetYaxis()->SetTitleOffset(1.2);
  grE100->SetTitle("");
  grE100->Draw("AP");
  
  TGraph * grE100_2sigma = new TGraph(NE_2sigma[1],aE_2sigma[1][0],aE_2sigma[1][1]);
  grE100_2sigma->SetMarkerStyle(6);
  grE100_2sigma->SetMarkerColor(kGreen);
  grE100_2sigma->Draw("Psame");

  TGraph * grE100_3sigma = new TGraph(NE_3sigma[1],aE_3sigma[1][0],aE_3sigma[1][1]);
  grE100_3sigma->SetMarkerStyle(6);
  grE100_3sigma->SetMarkerColor(kViolet);
  grE100_3sigma->Draw("Psame");

  
  // Make 125 GeV Plots
  TCanvas * cv125 = new TCanvas("cv125");
 
  cout << "====================================================" << endl;
  cout << "NE[2] = " << NE[2] << endl;
  cout << "NE_2sigma[2] = " << NE_2sigma[2] << endl;
  cout << "NE_3sigma[2] = " << NE_3sigma[2] << endl;
 
  TGraph * grE125 = new TGraph(NE[2],aE[2][0],aE[2][1]);
  grE125->SetMarkerColor(kRed);
  grE125->SetMarkerStyle(6);
  grE125->GetXaxis()->SetTitle("125 GeV CALSum (\"detector\" units)");
  grE125->GetXaxis()->CenterTitle();
  grE125->GetYaxis()->SetTitle("125 GeV Efrac (\"detector\" units)");
  grE125->GetYaxis()->CenterTitle();
  grE125->GetYaxis()->SetTitleOffset(1.2);
  grE125->SetTitle("");
  grE125->Draw("AP");
  
  TGraph * grE125_2sigma = new TGraph(NE_2sigma[2],aE_2sigma[2][0],aE_2sigma[2][1]);
  grE125_2sigma->SetMarkerStyle(6);
  grE125_2sigma->SetMarkerColor(kGreen);
  grE125_2sigma->Draw("Psame");

  TGraph * grE125_3sigma = new TGraph(NE_3sigma[2],aE_3sigma[2][0],aE_3sigma[2][1]);
  grE125_3sigma->SetMarkerStyle(6);
  grE125_3sigma->SetMarkerColor(kViolet);
  grE125_3sigma->Draw("Psame");
    
  // Make 150 GeV Plots
  TCanvas * cv150 = new TCanvas("cv150");
 
  cout << "====================================================" << endl;
  cout << "NE[3] = " << NE[3] << endl;
  cout << "NE_2sigma[3] = " << NE_2sigma[3] << endl;
  cout << "NE_3sigma[3] = " << NE_3sigma[3] << endl;
     
  TGraph * grE150 = new TGraph(NE[3],aE[3][0],aE[3][1]);
  grE150->SetMarkerColor(kRed);
  grE150->SetMarkerStyle(6);
  grE150->GetXaxis()->SetTitle("150 GeV CALSum (\"detector\" units)");
  grE150->GetXaxis()->CenterTitle();
  grE150->GetYaxis()->SetTitle("150 GeV Efrac (\"detector\" units)");
  grE150->GetYaxis()->CenterTitle();
  grE150->GetYaxis()->SetTitleOffset(1.2);
  grE150->SetTitle("");
  grE150->Draw("AP");
  
  TGraph * grE150_2sigma = new TGraph(NE_2sigma[3],aE_2sigma[3][0],aE_2sigma[3][1]);
  grE150_2sigma->SetMarkerStyle(6);
  grE150_2sigma->SetMarkerColor(kGreen);
  grE150_2sigma->Draw("Psame");

  TGraph * grE150_3sigma = new TGraph(NE_3sigma[3],aE_3sigma[3][0],aE_3sigma[3][1]);
  grE150_3sigma->SetMarkerStyle(6);
  grE150_3sigma->SetMarkerColor(kViolet);
  grE150_3sigma->Draw("Psame");
  
  // Make 175 GeV Plots
  TCanvas * cv175 = new TCanvas("cv175");
 
  cout << "====================================================" << endl;
  cout << "NE[4] = " << NE[4] << endl;
  cout << "NE_2sigma[4] = " << NE_2sigma[4] << endl;
  cout << "NE_3sigma[4] = " << NE_3sigma[4] << endl;
 
  TGraph * grE175 = new TGraph(NE[4],aE[4][0],aE[4][1]);
  grE175->SetMarkerColor(kRed);
  grE175->SetMarkerStyle(6);
  grE175->GetXaxis()->SetTitle("175 GeV CALSum (\"detector\" units)");
  grE175->GetXaxis()->CenterTitle();
  grE175->GetYaxis()->SetTitle("175 GeV Efrac (\"detector\" units)");
  grE175->GetYaxis()->CenterTitle();
  grE175->GetYaxis()->SetTitleOffset(1.2);
  grE175->SetTitle("");
  grE175->Draw("AP");
  
  TGraph * grE175_2sigma = new TGraph(NE_2sigma[4],aE_2sigma[4][0],aE_2sigma[4][1]);
  grE175_2sigma->SetMarkerStyle(6);
  grE175_2sigma->SetMarkerColor(kGreen);
  grE175_2sigma->Draw("Psame");

  TGraph * grE175_3sigma = new TGraph(NE_3sigma[4],aE_3sigma[4][0],aE_3sigma[4][1]);
  grE175_3sigma->SetMarkerStyle(6);
  grE175_3sigma->SetMarkerColor(kViolet);
  grE175_3sigma->Draw("Psame");
  
}

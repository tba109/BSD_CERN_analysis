#define cut_study_cxx
#include "cut_study.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

// This one counts non-exclusively. Good for seeing which cuts dominate.
void cut_study::Loop1()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  fChain->SetBranchStatus("*",1);
  
  Long64_t nbytes = 0, nb = 0;
  
  Double_t Ncut_Total = 0;
  Double_t Ncut_BSDLatePE = 0;
  Double_t Ncut_IsWithCal = 0;
  Double_t Ncut_CALSum = 0;
  Double_t Ncut_trig_layers = 0;
  Double_t Ncut_neg = 0;
  Double_t Nsurvive = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
      
      if( ptype != 1 ) continue;
      
      Ncut_Total++;
      
      if( !IsWithCal ) 
	{
	  Ncut_IsWithCal++;
	  // continue;
	}
  
      Bool_t neg = 0;
      for( int i = 1; i <= 10; i++ )
	if( fCHA[i] <= 0 )
	  {
	    neg = kTRUE;
	   }
      if( neg == kTRUE ) 
	{
	  Ncut_neg++;
	  // continue;
	}
          
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  // continue;
	}
      
      if( CALSum < 0 ) 
	{
	  Ncut_CALSum++;
	  // continue;
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
      if( Ntrig_layers < 6 )
	{
	  Ncut_trig_layers++;
	  // continue;
	}
      
      Nsurvive++;
    }
  
  cout << "Ncut_Total = " << Ncut_Total << endl;
  cout << "Ncut_IsWithCal = " << Ncut_IsWithCal << ", " << 100*Ncut_IsWithCal/Ncut_Total << endl;
  cout << "Ncut_neg = " << Ncut_neg << ", " << 100*Ncut_neg/Ncut_Total << endl;
  cout << "Ncut_BSDLatePE = " << Ncut_BSDLatePE << ", " << 100*Ncut_BSDLatePE/Ncut_Total << endl;
  cout << "Ncut_CALSum = " << Ncut_CALSum << ", " << 100*Ncut_CALSum/Ncut_Total << endl;
  cout << "Ncut_trig_layers = " << Ncut_trig_layers << ", " << 100*Ncut_trig_layers/Ncut_Total << endl;
  cout << "Nsurvive = " << Nsurvive << ", " << 100*Nsurvive/Ncut_Total << endl;
}


// This one counts exclusively.
void cut_study::Loop2()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  fChain->SetBranchStatus("*",1);
  
  Long64_t nbytes = 0, nb = 0;
  
  Double_t Ncut_Total = 0;
  Double_t Ncut_BSDLatePE = 0;
  Double_t Ncut_IsWithCal = 0;
  Double_t Ncut_CALSum = 0;
  Double_t Ncut_trig_layers = 0;
  Double_t Ncut_neg = 0;
  Double_t Nsurvive = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
      
      if( ptype != 2 ) continue;
      
      Ncut_Total++;
      
      if( !IsWithCal ) 
	{
	  Ncut_IsWithCal++;
	  continue;
	}
  
      Bool_t neg = 0;
      for( int i = 1; i <= 10; i++ )
	if( fCHA[i] <= 0 )
	  {
	    neg = kTRUE;
	   }
      if( neg == kTRUE ) 
	{
	  Ncut_neg++;
	  continue;
	}
          
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  continue;
	}
      
      if( CALSum < 0 ) 
	{
	  Ncut_CALSum++;
	  continue;
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
      if( Ntrig_layers < 6 )
	{
	  Ncut_trig_layers++;
	  continue;
	}
      
      Nsurvive++;
    }

  Int_t Nremain = 0;
  cout << "Ncut_Total = " << Ncut_Total << endl;
  Nremain = Ncut_Total;
  
  cout << "Ncut_IsWithCal = " << Ncut_IsWithCal << ", " << 100*Ncut_IsWithCal/Nremain << endl;
  Nremain -= Ncut_IsWithCal;

  cout << "Ncut_neg = " << Ncut_neg << ", " << 100*Ncut_neg/Nremain << endl;
  Nremain -= Ncut_neg;

  cout << "Ncut_BSDLatePE = " << Ncut_BSDLatePE << ", " << 100*Ncut_BSDLatePE/Nremain << endl;
  Nremain -= Ncut_BSDLatePE;

  cout << "Ncut_CALSum = " << Ncut_CALSum << ", " << 100*Ncut_CALSum/Nremain << endl;
  Nremain -= Ncut_CALSum;

  cout << "Ncut_trig_layers = " << Ncut_trig_layers << ", " << 100*Ncut_trig_layers/Nremain << endl;
  Nremain -= Ncut_trig_layers;

  cout << "Nsurvive = " << Nsurvive << ", " << 100*Nsurvive/Ncut_Total << endl;
  cout << "sanity: Nsurvive = " << Nsurvive << ", Nremain = " << Nremain << endl;
}

// This one counts exclusively with the revised cuts
void cut_study::Loop3()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  fChain->SetBranchStatus("*",1);
  
  Long64_t nbytes = 0, nb = 0;
  
  Double_t Ncut_Total = 0;
  Double_t Ncut_BSDLatePE = 0;
  Double_t Ncut_IsWithCal = 0;
  Double_t Ncut_CALSum = 0;
  Double_t Ncut_trig_layers = 0;
  Double_t Ncut_neg = 0;
  Double_t Nsurvive = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
      
      if( ptype != 2 ) continue;
      
      Ncut_Total++;
      
      if( !IsWithCal ) 
	{
	  Ncut_IsWithCal++;
	  continue;
	}
      
      /*  
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  continue;
	}
      */      
    
      /*
      if( CALSum < 0 ) 
	{
	  Ncut_CALSum++;
	  continue;
	}
      */    
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
      if( Ntrig_layers < 6 )
	{
	  Ncut_trig_layers++;
	  continue;
	}
        
      Bool_t neg = 0;
      for( int i = 1; i <= 10; i++ )
	if( fCHA[i] <= 0 )
	  {
	    neg = kTRUE;
	   }
      if( neg == kTRUE ) 
	{
	  Ncut_neg++;
	  // These are being adjusted now, so comment them out
	  // continue;
	}
      
      Nsurvive++;
    }

  Int_t Nremain = 0;
  cout << "Ncut_Total = " << Ncut_Total << endl;
  Nremain = Ncut_Total;
  
  cout << "Ncut_IsWithCal = " << Ncut_IsWithCal << ", " << 100*Ncut_IsWithCal/Nremain << endl;
  Nremain -= Ncut_IsWithCal;

  cout << "Ncut_neg = " << Ncut_neg << ", " << 100*Ncut_neg/Nremain << endl;
  // Nremain -= Ncut_neg;

  cout << "Ncut_BSDLatePE = " << Ncut_BSDLatePE << ", " << 100*Ncut_BSDLatePE/Nremain << endl;
  // Nremain -= Ncut_BSDLatePE;

  cout << "Ncut_CALSum = " << Ncut_CALSum << ", " << 100*Ncut_CALSum/Nremain << endl;
  // Nremain -= Ncut_CALSum;

  cout << "Ncut_trig_layers = " << Ncut_trig_layers << ", " << 100*Ncut_trig_layers/Nremain << endl;
  Nremain -= Ncut_trig_layers;

  cout << "Nsurvive = " << Nsurvive << ", " << 100*Nsurvive/Ncut_Total << endl;
  cout << "sanity: Nsurvive = " << Nsurvive << ", Nremain = " << Nremain << endl;
}


// This one counts exclusively with the revised cuts, but does it by energy, etc
void cut_study::Loop4()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  fChain->SetBranchStatus("*",1);
  
  Long64_t nbytes = 0, nb = 0;
  
  Double_t Ncut_Total[2][6];
  Double_t Ncut_BSDLatePE[2][6];
  Double_t Ncut_IsWithCal[2][6];
  Double_t Ncut_CALSum[2][6];
  Double_t Ncut_trig_layers[2][6];
  Double_t Ncut_neg[2][6];
  Double_t Nsurvive[2][6];
  
  for( int i = 0; i < 2; i++ )
    for( int j = 0; j < 6; j++ )
      {
	Ncut_Total[i][j] = 0;
	Ncut_BSDLatePE[i][j] = 0;
	Ncut_IsWithCal[i][j] = 0;
	Ncut_CALSum[i][j] = 0;
	Ncut_trig_layers[i][j] = 0;
	Ncut_neg[i][j] = 0;
	Nsurvive[i][j] = 0;
      }

  Int_t type_index = 0;
  Int_t energy_index = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
  
      type_index = 0;
      energy_index = 0;
      // First figure out type index and energy index 
      if( ptype == 1 )
	{
	  type_index = 0;
	  if( penergy == 50 )
	    energy_index = 0;
	  else if( penergy == 75 )
	    energy_index = 1;
	  else if( penergy == 100 )
	    energy_index = 2;
	  else if( penergy == 125 )
	    energy_index = 3;
	  else if( penergy == 150 )
	    energy_index = 4;
	  else if( penergy == 175 )
	    energy_index = 5;
	  else
	    continue;
	}
      else if( ptype == 2 )
	{
	  type_index = 1;
	  if( penergy == 250 )
	    energy_index = 0;
	  else if( penergy == 300 )
	    energy_index = 1;
	  else if( penergy == 350 )
	    energy_index = 2;
	  else
	    continue;
	}
      else
	continue;
     
      Ncut_Total[type_index][energy_index]++;
      
      if( !IsWithCal ) 
	{
	  Ncut_IsWithCal[type_index][energy_index]++;
	  continue;
	}
      
      /*  
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  continue;
	}
      */      
    
      /*
      if( CALSum < 0 ) 
	{
	  Ncut_CALSum++;
	  continue;
	}
      */    
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
      if( Ntrig_layers < 6 )
	{
	  Ncut_trig_layers[type_index][energy_index]++;
	  continue;
	}
        
      Bool_t neg = 0;
      for( int i = 1; i <= 10; i++ )
	if( fCHA[i] <= 0 )
	  {
	    neg = kTRUE;
	   }
      if( neg == kTRUE ) 
	{
	  Ncut_neg[type_index][energy_index]++;
	  // These are being adjusted now, so comment them out
	  // continue;
	}
      
      Nsurvive[type_index][energy_index]++;
    }

  for( int i = 0; i < 2; i++ )
    {
      for(int j = 0; j < 6; j++ )
	{
	  if( Ncut_Total[i][j] == 0 ) continue;
	  cout << "=========================================================================================" << endl;
	  cout << "type_index = " << i << ", " << "energy_index = " << j << endl;
	  Int_t Nremain = 0;
	  cout << "Ncut_Total = " << Ncut_Total[i][j] << endl;
	  Nremain = Ncut_Total[i][j];
	  
	  cout << "Survive IsWithCal " << Nremain - Ncut_IsWithCal[i][j] << ", Cuts " << 100*(Ncut_IsWithCal[i][j])/Nremain << "%" << endl;
	  Nremain -= Ncut_IsWithCal[i][j];
	  
	  cout << "Survive Ncut_trig_layers = " << Nremain - Ncut_trig_layers[i][j] << ", Cuts " << 100*(Ncut_trig_layers[i][j])/Nremain << "%" << endl;
	  Nremain -= Ncut_trig_layers[i][j];
	  
	  cout << "Nsurvive = " << Nsurvive[i][j] << ", " << 100*Nsurvive[i][j]/Ncut_Total[i][j] << endl;
	  cout << "sanity: Nsurvive = " << Nsurvive[i][j] << ", Nremain = " << Nremain << endl;
	  cout << "Conditioned Ncut_neg = " << Ncut_neg[i][j] << ", Conditioned " << 100*(Ncut_neg[i][j])/Nremain << "%" << endl;
	}
    }

  for( int i = 0; i < 2; i++ )
    {
      cout << "=========================================================================================" << endl;
      cout << "type_index = " << i << endl;
      Int_t Nremain = 0;
      cout << "Total Ncut_Total = " << Ncut_Total[i][0] + Ncut_Total[i][1] + Ncut_Total[i][2] + Ncut_Total[i][3] + Ncut_Total[i][4] + Ncut_Total[i][5] << endl;
      Nremain = Ncut_Total[i][0] + Ncut_Total[i][1] + Ncut_Total[i][2] + Ncut_Total[i][3] + Ncut_Total[i][4] + Ncut_Total[i][5];
      
      cout << "Survive IsWithCal " << "Cuts " << 100*(Ncut_IsWithCal[i][0] + Ncut_IsWithCal[i][1] + Ncut_IsWithCal[i][2] + Ncut_IsWithCal[i][3] +Ncut_IsWithCal[i][4] +Ncut_IsWithCal[i][5])/Nremain << "%" << endl;
      Nremain -= (Ncut_IsWithCal[i][0] + Ncut_IsWithCal[i][1] + Ncut_IsWithCal[i][2] + Ncut_IsWithCal[i][3] +Ncut_IsWithCal[i][4] +Ncut_IsWithCal[i][5]);
      
      cout << "Survive Ncut_trig_layers Cuts = " << 100*(Ncut_trig_layers[i][0] + Ncut_trig_layers[i][1] + Ncut_trig_layers[i][2] + Ncut_trig_layers[i][3] + Ncut_trig_layers[i][4] + Ncut_trig_layers[i][5])/Nremain << "%" << endl;
      Nremain -= Ncut_trig_layers[i][0] + Ncut_trig_layers[i][1] + Ncut_trig_layers[i][2] + Ncut_trig_layers[i][3] + Ncut_trig_layers[i][4] + Ncut_trig_layers[i][5];
      
      cout << "Nsurvive = " << (Nsurvive[i][0] + Nsurvive[i][1] + Nsurvive[i][2] + Nsurvive[i][3] + Nsurvive[i][4] + Nsurvive[i][5]) << endl;
      cout << "sanity: Nsurvive = " << Nsurvive[i][0] + Nsurvive[i][1] + Nsurvive[i][2] + Nsurvive[i][3] + Nsurvive[i][4] + Nsurvive[i][5] << ", Nremain = " << Nremain << endl;
      cout << "Conditioned Ncut_neg = " << Ncut_neg[i][0] + Ncut_neg[i][1] + Ncut_neg[i][2] + Ncut_neg[i][3] + Ncut_neg[i][4] + Ncut_neg[i][5] << ", Conditioned " << 100*(Ncut_neg[i][0] + Ncut_neg[i][1] + Ncut_neg[i][2] + Ncut_neg[i][3] + Ncut_neg[i][4] + Ncut_neg[i][5])/Nremain << "%" << endl;
    }
  
}

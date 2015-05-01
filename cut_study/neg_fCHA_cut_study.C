#define neg_fCHA_cut_study_cxx
#include "neg_fCHA_cut_study.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// Figure out where the fCHA[n] < 0 for n = 1 to 10 events are going for pions
void neg_fCHA_cut_study::Loop1()
{
    if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  TH1F * hA[10];
  char name[1000];
  for( int i = 0; i < 10; i++ )
    {
      sprintf(name,"hA_%d",i+1);
      hA[i] = new TH1F(name,"",200,-500,2200);
    }

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
      if( neg != kTRUE ) // pick out events where this is the case 
	{
	  Ncut_neg++;
	  continue;
	}
          
      /* Leave this one out since we want to see what negative events are doing
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  continue;
	}
      */

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

      for( int i = 0; i < 10; i++ )
	{
	  hA[i]->Fill(fCHA[i+1]);
	}
    }

  TCanvas * cv1 = new TCanvas("cv1");
  cv1->Divide(5,2);
  
  for( int i = 0; i < 10; i++ )
    {
      cv1->cd(i+1)->SetLogy();
      hA[i]->Draw();
    }
}

// Looks like fCHA[4] is the worse for saturating, so see what the others are doing when that one is saturated
void neg_fCHA_cut_study::Loop2()
{
    if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  TH1F * hA[10];
  char name[1000];
  for( int i = 0; i < 10; i++ )
    {
      sprintf(name,"hA_%d",i+1);
      hA[i] = new TH1F(name,"",200,-500,2200);
    }

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
      if( neg != kTRUE ) // pick out events where this is the case 
	{
	  Ncut_neg++;
	  continue;
	}
          
      if( fCHA[4] > 0 ) continue;

      /* Leave this one out since we want to see what negative events are doing
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  continue;
	}
      */

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

      for( int i = 0; i < 10; i++ )
	{
	  hA[i]->Fill(fCHA[i+1]);
	}
    }

  TCanvas * cv1 = new TCanvas("cv1");
  cv1->Divide(5,2);
  
  for( int i = 0; i < 10; i++ )
    {
      cv1->cd(i+1)->SetLogy();
      hA[i]->Draw();
    }
}

// Looks like fCHA[4] is saturated most often and fCHA[9] is saturated second most often
void neg_fCHA_cut_study::Loop3()
{
    if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  TH1F * hA[10];
  char name[1000];
  for( int i = 0; i < 10; i++ )
    {
      sprintf(name,"hA_%d",i+1);
      hA[i] = new TH1F(name,"",200,-500,2200);
    }

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
      if( neg != kTRUE ) // pick out events where this is the case 
	{
	  Ncut_neg++;
	  continue;
	}
          
      if( fCHA[4] > 0 && fCHA[9] > 0 ) continue;

      /* Leave this one out since we want to see what negative events are doing
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  continue;
	}
      */

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

      for( int i = 0; i < 10; i++ )
	{
	  hA[i]->Fill(fCHA[i+1]);
	}
    }

  TCanvas * cv1 = new TCanvas("cv1");
  cv1->Divide(5,2);
  
  for( int i = 0; i < 10; i++ )
    {
      cv1->cd(i+1)->SetLogy();
      hA[i]->Draw();
    }
}

// Okay, so now see what CALSum and Efrac look like for these saturated events.
void neg_fCHA_cut_study::Loop4()
{
    if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  TH1F * hCAL = new TH1F("hCAL","",100,-1000,30000);
  TH1F * hEfrac = new TH1F("hEfrac","",200,-0.5,2);

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
      if( neg != kTRUE ) // pick out events where this is the case 
	{
	  Ncut_neg++;
	  continue;
	}
          
      /* Leave this one out since we want to see what negative events are doing
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  continue;
	}
      */

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
      
      hCAL->Fill(CALSum);
      hEfrac->Fill(BSDLatePE/CALSum);
    }

  TCanvas * cv1 = new TCanvas("cv1");
  hCAL->Draw();
  
  TCanvas * cv2 = new TCanvas("cv2");
  hEfrac->Draw();

}

// Okay, so now let's choose that any events which have < 0 for fCHA end
// up getting the saturation value 2048
void neg_fCHA_cut_study::Loop5()
{
    if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  TH1F * hCAL = new TH1F("hCAL","",100,-1000,30000);
  TH1F * hEfrac = new TH1F("hEfrac","",200,-0.5,2);

  fChain->SetBranchStatus("*",1);
  
  Long64_t nbytes = 0, nb = 0;
  
  Double_t Ncut_Total = 0;
  Double_t Ncut_BSDLatePE = 0;
  Double_t Ncut_IsWithCal = 0;
  Double_t Ncut_CALSum = 0;
  Double_t Ncut_trig_layers = 0;
  Double_t Ncut_neg = 0;
  Double_t Nsurvive = 0;
  
  Double_t gainA[12] = {4.23,4.23,3.92,3.09,5.81,5.03,3.17,5.30,4.31,3.90,4.04,3.68};
  Double_t myfCHA[12];
  Double_t myBSDLatePE = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      /////////////////////////////////////////////////////////////////////////
      // These are just cuts for valid BSD/CAL data values
      // Need a BSD event which also includes the CAL
      ////////////////////////////////////////////////////////////////////////
      myBSDLatePE = 0;
      for( int i = 0; i < 12; i++ )
	myfCHA[i] = -1;
      
      if( ptype != 2 ) continue;
      
      Ncut_Total++;
      
      if( !IsWithCal ) 
	{
	  Ncut_IsWithCal++;
	  continue;
	}
  
      Bool_t neg = 0;
      for( int i = 1; i <= 10; i++ )
	{
	  if( fCHA[i] <= 0 )
	    {
	      neg = kTRUE;
	      myfCHA[i] = 2048;
	    }
	  else
	    myfCHA[i] = fCHA[i];
	}
      if( neg != kTRUE ) // pick out events where this is the case 
	{
	  Ncut_neg++;
	  continue;
	}

      for( int i = 1; i <= 10; i++ )
	myBSDLatePE += myfCHA[i]*0.25E-12/(gainA[i]*1.E6*1.6E-19);
          
      /* Leave this one out since we want to see what negative events are doing
      if( BSDLatePE < 0 ) 
	{
	  Ncut_BSDLatePE++;
	  continue;
	}
      */

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
      
      hCAL->Fill(CALSum);
      hEfrac->Fill(myBSDLatePE/CALSum);
    }

  TCanvas * cv1 = new TCanvas("cv1");
  hCAL->Draw();
  
  TCanvas * cv2 = new TCanvas("cv2");
  hEfrac->GetXaxis()->SetTitle("Efrac for Pions with Saturated Channels (detector units)");
  hEfrac->GetXaxis()->CenterTitle();
  hEfrac->Draw();

}

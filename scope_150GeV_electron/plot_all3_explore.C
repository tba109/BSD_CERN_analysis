// Tyler Anderson Tue Dec 11 16:43:54 EST 2012
// Find PMT 1B 150 GeV errors
void plot_all3_explore()
{
  gROOT->Reset();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat("emr");
  // gStyle->SetOptStat(1);
  
  //                 ksiourmen
  const Int_t Nfiles = 10;
  char ped_files[Nfiles][1000] = {
    "raw_files/F11p0_1p40_ped00000.txt",
    "raw_files/F11p40_1p80_ped00000.txt",
    "raw_files/F11p80_2p20_ped00000.txt",
    "raw_files/F12p20_2p60_ped00000.txt",
    "raw_files/F12p60_3p00_ped00000.txt",
    "raw_files/F13p00_3p40_ped00000.txt",
    "raw_files/F13p40_3p80_ped00000.txt",
    "raw_files/F13p80_4p20_ped00000.txt",
    "raw_files/F14p20_4p60_ped00000.txt",
    "raw_files/F14p60_5p00_ped00000.txt"};



  char event_files[Nfiles][1000] = {
    "raw_files/F11p00_1p4000000.txt",
    "raw_files/F11p40_1p8000000.txt",
    "raw_files/F11p80_2p2000000.txt",
    "raw_files/F12p20_2p6000000.txt",
    "raw_files/F12p60_3p0000000.txt",
    "raw_files/F13p00_3p4000000.txt",
    "raw_files/F13p40_3p8000000.txt",
    "raw_files/F13p80_4p2000000.txt",
    "raw_files/F14.20_4p6000000.txt",
    "raw_files/F14.60_5p0000000.txt"};
  
  Double_t start[Nfiles] = {1.0,1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6};
  Double_t end[Nfiles] = {1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6,5.0};
  Double_t mid[Nfiles];
  Double_t mean[Nfiles];
  Double_t mean_error[Nfiles];
  Double_t ped[Nfiles];
  Double_t ped_error[Nfiles];
  
  ifstream in;
  char dummy[1000];
  char name[1000];
  Double_t area;
  Double_t time;
  char comma;
  TH1F * hPed[Nfiles];
  TH1F * hEvent[Nfiles];
  char wait;
  
  /////////////////////////////////////////////////////////////////////////////////////
  // Find the ped data means and errors
  // Then fit the Peds to gaussian
  /////////////////////////////////////////////////////////////////////////////////////
  TCanvas * cPed = new TCanvas("cpeds");
  cPed->Divide(5,2);
  TF1 * fPed[12];
  char title[100];

  for( int i = 0; i < Nfiles; i++ )
    {
      cPed->cd(i+1);
      in.open(ped_files[i]);
      
      for( int j = 0; j < 12; j++ )
	{
	  in.getline(dummy,1000);
	  // cout << dummy << endl;
	}
      
      sprintf(name,"hPed_%d",i);
      hPed[i] = new TH1F(name,"",2000,-2000,10000);

      while( in.good() )
	{
	  in >> time >> comma >> area;
	  hPed[i]->Fill((-1*area*1.E12));
	}
      
      sprintf(title,"Pedestal %1.1f#mus to %1.1f#mus (pVs)",start[i]-0.1,end[i]-0.1);
      hPed[i]->GetXaxis()->SetTitle(title);
      hPed[i]->GetXaxis()->CenterTitle();
      hPed[i]->Draw();
      cout << start[i] << " " << end[i] << " " << hPed[i]->GetMean() << " " << hPed[i]->GetMeanError() << endl;
      
      ped[i] = hPed[i]->GetMean();
      ped_error[i] = hPed[i]->GetMeanError();

      sprintf(name,"fPed_%d",i);
      fPed[i] = new TF1(name,"gaus",-2000,10000);
      fPed[i]->SetNpx(100000);
      fPed[i]->SetParameter(1,ped[i]);
      fPed[i]->SetParameter(2,ped_error[i]);
      

      hPed[i]->Fit(fPed[i]);
      ped[i] = -1*fPed[i]->GetParameter(1);
      ped_error[i] = fPed[i]->GetParError(1);
            
      Double_t xmean = fPed[i]->GetParameter(1);
      Double_t xsigma = fPed[i]->GetParameter(2);
      
      hPed[i]->GetXaxis()->SetRangeUser(xmean-4*xsigma,xmean+4*xsigma);

      in.close();
    }

  //////////////////////////////////////////////////////////////////////////////////////
  // Find the event data means and errors
  /////////////////////////////////////////////////////////////////////////////////////
  TCanvas * cEvent = new TCanvas("cEvent");
  cEvent->Divide(5,2);

  TPaveStats * st[Nfiles];

  // TF1 fPeak[10];
  for( int i = 0; i < Nfiles; i++ )
    {
      cEvent->cd(i+1)->SetLogy(1);
      cEvent->cd(i+1)->SetRightMargin(0.1);
      cEvent->cd(i+1)->SetLeftMargin(0.15);      
      cEvent->cd(i+1)->SetTopMargin(0.01);      
      cEvent->cd(i+1);
      
      in.open(event_files[i]);

      for( int j = 0; j < 12; j++ )
	{
	  in.getline(dummy,1000);
	  // cout << dummy << endl;
	}
      
      sprintf(name,"h_%d",i);
      // hEvent[i] = new TH1F(name,"",600,-10000,2000);
      Float_t xmax = 2000;
      Float_t xmin = -100;
      //hEvent[i] = new TH1F(name,"",100,lowrange,100);
      hEvent[i] = new TH1F(name,"",200,xmin,xmax);

      while( in.good() )
	{
	  in >> time >> comma >> area;
	  
	  if( -1*1.E12*area < xmax)
	    hEvent[i]->Fill(-1*area*1.E12);
	}
      
      sprintf(title,"Signal %1.1f#mus to %1.1f#mus (pVs)",start[i],end[i]);
      hEvent[i]->GetXaxis()->SetTitle(title);
      hEvent[i]->GetXaxis()->CenterTitle();
      hEvent[i]->GetXaxis()->SetTitleSize(0.06);
      hEvent[i]->GetXaxis()->SetLabelSize(0.06);
      hEvent[i]->GetXaxis()->SetTitleOffset(0.81);
      hEvent[i]->GetYaxis()->SetTitle("Number of Events");
      hEvent[i]->GetYaxis()->CenterTitle();
      hEvent[i]->GetYaxis()->SetTitleSize(0.06);
      hEvent[i]->GetYaxis()->SetLabelSize(0.06);
      hEvent[i]->GetYaxis()->SetTitleOffset(1);
      hEvent[i]->Draw();
      
      cout << start[i] << " " << end[i] << " " << hEvent[i]->GetMean() << " " << hEvent[i]->GetMeanError() << endl;
      in.close();
      
      
      TF1 * f1t = new TF1("f1t","landau",xmin,xmax);
      hEvent[i]->Fit("f1t");
      TF1 * f2t = new TF1("f2t","f1t*x",xmin,xmax);
      
      mean[i] = -1*f2t->Integral(xmin,xmax)/f1t->Integral(xmin,xmax);
      // mean_error[i] = f1t->GetParameter(2);
      mean_error[i] = f1t->GetParError(1);
      cout << "mean = " << mean[i] << endl;

    }

  return;

  // Fitting a gaussian+landau. This didn't really work
  // Fitting just a landau always gave a better chisq/ndf
#if 0
  TCanvas * ct1 = new TCanvas("ct1");
  
  hEvent[0]->Scale(1./hEvent[0]->GetEntries());

  TF1 * fPeak1 = new TF1("fPeak1","landau",xmin,xmax);
  fPeak1->SetNpx(10000);
  fPeak1->SetLineColor(kBlack);
  hEvent[0]->Fit("fPeak1","R");
  

  Double_t lxmin = fPed[0]->GetParameter(1) - fPed[0]->GetParameter(2);
  Double_t lxmax = fPed[0]->GetParameter(1) + fPed[0]->GetParameter(2);
  TF1 * fPeak2 = new TF1("fPeak2","gaus",lxmin,lxmax);
  fPeak2->SetNpx(10000);
  fPeak2->SetLineColor(kGreen);
  fPeak2->SetParameter(1,fPed[0]->GetParameter(1));
  // cout << "parameter 1 set to " << fPed[0]->GetParameter(1) << endl;
  hEvent[0]->Fit("fPeak2","R+");
  
  Double_t par[6];

  fPeak1->GetParameters(&par[0]);
  fPeak2->GetParameters(&par[3]);
  
  TF1 * fPeak3 = new TF1("fPeak3","landau+gaus",xmin,xmax);
  fPeak3->SetParameters(par);
  hEvent[0]->Fit("fPeak3","R+");

#endif

  ////////////////////////////////////////////////////////////////////////////////////
  // Fit to a constant that will be subtracted off
  ///////////////////////////////////////////////////////////////////////////////////

  const Int_t Nproc = 10;
  
  for( int i = 0; i < Nproc; i++ ) 
    {
      mean[i] = -1*(mean[i]-ped[i]);
      start[i] = (start[i] - 0.1)*1000;
      end[i] = (end[i] - 0.1)*1000;
      mid[i] = (start[i] + end[i])*0.5; 
    }
  
  TCanvas * c2 = new TCanvas("cv2");
  TGraph * gr = new TGraph(Nproc,mid,mean);
  // TGraph * gr = new TGraph(14,end,mean);
  // TGraph * gr = new TGraph(14,mid,mean);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Time (ns)");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("");
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleOffset(1.2);
  gr->SetMarkerStyle(20);
  gr->Draw("AP");
  
  TF1 * f1 = new TF1("f1","expo+[2]",0,10000);
  f1->SetNpx(100000);
  gr->Fit("f1");
  
  cout << "Shower time constant is " << -1/f1->GetParameter(1) << " ns\n";
  
  // calculate error from having a wide center bin
  Double_t x2y_err[Nproc];
  for( int i = 0; i < Nproc; i++ )
    x2y_err[i] = f1->Eval(start[i]) - f1->Eval(end[i]);

  
  // Calculate the error by the derivative
  Double_t x2y_errb[Nproc];
  Double_t dy_exact = 0;
  Double_t dy_approx = 0;
  for( int i = 0; i < Nproc; i++ )
    {
      dy_exact = f1->Derivative(mid[i])*400;
      dy_approx = f1->Eval(end[i]) - f1->Eval(start[i]);
      x2y_errb[i] = fabs(dy_exact - dy_approx);
      cout << dy_exact << " " << dy_approx << " " << x2y_errb[i] << endl;
    }

  /////////////////////////////////////////////////////////////////////////////////////
  // This is greatly annoying, for some reason when I scale to PEs, the expo + [2] 
  // function stops fitting, but when I don't scale, it fits fine...                     
  // This is to get around that little problem                                          
  ////////////////////////////////////////////////////////////////////////////////////
  Double_t mean1[Nproc];
  
  for( int i = 0; i < Nproc; i++ )
    {
      // Convert to PEs/ns
      mean1[i] = (mean[i]-f1->GetParameter(2))*1.E-12/((50)*(1.6E-19)*3.31E6*400);
      
      // Add peak and ped mean errors in quadrature. 
      // These are the two largest sources of error. 
      // Peak error very much dominates.
      // Note that the error in the gain calibration is small enough to be neglected.
      // Also add error due the range in x value of the bin center
      // mean_error[i] = sqrt(mean_error[i]*mean_error[i] + ped_error[i]*ped_error[i] )*1.E-12/((50)*(1.6E-19)*3.31E6*400);
      // mean_error[i] = sqrt(mean_error[i]*mean_error[i] + ped_error[i]*ped_error[i] + x2y_err[i]*x2y_err[i])*1.E-12/((50)*(1.6E-19)*3.31E6*400);
      mean_error[i] = sqrt(mean_error[i]*mean_error[i] + ped_error[i]*ped_error[i] + x2y_errb[i]*x2y_errb[i])*1.E-12/((50)*(1.6E-19)*3.31E6*400);
    }
  
  TCanvas * c3 = new TCanvas("cv3");                                                      
  Double_t x_error[Nfiles];
  for( int i = 0; i < Nfiles; i++ )
    {
      x_error[i] = 0;
    }

  TGraphErrors * gr1 = new TGraphErrors(Nproc,mid,mean1,x_error,mean_error);
  gr1->SetTitle("");
  gr1->GetXaxis()->SetTitle("Time (ns)");
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->SetTitle("dN_{pe}/dt for 150 GeV Electron (counts/ns)");
  gr1->GetYaxis()->CenterTitle();
  gr1->GetYaxis()->SetTitleOffset(1.2);
  gr1->SetMarkerStyle(20);
  gr1->Draw("AP");
  
  TF1 * f11 = new TF1("f11","[2]*(expo)",0,10000);
  f11->SetNpx(100000);
  f11->SetParameter(0,f1->GetParameter(0));
  f11->SetParameter(1,f1->GetParameter(1));
  f11->FixParameter(2,1.E-12/((50)*(1.6E-19)*3.31E6*400));
  // f11->Draw("same");                                                                   
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gr1->Fit("f11");
  
  cout << "Time Constant = " << -1/(f11->GetParameter(1)) << " ns +/- " << f11->GetParError(1)/(f11->GetParameter(1)*f11->GetParameter(1)) << " ns" << endl;
  
  Double_t chisq = 0;
  for( int i = 0; i < Nproc; i++ )
    {
      chisq  += pow( (f11->Eval(mid[i])-mean1[i])/mean_error[i],2);
    }
  cout << "chisq/ndf = " << chisq << "/" << Nproc - 2 << endl;

  // Tyler Anderson Mon Jul 22 21:06:07 EDT 2013
  // cout << "here\n";
  // TF1 * fB = new TF1("fB","exp([0] + -1.*x/2700)",0,1000);
  // fB->SetNpx(100000);
  // gr1->Fit("fB");

}

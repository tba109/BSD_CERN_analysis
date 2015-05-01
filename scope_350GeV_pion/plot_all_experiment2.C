void plot_all_experiment()
{
  gROOT->Reset();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat("emr");
  const Int_t Nfiles = 17;
  char filenames[Nfiles][1000] = {
    "./raw_files/F11p25_1p65_pionlearningexample00000.txt",
    "./raw_files/F11p65_2p05_pionlearningexample00000.txt",
    "./raw_files/F12p05_2p45_pionlearningexample00000.txt",
    "./raw_files/F12p45_2p85_pionlearningexample00000.txt",
    "./raw_files/F12p85_3p25_pionlearningexample00000.txt",
    "./raw_files/F13p25_3p65_pionlearningexample00000.txt",
    "./raw_files/F13p65_4p05_pionlearningexample00000.txt",
    "./raw_files/F14p05_4p45_pionlearningexample00000.txt",
    "./raw_files/F14p45_4p85_pionlearningexample00000.txt",
    "./raw_files/F14p85_5p25_pionlearningexample00000.txt",
    "./raw_files/F15p25_5p65_pionlearningexample00000.txt",
    "./raw_files/F15p65_6p05_pionlearningexample00000.txt",
    "./raw_files/F16p05_6p45_pionlearningexample00000.txt",
    "./raw_files/F16p45_6p85_pionlearningexample00000.txt",
    "./raw_files/F16p85_7p25_pionlearningexample00000.txt",
    "./raw_files/F17p25_7p65_pionlearningexample00000.txt",
    "./raw_files/F17p65_8p05_pionlearningexample00000.txt"};
  
  Double_t start[Nfiles] = {1.25,1.65,2.05,2.45,2.85,3.25,3.65,4.05,4.45,4.85,5.25,5.65,6.05,6.45,6.85,7.25,7.65};
  Double_t end[Nfiles] = {1.65,2.05,2.45,2.85,3.25,3.65,4.05,4.45,4.85,5.25,5.65,6.05,6.45,6.85,7.25,7.65,8.05};
  Double_t mean[Nfiles];
  Double_t mean_error[Nfiles];
  Double_t mid[Nfiles];

  ifstream in;
  char dummy[1000];
  char name[1000];
  Double_t area;
  Double_t time;
  char comma;
  TH1F * h[Nfiles];
  char wait;
  TCanvas * c1 = new TCanvas("cv1");
  c1->Divide(5,3);
  
  TF1 * fPeak1[Nfiles];
  TF1 * fPeak2[Nfiles];
  TF1 * fPeak3[Nfiles];
  
  char name[100];
  Double_t par[6];

  Double_t x_landau_low[15] = {1000,1000,1000,1000,1000,800,600,600,600,600,600,600,600,600,600};
  Double_t x_gaus_high[15] =  { 550, 500, 450, 450, 450,450,450,450,450,450,450,450,450,450,450};
  
  for( int i = 0; i < 15; i++ )
    {
      c1->cd(i+1)->SetLogy(1);
      c1->cd(i+1)->SetRightMargin(0.1);
      c1->cd(i+1)->SetLeftMargin(0.15);      
      c1->cd(i+1)->SetTopMargin(0.01);      
      c1->cd(i+1);
      c1->cd(i+1)->SetLogy(1);
      // c1->cd(i+1);
      in.open(filenames[i]);
      
      for( int j = 0; j < 12; j++ )
	{
	  in.getline(dummy,1000);
	  // cout << dummy << endl;
	}
      
      sprintf(name,"h_%d",i);
      Float_t xmin = -2000;
      Float_t xmax = 10000;
      h[i] = new TH1F(name,"",200,xmin,xmax);

      while( in.good() )
	{
	  in >> time >> comma >> area;
	  
	  // if( area*1.E12 < -500 )
	  h[i]->Fill((-1*area*1.E12));
	}
      
      sprintf(name,"Signal %1.1f to %1.1f (pVs)",start[i]-0.1,end[i]-0.1);
      h[i]->GetXaxis()->SetTitle(name);
      h[i]->GetXaxis()->CenterTitle();
      h[i]->GetXaxis()->SetTitleSize(0.06);
      h[i]->GetXaxis()->SetLabelSize(0.06);
      h[i]->GetXaxis()->SetTitleOffset(0.81);
      h[i]->GetYaxis()->SetTitle("Number of Events");
      h[i]->GetYaxis()->CenterTitle();
      h[i]->GetYaxis()->SetTitleSize(0.06);
      h[i]->GetYaxis()->SetLabelSize(0.06);
      h[i]->GetYaxis()->SetTitleOffset(1);
      h[i]->Draw();

      cout << start[i] << " " << end[i] << " " << h[i]->GetMean() << " " << h[i]->GetMeanError() << endl;
      in.close();
      
      if( i < 9 )
	{
	  // Fit to a Landau + Gaussian
	  sprintf(name,"fPeak1_%d",i);
	  fPeak1[i] = new TF1(name,"landau",x_landau_low[i],10000);
	  fPeak1[i]->SetLineColor(kRed);
	  h[i]->Fit(fPeak1[i],"R");
	  
	  sprintf(name,"fPeak2_%d",i);
	  fPeak2[i] = new TF1(name,"gaus",0,x_gaus_high[i]);
	  fPeak2[i]->SetLineColor(kBlack);
	  h[i]->Fit(fPeak2[i],"R+");
	  
	  fPeak1[i]->GetParameters(&par[0]);
	  fPeak2[i]->GetParameters(&par[3]);
	  
	  sprintf(name,"fPeak3_%d",i);
	  fPeak3[i] = new TF1(name,"landau(0)+gaus(3)",-2000,10000);
	  fPeak3[i]->SetLineColor(kGreen);
	  fPeak3[i]->SetParameters(par);
	  h[i]->Fit(fPeak3[i],"R+");
	  
	}
      else
	{
	  // Fit to a Landau only
	  sprintf(name,"fPeak3_%d",i);
	  fPeak3[i] = new TF1(name,"landau",-2000,10000);
	  fPeak3[i]->SetLineColor(kGreen);
	  // fPeak3[i]->SetParameters(par);
	  h[i]->Fit(fPeak3[i],"R+");
	  par[0] = fPeak3[i]->GetParameter(0);
	  par[1] = fPeak3[i]->GetParameter(1);
	  par[2] = fPeak3[i]->GetParameter(2);
	}

      TF1 * fLandau = new TF1("fLandau","landau",-2000,10000);
      TF1 * fLandauX = new TF1("fLandauX","landau*x",-2000,10000);
      fLandau->SetParameters(&par[0]);
      fLandauX->SetParameters(&par[0]);      
      
      mean[i] = -1*fLandauX->Integral(-2000,10000)/fLandau->Integral(-2000,10000);
      // mean_error[i] = fPeak3[i]->GetParError(1);
      mean_error[i] = h[i]->GetMeanError();
      
      cout << "fLanda parameter error = " << fPeak3[i]->GetParError(1) << endl;

      // mean[i] = -1*h[i]->GetMean();
      // mean_error[i] = h[i]->GetMeanError();
      
      
    }
  
  return;
  const Int_t Nproc = 14;
  
  for( int i = 0; i < Nproc; i++ ) // Decided to avoid last 3 points b/c UMD changed beam energy
    {
      // Convert to dNpe/dt 
      // mean[i] = -1*(mean[i]+471.267)*1.E-12/((50)*(1.6E-19)*3.31E6*400);
      mean[i] = -1*(mean[i]+471.267);
      start[i] = (start[i] - 0.1)*1000;
      end[i] = (end[i] - 0.1)*1000;
      mid[i] = (start[i] + end[i])*0.5;
      // mean[i] = mean[i]/fcomp->Eval(mid[i]);
    }
  
  TCanvas * c2 = new TCanvas("cv2");
  TGraph * gr = new TGraph(14,start,mean);
  // TGraph * gr = new TGraph(14,end,mean);
  // TGraph * gr = new TGraph(14,mid,mean);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Time (ns)");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("dN_{pe}/dt for 350 GeV Pions (counts/ns)");
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleOffset(1.2);
  gr->SetMarkerStyle(20);
  gr->Draw("AP");
  
  TF1 * f1 = new TF1("f1","expo+[2]",0,10);
  f1->SetNpx(100000);
  gr->Fit("f1");
  
  // Calculate the error in y from the range of x values
  Double_t x2y_err[Nproc];
  for( int i = 0; i < Nproc; i++ )
    {
      x2y_err[i] = (f1->Eval(start[i]) - f1->Eval(end[i]));
      cout << "x2y_err["<<i<<"] = " << x2y_err[i] << " " << start[i] << " " << end[i] << endl;
    }

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
  cout << "Shower time constant is " << -1/f1->GetParameter(1) << " ns\n";
  
  /////////////////////////////////////////////////////////////////////////////////////
  // This is greatly annoying, for some reason when I scale to PEs, the expo + [2] 
  // function stops fitting, but when I don't scale, it fits fine...                     
  // This is to get around that little problem                                          
  ////////////////////////////////////////////////////////////////////////////////////
  Double_t mean1[Nproc];
  
  for( int i = 0; i < Nproc; i++ )
    {
      mean1[i] = (mean[i]-f1->GetParameter(2))*1.E-12/((50)*(1.6E-19)*3.31E6*400);
      // mean_error[i] = sqrt(mean_error[i]*mean_error[i])*1.E-12/((50)*(1.6E-19)*3.31E6*400);
      // mean_error[i] = sqrt(mean_error[i]*mean_error[i] + x2y_err[i]*x2y_err[i])*1.E-12/((50)*(1.6E-19)*3.31E6*400);
      mean_error[i] = sqrt(mean_error[i]*mean_error[i] + x2y_errb[i]*x2y_errb[i])*1.E-12/((50)*(1.6E-19)*3.31E6*400);
    }
  
  TCanvas * c2 = new TCanvas("cv2");                                                     
  
  Double_t x_errors[Nproc];
  for( int i = 0; i < Nproc; i++ )
    x_errors[i] = 0;
  
  TGraphErrors * gr1 = new TGraphErrors(Nproc, mid, mean1, x_errors, mean_error);
  gr1->SetTitle("");
  gr1->GetXaxis()->SetTitle("Time (ns)");
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->SetTitle("dN_{pe}/dt for 350 GeV Pions (counts/ns)");
  gr1->GetYaxis()->CenterTitle();
  gr1->GetYaxis()->SetTitleOffset(1.2);
  gr1->SetMarkerStyle(20);
  gr1->Draw("AP");
  
  TF1 * f11 = new TF1("f11","[2]*(expo)",0,10000);
  f11->SetNpx(100000);
  f11->SetParameter(0,f1->GetParameter(0));
  f11->SetParameter(1,f1->GetParameter(1));
  f11->FixParameter(2,1.E-12/((50)*(1.6E-19)*3.31E6*400));
  // f11->SetParameter(2,1.E-12/((50)*(1.6E-19)*3.31E6*400));
  // f11->Draw("same");                                                                   
  gStyle->SetOptFit(0);
  gr1->Fit("f11");
  
  cout << "Final time constant = " << -1/f11->GetParameter(1) << " +/- " << f11->GetParError(1)/(f11->GetParameter(1)*f11->GetParameter(1)) << endl;

  Double_t chisq = 0;
  Double_t chisq2 = 0;
  for( int i = 0; i < Nproc; i++ )
    {
      chisq  += pow( (f11->Eval(mid[i])-mean1[i])/mean_error[i],2);
      chisq2 += pow( (f11->Eval(mid[i])-mean1[i])/f11->Eval(mid[i]),2);
    }
  
  cout << "chisq/ndf = " << chisq << "/" << Nproc - 2 << endl;
  cout << "chisq2/ndf = " << chisq2 << "/" << Nproc - 2 << endl;
}

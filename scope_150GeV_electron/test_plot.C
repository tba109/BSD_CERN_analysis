{
  Double_t mid[10] = {1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8};
  Double_t ped[10] = {-166,-200,-172,-160,-165,-174,-167,-160,-169,-169};
  Double_t peak[10] = {-169,-189,-159,-143,-135,-150,-141,-136,-152,-147};

  TCanvas * cv1 = new TCanvas("cv1");
  TGraph * grPed = new TGraph(10,mid,ped);
  grPed->SetMarkerStyle(20);
  grPed->Draw("APL");

  TCanvas * cv2 = new TCanvas("cv2");
  TGraph * grPeak = new TGraph(10,mid,peak);
  grPeak->SetMarkerStyle(20);
  grPeak->Draw("APL");
  
  Double_t pp[10];
  for( int i = 0; i < 10; i++ )
    pp[i] = peak[i] - ped[i];

  TCanvas * cv3 = new TCanvas("cv3");
  TGraph * grPP = new TGraph(10,mid,pp);
  grPP->SetMarkerStyle(20);
  grPP->Draw("APL");

  Double_t start[6] = {-1000,-900,-800,-700,-600,-500};
  Double_t tc[6] = {1.04,1.66,1.70,1.92,2.26,2.56};

  TGraph * gr4 = new TGraph(6,start,tc);
  gr4->SetMarkerStyle(20);
  TCanvas * cv4 = new TCanvas("cv4");
  gr4->Draw("APL");

}

{
  // Double_t mid[9] = {1.45,1.85,2.25,2.65,3.05,3.45,3.85,4.25,4.65};
  // Double_t height[9] = {-1430,-1023,-922,-694,-729,-576,-608,-547,-490};
  
  Double_t mid[9] = {2.25,2.65,3.05,3.45,3.85,4.25,4.65};
  Double_t height[9] = {-922,-694,-729,-576,-608,-547,-490};
    
  TGraph * gr1 = new TGraph(7,mid,height);
  gr1->SetMarkerStyle(20);
  gr1->Draw("APL");
  TF1 * f1 = new TF1("f1","-expo +[2]",1,5);
  gr1->Fit("f1");
}

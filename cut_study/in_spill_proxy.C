{
  TChain * ch = new TChain("t");
  ch->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-24-15.root");
  
  ch->Draw("(0x7FFFFF & feventno):fcurtime");
  TGraph * gr = (TGraph*)gPad->GetPrimitive("Graph")->Clone();
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Time");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("Event Number");
  gr->GetYaxis()->CenterTitle();

  ch->Draw("(0x7FFFFF & feventno):fcurtime","IsWithCal");
  TGraph * grCut = (TGraph*)gPad->GetPrimitive("Graph")->Clone();
  grCut->SetMarkerColor(kGreen);
  
  gr->Draw("AP");
  grCut->Draw("Psame");
}

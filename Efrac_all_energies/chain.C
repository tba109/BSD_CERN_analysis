{
  gROOT->Reset();
  gStyle->SetOptStat(1);
  TChain * chPi350 = new TChain("t");
  chPi350->Add("/home/tyler/CERN2012/BSDCALCombined/*00-00-18*.root");
  chPi350->Add("/home/tyler/CERN2012/BSDCALCombined/*00-08-38*.root");
  chPi350->Add("/home/tyler/CERN2012/BSDCALCombined/*00-13-56*.root");
  chPi350->Add("/home/tyler/CERN2012/BSDCALCombined/*00-18-29*.root");
  chPi350->Add("/home/tyler/CERN2012/BSDCALCombined/*00-24-21*.root");
  chPi350->Add("/home/tyler/CERN2012/BSDCALCombined/*00-28-39*.root");
  chPi350->Add("/home/tyler/CERN2012/BSDCALCombined/*00-33-04*.root");
  

  // These are from a different day, so check b/c the BSD has moved
  // And the BCD/TCD is in place
  TChain * chPi350b = new TChain("t");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*18-24-46*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*18-33-03*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*18-43-53*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*18-51-40*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*18-59-26*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*19-07-23*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*19-20-32*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*19-28-15*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*19-36-26*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*19-55-23*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*20-03-16*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*20-13-18*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*20-28-18*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*20-51-09*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*20-36-11*.root");
  chPi350b->Add("/home/tyler/CERN2012/BSDCALCombined/*21-02-29*.root");
  
  TChain * chE150 = new TChain("t");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-24-15.root");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-24-15.root");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-10-16-35.root");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-10-24-38.root");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-10-30-54.root");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-46-34.root");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-37-37.root");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-19-36-38.root");
  chE150->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-19-49-55.root");

  TChain * chE150b = new TChain("t");
  chE150b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-13-01-07.root");
  chE150b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-13-09-10.root");
  chE150b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-13-17-02.root");
  
  TChain * chE100b = new TChain("t");
  chE100b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-13-32-47.root");
  chE100b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-13-37-35.root");
  
  TChain * chE175b = new TChain("t");
  chE175b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-15-09-26.root");

  TChain * chE125b = new TChain("t");
  chE125b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-15-31-00.root");
  chE125b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-15-35-09.root");
  
  TChain * chE50b = new TChain("t");
  chE50b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-15-44-33.root");
  chE50b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-15-48-40.root");

  TChain * chE75b = new TChain("t");
  chE75b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-15-53-21.root");
  chE75b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-15-57-47.root");

  TChain * chPi300b = new TChain("t");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-17-12-59.root");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-17-21-24.root");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-17-27-39.root");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-17-31-19.root");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-17-38-59.root");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-17-45-45.root");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-17-51-53.root");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-17-59-01.root");
  chPi300b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-18-06-43.root");

  TChain * chPi250b = new TChain("t");
  chPi250b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-21-12-45.root");
  chPi250b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-21-18-37.root");
  chPi250b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-21-25-05.root");
  chPi250b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-21-31-44.root");
  chPi250b->Add("/home/tyler/CERN2012/BSDCALCombined/*11-30-2012-21-38-17.root");
}

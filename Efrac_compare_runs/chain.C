{
  gROOT->Reset();
  gStyle->SetOptStat(1);
  TChain * chPi = new TChain("t");
  chPi->Add("/home/tyler/CERN2012/BSDCALCombined/*00-00-18*.root");
  chPi->Add("/home/tyler/CERN2012/BSDCALCombined/*00-08-38*.root");
  chPi->Add("/home/tyler/CERN2012/BSDCALCombined/*00-13-56*.root");
  chPi->Add("/home/tyler/CERN2012/BSDCALCombined/*00-18-29*.root");
  chPi->Add("/home/tyler/CERN2012/BSDCALCombined/*00-24-21*.root");
  chPi->Add("/home/tyler/CERN2012/BSDCALCombined/*00-28-39*.root");
  chPi->Add("/home/tyler/CERN2012/BSDCALCombined/*00-33-04*.root");
  

  // These are from a different day, so check b/c the BSD has moved
  // And the BCD/TCD is in place
  TChain * chPi2 = new TChain("t");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*18-24-46*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*18-33-03*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*18-43-53*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*18-51-40*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*18-59-26*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*19-07-23*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*19-20-32*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*19-28-15*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*19-36-26*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*19-55-23*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*20-03-16*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*20-13-18*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*20-28-18*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*20-51-09*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*20-36-11*.root");
  chPi2->Add("/home/tyler/CERN2012/BSDCALCombined/*21-02-29*.root");
  

  TChain * chE = new TChain("t");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-24-15.root");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-24-15.root");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-10-16-35.root");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-10-24-38.root");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-10-30-54.root");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-46-34.root");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-13-37-37.root");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-19-36-38.root");
  chE->Add("/home/tyler/CERN2012/BSDCALCombined/*11-27-2012-19-49-55.root");

  
}

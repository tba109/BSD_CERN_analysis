{
  gROOT->Reset();
  gStyle->SetOptStat(1);
  TChain * ch = new TChain("t");
  TChain * chb = new TChain("t");
  // 150 GeV electrons
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*00-00-18*.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*00-08-38*.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*00-13-56*.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*00-18-29*.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*00-24-21*.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*00-28-39*.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*00-33-04*.root");
  
  // 350 GeV pions
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*18-24-46*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*18-33-03*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*18-43-53*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*18-51-40*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*18-59-26*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*19-07-23*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*19-20-32*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*19-28-15*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*19-36-26*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*19-55-23*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*20-03-16*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*20-13-18*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*20-28-18*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*20-51-09*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*20-36-11*.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*21-02-29*.root");
  
  // 350 GeV pions
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-13-24-15.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-13-24-15.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-10-16-35.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-10-24-38.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-10-30-54.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-13-46-34.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-13-37-37.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-19-36-38.root");
  ch->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-27-2012-19-49-55.root");

  // 150 GeV electrons
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-13-01-07.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-13-09-10.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-13-17-02.root");
  
  // 100 GeV electrons
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-13-32-47.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-13-37-35.root");
  
  // 175 GeV electrons
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-15-09-26.root");

  // 125 GeV electrons
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-15-31-00.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-15-35-09.root");
  
  // 50 GeV electrons
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-15-44-33.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-15-48-40.root");

  // 75 GeV electrons
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-15-53-21.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-15-57-47.root");

  // 300 GeV pions
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-17-12-59.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-17-21-24.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-17-27-39.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-17-31-19.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-17-38-59.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-17-45-45.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-17-51-53.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-17-59-01.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-18-06-43.root");

  // 250 GeV pions
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-21-12-45.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-21-18-37.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-21-25-05.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-21-31-44.root");
  chb->Add("/home/tyler/BSD/CERN2012/BSDCALCombined/*11-30-2012-21-38-17.root");
}

{

  gInterpreter->AddIncludePath("/afs/cern.ch/work/c/chanon/MG5_aMC_v2_5_5/Delphes/");
  //gInterpreter->AddIncludePath("/afs/cern.ch/work/c/chanon/MG5_aMC_v2_5_5/ExRootAnalysis/");
  gSystem->Load("libDelphes");
  //gSystem->Load("libExRootAnalysis");
  //gSystem->Load("libPhysics");
  gROOT->ProcessLine(".L TestExRootAnalysis.C+");
  TestExRootAnalysis();
}

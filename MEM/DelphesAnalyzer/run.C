
{
//  gInterpreter->AddIncludePath("/afs/cern.ch/work/c/chanon/MG5_aMC_v2_5_5/Delphes/");
  gInterpreter->AddIncludePath("/home/junho/MG5_aMC_v2_6_1/Delphes/");
  gSystem->Load("libDelphes");
//  gROOT->ProcessLine(".L TestExRootAnalysis.C+");
  gROOT->ProcessLine(".L TestExRootAnalysis.C");
  TestExRootAnalysis();
}

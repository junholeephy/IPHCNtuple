



void ReadEntries(string file){

  TFile* f = TFile::Open(file.c_str());
  TTree* t = (TTree*) f->Get("Tree");
  cout << t->GetEntries() << endl;
  return;
}

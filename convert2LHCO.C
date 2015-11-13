// convert IPHCNtupleTTH-derived TTrees to LHCO format
// infile: input file name
// outfile: output file name
// nLep: number of saved leptons
// nJet: number of saved jets
// chan(should be always -1 if not running on ttH sample): -1 (all events), 0 (H->WW), 1 (H->ZZ), 2 (H->tautau)
// selec(pre-selection pass flag, https://twiki.cern.ch/twiki/bin/view/CMS/TTHmultileptonsFor13TeV): set to -1 for now

void convert2LHCO(std::string infile,std::string outfile,int nLep,int nJet,int chan,int selec)
{
   gROOT->SetBatch(1);   

   if( nLep > 4 )
     {
	std::cout << "Max number of leptons is 4" << std::endl;
	gApplication->Terminate();
     }
   if( nJet > 6 )
     {
	std::cout << "Max number of jets is 6" << std::endl;
	gApplication->Terminate();
     }   
   
   std::string fileName = infile;
   std::ofstream fout(outfile.c_str());

   int channel;
   int sel;
   double weight;
   
   double l1_pt, l1_eta, l1_phi, l1_m;
   int l1_type, l1_charge;
   double l2_pt, l2_eta, l2_phi, l2_m;
   int l2_type, l2_charge;
   double l3_pt, l3_eta, l3_phi, l3_m;
   int l3_type, l3_charge;
   double l4_pt, l4_eta, l4_phi, l4_m;
   int l4_type, l4_charge;
   double j1_pt, j1_eta, j1_phi, j1_m, j1_btag;
   int j1_ntrk;
   double j2_pt, j2_eta, j2_phi, j2_m, j2_btag;
   int j2_ntrk;
   double j3_pt, j3_eta, j3_phi, j3_m, j3_btag;
   int j3_ntrk;
   double j4_pt, j4_eta, j4_phi, j4_m, j4_btag;
   int j4_ntrk;
   double j5_pt, j5_eta, j5_phi, j5_m, j5_btag;
   int j5_ntrk;
   double j6_pt, j6_eta, j6_phi, j6_m, j6_btag;
   int j6_ntrk;
   double met, met_phi;
   
   TChain ch("tr");
   
   ch.SetBranchAddress("channel",&channel);
   ch.SetBranchAddress("sel",&sel);
   ch.SetBranchAddress("weight",&weight);
   
   ch.SetBranchAddress("l1_pt",&l1_pt);
   ch.SetBranchAddress("l1_eta",&l1_eta);
   ch.SetBranchAddress("l1_phi",&l1_phi);
   ch.SetBranchAddress("l1_m",&l1_m);
   ch.SetBranchAddress("l1_type",&l1_type);
   ch.SetBranchAddress("l1_charge",&l1_charge);

   ch.SetBranchAddress("l2_pt",&l2_pt);
   ch.SetBranchAddress("l2_eta",&l2_eta);
   ch.SetBranchAddress("l2_phi",&l2_phi);
   ch.SetBranchAddress("l2_m",&l2_m);
   ch.SetBranchAddress("l2_type",&l2_type);
   ch.SetBranchAddress("l2_charge",&l2_charge);

   ch.SetBranchAddress("l3_pt",&l3_pt);
   ch.SetBranchAddress("l3_eta",&l3_eta);
   ch.SetBranchAddress("l3_phi",&l3_phi);
   ch.SetBranchAddress("l3_m",&l3_m);
   ch.SetBranchAddress("l3_type",&l3_type);
   ch.SetBranchAddress("l3_charge",&l3_charge);

   ch.SetBranchAddress("l4_pt",&l4_pt);
   ch.SetBranchAddress("l4_eta",&l4_eta);
   ch.SetBranchAddress("l4_phi",&l4_phi);
   ch.SetBranchAddress("l4_m",&l4_m);
   ch.SetBranchAddress("l4_type",&l4_type);
   ch.SetBranchAddress("l4_charge",&l4_charge);
   
   ch.SetBranchAddress("j1_pt",&j1_pt);
   ch.SetBranchAddress("j1_eta",&j1_eta);
   ch.SetBranchAddress("j1_phi",&j1_phi);
   ch.SetBranchAddress("j1_m",&j1_m);
   ch.SetBranchAddress("j1_ntrk",&j1_ntrk);
   ch.SetBranchAddress("j1_btag",&j1_btag);

   ch.SetBranchAddress("j2_pt",&j2_pt);
   ch.SetBranchAddress("j2_eta",&j2_eta);
   ch.SetBranchAddress("j2_phi",&j2_phi);
   ch.SetBranchAddress("j2_m",&j2_m);
   ch.SetBranchAddress("j2_ntrk",&j2_ntrk);
   ch.SetBranchAddress("j2_btag",&j2_btag);

   ch.SetBranchAddress("j3_pt",&j3_pt);
   ch.SetBranchAddress("j3_eta",&j3_eta);
   ch.SetBranchAddress("j3_phi",&j3_phi);
   ch.SetBranchAddress("j3_m",&j3_m);
   ch.SetBranchAddress("j3_ntrk",&j3_ntrk);
   ch.SetBranchAddress("j3_btag",&j3_btag);

   ch.SetBranchAddress("j4_pt",&j4_pt);
   ch.SetBranchAddress("j4_eta",&j4_eta);
   ch.SetBranchAddress("j4_phi",&j4_phi);
   ch.SetBranchAddress("j4_m",&j4_m);
   ch.SetBranchAddress("j4_ntrk",&j4_ntrk);
   ch.SetBranchAddress("j4_btag",&j4_btag);

   ch.SetBranchAddress("j5_pt",&j5_pt);
   ch.SetBranchAddress("j5_eta",&j5_eta);
   ch.SetBranchAddress("j5_phi",&j5_phi);
   ch.SetBranchAddress("j5_m",&j5_m);
   ch.SetBranchAddress("j5_ntrk",&j5_ntrk);
   ch.SetBranchAddress("j5_btag",&j5_btag);

   ch.SetBranchAddress("j6_pt",&j6_pt);
   ch.SetBranchAddress("j6_eta",&j6_eta);
   ch.SetBranchAddress("j6_phi",&j6_phi);
   ch.SetBranchAddress("j6_m",&j6_m);
   ch.SetBranchAddress("j6_ntrk",&j6_ntrk);
   ch.SetBranchAddress("j6_btag",&j6_btag);
   
   ch.SetBranchAddress("met",&met);
   ch.SetBranchAddress("met_phi",&met_phi);
   
   ch.Add(fileName.c_str());
   
   int  l1_lhco_id = -666.;
   int  l2_lhco_id = -666.;
   int  l3_lhco_id = -666.;
   int  l4_lhco_id = -666.;
   if ()
  
  
   std::string fline00 = "#   typ     eta    phi       pt  jmass  ntrk  btag   had/em  dummy dummy";
   
   std::string del = "     ";
   std::string trig = "5";
   
   int nEvents = ch.GetEntries();
   for(int i=0;i<nEvents;i++)
     {
	ch.GetEntry(i);

	if( channel != chan && chan >= 0 ) continue;
	if( sel != selec && selec >= 0 ) continue;
	
	std::string fline0 = "0      " + std::string(Form("%d",i)) + "     " + trig;

	int nobj = 1;
	int nlep = 0;
	int njet = 0;
	
	// l1
	if( l1_m < 0 ) l1_m = 0.;
	std::string l1_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,l1_lhco_id,l1_eta,Phi_0_2Pi(l1_phi),l1_pt,0.0,l1_charge,0,0,0,0));
	if( l1_pt != -666 ) {nobj++;nlep++;}

	// l2
	if( l2_m < 0 ) l2_m = 0.;
	std::string l2_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,l2_lhco_id,l2_eta,Phi_0_2Pi(l2_phi),l2_pt,0.0,l2_charge,0,0,0,0));
	if( l2_pt != -666 ) {nobj++;nlep++;}

	// l3
	if( l3_m < 0 ) l3_m = 0.;
	std::string l3_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,l3_lhco_id,l3_eta,Phi_0_2Pi(l3_phi),l3_pt,0.0,l3_charge,0,0,0,0));
	if( l3_pt != -666 && nLep > 2 ) {nobj++;nlep++;}

	// l4
	if( l4_m < 0 ) l4_m = 0.;
	std::string l4_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,l4_lhco_id,l4_eta,Phi_0_2Pi(l4_phi),l4_pt,0.0,l4_charge,0,0,0,0));
	if( l4_pt != -666 && nLep > 3 ) {nobj++;nlep++;}

	// j1
	if( j1_m < 0 ) j1_m = 0.;
	int j1_btag_res = (j1_btag > 0.898) ? 2 : 0;
	std::string j1_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,4,j1_eta,Phi_0_2Pi(j1_phi),j1_pt,j1_m,1,j1_btag_res,0,0,0));
	if( j1_pt != -666 && nJet > 0 ) {nobj++;njet++;}

	// j2
	if( j2_m < 0 ) j2_m = 0.;
	int j2_btag_res = (j2_btag > 0.898) ? 2 : 0;
	std::string j2_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,4,j2_eta,Phi_0_2Pi(j2_phi),j2_pt,j2_m,1,j2_btag_res,0,0,0));
	if( j2_pt != -666 && nJet > 1 ) {nobj++;njet++;}

	// j3
	if( j3_m < 0 ) j3_m = 0.;
	int j3_btag_res = (j3_btag > 0.898) ? 2 : 0;
	std::string j3_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,4,j3_eta,Phi_0_2Pi(j3_phi),j3_pt,j3_m,1,j3_btag_res,0,0,0));
	if( j3_pt != -666 && nJet > 2 ) {nobj++;njet++;}

	// j4
	if( j4_m < 0 ) j4_m = 0.;
	int j4_btag_res = (j4_btag > 0.898) ? 2 : 0;
	std::string j4_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,4,j4_eta,Phi_0_2Pi(j4_phi),j4_pt,j4_m,1,j4_btag_res,0,0,0));
	if( j4_pt != -666 && nJet > 3 ) {nobj++;njet++;}

	// j5
	if( j5_m < 0 ) j5_m = 0.;
	int j5_btag_res = (j5_btag > 0.898) ? 2 : 0;
	std::string j5_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,4,j5_eta,Phi_0_2Pi(j5_phi),j5_pt,j5_m,1,j5_btag_res,0,0,0));
	if( j5_pt != -666 && nJet > 4 ) {nobj++;njet++;}

	// j6
	if( j6_m < 0 ) j6_m = 0.;
	int j6_btag_res = (j6_btag > 0.898) ? 2 : 0;
	std::string j6_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,4,j6_eta,Phi_0_2Pi(j6_phi),j6_pt,j6_m,1,j6_btag_res,0,0,0));
	if( j6_pt != -666 && nJet > 5 ) {nobj++;njet++;}

	// met
	std::string met_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,6,0.,met_phi,met,0.,0,0,0,0,0));
	if( met != -666 ) nobj++;
		
	if( nlep < nLep ) continue;
	if( njet < nJet ) continue;
	
	fout << fline00 << std::endl;
	fout << fline0 << std::endl;
	if( l1_pt != -666 ) fout << l1_fline << std::endl;
	if( l2_pt != -666 ) fout << l2_fline << std::endl;
	if( l3_pt != -666 && nLep > 2 ) fout << l3_fline << std::endl;
	if( l4_pt != -666 && nLep > 3 ) fout << l4_fline << std::endl;
	if( j1_pt != -666 && nJet > 0 ) fout << j1_fline << std::endl;
	if( j2_pt != -666 && nJet > 1 ) fout << j2_fline << std::endl;
	if( j3_pt != -666 && nJet > 2 ) fout << j3_fline << std::endl;
	if( j4_pt != -666 && nJet > 3 ) fout << j4_fline << std::endl;
	if( j5_pt != -666 && nJet > 4 ) fout << j5_fline << std::endl;
	if( j6_pt != -666 && nJet > 5 ) fout << j6_fline << std::endl;
	if( met != -666 ) fout << met_fline << std::endl;
     }   
   
   gApplication->Terminate();
}


float Phi_0_2Pi(float phi)
{
 float phi_0_2pi = phi;
 if (phi>= TMath::Pi()) phi_0_2pi  -= 2.*TMath::Pi();
 if (phi<0.)            phi_0_2pi  += 2.*TMath::Pi();
 return phi_0_2pi;
}

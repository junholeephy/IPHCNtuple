// convert IPHCNtupleTTH-derived TTrees to LHCO format
// infile: input file name
// outfile: output file name
// proc: characterizes 

//void convert2LHCO(std::string infile,std::string outfile,int nLep,int nJet,int chan,int selec)
void convert2LHCO_MC(std::string infile,std::string outfile,int proc)
{
   gROOT->SetBatch(1);   

   if( proc<-1 || proc > 6 )
     {
	std::cout << "proc can only take following values: -1,1,2,3,4,5,6" << std::endl;
	std::cout << "3l final state specific" << std::endl;	
	std::cout << "1,2,3,4: specific to the ttH final state, cf patches provided to madweight" << std::endl;
	std::cout << "5: specific to the ttZ with l+l-l+ final state" << std::endl;
	std::cout << "6: specific to the ttZ with l+l-l- final state" << std::endl;
	std::cout << "-1: no selection on the final state applied" << std::endl;
	
	gApplication->Terminate();
     }
  
   std::string fileName = infile;
   std::ofstream fout(outfile.c_str());

  int proc;
  
  Int_t mc_truth_h0_id;
  TLorentzVector* mc_truth_h0_p4 = 0;
  Int_t mc_truth_h0Wl1_id = -555;
  TLorentzVector* mc_truth_h0Wl1_p4 = 0;
  Int_t mc_truth_h0Wl2_id = -555;
  TLorentzVector* mc_truth_h0Wl2_p4 = 0;
  TLorentzVector* mc_truth_h0Wnu1_p4 = 0;
  TLorentzVector* mc_truth_h0Wnu2_p4 = 0;
  Int_t mc_truth_t1_id = -555;
  TLorentzVector* mc_truth_t1_p4 = 0;
  Int_t mc_truth_t2_id = -555;
  TLorentzVector* mc_truth_t2_p4 = 0;
  Int_t mc_truth_tb1_id = -555;
  TLorentzVector* mc_truth_tb1_p4 = 0;
  Int_t mc_truth_tb2_id = -555;
  TLorentzVector* mc_truth_tb2_p4 = 0;
  Int_t mc_truth_tWl1_id = -555; 
  TLorentzVector* mc_truth_tWl1_p4 = 0;
  TLorentzVector* mc_truth_tWnu1_p4 = 0;
  Int_t mc_truth_tWl2_id = -555;
  TLorentzVector* mc_truth_tWl2_p4 = 0;
  TLorentzVector* mc_truth_tWnu2_p4 = 0;
  Int_t mc_truth_tWq11_id = -555;
  TLorentzVector* mc_truth_tWq11_p4 = 0;
  Int_t mc_truth_tWq21_id = -555;
  TLorentzVector* mc_truth_tWq21_p4 = 0;
  Int_t mc_truth_tWq12_id = -555;
  TLorentzVector* mc_truth_tWq12_p4 = 0;
  Int_t mc_truth_tWq22_id = -555;
  TLorentzVector* mc_truth_tWq22_p4 = 0;
  Int_t mc_truth_Z_id = -555;
  TLorentzVector* mc_truth_Z_p4 = 0;
  Int_t mc_truth_Zl1_id = -555;
  TLorentzVector* mc_truth_Zl1_p4 = 0;
  Int_t mc_truth_Zl2_id = -555;
  TLorentzVector* mc_truth_Zl2_p4 = 0;
  
  Int_t mc_truth_tWtau1_id = -555; 
  TLorentzVector* mc_truth_tWtau1_p4 = 0;
  Int_t mc_truth_tWtau2_id = -555;
  TLorentzVector* mc_truth_tWtau2_p4 = 0;
 
  TBranch* b_mc_truth_h0_id;
  TBranch* b_mc_truth_h0_p4;
  TBranch* b_mc_truth_h0Wl1_id;
  TBranch* b_mc_truth_h0Wl1_p4;
  TBranch* b_mc_truth_h0Wl2_id;
  TBranch* b_mc_truth_h0Wl2_p4;
  TBranch* b_mc_truth_h0Wnu1_p4;
  TBranch* b_mc_truth_h0Wnu2_p4;
  TBranch* b_mc_truth_t1_id;
  TBranch* b_mc_truth_t1_p4;
  TBranch* b_mc_truth_t2_id;
  TBranch* b_mc_truth_t2_p4;
  TBranch* b_mc_truth_tb1_id;
  TBranch* b_mc_truth_tb1_p4;
  TBranch* b_mc_truth_tb2_id;
  TBranch* b_mc_truth_tb2_p4;
  TBranch* b_mc_truth_tWl1_id;
  TBranch* b_mc_truth_tWl1_p4;
  TBranch* b_mc_truth_tWl2_id;
  TBranch* b_mc_truth_tWl2_p4;
  TBranch* b_mc_truth_tWnu1_p4;
  TBranch* b_mc_truth_tWnu2_p4;
  TBranch* b_mc_truth_tWq11_id;
  TBranch* b_mc_truth_tWq11_p4;
  TBranch* b_mc_truth_tWq21_id;
  TBranch* b_mc_truth_tWq21_p4;
  TBranch* b_mc_truth_tWq12_id;
  TBranch* b_mc_truth_tWq12_p4;
  TBranch* b_mc_truth_tWq22_id;
  TBranch* b_mc_truth_tWq22_p4;
  TBranch* b_mc_truth_Z_id;
  TBranch* b_mc_truth_Z_p4;
  TBranch* b_mc_truth_Zl1_id;
  TBranch* b_mc_truth_Zl1_p4;
  TBranch* b_mc_truth_Zl2_id;
  TBranch* b_mc_truth_Zl2_p4;
  
  TBranch* b_mc_truth_tWtau1_id;
  TBranch* b_mc_truth_tWtau1_p4;
  TBranch* b_mc_truth_tWtau2_id;
  TBranch* b_mc_truth_tWtau2_p4;

  TChain ch("FlatTree/tree");
   
  ch.SetBranchAddress("mc_truth_h0_id",&mc_truth_h0_id,&b_mc_truth_h0_id);
  ch.SetBranchAddress("mc_truth_h0_p4",&mc_truth_h0_p4,&b_mc_truth_h0_p4);
  ch.SetBranchAddress("mc_truth_h0Wl1_id",&mc_truth_h0Wl1_id,&b_mc_truth_h0Wl1_id);
  ch.SetBranchAddress("mc_truth_h0Wl1_p4",&mc_truth_h0Wl1_p4,&b_mc_truth_h0Wl1_p4);
  ch.SetBranchAddress("mc_truth_h0Wl2_id",&mc_truth_h0Wl2_id,&b_mc_truth_h0Wl2_id);
  ch.SetBranchAddress("mc_truth_h0Wl2_p4",&mc_truth_h0Wl2_p4,&b_mc_truth_h0Wl2_p4);
  ch.SetBranchAddress("mc_truth_h0Wnu1_p4",&mc_truth_h0Wnu1_p4,&b_mc_truth_h0Wnu1_p4);
  ch.SetBranchAddress("mc_truth_h0Wnu2_p4",&mc_truth_h0Wnu2_p4,&b_mc_truth_h0Wnu2_p4);
  ch.SetBranchAddress("mc_truth_t1_id",&mc_truth_t1_id,&b_mc_truth_t1_id);
  ch.SetBranchAddress("mc_truth_t1_p4",&mc_truth_t1_p4,&b_mc_truth_t1_p4);
  ch.SetBranchAddress("mc_truth_t2_id",&mc_truth_t2_id,&b_mc_truth_t2_id);
  ch.SetBranchAddress("mc_truth_t2_p4",&mc_truth_t2_p4,&b_mc_truth_t2_p4);
  ch.SetBranchAddress("mc_truth_tb1_id",&mc_truth_tb1_id,&b_mc_truth_tb1_id);
  ch.SetBranchAddress("mc_truth_tb1_p4",&mc_truth_tb1_p4,&b_mc_truth_tb1_p4);
  ch.SetBranchAddress("mc_truth_tb2_id",&mc_truth_tb2_id,&b_mc_truth_tb2_id);
  ch.SetBranchAddress("mc_truth_tb2_p4",&mc_truth_tb2_p4,&b_mc_truth_tb2_p4);
  ch.SetBranchAddress("mc_truth_tWl1_id",&mc_truth_tWl1_id,&b_mc_truth_tWl1_id);
  ch.SetBranchAddress("mc_truth_tWl1_p4",&mc_truth_tWl1_p4,&b_mc_truth_tWl1_p4);
  ch.SetBranchAddress("mc_truth_tWl2_id",&mc_truth_tWl2_id,&b_mc_truth_tWl2_id);
  ch.SetBranchAddress("mc_truth_tWl2_p4",&mc_truth_tWl2_p4,&b_mc_truth_tWl2_p4);
  ch.SetBranchAddress("mc_truth_tWnu1_p4",&mc_truth_tWnu1_p4,&b_mc_truth_tWnu1_p4);
  ch.SetBranchAddress("mc_truth_tWnu2_p4",&mc_truth_tWnu2_p4,&b_mc_truth_tWnu2_p4);
  ch.SetBranchAddress("mc_truth_tWq11_id",&mc_truth_tWq11_id,&b_mc_truth_tWq11_id);
  ch.SetBranchAddress("mc_truth_tWq11_p4",&mc_truth_tWq11_p4,&b_mc_truth_tWq11_p4);
  ch.SetBranchAddress("mc_truth_tWq21_id",&mc_truth_tWq21_id,&b_mc_truth_tWq21_id);
  ch.SetBranchAddress("mc_truth_tWq21_p4",&mc_truth_tWq21_p4,&b_mc_truth_tWq21_p4);
  ch.SetBranchAddress("mc_truth_tWq12_id",&mc_truth_tWq12_id,&b_mc_truth_tWq12_id);
  ch.SetBranchAddress("mc_truth_tWq12_p4",&mc_truth_tWq12_p4,&b_mc_truth_tWq12_p4);
  ch.SetBranchAddress("mc_truth_tWq22_id",&mc_truth_tWq22_id,&b_mc_truth_tWq22_id);
  ch.SetBranchAddress("mc_truth_tWq22_p4",&mc_truth_tWq22_p4,&b_mc_truth_tWq22_p4);
  ch.SetBranchAddress("mc_truth_Z_id",&mc_truth_Z_id,&b_mc_truth_Z_id);
  ch.SetBranchAddress("mc_truth_Z_p4",&mc_truth_Z_p4,&b_mc_truth_Z_p4);
  ch.SetBranchAddress("mc_truth_Zl1_id",&mc_truth_Zl1_id,&b_mc_truth_Zl1_id);
  ch.SetBranchAddress("mc_truth_Zl1_p4",&mc_truth_Zl1_p4,&b_mc_truth_Zl1_p4);
  ch.SetBranchAddress("mc_truth_Zl2_id",&mc_truth_Zl2_id,&b_mc_truth_Zl2_id);
  ch.SetBranchAddress("mc_truth_Zl2_p4",&mc_truth_Zl2_p4,&b_mc_truth_Zl2_p4);

  ch.SetBranchAddress("mc_truth_tWtau1_id",&mc_truth_tWtau1_id,&b_mc_truth_tWtau1_id);
  ch.SetBranchAddress("mc_truth_tWtau1_p4",&mc_truth_tWtau1_p4,&b_mc_truth_tWtau1_p4);
  ch.SetBranchAddress("mc_truth_tWtau2_id",&mc_truth_tWtau2_id,&b_mc_truth_tWtau2_id);
  ch.SetBranchAddress("mc_truth_tWtau2_p4",&mc_truth_tWtau2_p4,&b_mc_truth_tWtau2_p4);

  ch.Add(fileName.c_str());
  
  // needed for LHCO  definition
  //
  float  mc_truth_h0Wl1_phi = -555.;
  float  mc_truth_h0Wl2_phi = -555.;
  float  mc_truth_tWl1_phi  = -555.; 
  float  mc_truth_tWl2_phi  = -555.;
  float  mc_truth_tWq11_phi = -555.;
  float  mc_truth_tWq21_phi = -555.;
  float  mc_truth_tWq12_phi = -555.;
  float  mc_truth_tWq22_phi = -555.;
  float  mc_truth_tb1_phi   = -555.;
  float  mc_truth_tb2_phi   = -555.;
  float  mc_truth_Zl1_phi   = -555.;
  float  mc_truth_Zl2_phi   = -555.;
  
  int  mc_truth_h0Wl1_lhco_id = -666.;
  int  mc_truth_h0Wl2_lhco_id = -666.;
  int  mc_truth_tWl1_lhco_id  = -666.; 
  int  mc_truth_tWl2_lhco_id  = -666.;
  int  mc_truth_tWl1_lhco_id  = -666.; 
  int  mc_truth_tWl2_lhco_id  = -666.;
  int  mc_truth_Zl1_lhco_id   = -666.; 
  int  mc_truth_Zl2_lhco_id   = -666.;
  
  //
  std::string fline00 = "#   typ     eta    phi       pt  jmass  ntrk  btag   had/em  dummy dummy";
  
  std::string del = "	  ";
  std::string trig = "8";
  
  int nEvents = ch.GetEntries();
  for(int i=0;i<nEvents;i++)
    {
       ch.GetEntry(i);

	//if( channel != chan && chan >= 0 ) continue;
	//if( sel != selec && selec >= 0 ) continue;
	
	std::string fline0 = "0      " + std::string(Form("%d",i)) + "     " + trig;

	int nobj = 1;
        
	//
        // LHCO lepton ID convention
	//
	if(abs(mc_truth_h0Wl1_id)==11)      mc_truth_h0Wl1_lhco_id = 1;
	else if(abs(mc_truth_h0Wl1_id)==13) mc_truth_h0Wl1_lhco_id = 2;
	if(abs(mc_truth_h0Wl2_id)==11)      mc_truth_h0Wl2_lhco_id = 1;
	else if(abs(mc_truth_h0Wl2_id)==13) mc_truth_h0Wl2_lhco_id = 2;
	if(abs(mc_truth_tWl1_id)==11)       mc_truth_tWl1_lhco_id  = 1;
	else if(abs(mc_truth_tWl1_id)==13)  mc_truth_tWl1_lhco_id  = 2;
	if(abs(mc_truth_tWl2_id)==11)       mc_truth_tWl2_lhco_id  = 1;
	else if(abs(mc_truth_tWl2_id)==13)  mc_truth_tWl2_lhco_id  = 2;
	if(abs(mc_truth_Zl1_id)==11)        mc_truth_Zl1_lhco_id   = 1;
	else if(abs(mc_truth_Zl1_id)==13)   mc_truth_Zl1_lhco_id   = 2;
	if(abs(mc_truth_Zl2_id)==11)        mc_truth_Zl2_lhco_id   = 1;
	else if(abs(mc_truth_Zl2_id)==13)   mc_truth_Zl2_lhco_id   = 2;
	
	//
	// LHCO phi convention (really needed ??)
	//
	mc_truth_h0Wl1_phi = Phi_0_2Pi(mc_truth_h0Wl1_p4->Phi());
        mc_truth_h0Wl2_phi = Phi_0_2Pi(mc_truth_h0Wl2_p4->Phi());
        mc_truth_tWl1_phi  = Phi_0_2Pi(mc_truth_tWl1_p4->Phi());
        mc_truth_tWl2_phi  = Phi_0_2Pi(mc_truth_tWl2_p4->Phi());
        mc_truth_tWq11_phi = Phi_0_2Pi(mc_truth_tWq11_p4->Phi());
        mc_truth_tWq21_phi = Phi_0_2Pi(mc_truth_tWq21_p4->Phi());
        mc_truth_tWq12_phi = Phi_0_2Pi(mc_truth_tWq12_p4->Phi());
        mc_truth_tWq22_phi = Phi_0_2Pi(mc_truth_tWq22_p4->Phi());
	mc_truth_tb1_phi   = Phi_0_2Pi(mc_truth_tb1_p4->Phi());
        mc_truth_tb2_phi   = Phi_0_2Pi(mc_truth_tb2_p4->Phi());
        mc_truth_Zl1_phi   = Phi_0_2Pi(mc_truth_Zl1_p4->Phi());
        mc_truth_Zl2_phi   = Phi_0_2Pi(mc_truth_Zl2_p4->Phi());
       
	
	//std::cout <<"---------------------------"<<std::endl;
	//std::cout <<mc_truth_h0Wl1_id<<" "<<mc_truth_h0Wl2_id<<" "<<mc_truth_tWl1_id<<" "<<mc_truth_tWl2_id<<std::endl;
	//std::cout <<mc_truth_tWq11_p4->Eta()<<" "<<mc_truth_tWq21_p4->Eta()<<" "<<mc_truth_tWq12_p4->Eta()<<" "<<mc_truth_tWq22_p4->Eta()<<std::endl;
	//std::cout <<mc_truth_Zl1_id<<" "<<mc_truth_Zl2_id<<std::endl;
	
	
	// l1
	std::string l1_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,mc_truth_h0Wl1_lhco_id,mc_truth_h0Wl1_p4->Eta(),mc_truth_h0Wl1_phi,mc_truth_h0Wl1_p4->Pt(),0.0,mc_truth_h0Wl1_id/abs(mc_truth_h0Wl1_id),0,0,0,0));
	if( mc_truth_h0Wl1_id != -666 ) {nobj++;}

	// l2
	std::string l2_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,mc_truth_h0Wl2_lhco_id,mc_truth_h0Wl2_p4->Eta(),mc_truth_h0Wl2_phi,mc_truth_h0Wl2_p4->Pt(),0.0,mc_truth_h0Wl2_id/abs(mc_truth_h0Wl2_id),0,0,0,0));
	if( mc_truth_h0Wl2_id != -666 ) {nobj++;}

	// l3
	std::string l3_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,mc_truth_tWl1_lhco_id,mc_truth_tWl1_p4->Eta(),mc_truth_tWl1_phi,mc_truth_tWl1_p4->Pt(),0.0,mc_truth_tWl1_id/abs(mc_truth_tWl1_id),0,0,0,0));
	if( mc_truth_tWl1_id != -666 ) {nobj++;}

	// l4
	std::string l4_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,mc_truth_tWl2_lhco_id,mc_truth_tWl2_p4->Eta(),mc_truth_tWl2_phi,mc_truth_tWl2_p4->Pt(),0.0,mc_truth_tWl2_id/abs(mc_truth_tWl2_id),0,0,0,0));
	if( mc_truth_tWl2_id != -666 ) {nobj++;}

        // l5
	std::string l5_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,mc_truth_Zl1_lhco_id,mc_truth_Zl1_p4->Eta(),mc_truth_Zl1_phi,mc_truth_Zl1_p4->Pt(),0.0,mc_truth_Zl1_id/abs(mc_truth_Zl1_id),0,0,0,0));
	if( mc_truth_Zl1_id != -666 ) {nobj++;}

        // l6
	std::string l6_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,mc_truth_Zl2_lhco_id,mc_truth_Zl2_p4->Eta(),mc_truth_Zl2_phi,mc_truth_Zl2_p4->Pt(),0.0,mc_truth_Zl2_id/abs(mc_truth_Zl2_id),0,0,0,0));
	if( mc_truth_Zl2_id != -666 ) {nobj++;}

        // j1
        std::string j1_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,4,mc_truth_tWq11_p4->Eta(),mc_truth_tWq11_phi,mc_truth_tWq11_p4->Pt(),0.0,1,0,0,0,0));
        if( mc_truth_tWq11_id != -666 ) {nobj++;}

        // j2
        std::string j2_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,4,mc_truth_tWq21_p4->Eta(),mc_truth_tWq21_phi,mc_truth_tWq21_p4->Pt(),0.0,1,0,0,0,0));
        if( mc_truth_tWq21_id != -666 ) {nobj++;}

        // j3
        std::string j3_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,4,mc_truth_tWq12_p4->Eta(),mc_truth_tWq12_phi,mc_truth_tWq12_p4->Pt(),0.0,1,0,0,0,0));
        if( mc_truth_tWq12_id != -666 ) {nobj++;}

        // j4
        std::string j4_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,4,mc_truth_tWq22_p4->Eta(),mc_truth_tWq22_phi,mc_truth_tWq22_p4->Pt(),0.0,1,0,0,0,0));
        if( mc_truth_tWq22_id != -666 ) {nobj++;}

        // bj1
	std::string b1_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					      nobj,4,mc_truth_tb1_p4->Eta(),mc_truth_tb1_phi,mc_truth_tb1_p4->Pt(),mc_truth_tb1_p4->M(),1,2,0,0,0));
	nobj++;
	 
	// bj2 
	std::string b2_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
				              nobj,4,mc_truth_tb2_p4->Eta(),mc_truth_tb2_phi,mc_truth_tb2_p4->Pt(),mc_truth_tb2_p4->M(),1,2,0,0,0));
	nobj++;
	
	// met
	std::string met_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
					       nobj,6,0.,0.,0.,0.,0.,0.,0.,0.,0.));
	nobj++;
	
	
	
	bool print = false;
	
	if( ( proc == 1 || proc == 2 ) && 
	    ( mc_truth_h0Wl1_id != -666 && mc_truth_h0Wl2_id != -666 ) && 
	    ( mc_truth_tWl1_id != -666 || mc_truth_tWl2_id != -666 ) && 
	    !(mc_truth_tWl1_id != -666 && mc_truth_tWl2_id != -666) ) print = true; 
	
	if( ( proc == 3 || proc == 4 ) && 
	    ( mc_truth_h0Wl1_id != -666 || mc_truth_h0Wl2_id != -666 ) && 
	    ( mc_truth_tWl1_id != -666 && mc_truth_tWl2_id != -666 ) && 
	    !(mc_truth_h0Wl1_id != -666 && mc_truth_h0Wl2_id != -666) ) print = true; 
	
	if( proc == 5 && mc_truth_Zl1_id != -666 && mc_truth_Zl2_id != -666 &&
	    ((mc_truth_tWl1_id != -666 && mc_truth_tWl1_id>0 && mc_truth_tWl2_id == -666 && mc_truth_tWtau2_id == -666) || 
	     (mc_truth_tWl2_id != -666 && mc_truth_tWl2_id>0 && mc_truth_tWl1_id == -666 && mc_truth_tWtau1_id == -666)) ) print = true; 
	  
	if( proc == 6 && mc_truth_Zl1_id != -666 && mc_truth_Zl2_id != -666 &&
	    ((mc_truth_tWl1_id != -666 && mc_truth_tWl1_id<0 && mc_truth_tWl2_id == -666 && mc_truth_tWtau2_id == -666) || 
	     (mc_truth_tWl2_id != -666 && mc_truth_tWl2_id<0 && mc_truth_tWl1_id == -666 && mc_truth_tWtau1_id == -666)) ) print = true; 
	  	  
	if( proc == -1 ) print = true;
	
 
	
	if ( print == true )
	 {
	
	  fout << fline00 << std::endl;
	  fout << fline0 << std::endl;
	
	  if( mc_truth_h0Wl1_id != -666 )   fout << l1_fline << std::endl;
	  if( mc_truth_h0Wl2_id != -666 )   fout << l2_fline << std::endl;
	  if( mc_truth_tWl1_id  != -666 )   fout << l3_fline << std::endl;
	  if( mc_truth_tWl2_id  != -666 )   fout << l4_fline << std::endl;
	  if( mc_truth_Zl1_id   != -666 )   fout << l5_fline << std::endl;
	  if( mc_truth_Zl2_id   != -666 )   fout << l6_fline << std::endl;	     
	  if( mc_truth_tWq11_id != -666 )   fout << j1_fline << std::endl;
	  if( mc_truth_tWq21_id != -666 )   fout << j2_fline << std::endl;
	  if( mc_truth_tWq12_id != -666 )   fout << j3_fline << std::endl;
	  if( mc_truth_tWq22_id != -666 )   fout << j4_fline << std::endl;
	  if( mc_truth_tb1_id   != -666 )   fout << b1_fline << std::endl;
	  if( mc_truth_tb2_id   != -666 )   fout << b2_fline << std::endl;
          
	  fout << met_fline << std::endl;	
         }
	
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

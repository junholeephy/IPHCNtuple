#include "../include/Hist.h"
#include "../include/Ranges.h"

#include <boost/bind.hpp>                                                                                                   
#include <boost/function.hpp>                                                                                               
#include <boost/type_traits.hpp>

Hist::Hist()
{
//   std::string foutlog = "log/"+std::string(flog);
   std::string foutlog = "log/list.txt";
   _fevc = fopen(foutlog.c_str(),"w");
//   std::string foutlogVal = "log/"+std::string(flog)+".val";
   std::string foutlogVal = "log/list.val";
   _fevcVal.open(foutlogVal.c_str());
   
   _v_ElectronTight = new std::vector<Electron>;
   _v_MuonTight = new std::vector<Muon>;
   _v_JetTight = new std::vector<Jet>;
}

Hist::~Hist()
{
   delete _v_ElectronTight;
   delete _v_MuonTight;
   delete _v_JetTight;
}

void Hist::init()
{
   rnd = new TRandom3(666);
   
   _fout = new TFile("hist/output.root","RECREATE");

   _trout = new TTree("tr","tree");
   
   _trout->Branch("weight",&m_weight,"weight/D");
   _trout->Branch("channel",&m_channel,"channel/I");
   _trout->Branch("sel",&m_sel,"sel/I");
   
   _trout->Branch("l1_type",&m_l1_type,"l1_type/I");
   _trout->Branch("l1_charge",&m_l1_charge,"l1_charge/I");                                                  
   _trout->Branch("l1_pt",&m_l1_pt,"l1_pt/D");
   _trout->Branch("l1_eta",&m_l1_eta,"l1_eta/D");
   _trout->Branch("l1_phi",&m_l1_phi,"l1_phi/D");
   _trout->Branch("l1_m",&m_l1_m,"l1_m/D");
   
   _trout->Branch("l2_type",&m_l2_type,"l2_type/I");
   _trout->Branch("l2_charge",&m_l2_charge,"l2_charge/I");
   _trout->Branch("l2_pt",&m_l2_pt,"l2_pt/D");
   _trout->Branch("l2_eta",&m_l2_eta,"l2_eta/D");
   _trout->Branch("l2_phi",&m_l2_phi,"l2_phi/D");
   _trout->Branch("l2_m",&m_l2_m,"l2_m/D");

   _trout->Branch("l3_type",&m_l3_type,"l3_type/I");
   _trout->Branch("l3_charge",&m_l3_charge,"l3_charge/I");
   _trout->Branch("l3_pt",&m_l3_pt,"l3_pt/D");
   _trout->Branch("l3_eta",&m_l3_eta,"l3_eta/D");
   _trout->Branch("l3_phi",&m_l3_phi,"l3_phi/D");
   _trout->Branch("l3_m",&m_l3_m,"l3_m/D");

   _trout->Branch("l4_type",&m_l4_type,"l4_type/I");
   _trout->Branch("l4_charge",&m_l4_charge,"l4_charge/I");
   _trout->Branch("l4_pt",&m_l4_pt,"l4_pt/D");
   _trout->Branch("l4_eta",&m_l4_eta,"l4_eta/D");
   _trout->Branch("l4_phi",&m_l4_phi,"l4_phi/D");
   _trout->Branch("l4_m",&m_l4_m,"l4_m/D");
   
   _trout->Branch("j1_ntrk",&m_j1_ntrk,"j1_ntrk/I");
   _trout->Branch("j1_btag",&m_j1_btag,"j1_btag/D");
   _trout->Branch("j1_pt",&m_j1_pt,"j1_pt/D");
   _trout->Branch("j1_eta",&m_j1_eta,"j1_eta/D");
   _trout->Branch("j1_phi",&m_j1_phi,"j1_phi/D");
   _trout->Branch("j1_m",&m_j1_m,"j1_m/D");

   _trout->Branch("j2_ntrk",&m_j2_ntrk,"j2_ntrk/I");                                                        
   _trout->Branch("j2_btag",&m_j2_btag,"j2_btag/D");                                                        
   _trout->Branch("j2_pt",&m_j2_pt,"j2_pt/D");                                                              
   _trout->Branch("j2_eta",&m_j2_eta,"j2_eta/D");                                                           
   _trout->Branch("j2_phi",&m_j2_phi,"j2_phi/D");                                                           
   _trout->Branch("j2_m",&m_j2_m,"j2_m/D");

   _trout->Branch("j3_ntrk",&m_j3_ntrk,"j3_ntrk/I");                                                        
   _trout->Branch("j3_btag",&m_j3_btag,"j3_btag/D");
   _trout->Branch("j3_pt",&m_j3_pt,"j3_pt/D");                                                              
   _trout->Branch("j3_eta",&m_j3_eta,"j3_eta/D");                                                           
   _trout->Branch("j3_phi",&m_j3_phi,"j3_phi/D");                                                           
   _trout->Branch("j3_m",&m_j3_m,"j3_m/D");

   _trout->Branch("j4_ntrk",&m_j4_ntrk,"j4_ntrk/I");
   _trout->Branch("j4_btag",&m_j4_btag,"j4_btag/D");
   _trout->Branch("j4_pt",&m_j4_pt,"j4_pt/D");                                                              
   _trout->Branch("j4_eta",&m_j4_eta,"j4_eta/D");                                                           
   _trout->Branch("j4_phi",&m_j4_phi,"j4_phi/D");                                                           
   _trout->Branch("j4_m",&m_j4_m,"j4_m/D");

   _trout->Branch("j5_ntrk",&m_j5_ntrk,"j5_ntrk/I");
   _trout->Branch("j5_btag",&m_j5_btag,"j5_btag/D");
   _trout->Branch("j5_pt",&m_j5_pt,"j5_pt/D");                                                              
   _trout->Branch("j5_eta",&m_j5_eta,"j5_eta/D");                                                           
   _trout->Branch("j5_phi",&m_j5_phi,"j5_phi/D");                                                           
   _trout->Branch("j5_m",&m_j5_m,"j5_m/D");

   _trout->Branch("j6_ntrk",&m_j6_ntrk,"j6_ntrk/I");
   _trout->Branch("j6_btag",&m_j6_btag,"j6_btag/D");
   _trout->Branch("j6_pt",&m_j6_pt,"j6_pt/D");                                                              
   _trout->Branch("j6_eta",&m_j6_eta,"j6_eta/D");                                                           
   _trout->Branch("j6_phi",&m_j6_phi,"j6_phi/D");                                                           
   _trout->Branch("j6_m",&m_j6_m,"j6_m/D");
   
   _trout->Branch("met_phi",&m_met_phi,"met_phi/D");
   _trout->Branch("met",&m_met,"met/D");
   
/*   hname.clear(); 

   int histname_dilep_n = 1;
   histname_dilep[0] = "h_mll_";

   int histname_lep_n = 0;
//   histname_lep[0] = "h_l_d0s_";
   
   int histname_lep_2d_n = 0;
//   histname_lep_2d[0] = "h_l1_origin_vs_type_";

   type_n = 3;                                                                                                              
   type[0] = "nonQCD";                                                                                                      
   type[1] = "QCD";                                                                                                         
   type[2] = "ALL";

   sel_n = 1;
   sel[0] = "nosel";

   int dilchan_n = 3;
   dilchan[0] = "ee";
   dilchan[1] = "mm";
   dilchan[2] = "em";
   
   int dilchar_n = 2;
   dilchar[0] = "ss";
   dilchar[1] = "os";
   
   int jets_n = 1;
   jets[0]  = "jge0";
   
   sys_low_n = 1;
   sys_low[0]   = "";
//   sys_low[1]   = "_pu_low";

   sys_up_n = sys_low_n-1;
//   sys_up[0]   = "_pu_up";

   sys_n = sys_low_n + sys_up_n;
   
   for(int is1=0;is1<sys_low_n;is1++)
     {
	sys[is1] = sys_low[is1];
     }   
   for(int is2=0;is2<sys_up_n;is2++)
     {
	sys[sys_low_n+is2] = sys_up[is2];
     }   

   _s_Dilepton = new std::vector<std::pair<std::vector<std::string>,double*> >();
   _m1d_Dilepton = new std::map<std::string, TH1D*>();

   std::vector<double*> set_hist_jet;
   set_hist_jet.clear();

   std::vector<double*> set_hist_jet_2d;
   set_hist_jet_2d.clear();

   std::vector<double*> set_hist_jet_3d;
   set_hist_jet_3d.clear();

   set_hist_jet.push_back(RANGE::set_j_pt);
   set_hist_jet.push_back(RANGE::set_j_JP);
//   set_hist_jet.push_back(RANGE::set_j_ntrk);
//   set_hist_jet.push_back(RANGE::set_j_njet);
//   set_hist_jet.push_back(RANGE::set_j_eta);
   
////   set_hist_jet.push_back(RANGE::set_j_mupt);
//   set_hist_jet.push_back(RANGE::set_j_eta);
//   set_hist_jet.push_back(RANGE::set_j_phi);
//   set_hist_jet.push_back(RANGE::set_j_mass);
//   set_hist_jet.push_back(set_j_mueta);
//   set_hist_jet.push_back(set_j_muphi);
//   set_hist_jet.push_back(set_j_muptrel);
//   set_hist_jet.push_back(RANGE::set_j_svntrk);
//   set_hist_jet.push_back(RANGE::set_j_svmass);
//   set_hist_jet.push_back(set_j_npv);
//   set_hist_jet.push_back(set_j_npu);
//   set_hist_jet.push_back(RANGE::set_j_nmuon);
//   set_hist_jet.push_back(RANGE::set_j_nsv);
   set_hist_jet.push_back(RANGE::set_j_nseltrk);
   set_hist_jet.push_back(RANGE::set_j_eta);
   set_hist_jet.push_back(RANGE::set_j_ntrk);
   set_hist_jet.push_back(RANGE::set_j_njet);
   set_hist_jet.push_back(RANGE::set_j_nsv);
   
//   set_hist_jet_2d.push_back(set_j_2d_ntrk_vs_pt);
//   set_hist_jet_2d.push_back(set_j_2d_ntrk_vs_mupt);

   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_CSV);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_nseltrk);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_nseltrk_vs_ntrkgen);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_nseltrk_vs_njet);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_nseltrk_vs_jeta);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_nseltrk_vs_mupt);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_nseltrk_vs_njet);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_nseltrk_vs_mupt);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_nseltrk_vs_jeta);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_njet_vs_mupt);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_njet_vs_jeta);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_mupt_vs_jeta);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_nseltrk_vs_nsv);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_njet_vs_nsv);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_mupt_vs_nsv);
  // set_hist_jet_2d.push_back(RANGE::set_j_2d_jeta_vs_nsv);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_nseltrk);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_njet);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_mupt);
   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_jeta);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_nsv);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_pt);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_pt_vs_nseltrk);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_pt_vs_njet);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_pt_vs_mupt);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_pt_vs_jeta);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_pt_vs_nsv);
   
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_pt);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_mupt);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_JP_vs_ntrk);
//   set_hist_jet_2d.push_back(RANGE::set_j_2d_ntrk_vs_mupt);
   
//   set_hist_jet_3d.push_back(set_j_3d_ntrk_vs_pt_vs_mupt);
   
//   set_hist_event.push_back(set_event_mjjMax);

   std::string titl;
   
   for(int j=0;j<flav_n;j++)
     {	
	for(int ss=0;ss<sel_n;ss++)
	  {
	     for(int se=0;se<eta_n;se++)
	       {
		  for(int ssb=0;ssb<selb_n;ssb++)
		    {
		       for(int ssa=0;ssa<seladd_n;ssa++)
			 {
			    for(int h=0;h<histname_jet_n;h++)
			      {
				 std::string hn = histname_jet[h]+flav[j]+"_"+sel[ss]+"_"+eta[se]+"_"+selb[ssb]+seladd[ssa];
				 hname.push_back(hn);
				 
				 for(int s=0;s<sys_n;s++)
				   {
				      titl = hn+sys[s];
				      std::vector<std::string> names;
				      names.clear();
				      names.push_back(titl);
				      names.push_back(sys[s]);
				      if( h == 7 )
					_s_Jet->push_back(std::make_pair(names,RANGE::set_j_mupt[ss]));
//				      else if( h == 1 )
//				      else if( h == 4 )
//					_s_Jet->push_back(std::make_pair(names,RANGE::set_j_JP[ss]));
				      else
					_s_Jet->push_back(std::make_pair(names,set_hist_jet.at(h)));
				   }				 
			      }

			    for(int h=0;h<histname_jet_2d_n;h++)
			      {
				 std::string hn = histname_jet_2d[h]+flav[j]+"_"+sel[ss]+"_"+eta[se]+"_"+selb[ssb]+seladd[ssa];
				 hname.push_back(hn);
				 
				 for(int s=0;s<sys_n;s++)
				   {
				      titl = hn+sys[s];
				      std::vector<std::string> names;
				      names.clear();
				      names.push_back(titl);
				      names.push_back(sys[s]);

				      _s2_Jet->push_back(std::make_pair(names,set_hist_jet_2d.at(h)));
				   }				 
			      }

			    for(int h=0;h<histname_jet_3d_n;h++)
			      {
				 std::string hn = histname_jet_3d[h]+flav[j]+"_"+sel[ss]+"_"+eta[se]+"_"+selb[ssb]+seladd[ssa];
				 hname.push_back(hn);
				 
				 for(int s=0;s<sys_n;s++)
				   {
				      titl = hn+sys[s];
				      std::vector<std::string> names;
				      names.clear();
				      names.push_back(titl);
				      names.push_back(sys[s]);

				      _s3_Jet->push_back(std::make_pair(names,set_hist_jet_3d.at(h)));
				   }				 
			      }
			 }
		    }		  
	       }
	  }
     }   
   
   for(int i=0;i<_s_Jet->size();i++)
     {
	TH1D *_h1d = new TH1D(_s_Jet->at(i).first.at(0).c_str(),
			      _s_Jet->at(i).first.at(0).c_str(),
			      _s_Jet->at(i).second[0],
			      _s_Jet->at(i).second[1],
			      _s_Jet->at(i).second[2]);
	
	_h1d->Sumw2();
	
	_m1d_Jet->insert(std::pair<std::string,TH1D*>(_s_Jet->at(i).first.at(0),_h1d));
     }
	
   for(int i=0;i<_s2_Jet->size();i++)
     {
	TH2D *_h2d = new TH2D(_s2_Jet->at(i).first.at(0).c_str(),
			      _s2_Jet->at(i).first.at(0).c_str(),
			      _s2_Jet->at(i).second[0],
			      _s2_Jet->at(i).second[1],
			      _s2_Jet->at(i).second[2],
			      _s2_Jet->at(i).second[3],
			      _s2_Jet->at(i).second[4],
			      _s2_Jet->at(i).second[5]);
	
	_h2d->Sumw2();
	
	_m2d_Jet->insert(std::pair<std::string,TH2D*>(_s2_Jet->at(i).first.at(0),_h2d));
     }

   for(int i=0;i<_s3_Jet->size();i++)
     {
	TH3D *_h3d = new TH3D(_s3_Jet->at(i).first.at(0).c_str(),
			      _s3_Jet->at(i).first.at(0).c_str(),
			      _s3_Jet->at(i).second[0],
			      _s3_Jet->at(i).second[1],
			      _s3_Jet->at(i).second[2],
			      _s3_Jet->at(i).second[3],
			      _s3_Jet->at(i).second[4],
			      _s3_Jet->at(i).second[5],
			      _s3_Jet->at(i).second[6],
			      _s3_Jet->at(i).second[7],
			      _s3_Jet->at(i).second[8]);
	
	_h3d->Sumw2();
	
	_m3d_Jet->insert(std::pair<std::string,TH3D*>(_s3_Jet->at(i).first.at(0),_h3d));
     }
   
   if( doRW )
     {       
	std::cout << "Read corrections" << std::endl;
	
	TFile fcor("corrections.root");
	
	if( !fcor.IsOpen() )
	  {
	     std::cout << "Input file with RW factors can not be opened !" << std::endl;
	     exit(1);
	  }
	if( doRW == 1 )
	  {	     
	     for(int isf=0;isf<rw_n;isf++)
	       {
		  std::vector<std::vector<std::pair<double,double> > > v_ptBin;
		  std::vector<std::vector<std::pair<double,double> > > v_ptSf;
		  v_ptBin.clear();
		  v_ptSf.clear();
		  
		  for(int ipt=0;ipt<sel_n;ipt++)
		    {		       
		       std::string sfname = "sf_"+sel[ipt]+"_nosel_"+name_rw[isf];
		       TH1D *hSF = (TH1D*)fcor.Get(sfname.c_str());
		       std::cout << sfname << std::endl;
		       
		       std::vector<std::pair<double,double> > v_valBin;
		       std::vector<std::pair<double,double> > v_valSf;
		       v_valBin.clear();
		       v_valSf.clear();
		       
		       int nb = hSF->GetXaxis()->GetNbins();
		       for(int ib=1;ib<=nb;ib++)
			 {		  
			    float pt1 = hSF->GetXaxis()->GetBinLowEdge(ib);
			    float pt2 = hSF->GetXaxis()->GetBinUpEdge(ib);
			    float sf = hSF->GetBinContent(ib);
			    float sferr = hSF->GetBinError(ib);
			    
			    std::pair<double,double> binfo = std::make_pair(pt1,pt2);
			    std::pair<double,double> sfinfo = std::make_pair(sf,sferr);
			    v_valBin.push_back(binfo);
			    v_valSf.push_back(sfinfo);
			 }
		       v_ptBin.push_back(v_valBin);
		       v_ptSf.push_back(v_valSf);
		    }	     
		  rwBin.push_back(v_ptBin);
		  rwSf.push_back(v_ptSf);
		  //
		  v_ptBin.clear();
		  v_ptSf.clear();
		  
		  for(int ipt=0;ipt<sel_n;ipt++)
		    {		       
		       std::string sfname = "sf_"+sel[ipt]+"_"+selb[1]+"_"+name_rw[isf];
		       TH1D *hSF = (TH1D*)fcor.Get(sfname.c_str());
		       std::cout << sfname << std::endl;
		       
		       std::vector<std::pair<double,double> > v_valBin;
		       std::vector<std::pair<double,double> > v_valSf;
		       v_valBin.clear();
		       v_valSf.clear();
		       
		       int nb = hSF->GetXaxis()->GetNbins();
		       for(int ib=1;ib<=nb;ib++)
			 {		  
			    float pt1 = hSF->GetXaxis()->GetBinLowEdge(ib);
			    float pt2 = hSF->GetXaxis()->GetBinUpEdge(ib);
			    float sf = hSF->GetBinContent(ib);
			    float sferr = hSF->GetBinError(ib);
			    
			    std::pair<double,double> binfo = std::make_pair(pt1,pt2);
			    std::pair<double,double> sfinfo = std::make_pair(sf,sferr);
			    v_valBin.push_back(binfo);
			    v_valSf.push_back(sfinfo);
			 }
		       v_ptBin.push_back(v_valBin);
		       v_ptSf.push_back(v_valSf);
		    }	     
		  rwBin_btag.push_back(v_ptBin);
		  rwSf_btag.push_back(v_ptSf);
	       }		
	  }
	if( doRW == 2 )
	  {
	     std::vector<std::vector<std::pair<double,double> > > v_pt2DBinX;
	     std::vector<std::vector<std::pair<double,double> > > v_pt2DBinY;
	     std::vector<std::vector<std::pair<double,double> > > v_pt2DSf;
	     v_pt2DBinX.clear();
	     v_pt2DBinY.clear();
	     v_pt2DSf.clear();
		  
	     for(int ipt=0;ipt<sel_n;ipt++)
	       {		       
		  std::string sfname = "sf_"+sel[ipt]+"_nosel_"+name_rw_2d;
		  TH2D *hSF = (TH2D*)fcor.Get(sfname.c_str());
		  std::cout << sfname << std::endl;
		       
		  std::vector<std::pair<double,double> > v_val2DBinX;
		  std::vector<std::pair<double,double> > v_val2DBinY;
		  std::vector<std::pair<double,double> > v_val2DSf;
		  v_val2DBinX.clear();
		  v_val2DBinY.clear();
		  v_val2DSf.clear();
		       
		  int nbX = hSF->GetXaxis()->GetNbins();
		  for(int ib=1;ib<=nbX;ib++)
		    {		  
		       float pt1 = hSF->GetXaxis()->GetBinLowEdge(ib);
		       float pt2 = hSF->GetXaxis()->GetBinUpEdge(ib);
		       
		       std::pair<double,double> binfo = std::make_pair(pt1,pt2);
		       v_val2DBinX.push_back(binfo);
		    }
		  v_pt2DBinX.push_back(v_val2DBinX);

		  int nbY = hSF->GetYaxis()->GetNbins();
		  for(int ib=1;ib<=nbY;ib++)
		    {		  
		       float pt1 = hSF->GetYaxis()->GetBinLowEdge(ib);
		       float pt2 = hSF->GetYaxis()->GetBinUpEdge(ib);
		       
		       std::pair<double,double> binfo = std::make_pair(pt1,pt2);
		       v_val2DBinY.push_back(binfo);
		    }
		  v_pt2DBinY.push_back(v_val2DBinY);

		  for(int ibx=1;ibx<=nbX;ibx++)
		    {		  
		       for(int iby=1;iby<=nbY;iby++)
			 {		  
			    float sf = hSF->GetBinContent(ibx,iby);
			    float sferr = hSF->GetBinError(ibx,iby);
			    
			    std::pair<double,double> sfinfo = std::make_pair(sf,sferr);
			    v_val2DSf.push_back(sfinfo);
			 }		       
		    }
		  
		  v_pt2DSf.push_back(v_val2DSf);
	       }
	     
	     rw2DBinX.push_back(v_pt2DBinX);
	     rw2DBinY.push_back(v_pt2DBinY);
	     
	     rw2DSf.push_back(v_pt2DSf);
	  }	
	if( doRW == 3 )
	  {
	     std::vector<std::vector<std::pair<double,double> > > v_pt3DBinX;
	     std::vector<std::vector<std::pair<double,double> > > v_pt3DBinY;
	     std::vector<std::vector<std::pair<double,double> > > v_pt3DBinZ;
	     std::vector<std::vector<std::pair<double,double> > > v_pt3DSf;
	     v_pt3DBinX.clear();
	     v_pt3DBinY.clear();
	     v_pt3DBinZ.clear();
	     v_pt3DSf.clear();
		  
	     for(int ipt=0;ipt<sel_n;ipt++)
	       {		       
		  std::string sfname = "sf_"+sel[ipt]+"_nosel_"+name_rw_3d;
		  TH3D *hSF = (TH3D*)fcor.Get(sfname.c_str());
		  std::cout << sfname << std::endl;
		       
		  std::vector<std::pair<double,double> > v_val3DBinX;
		  std::vector<std::pair<double,double> > v_val3DBinY;
		  std::vector<std::pair<double,double> > v_val3DBinZ;
		  std::vector<std::pair<double,double> > v_val3DSf;
		  v_val3DBinX.clear();
		  v_val3DBinY.clear();
		  v_val3DBinZ.clear();
		  v_val3DSf.clear();
		       
		  int nbX = hSF->GetXaxis()->GetNbins();
		  for(int ib=1;ib<=nbX;ib++)
		    {		  
		       float pt1 = hSF->GetXaxis()->GetBinLowEdge(ib);
		       float pt2 = hSF->GetXaxis()->GetBinUpEdge(ib);
		       
		       std::pair<double,double> binfo = std::make_pair(pt1,pt2);
		       v_val3DBinX.push_back(binfo);
		    }
		  v_pt3DBinX.push_back(v_val3DBinX);

		  int nbY = hSF->GetYaxis()->GetNbins();
		  for(int ib=1;ib<=nbY;ib++)
		    {		  
		       float pt1 = hSF->GetYaxis()->GetBinLowEdge(ib);
		       float pt2 = hSF->GetYaxis()->GetBinUpEdge(ib);
		       
		       std::pair<double,double> binfo = std::make_pair(pt1,pt2);
		       v_val3DBinY.push_back(binfo);
		    }
		  v_pt3DBinY.push_back(v_val3DBinY);

		  int nbZ = hSF->GetZaxis()->GetNbins();
		  for(int ib=1;ib<=nbZ;ib++)
		    {		  
		       float pt1 = hSF->GetZaxis()->GetBinLowEdge(ib);
		       float pt2 = hSF->GetZaxis()->GetBinUpEdge(ib);
		       
		       std::pair<double,double> binfo = std::make_pair(pt1,pt2);
		       v_val3DBinZ.push_back(binfo);
		    }
		  v_pt3DBinZ.push_back(v_val3DBinZ);
		  
		  for(int ibx=1;ibx<=nbX;ibx++)
		    {		  
		       for(int iby=1;iby<=nbY;iby++)
			 {		  
			    for(int ibz=1;ibz<=nbZ;ibz++)
			      {		  
				 float sf = hSF->GetBinContent(ibx,iby,ibz);
				 float sferr = hSF->GetBinError(ibx,iby,ibz);
				 
				 std::pair<double,double> sfinfo = std::make_pair(sf,sferr);
				 v_val3DSf.push_back(sfinfo);
			      }
			 }		       
		    }
		  
		  v_pt3DSf.push_back(v_val3DSf);
	       }
	     
	     rw3DBinX.push_back(v_pt3DBinX);
	     rw3DBinY.push_back(v_pt3DBinY);
	     rw3DBinZ.push_back(v_pt3DBinZ);
	     
	     rw3DSf.push_back(v_pt3DSf);
	  }	
	fcor.Close();
     }
   
   jesTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters("Summer13_V5_DATA_UncertaintySources_AK5PF.txt", "Total")));      

   // JER taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
   
   cJER[0] = 1.052; // 0.0-0.5
   cJER[1] = 1.057; // 0.5-1.1
   cJER[2] = 1.096; // 1.1-1.7
   cJER[3] = 1.134; // 1.7-2.3
   cJER[4] = 1.288; // 2.3-5.0
   
   cJER_down[0] = 0.990;
   cJER_down[1] = 1.001;
   cJER_down[2] = 1.032;
   cJER_down[3] = 1.042;
   cJER_down[4] = 1.089;

   cJER_up[0] = 1.115;
   cJER_up[1] = 1.114;
   cJER_up[2] = 1.161;
   cJER_up[3] = 1.228;
   cJER_up[4] = 1.488;
   
   for(int j=0;j<flav_n;j++)
     {	
	for(int ss=0;ss<sel_n;ss++)
	  {
	     for(int se=0;se<eta_n;se++)
	       {   
		  for(int ssb=0;ssb<selb_n;ssb++)
		    {
		       for(int ssa=0;ssa<seladd_n;ssa++)
			 {
			    for(int s=0;s<sys_n;s++)
			      {
				 for(int h=0;h<histname_jet_n;h++)
				   {
				      histNAMES[j][ss][se][ssb][ssa][h][s] =
					histname_jet[h]+flav[j]+"_"+sel[ss]+"_"+eta[se]+"_"+selb[ssb]+seladd[ssa]+sys[s];
				   }

				 for(int h=0;h<histname_jet_2d_n;h++)
				   {
				      histNAMES_2d[j][ss][se][ssb][ssa][h][s] =
					histname_jet_2d[h]+flav[j]+"_"+sel[ss]+"_"+eta[se]+"_"+selb[ssb]+seladd[ssa]+sys[s];
				   }

				 for(int h=0;h<histname_jet_3d_n;h++)
				   {
				      histNAMES_3d[j][ss][se][ssb][ssa][h][s] =
					histname_jet_3d[h]+flav[j]+"_"+sel[ss]+"_"+eta[se]+"_"+selb[ssb]+seladd[ssa]+sys[s];
				   }
			      }			    
			 }		       
		    }		  
	       }
	  }
     }

   if( !isdata )
     {	
	nPtRelPtBins = 15;
	TString PtRelPtBin[] = {
	   "Pt2030", "Pt3040", "Pt4050", "Pt5060", "Pt6070"
	   , "Pt7080", "Pt80100", "Pt100120", "Pt120160", "Pt160210"
	   , "Pt210260", "Pt260320", "Pt320400", "Pt400500", "Pt500" 
	};
	
	for( int ptb=0;ptb<nPtRelPtBins;ptb++ )
	  {      
	     for(int ib=0;ib<100;ib++ )
	       {
		  BTemplateCorrections[ib][ptb][0] = 1.;
		  BTemplateCorrections[ib][ptb][1] = 1.;
	       }
	     
	     std::ifstream MnusCorrectionsFile; 
	     MnusCorrectionsFile.open("./PtRelFall12/EnergyFraction_" + PtRelPtBin[ptb] + "_m5.txt");
	     while( MnusCorrectionsFile )
	       {
		  float xBin, efcorr;
		  MnusCorrectionsFile >> xBin >> efcorr;
		  if ( efcorr > 4. ) efcorr = 1.;
		  int ib = int(xBin/0.02);
		  BTemplateCorrections[ib][ptb][0] = efcorr;
	       }
	     
	     std::ifstream PlusCorrectionsFile; 
	     PlusCorrectionsFile.open("./PtRelFall12/EnergyFraction_" + PtRelPtBin[ptb] + "_p5.txt");
	     while( PlusCorrectionsFile )
	       {
		  float xBin, efcorr;
		  PlusCorrectionsFile >> xBin >> efcorr;
		  if ( efcorr > 4. ) efcorr = 1.;
		  int ib = int(xBin/0.02);
		  BTemplateCorrections[ib][ptb][1] = efcorr;
	       }
	  }
     }   
*/	
   std::cout << "Initialisation done" << std::endl;
}

void Hist::fill()
{/*
//   double w = noe/(xsec*eff);
//   if( !isdata ) w = (9238642/(1024000.0*0.0395))/w;

   double ran = rnd->Rndm();
//   if( ran >= 0.5 ) return;
////   if( ran < 0.5 ) return;
   
   double w = 1.;
   
   if( !isdata )
     {	
	double w015 = 1722681/(7.022E8*0.0039);
	double w020 = 8486904/(2.87E8*0.0065);
	double w030 = 9560265/(6.609E7*0.0122);
	double w050 =10365230/(8082000.0*0.0218);
	double w080 = 9238642/(1024000.0*0.0395);
	double w120 = 8501935/(157800.0*0.0473);
	double w170 = 7670288/(34020.0*0.0676);
	double w300 = 7839607/(1757.0*0.0864);
	double w470 = 3783069/(115.2*0.1024);
	double w600 = 4119000/(27.01*0.0996);
	double w800 = 4107853/(3.57*0.1033);
	double w1000= 3873970/(0.774*0.1097);

	if( year2011 )
	  {
	     w015 = 1901694/(580700000*0.0035);
	     w020 =12133278/(235400000*0.0045);
	     w030 =11610144/(53470000*0.01);
	     w050 = 9887094/(6317000*0.0245);
	     w080 =10497146/(781100*0.0345);
	     w120 = 7913447/(116000*0.0545);
	     w170 = 8116425/(24690*0.0625);
	     w300 = 7870022/(1153*0.0923);
	     w470 = 3812536/(71.81*0.086);
	     w600 = 4149919/(15.15*0.094);
	     w800 = 3945269/(1.863*0.1033);
	     w1000= 4133903/(0.3392*0.094);
	  }	
	
	w015 = w080 / w015;
	w020 = w080 / w020;
	w030 = w080 / w030;
	w050 = w080 / w050;
	w120 = w080 / w120;
	w170 = w080 / w170;
	w300 = w080 / w300;
	w470 = w080 / w470;
	w600 = w080 / w600;
	w800 = w080 / w800;
	w1000= w080 / w1000;
	w080 = 1.;
	
	if( pthat >=  15. && pthat <  20. ) w = w015;
	else if ( pthat >=  20. && pthat <  30. ) w = w020;
	else if ( pthat >=  30. && pthat <  50. ) w = w030;
	else if ( pthat >=  50. && pthat <  80. ) w = w050;
	else if ( pthat >=  80. && pthat < 120. ) w = w080;
	else if ( pthat >= 120. && pthat < 170. ) w = w120;
	else if ( pthat >= 170. && pthat < 300. ) w = w170;
	else if ( pthat >= 300. && pthat < 470. ) w = w300;
	else if ( pthat >= 470. && pthat < 600. ) w = w470;
	else if ( pthat >= 600. && pthat < 800. ) w = w600;
	else if ( pthat >= 800. && pthat < 1000.) w = w800;
	else if ( pthat >= 1000. )                w = w1000;	
     }   
   
   float mc_weight = ntP->mcweight;
   if( !isdata ) w *= mc_weight;
   
   int npu = ntP->nPU;
   int npv = ntP->nPV;
   if( year2011 )
     {
	if( npu > 34 ) npu = 34;
     }   
   else
     {
	if( ntP->Run < 0 ) npu = int(ntP->nPUtrue);
	if( npu > 49 ) npu = 49;
     }   
//   float wPU = WeightPU[npu];
//   if( !isdata ) w *= wPU;

//   int jidx = 0; // highest-Pt jet
   int jidx = -1; // all jets

   // Event
   
   if( ntP->nMuon == 0 ) return;

   int ntNJET = ntP->nJet;

   // 'overcut' for trigger selection
   int nJet20 = 0;
   for( int ij=0;ij<ntNJET;ij++ )
     {
	if( ntP->Jet_pt[ij] > 20. ) nJet20++;
     }	
   if( nJet20 < 2 ) return;
   
   bool filterPTHAT = false;
   bool filterNJET = false;
   bool filterPT = false;
   if( pthat > 0 && pthat < 170. )
     {
	for( int ij=0;ij<ntNJET;ij++ )
	  {
	     if ( ntP->Jet_pt[ij] > pthat * 7. ) filterPTHAT = true;
	  }	
     }   
   if( ntNJET > 15 ) filterNJET = true;
   for( int ij=0;ij<ntNJET;ij++ )
     {
	if( ntP->Jet_pt[ij] > 1600. ) filterPT = true;
     }   
   if( filterPTHAT || filterNJET || filterPT ) return;
   
   // Trigger
   
   bool passTrigJet20 = passTrig(35);
   bool passTrigJet40 = passTrig(41);
   bool passTrigJet70 = passTrig(44);
   bool passTrigJet110 = passTrig(47);
   bool passTrigJet300 = passTrig(50);
   bool passTrigJet300_HLT = passTrig(18);

   if( year2011 ) passTrigJet300 = passTrigJet300_HLT;
   
   int bitIdx;
   int trigIdx;
   
   int njet20 = 0, njet40 = 0, njet70 = 0,
     njet110 = 0, njet300 = 0;
   
   for(int ij=0;ij<ntNJET;ij++)
     {
	float jpt = ntP->Jet_pt[ij];
	float jeta = fabs(ntP->Jet_eta[ij]);

	if( jeta > 2.4 ) continue;
	
	if( jpt > 20. )  njet20++;
	if( jpt > 50. )  njet40++;
	if( jpt > 80. )  njet70++;
	if( jpt > 120. ) njet110++;
	if( jpt > 320. ) njet300++;
     }      
   if( njet20 == 0 ) passTrigJet20 = false;
   if( njet40 == 0 ) passTrigJet40 = false;
   if( njet70 == 0 ) passTrigJet70 = false;
   if( njet110 == 0 ) passTrigJet110 = false;
   if( njet300 == 0 ) passTrigJet300 = false;
   
   if( !passTrigJet20 &&
       !passTrigJet40 &&
       !passTrigJet70 &&
       !passTrigJet110 && 
       !passTrigJet300 ) return;
   
//   if( njet300 == 0 && year2011 && sp2011 ) return;
//   if( njet300 > 0 && year2011 && !sp2011 ) return;
   
   // Loop on jets
   
   int nJetMax = (jidx >= 0 && ntNJET >= jidx) ? jidx+1 : ntNJET;

   for(int ij=0;ij<nJetMax;ij++)
     {
////	if( ntP->Jet_Proba[ij] == 0 ) continue;
	
	ntrkgen = 0;
	
//	std::cout << "jet#" << ij << " " << ntP->nTrkInc << std::endl;
	
	v_jet = new TLorentzVector(0,0,0,0);
	v_jet_sys_jesTotalUp = new TLorentzVector(0,0,0,0);
	v_jet_sys_jesTotalLow = new TLorentzVector(0,0,0,0);
	v_jet_sys_jerTotalUp = new TLorentzVector(0,0,0,0);
	v_jet_sys_jerTotalLow = new TLorentzVector(0,0,0,0);
	v_mu = new TLorentzVector(0,0,0,0);

	bool passSel = true;
	
	float jeta = ntP->Jet_eta[ij];
	float jphi = ntP->Jet_phi[ij];
	float jm = ntP->Jet_mass[ij];
	float residual = 1.;
	if( ntP->Run > 0 ) residual = ntP->Jet_residual[ij];
	float jpt = ntP->Jet_pt[ij] * residual;

	if( isdata && year2011 && sp2011 && jpt < 320. ) continue;
	if( isdata && year2011 && !sp2011 && jpt >= 320. ) continue;
	
	// JER

	v_jet_sys_jerTotalUp->SetPtEtaPhiM(jpt,jeta,jphi,jm);
	v_jet_sys_jerTotalLow->SetPtEtaPhiM(jpt,jeta,jphi,jm);
	
	if( !isdata )
	  {	    
	     int etaIdx = -1;
	     if( fabs(jeta) >= 0. && fabs(jeta) < 0.5 ) etaIdx = 0;
	     if( fabs(jeta) >= 0.5 && fabs(jeta) < 1.1 ) etaIdx = 1;
	     if( fabs(jeta) >= 1.1 && fabs(jeta) < 1.7 ) etaIdx = 2;
	     if( fabs(jeta) >= 1.7 && fabs(jeta) < 2.3 ) etaIdx = 3;
	     if( fabs(jeta) >= 2.3 && fabs(jeta) < 5.0 ) etaIdx = 4;
	
	     double jpt_c = jpt;
	     if( etaIdx >= 0 )
	       {		  
		  double genpt = ntP->Jet_genpt[ij];
		  if( genpt >= 0. )
		    {
		       jpt_c = genpt + cJER[etaIdx]*(jpt-genpt);
		       
		       double jpt_c_down = genpt + cJER_down[etaIdx]*(jpt-genpt);
		       double jpt_c_up = genpt + cJER_up[etaIdx]*(jpt-genpt);
		       
		       if( jpt_c < 0. ) jpt_c = 0.;
		       if( jpt_c_down < 0. ) jpt_c_down = 0.;
		       if( jpt_c_up < 0. ) jpt_c_up = 0.;
		       
		       v_jet_sys_jerTotalUp->SetPtEtaPhiM(jpt_c_up,jeta,jphi,jm);
		       v_jet_sys_jerTotalLow->SetPtEtaPhiM(jpt_c_down,jeta,jphi,jm);
		  
//		  std::cout << jpt_c << "-" << jpt_c_down << "+" << jpt_c_up << std::endl;
		    }
	       }	     
	     
	     jpt = jpt_c;
	  }	
	     
	v_jet->SetPtEtaPhiM(jpt,jeta,jphi,jm);

	if( jpt < 20. || jpt > 1500. ) passSel = false;
	if( fabs(jeta) > 2.4 ) passSel = false;

//	if( ntP->Evt == 983779189 && ij == 1 )
//	  std::cout << "jet#" << ij << " pt=" << ntP->Jet_pt[ij] <<
//	  " eta=" << ntP->Jet_eta[ij] <<
//	  " passSel=" << passSel << std::endl;
	
	int njet = ntNJET;
	
	int nsv = ntP->nSV;
	int nsvj = ntP->Jet_nLastSV[ij];
	
	int ntrk = ntP->Jet_ntracks[ij];
	int nseltrk = ntP->Jet_nseltracks[ij];
////	if( ntrk == 0 ) passSel = false;
	
	v_jet_sys_jesTotalUp->SetPtEtaPhiM(jpt,jeta,jphi,jm);
	v_jet_sys_jesTotalLow->SetPtEtaPhiM(jpt,jeta,jphi,jm);

	// JES
	jesTotal->setJetPt(jpt);
	jesTotal->setJetEta(jeta);
	double uncert = jesTotal->getUncertainty(true);
	
	v_jet_sys_jesTotalUp->SetPtEtaPhiE(jpt*(1.+uncert),
					   jeta,
					   v_jet->Phi(),
					   v_jet->E()*(1.+uncert));
	
	v_jet_sys_jesTotalLow->SetPtEtaPhiE(jpt*(1.-uncert),
					    jeta,
					    v_jet->Phi(),
					    v_jet->E()*(1.-uncert));
	
	if( !sp2011 )
	  {	     
	     if( isdata )
	       {
		  if( jpt <= 50 && !passTrigJet20 ) passSel = false;
		  if( jpt > 50 && jpt <= 80 && !passTrigJet40 && !passTrigJet20 ) passSel = false;
		  if( jpt > 80 && jpt <= 120 && !passTrigJet70 && !passTrigJet40 ) passSel = false;
		  if( jpt > 120 && jpt <= 320 && !passTrigJet110 && !passTrigJet70 ) passSel = false;
		  if( jpt > 320 ) 
		    {
		       if( ntP->Run < 190000 && !passTrigJet110 && !passTrigJet70 ) passSel = false;
		       if( ntP->Run > 190000 && !passTrigJet300 ) passSel = false;
		    }	
	       }
	     else
	       {   
		  if( jpt <= 50 && !passTrigJet20 ) passSel = false;
		  if( jpt > 50 && jpt <= 80 && !passTrigJet40 ) passSel = false;
		  if( jpt > 80 && jpt <= 120 && !passTrigJet70 ) passSel = false;
		  if( jpt > 120 && !passTrigJet110 ) passSel = false;
	       }
	  }
	
	if( year2011 )
	  {
	     if( jpt >= 320 && !passTrigJet300 ) passSel = false; 
	  }	

	// muon-in-jet
//	std::vector<int> muidx;
	muidx.clear();

	for(int imu=0;imu<ntP->nMuon;imu++)
	  {
	     if( ntP->Muon_IdxJet[imu] == ij )
	       {	     
		  if( ntP->Muon_isGlobal[imu] == 0 ) continue;
		  if( ntP->Muon_pt[imu] < muptcut ) continue;
		  if( ntP->Muon_nMuHit[imu] == 0 ) continue;
		  if( ntP->Muon_nTkHit[imu] < 11 ) continue;
		  if( ntP->Muon_nPixHit[imu] < 2 ) continue;
		  if( ntP->Muon_nOutHit[imu] > 2 ) continue;
		  if( ntP->Muon_nMatched[imu] < 2 ) continue;
		  if( ntP->Muon_chi2[imu] > 10. ) continue;
		  if( ntP->Muon_chi2Tk[imu] > 10. ) continue;
		  if( fabs(ntP->Muon_eta[imu]) > 2.4 ) continue;
		  muidx.push_back(imu);
	       }
	  }
	bool passMuonInJet = (muidx.size() > 0);
	
	int jFlavour = ntP->Jet_flavour[ij];
	if( jFlavour >= 1 && jFlavour <= 3 ) jFlavour = 1;

	std::string flavch = "NOTFOUND";
	int flavchI = 0;
	if( jFlavour == 5 )  {flavch = "bjet"; flavchI = 1;}
	if( jFlavour == 4 )  {flavch = "cjet"; flavchI = 2;}
	if( jFlavour == 1 || jFlavour == 21 )  {flavch = "ljet"; flavchI = 3;}

	std::string selch = "NOTFOUND";
	std::string etach = "NOTFOUND";
	std::string selbch = "NOTFOUND";
   
	float drjmu = -666;
	if( passMuonInJet )
	  {
	     for(int im=0;im<muidx.size();im++)
	       {		  
		  float mupt = ntP->Muon_pt[muidx[im]];
		  float mueta = ntP->Muon_eta[muidx[im]];
		  float muphi = ntP->Muon_phi[muidx[im]];
		  float mum = 0.10566;
		  v_mu->SetPtEtaPhiM(mupt,mueta,muphi,mum);
		  drjmu = v_jet->DeltaR(*v_mu);
		  if( drjmu < 0.4 ) break;
	       }	     
	  }	
	
	if( drjmu > 0.4 || drjmu < 0 ) passSel = false;
	
	// END SELECTION
	
	if( passSel )
	  {
//	     if( fabs(jFlavour) == 5 && ntP->Jet_CombSvx[ij] >= 0.679 )
//	       _fevc << ntP->Evt << " " << ntP->Jet_Proba[ij]*10. << " " << jpt << " " << w*WeightPU[npu] << std::endl;
	     
	     bool selbch_res[selb_n];
	     for( int ibt=0;ibt<selb_n;ibt++ ) selbch_res[ibt] = false;
	     
	     float jbtCombSvx = ntP->Jet_CombSvx[ij];
	     float jbtTCHP = ntP->Jet_Ip3P[ij];
	     
	     selbch_res[0] = true;
	     if( jbtCombSvx >= 0.244 ) selbch_res[1] = true;
	     if( jbtCombSvx >= 0.679 ) selbch_res[2] = true;
	     if( jbtCombSvx >= 0.898 ) selbch_res[3] = true;
	     if( jbtTCHP >= 1.19 ) selbch_res[4] = true;
	     if( jbtTCHP >= 1.93 ) selbch_res[5] = true;
	     if( jbtTCHP >= 3.41 ) selbch_res[6] = true;
	     if( jbtCombSvx < 0.898 ) selbch_res[7] = true;
//	     if( jbtTCHP < 3.41 ) selbch_res[4] = true;
	     
//	     if( ntP->Jet_Proba[ij] == 0. && ntP->Jet_CombSvx[ij] > 0.244 )
//	       std::cout << ntP->Jet_CombSvx[ij] << std::endl;

	     bool selach_res[seladd_n];
	     for( int ibt=0;ibt<seladd_n;ibt++ ) selach_res[ibt] = false;

	     selach_res[0] = true;
//	     if( ntrk >= 1 ) selach_res[1] = true;
//	     if( ntrk >= 2 ) selach_res[2] = true;
//	     if( ntrk >= 4 ) selach_res[3] = true;
//	     if( ntrk >= 7 ) selach_res[4] = true;
//	     if( ntrk >= 10 ) selach_res[5] = true;
	     
	     for(int isys=0;isys<sys_n;isys++)
	       {	
		  float jPT = getPt(sys[isys]);
		  
		  double sfj = w;
		  
		  // l-jets are underestimated in MC by 1.27
		  if( flavchI == 3 ) sfj *= 1.27;
		  
		  // additional scaling to study the effect on SF
//		  if( flavchI == 1 ) sfj *= 0.80;
//		  if( flavchI == 1 ) sfj *= 1.20;
//		  if( flavchI == 2 ) sfj *= 0.80;
//		  if( flavchI == 2 ) sfj *= 1.20;
//		  if( flavchI == 3 ) sfj *= 0.80;
//		  if( flavchI == 3 ) sfj *= 1.20;
		  
		  // pileup reweighting sys
		  float wPU = WeightPU[npu];
		  if( sys[isys] == "_pu_up" ) wPU = WeightPUmax[npu];
		  if( sys[isys] == "_pu_low" ) wPU = WeightPUmin[npu];
		  if( !isdata ) sfj *= wPU;

		  // gluon splitting sys
		  if( (sys[isys] == "_gsplit_low" || sys[isys] == "_gsplit_up") &&
		      !isdata )
		    {		       
		       bool GSPc = false, GSPb = false;
		       
		       if( flavch == "cjet" )
			 {
			    for( int k=0;k<ntP->ncQuarks;k++ )
			      {
				 if( ntP->cQuark_status[k] != 2 ) continue;
				 double dRqj = DeltaR(ntP->Jet_eta[ij],
						      ntP->Jet_phi[ij], 
						      ntP->cQuark_eta[k], 
						      ntP->cQuark_phi[k]);
				 if( k == ntP->ncQuarks-1 ) continue;
				 if( dRqj > 0.5 ) continue;
				 for( int l=k+1;l<ntP->ncQuarks;l++ )
				   {
				      if( ntP->cQuark_status[l] != 2 ) continue;
				      double dRqj2 = DeltaR(ntP->Jet_eta[ij], 
							    ntP->Jet_phi[ij], 
							    ntP->cQuark_eta[l], 
							    ntP->cQuark_phi[l]);
				      if( ntP->cQuark_pdgID[k] * ntP->cQuark_pdgID[l] > 0 ) continue;
				      if( dRqj2 < 0.5 ) GSPc = true;
				   }			    
			      }		       
			 }		  
		       
		       if( flavch == "bjet" )
			 {
			    for( int k=0;k<ntP->nbQuarks;k++ )
			      {
				 if( ntP->bQuark_status[k] != 2 ) continue;
				 double dRqj = DeltaR(ntP->Jet_eta[ij],
						      ntP->Jet_phi[ij], 
						      ntP->bQuark_eta[k], 
						      ntP->bQuark_phi[k]);
				 if( k == ntP->nbQuarks-1 ) continue;
				 if( dRqj > 0.5 ) continue;
				 for( int l=k+1;l<ntP->nbQuarks;l++ )
				   {
				      if( ntP->bQuark_status[l] != 2 ) continue;
				      double dRqj2 = DeltaR(ntP->Jet_eta[ij], 
							    ntP->Jet_phi[ij], 
							    ntP->bQuark_eta[l], 
							    ntP->bQuark_phi[l]);
				      if( ntP->bQuark_pdgID[k] * ntP->bQuark_pdgID[l] > 0 ) continue;
				      if( dRqj2 < 0.5 ) GSPb = true;
				   }			    
			      }		       
			 }		  
		  
		       if( 
			  sys[isys] == "_gsplit_low" && 
			  ((GSPc && flavch == "cjet") || (GSPb && flavch == "bjet"))
			 ) sfj *= 0.5;
		       if( 
			  sys[isys] == "_gsplit_up" &&
			  ((GSPc && flavch == "cjet") || (GSPb && flavch == "bjet"))
			 ) sfj *= 1.5;
		    }
		  
		  // b fragmentation sys
		  if( (sys[isys] == "_bfrag_low" || sys[isys] == "_bfrag_up") &&
		      !isdata )
		    {		       
		       float WeightBFrag = 1.;
		       float EnergyFraction = 0.; 
		       int iB = -1, iptBin = 0, efbin = -1;
		       if( flavch == "bjet" && ( sys[isys] == "_bfrag_low" || sys[isys] == "_bfrag_up" ) )
			 {
			    if( jPT > 500 ) iptBin = 14;
			    else if( jPT > 400 ) iptBin = 13;
			    else if( jPT > 320 ) iptBin = 12;  
			    else if( jPT > 260 ) iptBin = 11;
			    else if( jPT > 210 ) iptBin = 10;
			    else if( jPT > 160 ) iptBin =  9;
			    else if( jPT > 120 ) iptBin =  8;  
			    else if( jPT > 100 ) iptBin =  7;  
			    else if( jPT >  80 ) iptBin =  6;  
			    else if( jPT >  70 ) iptBin =  5;  
			    else if( jPT >  60 ) iptBin =  4;  
			    else if( jPT >  50 ) iptBin =  3;  
			    else if( jPT >  40 ) iptBin =  2;  
			    else if( jPT >  30 ) iptBin =  1;  
			    else                 iptBin =  0;
			    
			    float B_Mass = 0.;
			    for( int ib=0;ib<ntP->nBHadrons;ib++ )
			      {
				 if( DeltaR(ntP->Jet_eta[ij],
					    ntP->Jet_phi[ij], 
					    ntP->BHadron_eta[ib], 
					    ntP->BHadron_phi[ib]) < 0.5 )
				   {
				      if( ntP->BHadron_mass[ib] > B_Mass ) 
					{
					   B_Mass = ntP->BHadron_mass[ib];
					   iB = ib;
					}				 
				   }			    
			      }
			    
			    if( iB >= 0 ) 
			      {
				 EnergyFraction = ntP->BHadron_pT[iB]/ntP->Jet_genpt[ij];
				 efbin = int( EnergyFraction / 0.02 );
				 if( efbin >= 0 && efbin < 100 ) 
				   {
				      if( sys[isys] == "_bfrag_low" ) WeightBFrag = BTemplateCorrections[efbin][iptBin][0];
				      if( sys[isys] == "_bfrag_up" ) WeightBFrag = BTemplateCorrections[efbin][iptBin][1];
				   }			    
			      }
			    
			    sfj *= WeightBFrag;
			 }
		    }		  
		  
		  // c->D fragmentation sys
		  if( (sys[isys] == "_cdfrag_low" || sys[isys] == "_cdfrag_up") &&
		      !isdata )
		    {		       
		       if( flavch == "cjet" || flavch == "bjet" )
			 {
			    bool isDplusMu = false, isDzeroMu = false, isDsubsMu = false;
			    
			    int ndaughters = 0;
			    for( int k=0;k<ntP->nDHadrons;k++ )
			      {
				 double dR = DeltaR(ntP->DHadron_eta[k], 
						    ntP->DHadron_phi[k], 
						    ntP->Jet_eta[ij], 
						    ntP->Jet_phi[ij]);
				 if( dR > 0.4 ) continue;
				 bool isSemiMu = false;
				 int nd = ntP->DHadron_nDaughters[k];
				 for( int kk=0;kk<nd;kk++ )
				   {
				      if( abs(ntP->DHadron_DaughtersPdgID[kk+ndaughters]) == 13 ) isSemiMu = true;
				   }
			    
				 ndaughters += nd;
				 
				 if( !isSemiMu ) continue;
				 if( abs(ntP->DHadron_pdgID[k]) == 411 ) isDplusMu = true;
				 if( abs(ntP->DHadron_pdgID[k]) == 421 ) isDzeroMu = true;
				 if( abs(ntP->DHadron_pdgID[k]) == 431 ) isDsubsMu = true;
			      }
			    
			    // weight for D->mu decay: Pythia vs PDG2013
			    if( sys[isys] == "_cdfrag_low" )
			      {			    
				 if( isDplusMu ) sfj *= 0.176 / 0.172;
				 if( isDzeroMu ) sfj *= 0.067 / 0.077;
				 if( isDsubsMu ) sfj *= 0.067 / 0.080;
			      }		       
			 }
		    }		  
		  
		  // c fragmentation sys
		  int nChargedGen1_c = -1, nChargedGen2_c = -1;
		  if( !isdata )
		    {
		       if( flavch == "cjet" )
			 {
			    bool hasCquark = 0;
			    for( int c=0;c<ntP->ncQuarks;c++ )
			      {
				 double dRc = DeltaR(ntP->cQuark_eta[c], 
						     ntP->cQuark_phi[c], 
						     ntP->Jet_eta[ij], 
						     ntP->Jet_phi[ij]);
				 if( dRc < 0.3 ) hasCquark = 1;
			      }
			    
			    if( hasCquark )
			      {				 
				 bool isDplus = false, isDzero = false, isDsubs = false, isDbary = false;
				 for( int k=0;k<ntP->nDHadrons;k++ )
				   {
				      double dR = DeltaR(ntP->DHadron_eta[k], 
							 ntP->DHadron_phi[k], 
							 ntP->Jet_eta[ij], 
							 ntP->Jet_phi[ij]);
				      if( dR > 0.4 ) continue;
				 
				      if( nChargedGen1_c >= 0 && nChargedGen2_c < 0 ) nChargedGen2_c = ntP->DHadron_nCharged[k];
				      if( nChargedGen1_c < 0 )  nChargedGen1_c = ntP->DHadron_nCharged[k];
				      
				      if( abs(ntP->DHadron_pdgID[k]) == 411 ) isDplus = true;
				      if( abs(ntP->DHadron_pdgID[k]) == 421 ) isDzero = true;
				      if( abs(ntP->DHadron_pdgID[k]) == 431 ) isDsubs = true;
				      if((abs(ntP->DHadron_pdgID[k])/1000)%10 == 4 ) isDbary = true;
				   }		       
				 
				 // weight for c->D fragmentation rate: Pythia vs PDG2008
				 if( sys[isys] == "_cfrag_low" )
				   {
				      // DB
				      if( isDplus ) sfj *= 1.37; // PDG2008(0.246+-0.020)
				      if( isDzero ) sfj *= 0.91; // PDG2008(0.565+-0.032)
				      if( isDsubs ) sfj *= 0.67; // PDG2008(0.080+-0.017)
				      // 0.185072, 0.58923, 0.115961
				   }		       
			      }			    
			 }
		    }		  
		  
		  // K0s Lambda sys
		  if( (sys[isys] == "_ksl_low" || sys[isys] == "_ksl_up") &&
		      !isdata )
		    {		       
		       int nK0s = 0, nLambda = 0;
		       if( flavch == "ljet" )
			 {
			    for( int k=0;k<ntP->nGenV0;k++ )
			      {
				 double dR = DeltaR(ntP->GenV0_eta[k], 
						    ntP->GenV0_phi[k], 
						    ntP->Jet_eta[ij], 
						    ntP->Jet_phi[ij]);
				 if( dR > 0.3 ) continue;
				 int pdgid = abs(ntP->GenV0_pdgID[k]);
				 if( pdgid == 310 )  nK0s++;
				 if( pdgid == 3122 ) nLambda++;
			      }
			    if( sys[isys] == "_ksl_up" )
			      {			    
				 if( nK0s > 0 )    sfj *= 1.3;
				 if( nLambda > 0 ) sfj *= 1.5;
			      }
			 }
		    }		  
		  
		  // generated ntrk in b- and c- hadron decays sys
		  if( !isdata )
		    {		       		  
		       // Weight in number of generated tracks from c-hadrons (+_3%)
		       float WeightCtrkMin[] = 
			 {
			    1.83175, 1.09414, 0.969433, 0.957613, 0.944791, 0.93478, 0.951504, 0.933117, 0.933552, 1,
			    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
			 };
		       float WeightCtrkMax[] = 
			 {
			    0.941931, 0.935278, 0.959889, 1.05731, 1.11384, 1.28445, 1, 1, 1, 1,
			    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
			 };
		       // Weight in number of generated tracks from b-hadrons (+_3%)
		       float WeightBtrkMin[] = 
			 {
			    3.4573, 1.70516, 1.12162, 1.05566, 0.945077, 0.996961, 0.925635, 0.987391, 0.917369, 0.969089,
			    0.904673, 0.947409, 0.894847, 0.939174, 0.887052, 1, 1, 1, 1, 1
			 };
		       float WeightBtrkMax[] = 
			 {
			    1.06921, 0.876405, 0.89011, 0.937056, 0.960934, 1.09521, 1.00306, 1.17408, 1.01393, 1.22725,
			    1.04056, 1.35842, 1.08842, 1.55228, 1.11446, 1, 1, 1, 1, 1
			 };

		       int nChargedGen1_b = -1, nChargedGen2_b = -1;
		       
		       if( flavch == "bjet" )
			 {
			    for( int k=0;k<ntP->nBHadrons;k++ )
			      {
				 if( ntP->BHadron_hasBdaughter[k] == 1 ) continue;
				 double dR = DeltaR(ntP->BHadron_eta[k], 
						    ntP->BHadron_phi[k], 
						    ntP->Jet_eta[ij], 
						    ntP->Jet_phi[ij]);
				 if( dR > 0.4 ) continue;
				 int k1 = ntP->BHadron_DHadron1[k];
				 int k2 = ntP->BHadron_DHadron2[k];
				 if( nChargedGen1_b >= 0 && nChargedGen2_b < 0 ) 
				   {
				      nChargedGen2_b  = ntP->BHadron_nCharged[k];
				      if( k1 >= 0 ) nChargedGen2_b += ntP->DHadron_nCharged[k1];
				      if( k2 >= 0 ) nChargedGen2_b += ntP->DHadron_nCharged[k2];
				   }
				 if( nChargedGen1_b < 0 ) 
				   {
				      nChargedGen1_b  = ntP->BHadron_nCharged[k];
				      if( k1 >= 0 ) nChargedGen1_b += ntP->DHadron_nCharged[k1];
				      if( k2 >= 0 ) nChargedGen1_b += ntP->DHadron_nCharged[k2];
				   }
			      }
			 }
		       
		       // Weight in track multiplicity from c-hadron and b-hadron decays
		       int nchgen1_b = nChargedGen1_b, nchgen2_b = nChargedGen2_b;
		       if( nchgen1_b > 19 ) nchgen1_b = 19;
		       if( nchgen2_b > 19 ) nchgen2_b = 19;

		       int nchgen1_c = nChargedGen1_c, nchgen2_c = nChargedGen2_c;
		       if( nchgen1_c > 19 ) nchgen1_c = 19;
		       if( nchgen2_c > 19 ) nchgen2_c = 19;
		       
		       if( sys[isys] == "_ntrkgen_low" )
			 {
			    if( flavch == "cjet" )
			      {
				 if( nchgen1_c >= 0 ) sfj *= WeightCtrkMin[nchgen1_c];
				 if( nchgen2_c >= 0 ) sfj *= WeightCtrkMin[nchgen2_c];
			      }
			    if( flavch == "bjet" )
			      {
				 if( nchgen1_b >= 0 ) sfj *= WeightBtrkMin[nchgen1_b];
				 if( nchgen2_b >= 0 ) sfj *= WeightBtrkMin[nchgen2_b];
			      }
			 }
		       if( sys[isys] == "_ntrkgen_up" )
			 {
			    if( flavch == "cjet" )
			      {			    
				 if( nchgen1_c >= 0 ) sfj *= WeightCtrkMax[nchgen1_c];
				 if( nchgen2_c >= 0 ) sfj *= WeightCtrkMax[nchgen2_c];
			      }		       
			    if( flavch == "bjet" )
			      {			    
				 if( nchgen1_b >= 0 ) sfj *= WeightBtrkMax[nchgen1_b];
				 if( nchgen2_b >= 0 ) sfj *= WeightBtrkMax[nchgen2_b];
			      }
			 }
		       
		       if( flavch == "bjet" ) ntrkgen = nChargedGen1_b;
		       else if( flavch == "cjet" ) ntrkgen = nChargedGen1_c;		       
		    } 
		  // end sys
		  
		  int jPtBin = 0;   
		  if( jPT >= jPtMin[0] && jPT <= jPtMax[0] )   {selch = "pt20t30"; jPtBin = 1;}
		  if( jPT > jPtMin[1] && jPT <= jPtMax[1] )   {selch = "pt30t40"; jPtBin = 2;}
		  if( jPT > jPtMin[2] && jPT <= jPtMax[2] )   {selch = "pt40t50"; jPtBin = 3;}
		  if( jPT > jPtMin[3] && jPT <= jPtMax[3] )   {selch = "pt50t60"; jPtBin = 4;}
		  if( jPT > jPtMin[4] && jPT <= jPtMax[4] )   {selch = "pt60t70"; jPtBin = 5;}
		  if( jPT > jPtMin[5] && jPT <= jPtMax[5] )   {selch = "pt70t80"; jPtBin = 6;}
		  if( jPT > jPtMin[6] && jPT <= jPtMax[6] )   {selch = "pt80t100"; jPtBin = 7;}
		  if( jPT > jPtMin[7] && jPT <= jPtMax[7] )  {selch = "pt100t120"; jPtBin = 8;}
		  if( jPT > jPtMin[8] && jPT <= jPtMax[8] )  {selch = "pt120t160"; jPtBin = 9;}
		  if( jPT > jPtMin[9] && jPT <= jPtMax[9] )  {selch = "pt160t210"; jPtBin = 10;}
		  if( jPT > jPtMin[10] && jPT <= jPtMax[10] )  {selch = "pt210t260"; jPtBin = 11;}
		  if( jPT > jPtMin[11] && jPT <= jPtMax[11] )  {selch = "pt260t320"; jPtBin = 12;}
		  if( jPT > jPtMin[12] && jPT <= jPtMax[12] )  {selch = "pt320t400"; jPtBin = 13;}
		  if( jPT > jPtMin[13] && jPT <= jPtMax[13] )  {selch = "pt400t500"; jPtBin = 14;}
		  if( jPT > jPtMin[14] && jPT <= jPtMax[14] )  {selch = "pt500t600"; jPtBin = 15;}
		  if( jPT > jPtMin[15] && jPT <= jPtMax[15] )  {selch = "pt600t800"; jPtBin = 16;}
		  if( jPT > jPtMin[16] && jPT <= jPtMax[16] )  {selch = "pt800t1000"; jPtBin = 17;}
		  if( jPT > jPtMin[17] && jPT <= jPtMax[17] )  {selch = "pt1000t1400"; jPtBin = 18;}
//		  if( jPT > jPtMin[17] && jPT <= jPtMax[17] )  {selch = "pt20t25"; jPtBin = 18;}
//		  if( jPT > jPtMin[18] && jPT <= jPtMax[18] )  {selch = "pt25t30"; jPtBin = 19;}
		  double jPtMinCur = 0.;
		  double jPtMaxCur = 1250.;
		  
		  if( jPtBin > 0 )
		    {		       
		       jPtMinCur = jPtMin[jPtBin-1];
		       jPtMaxCur = jPtMax[jPtBin-1];
		    }		  

		  double mupt = -1.;
		  if( muidx.size() > 0 )
		    mupt = ntP->Muon_pt[muidx[0]];
		  
		  bool jeta_res[eta_n];
		  for( int ibt=0;ibt<eta_n;ibt++ ) jeta_res[ibt] = false;
		  		  
		  jeta_res[0] = true;
//		  if( fabs(jeta) >= 0.0 && fabs(jeta) <= 0.6 )  jetaI = 1;
//		  if( fabs(jeta) > 0.6 && fabs(jeta) <= 1.2 )   jetaI = 2;
//		  if( fabs(jeta) > 1.2 && fabs(jeta) <= 1.8 )   jetaI = 3;
//		  if( fabs(jeta) > 1.8 && fabs(jeta) <= 2.5 )   jetaI = 4;

		  // apply RW
		  double rwSF = 1.;
		  double rwSF_btag = rwSF;
		  double sfj_btag = sfj;
		  if( doRW == 1 && !isdata )
		    {
		       for(int iv=0;iv<rw_n;iv++)
			 {			    
			    std::vector<std::pair<double,double> > rwBinFound =
			      rwBin[iv][jPtBin];
			    std::vector<std::pair<double,double> > rwSfFound =
			      rwSf[iv][jPtBin];

			    std::vector<std::pair<double,double> > rwBinFound_btag =
			      rwBin_btag[iv][jPtBin];
			    std::vector<std::pair<double,double> > rwSfFound_btag =
			      rwSf_btag[iv][jPtBin];
			    
			    if( strcmp(name_rw[iv].c_str(),"1d_ntrk") == 0 )
			      {		
				 rwSF *= getSF(rwBinFound,rwSfFound,ntrk);
				 rwSF_btag *= getSF(rwBinFound_btag,rwSfFound_btag,ntrk);
			      }
			    if( strcmp(name_rw[iv].c_str(),"1d_nseltrk") == 0 )
			      {		
				 rwSF *= getSF(rwBinFound,rwSfFound,nseltrk);
				 rwSF_btag *= getSF(rwBinFound_btag,rwSfFound_btag,nseltrk);
			      }
			    if( strcmp(name_rw[iv].c_str(),"1d_njet") == 0 )
			      {		
				 rwSF *= getSF(rwBinFound,rwSfFound,njet);
				 rwSF_btag *= getSF(rwBinFound_btag,rwSfFound_btag,njet);
			      }
			    if( strcmp(name_rw[iv].c_str(),"1d_nsv") == 0 )
			      {		
				 rwSF *= getSF(rwBinFound,rwSfFound,nsv);
				 rwSF_btag *= getSF(rwBinFound_btag,rwSfFound_btag,nsv);
			      }
			    if( strcmp(name_rw[iv].c_str(),"1d_mupt") == 0 )
			      {
				 if( muidx.size() > 0 )
				   {				      
				      rwSF *= getSF(rwBinFound,rwSfFound,mupt);
				      rwSF_btag *= getSF(rwBinFound_btag,rwSfFound_btag,mupt);
				   }				 
			      }
			    if( strcmp(name_rw[iv].c_str(),"1d_eta") == 0 )
			      {		
				 double jeta = ntP->Jet_eta[ij];
				 rwSF *= getSF(rwBinFound,rwSfFound,jeta);
				 rwSF_btag *= getSF(rwBinFound_btag,rwSfFound_btag,jeta);
			      }
			    if( strcmp(name_rw[iv].c_str(),"1d_pt") == 0 )
			      {		
				 double jptmax = 1250.;
				 double jptmin = 0.;
				 if( jPtBin > 0 ) 
				   {
				      jptmax = jPtMax[jPtBin-1];
				      jptmin = jPtMin[jPtBin-1];
				   }	     
				 double vPT = (jPT-jptmin)/(jptmax-jptmin);
				 
				 rwSF *= getSF(rwBinFound,rwSfFound,vPT);
				 rwSF_btag *= getSF(rwBinFound_btag,rwSfFound_btag,vPT);
			      }
			 }		       
		    }
		  if( doRW == 2 && !isdata )
		    {
		       std::vector<std::pair<double,double> > rw2DBinXFound =
			 rw2DBinX[0][jPtBin];
		       std::vector<std::pair<double,double> > rw2DBinYFound =
			 rw2DBinY[0][jPtBin];
		       std::vector<std::pair<double,double> > rw2DSfFound =
			 rw2DSf[0][jPtBin];

		       if( 
			  strcmp(name_rw_2d.c_str(),"2d_ntrk_vs_mupt") == 0
			 )
			 {
			    rwSF *= get2DSF(rw2DBinXFound,rw2DBinYFound,
					    rw2DSfFound,
					    mupt,
					    ntrk);
			 }
		       else if(
			       strcmp(name_rw_2d.c_str(),"2d_nseltrk_vs_njet") == 0
			      )
			 {
			    rwSF *= get2DSF(rw2DBinXFound,rw2DBinYFound,
					    rw2DSfFound,
					    njet,
					    nseltrk);
			 }
		    }
		  if( doRW == 3 && !isdata )
		    {
		       std::vector<std::pair<double,double> > rw3DBinXFound =
			 rw3DBinX[0][jPtBin];
		       std::vector<std::pair<double,double> > rw3DBinYFound =
			 rw3DBinY[0][jPtBin];
		       std::vector<std::pair<double,double> > rw3DBinZFound =
			 rw3DBinZ[0][jPtBin];
		       std::vector<std::pair<double,double> > rw3DSfFound =
			 rw3DSf[0][jPtBin];

		       if( 
			  strcmp(name_rw_3d.c_str(),"3d_ntrk_pt_mupt") == 0
			 )
			 {
			    rwSF *= get3DSF(rw3DBinXFound,rw3DBinYFound,rw3DBinZFound,
					    rw3DSfFound,
					    mupt,
					    (jPT-jPtMinCur)/(jPtMaxCur-jPtMinCur),
					    ntrk);
			 }
		    }
		  sfj *= rwSF;
		  sfj_btag *= rwSF_btag;
		  
		  // 1d
		  std::vector<std::string> histNAMESSEL;
		  std::vector<int> histSYS;
		  std::vector<int> histVAR;
		  std::vector<double> sfarr;
		  sfarr.clear();
		  for(int ih=0;ih<histname_jet_n;ih++)
		    {
		       for(int ihb=0;ihb<selb_n;ihb++)
			 {
			    if( !selbch_res[ihb] ) continue;

			    for(int iha=0;iha<seladd_n;iha++)
			      {
				 if( !selach_res[iha] ) continue;
			    
				 for(int ihe=0;ihe<eta_n;ihe++)
				   {
				      if( !jeta_res[ihe] ) continue;
				      
				      histNAMESSEL.push_back(histNAMES[flavchI][jPtBin][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL.push_back(histNAMES[0][jPtBin][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL.push_back(histNAMES[flavchI][0][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL.push_back(histNAMES[0][0][ihe][ihb][iha][ih][isys]);
				      histSYS.push_back(isys);
				      histVAR.push_back(ih);
				   }
			      }
			 }		       
		    }
		  
		  int nHISTSEL = histNAMESSEL.size();
		  for(int ih=0;ih<nHISTSEL;ih++)
		    {		       
		       TH1D *h = _m1d_Jet->find(histNAMESSEL.at(ih))->second;

		       int hidx = int(ih/4);
		       fillThis = true;
		       float var;
		       std::vector<float> varv;
		       bool single = 1;
		       std::string varName = histname_jet[histVAR[hidx]];
		       if( strcmp(varName.c_str(),"h_j1_svntrk_") == 0 ||
			   strcmp(varName.c_str(),"h_j1_svmass_") == 0 )
			 {			    
			    varv = getVarVec(sys[histSYS[hidx]],ij,varName,jPtBin);
			    single = 0;
			 }		       
		       else
			 var = getVar(sys[histSYS[hidx]],ij,varName,jPtBin);
		       
		       if( fillThis )
			 {			    
			    //			 h->Fill(var,sfarr[ih]);
			    if( single ) 
			      {
				 h->Fill(var,sfj);
			      }
			    else 
			      {
				 int varvs = varv.size();
				 for(int iv=0;iv<varvs;iv++)
				   {
				      h->Fill(varv[iv],sfj);
				   }				 
			      }			    
			 }		       
		    }
		  histNAMESSEL.clear();

		  // 2d
		  std::vector<std::string> histNAMESSEL_2d;
		  std::vector<int> histSYS_2d;
		  std::vector<int> histVAR_2d;
		  for(int ih=0;ih<histname_jet_2d_n;ih++)
		    {
		       for(int ihb=0;ihb<selb_n;ihb++)
			 {
			    if( !selbch_res[ihb] ) continue;

			    for(int iha=0;iha<seladd_n;iha++)
			      {
				 if( !selach_res[iha] ) continue;
			    
				 for(int ihe=0;ihe<eta_n;ihe++)
				   {
				      if( !jeta_res[ihe] ) continue;
				      
				      histNAMESSEL_2d.push_back(histNAMES_2d[flavchI][jPtBin][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL_2d.push_back(histNAMES_2d[0][jPtBin][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL_2d.push_back(histNAMES_2d[flavchI][0][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL_2d.push_back(histNAMES_2d[0][0][ihe][ihb][iha][ih][isys]);
				      histSYS_2d.push_back(isys);
				      histVAR_2d.push_back(ih);
				   }				 
			      }
			 }		       
		    }

		  int nHISTSEL_2d = histNAMESSEL_2d.size();
		  for(int ih=0;ih<nHISTSEL_2d;ih++)
		    {		       
		       TH2D *h_2d = _m2d_Jet->find(histNAMESSEL_2d.at(ih))->second;
		       
		       int hidx = int(ih/4);
		       fillThis = true;
		       std::pair<float,float> var = getVar2d(sys[histSYS_2d[hidx]],ij,histname_jet_2d[histVAR_2d[hidx]],jPtBin);
		       float varX = var.first;
		       float varY = var.second;
		       
		       if( fillThis )
			 h_2d->Fill(varX,varY,sfj);
		    }		  
		  histNAMESSEL_2d.clear();
		  
		  // 3d
		  std::vector<std::string> histNAMESSEL_3d;
		  std::vector<int> histSYS_3d;
		  std::vector<int> histVAR_3d;
		  for(int ih=0;ih<histname_jet_3d_n;ih++)
		    {
		       for(int ihb=0;ihb<selb_n;ihb++)
			 {
			    if( !selbch_res[ihb] ) continue;

			    for(int iha=0;iha<seladd_n;iha++)
			      {
				 if( !selach_res[iha] ) continue;
			    
				 for(int ihe=0;ihe<eta_n;ihe++)
				   {
				      if( !jeta_res[ihe] ) continue;
				      
				      histNAMESSEL_3d.push_back(histNAMES_3d[flavchI][jPtBin][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL_3d.push_back(histNAMES_3d[0][jPtBin][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL_3d.push_back(histNAMES_3d[flavchI][0][ihe][ihb][iha][ih][isys]);
				      histNAMESSEL_3d.push_back(histNAMES_3d[0][0][ihe][ihb][iha][ih][isys]);
				      histSYS_3d.push_back(isys);
				      histVAR_3d.push_back(ih);
				   }				 
			      }
			 }		       
		    }

		  int nHISTSEL_3d = histNAMESSEL_3d.size();
		  for(int ih=0;ih<nHISTSEL_3d;ih++)
		    {		       
		       TH3D *h_3d = _m3d_Jet->find(histNAMESSEL_3d.at(ih))->second;
		       
		       int hidx = int(ih/4);
		       fillThis = true;
		       std::vector<float> var = getVar3d(sys[histSYS_3d[hidx]],ij,histname_jet_3d[histVAR_3d[hidx]],jPtBin);
		       float varX = var[0];
		       float varY = var[1];
		       float varZ = var[2];
		       
		       if( fillThis )
			 h_3d->Fill(varX,varY,varZ,sfj);
		    }		  
		  histNAMESSEL_3d.clear();		  
	       }
	  }
	
	delete v_mu;
	delete v_jet;
	delete v_jet_sys_jesTotalUp;
	delete v_jet_sys_jesTotalLow;	
	delete v_jet_sys_jerTotalUp;
	delete v_jet_sys_jerTotalLow;	
     }
*/   
}

void Hist::close()
{
   _fout->Write();
   _fout->Close();
//   _fevc.close();
   _fevcVal.close();
   
//   delete rnd;
}

bool Hist::printout(bool doPrint)
{
   m_weight = -666;
   
   m_channel = -666;
   m_sel = -666;
   
   m_l1_type = -666;
   m_l1_charge = -666;
   m_l1_pt = -666;
   m_l1_eta = -666;
   m_l1_phi = -666;
   m_l1_m = -666;
   
   m_l2_type = -666;
   m_l2_charge = -666;
   m_l2_pt = -666;
   m_l2_eta = -666;
   m_l2_phi = -666;
   m_l2_m = -666;

   m_l3_type = -666;
   m_l3_charge = -666;
   m_l3_pt = -666;
   m_l3_eta = -666;
   m_l3_phi = -666;
   m_l3_m = -666;

   m_l4_type = -666;
   m_l4_charge = -666;
   m_l4_pt = -666;
   m_l4_eta = -666;
   m_l4_phi = -666;
   m_l4_m = -666;
   
   m_j1_ntrk = -666;
   m_j1_btag = -666;
   m_j1_pt = -666;
   m_j1_eta = -666;
   m_j1_phi = -666;
   m_j1_m = -666;

   m_j2_ntrk = -666;
   m_j2_btag = -666;
   m_j2_pt = -666;
   m_j2_eta = -666;
   m_j2_phi = -666;
   m_j2_m = -666;

   m_j3_ntrk = -666;
   m_j3_btag = -666;
   m_j3_pt = -666;
   m_j3_eta = -666;
   m_j3_phi = -666;
   m_j3_m = -666;

   m_j4_ntrk = -666;
   m_j4_btag = -666;
   m_j4_pt = -666;
   m_j4_eta = -666;
   m_j4_phi = -666;
   m_j4_m = -666;

   m_j5_ntrk = -666;
   m_j5_btag = -666;
   m_j5_pt = -666;
   m_j5_eta = -666;
   m_j5_phi = -666;
   m_j5_m = -666;

   m_j6_ntrk = -666;
   m_j6_btag = -666;
   m_j6_pt = -666;
   m_j6_eta = -666;
   m_j6_phi = -666;
   m_j6_m = -666;
   
   m_met_phi = -666;
   m_met = -666;
   
   _v_ElectronTight->clear();
   _v_MuonTight->clear();
   _v_JetTight->clear();
   
   for(int i=0;i<_v_Electron->size();i++)
     {
	if( !_v_Electron->at(i).passPtEta() ) continue;
	if( _v_Electron->at(i).isLoose() )
	  _v_ElectronTight->push_back(_v_Electron->at(i));
     }       

   for(int i=0;i<_v_Muon->size();i++)
     {
	if( !_v_Muon->at(i).passPtEta() ) continue;
	if( _v_Muon->at(i).isLoose() )
	  _v_MuonTight->push_back(_v_Muon->at(i));
     }      

   for(int i=0;i<_v_Jet->size();i++)
     {
	if( _v_Jet->at(i).isTight() && _v_Jet->at(i).pt() > 25. &&
	    fabs(_v_Jet->at(i).eta()) < 2.4 )
	  _v_JetTight->push_back(_v_Jet->at(i));
     }      
   
   std::vector<int> res = filterPt(0,0,0,
				   _v_ElectronTight,
				   _v_MuonTight);

   if( res[0] >= 0 )
     {
	int id = _v_Event->at(0).id();
	int run = _v_Event->at(0).run();
	int lumi = _v_Event->at(0).lumi();
	
	float metpt = _v_Event->at(0).metpt();
	float metphi = _v_Event->at(0).metphi();

	m_channel = _v_Event->at(0).tth_channel();
	
	m_met = metpt;
	m_met_phi = metphi;
	
	m_weight = _v_Event->at(0).mc_weight();
	
	int njets = _v_JetTight->size();
	int nbjets = 0;
	for(int ij=0;ij<njets;ij++)
	  {
	     if( _v_JetTight->at(ij).CSVv2() >= 0.244 ) nbjets++;
	  }	
	
	int lidx1 = res[2];
	int lidx2 = res[3];
	int l1ise = res[4];
	int l2ise = res[5];

	int lidx3 = res[6];
	int lidx4 = res[7];
	int l3ise = res[8];
	int l4ise = res[9];
	
	int l1id = 0;
	int l1charge = 0;
	float l1pt = 0., l1eta = 0., l1phi = 0., l1m = 0.;
 	int l2id = 0;
	int l2charge = 0;
	float l2pt = 0., l2eta = 0., l2phi = 0., l2m = 0.;
	bool l1_isLoose = 0;
	bool l2_isLoose = 0;
	bool l1_isTight = 0;
	bool l2_isTight = 0;
	bool l1_isTightMVA = 0;
	bool l2_isTightMVA = 0;
	bool l1_passCF = 0;
	bool l2_passCF = 0;

 	int l3id = 0;
	int l3charge = 0;
	float l3pt = 0., l3eta = 0., l3phi = 0., l3m = 0.;

 	int l4id = 0;
	int l4charge = 0;
	float l4pt = 0., l4eta = 0., l4phi = 0., l4m = 0.;
	
	if( l1ise == 1 ) 
	  {
	     l1id = _v_ElectronTight->at(lidx1).id();
	     l1charge = _v_ElectronTight->at(lidx1).charge();
	     l1pt = _v_ElectronTight->at(lidx1).pt();
	     l1eta = _v_ElectronTight->at(lidx1).eta();
	     l1phi = _v_ElectronTight->at(lidx1).phi();
	     l1m = _v_ElectronTight->at(lidx1).m();
	     l1_isLoose = _v_ElectronTight->at(lidx1).isLoose();
	     l1_isTight = _v_ElectronTight->at(lidx1).isTight();
	     l1_isTightMVA = _v_ElectronTight->at(lidx1).isTightMVA();
	     l1_passCF = _v_ElectronTight->at(lidx1).passChargeFlip();
	     if( l2ise == 1 ) 
	       {
		  l2id = _v_ElectronTight->at(lidx2).id();
		  l2charge = _v_ElectronTight->at(lidx2).charge();
		  l2pt = _v_ElectronTight->at(lidx2).pt();
		  l2eta = _v_ElectronTight->at(lidx2).eta();
		  l2phi = _v_ElectronTight->at(lidx2).phi();
		  l2m = _v_ElectronTight->at(lidx2).m();
		  l2_isLoose = _v_ElectronTight->at(lidx2).isLoose();
		  l2_isTight = _v_ElectronTight->at(lidx2).isTight();
		  l2_isTightMVA = _v_ElectronTight->at(lidx2).isTightMVA();
		  l2_passCF = _v_ElectronTight->at(lidx2).passChargeFlip();
	       }
	     else
	       {
		  l2id = _v_MuonTight->at(lidx2).id();
		  l2charge = _v_MuonTight->at(lidx2).charge();
		  l2pt = _v_MuonTight->at(lidx2).pt();
		  l2eta = _v_MuonTight->at(lidx2).eta();
		  l2phi = _v_MuonTight->at(lidx2).phi();
		  l2m = _v_MuonTight->at(lidx2).m();
		  l2_isLoose = _v_MuonTight->at(lidx2).isLoose();
		  l2_isTight = _v_MuonTight->at(lidx2).isTight();
		  l2_isTightMVA = _v_MuonTight->at(lidx2).isTightMVA();
		  l2_passCF = _v_MuonTight->at(lidx2).passChargeFlip();
	       }	     
	  }	
	else 
	  {
	     l1id = _v_MuonTight->at(lidx1).id();
	     l1charge = _v_MuonTight->at(lidx1).charge();
	     l1pt = _v_MuonTight->at(lidx1).pt();
	     l1eta = _v_MuonTight->at(lidx1).eta();
	     l1phi = _v_MuonTight->at(lidx1).phi();
	     l1m = _v_MuonTight->at(lidx1).m();
	     l1_isLoose = _v_MuonTight->at(lidx1).isLoose();
	     l1_isTight = _v_MuonTight->at(lidx1).isTight();
	     l1_isTightMVA = _v_MuonTight->at(lidx1).isTightMVA();
	     l1_passCF = _v_MuonTight->at(lidx1).passChargeFlip();
	     if( l2ise == 1 ) 
	       {
		  l2id = _v_ElectronTight->at(lidx2).id();
		  l2charge = _v_ElectronTight->at(lidx2).charge();
		  l2pt = _v_ElectronTight->at(lidx2).pt();
		  l2eta = _v_ElectronTight->at(lidx2).eta();
		  l2phi = _v_ElectronTight->at(lidx2).phi();
		  l2m = _v_ElectronTight->at(lidx2).m();
		  l2_isLoose = _v_ElectronTight->at(lidx2).isLoose();
		  l2_isTight = _v_ElectronTight->at(lidx2).isTight();
		  l2_isTightMVA = _v_ElectronTight->at(lidx2).isTightMVA();
		  l2_passCF = _v_ElectronTight->at(lidx2).passChargeFlip();
	       }
	     else
	       {
		  l2id = _v_MuonTight->at(lidx2).id();
		  l2charge = _v_MuonTight->at(lidx2).charge();
		  l2pt = _v_MuonTight->at(lidx2).pt();
		  l2eta = _v_MuonTight->at(lidx2).eta();
		  l2phi = _v_MuonTight->at(lidx2).phi();
		  l2m = _v_MuonTight->at(lidx2).m();
		  l2_isLoose = _v_MuonTight->at(lidx2).isLoose();
		  l2_isTight = _v_MuonTight->at(lidx2).isTight();
		  l2_isTightMVA = _v_MuonTight->at(lidx2).isTightMVA();
		  l2_passCF = _v_MuonTight->at(lidx2).passChargeFlip();
	       }
	  }

	if( lidx3 != -1 && l3ise == 1 )
	  {
	     l3id = _v_ElectronTight->at(lidx3).id();
	     l3charge = _v_ElectronTight->at(lidx3).charge();
	     l3pt = _v_ElectronTight->at(lidx3).pt();
	     l3eta = _v_ElectronTight->at(lidx3).eta();
	     l3phi = _v_ElectronTight->at(lidx3).phi();
	     l3m = _v_ElectronTight->at(lidx3).m();
	  }
	if( lidx3 != -1 && l3ise == 0 )
	  {
	     l3id = _v_MuonTight->at(lidx3).id();
	     l3charge = _v_MuonTight->at(lidx3).charge();
	     l3pt = _v_MuonTight->at(lidx3).pt();
	     l3eta = _v_MuonTight->at(lidx3).eta();
	     l3phi = _v_MuonTight->at(lidx3).phi();
	     l3m = _v_MuonTight->at(lidx3).m();
	  }	

	if( lidx4 != -1 && l4ise == 1 )
	  {
	     l4id = _v_ElectronTight->at(lidx4).id();
	     l4charge = _v_ElectronTight->at(lidx4).charge();
	     l4pt = _v_ElectronTight->at(lidx4).pt();
	     l4eta = _v_ElectronTight->at(lidx4).eta();
	     l4phi = _v_ElectronTight->at(lidx4).phi();
	     l4m = _v_ElectronTight->at(lidx4).m();
	  }
	if( lidx4 != -1 && l4ise == 0 )
	  {
	     l4id = _v_MuonTight->at(lidx4).id();
	     l4charge = _v_MuonTight->at(lidx4).charge();
	     l4pt = _v_MuonTight->at(lidx4).pt();
	     l4eta = _v_MuonTight->at(lidx4).eta();
	     l4phi = _v_MuonTight->at(lidx4).phi();
	     l4m = _v_MuonTight->at(lidx4).m();
	  }	
	
	// leptons
	if( l1id != 0 )
	  {
	     m_l1_type = (abs(l1id) == 11) ? 1 : 2;
	     m_l1_charge = l1charge;
	     m_l1_pt = l1pt;
	     m_l1_eta = l1eta;
	     m_l1_phi = l1phi;
	     m_l1_m = l1m;
	  }	
	if( l2id != 0 )
	  {
	     m_l2_type = (abs(l2id) == 11) ? 1 : 2;
	     m_l2_charge = l2charge;
	     m_l2_pt = l2pt;
	     m_l2_eta = l2eta;
	     m_l2_phi = l2phi;
	     m_l2_m = l2m;
	  }
	if( l3id != 0 )
	  {
	     m_l3_type = (abs(l3id) == 11) ? 1 : 2;
	     m_l3_charge = l3charge;
	     m_l3_pt = l3pt;
	     m_l3_eta = l3eta;
	     m_l3_phi = l3phi;
	     m_l3_m = l3m;
	  }
	if( l4id != 0 )
	  {
	     m_l4_type = (abs(l4id) == 11) ? 1 : 2;
	     m_l4_charge = l4charge;
	     m_l4_pt = l4pt;
	     m_l4_eta = l4eta;
	     m_l4_phi = l4phi;
	     m_l4_m = l4m;
	  }
	
	// jets
	for(int ij=0;ij<njets;ij++)
	  {
	     if( ij == 0 )
	       {		  
		  m_j1_ntrk = _v_JetTight->at(ij).ntrk();
		  m_j1_btag = _v_JetTight->at(ij).CSVv2();
		  m_j1_pt = _v_JetTight->at(ij).pt();
		  m_j1_eta = _v_JetTight->at(ij).eta();
		  m_j1_phi = _v_JetTight->at(ij).phi();
		  m_j1_m = _v_JetTight->at(ij).m();
	       }	     
	     if( ij == 1 )
	       {		  
		  m_j2_ntrk = _v_JetTight->at(ij).ntrk();
		  m_j2_btag = _v_JetTight->at(ij).CSVv2();
		  m_j2_pt = _v_JetTight->at(ij).pt();
		  m_j2_eta = _v_JetTight->at(ij).eta();
		  m_j2_phi = _v_JetTight->at(ij).phi();
		  m_j2_m = _v_JetTight->at(ij).m();
	       }	     
	     if( ij == 2 )
	       {		  
		  m_j3_ntrk = _v_JetTight->at(ij).ntrk();
		  m_j3_btag = _v_JetTight->at(ij).CSVv2();
		  m_j3_pt = _v_JetTight->at(ij).pt();
		  m_j3_eta = _v_JetTight->at(ij).eta();
		  m_j3_phi = _v_JetTight->at(ij).phi();
		  m_j3_m = _v_JetTight->at(ij).m();
	       }	     
	     if( ij == 3 )
	       {		  
		  m_j4_ntrk = _v_JetTight->at(ij).ntrk();
		  m_j4_btag = _v_JetTight->at(ij).CSVv2();
		  m_j4_pt = _v_JetTight->at(ij).pt();
		  m_j4_eta = _v_JetTight->at(ij).eta();
		  m_j4_phi = _v_JetTight->at(ij).phi();
		  m_j4_m = _v_JetTight->at(ij).m();
	       }	     
	     if( ij == 4 )
	       {		  
		  m_j5_ntrk = _v_JetTight->at(ij).ntrk();
		  m_j5_btag = _v_JetTight->at(ij).CSVv2();
		  m_j5_pt = _v_JetTight->at(ij).pt();
		  m_j5_eta = _v_JetTight->at(ij).eta();
		  m_j5_phi = _v_JetTight->at(ij).phi();
		  m_j5_m = _v_JetTight->at(ij).m();
	       }	     
	     if( ij == 5 )
	       {		  
		  m_j6_ntrk = _v_JetTight->at(ij).ntrk();
		  m_j6_btag = _v_JetTight->at(ij).CSVv2();
		  m_j6_pt = _v_JetTight->at(ij).pt();
		  m_j6_eta = _v_JetTight->at(ij).eta();
		  m_j6_phi = _v_JetTight->at(ij).phi();
		  m_j6_m = _v_JetTight->at(ij).m();
	       }	     
	  }       
	
	passSel = 0x0;
	passSel |= 1   << 0; 
//	passSel |= (l1_isLoose && l2_isLoose) << 0;
	bool isEESS = (res[1] == 0 && l1ise && l2ise && l1_passCF && l2_passCF);
	bool isMMSS = (res[1] == 0 && !l1ise && !l2ise && l1_passCF && l2_passCF);
	passSel |= isEESS   << 1;
	passSel |= isMMSS   << 2;
	bool pass_lepPt = (l1pt >= 20. && l2pt >= 20.);
	passSel |= pass_lepPt   << 3;
	bool pass_isTight = (l1_isTight && l2_isTight);
	passSel |= pass_isTight   << 4;
	bool pass_tightMVA = (l1_isTightMVA && l2_isTightMVA);
	passSel |= pass_tightMVA   << 5;
	bool pass_nbjets = (nbjets >= 2);
	passSel |= pass_nbjets   << 6;
	bool pass_njets = (njets >= 4);
	passSel |= pass_njets   << 7;

	if( CHECK_BIT(passSel,0) ) m_sel = 0; // >=2lep
	if( CHECK_BIT(passSel,1) )
	  {	     
	     m_sel = 10; // eeSS
	     if( CHECK_BIT(passSel,3) )
	       {
		  m_sel = 11; // eeSS_pt20
		  if( CHECK_BIT(passSel,6) )
		    {
		       m_sel = 12; // eeSS_pt20_bj2
		       if( CHECK_BIT(passSel,7) )
			 m_sel = 13; // eeSS_pt20_bj2_j4
		    }		  
	       }
	  }
	if( CHECK_BIT(passSel,2) )
	  {	     
	     m_sel = 20; // mmSS
	     if( CHECK_BIT(passSel,3) )
	       {
		  m_sel = 21; // mmSS_pt20
		  if( CHECK_BIT(passSel,6) )
		    {
		       m_sel = 22; // mmSS_pt20_bj2
		       if( CHECK_BIT(passSel,7) )
			 m_sel = 23; // mmSS_pt20_bj2_j4
		    }		  
	       }
	  }
	
	if( CHECK_BIT(passSel,0) &&
	    CHECK_BIT(passSel,2) &&
	    CHECK_BIT(passSel,3) &&
//	    CHECK_BIT(passSel,4) &&
	    CHECK_BIT(passSel,5)
////	    CHECK_BIT(passSel,5) &&
//	    CHECK_BIT(passSel,6) &&
//	    CHECK_BIT(passSel,7)
	  )
	  {	     
	     fprintf(_fevc,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %+2d  %6.2f %+4.2f %+4.2f    %6.1f  %+4.2f    %d \n",
		     run, lumi, id,
		     l1id, l1pt, l1eta, l1phi,
		     l2id, l2pt, l2eta, l2phi,
		     metpt, metphi,
		     njets);
	     
	     _trout->Fill();
	  }	
     }   
   
   return 1;
}

std::vector<int> Hist::filterPt(double Pt1, double Pt2, double Pt3,
				std::vector<Electron>* ntElectron,
				std::vector<Muon>* ntMuon)
{
   int llChan = -1;
   int llChar = -1;
   
   bool foundPt1 = false;

   std::vector<std::pair<int,double> > el_pt;
   std::vector<std::pair<int,double> > mu_pt;
   std::vector<std::pair<int,double> > ll_pt;
   
   // descending
   std::sort(ntElectron->begin(),ntElectron->end(),Electron::sortPtPredicate);
   std::sort(ntMuon->begin(),ntMuon->end(),Muon::sortPtPredicate);
   
   for(int i=0;i<ntElectron->size();i++)
     {	
	el_pt.push_back(std::make_pair(i,ntElectron->at(i).pt()));
     }   
   for(int i=0;i<ntMuon->size();i++)
     {
	mu_pt.push_back(std::make_pair(i+1000,ntMuon->at(i).pt()));
     }
   
   ll_pt.reserve( el_pt.size() + mu_pt.size() );
   ll_pt.insert( ll_pt.end(), el_pt.begin(), el_pt.end() );
   ll_pt.insert( ll_pt.end(), mu_pt.begin(), mu_pt.end() );

   // descending
   std::sort( ll_pt.begin(), ll_pt.end(),
	      boost::bind(&std::pair<int, double>::second, _1) >
	      boost::bind(&std::pair<int, double>::second, _2));

   bool found_pt1 = false;
   bool found_pt2 = false;

   std::vector<int> id_remove;

   for(int i=0;i<ll_pt.size();i++)
     {
	if( found_pt1 && found_pt2 )
	  {
	     if( ll_pt.at(i).second < Pt3 )
	       id_remove.push_back(i);
	  }	

	if( ! found_pt2 )
	  {	     
	     if( ll_pt.at(i).second > Pt1 )
	       found_pt1 = true;

	     if( (found_pt1 && ll_pt.at(i).second < Pt2) || !found_pt1 )
	       id_remove.push_back(i);
	     else
	       found_pt2 = true;
	  }	
     }   

   for(int j=id_remove.size()-1;j>=0;j--)
     {
	int index = ll_pt.at(id_remove.at(j)).first;
	
	bool is_electron = true;
	     
	if( index >= 1000 ) 
	  {
	     index = index - 1000;
	     is_electron = false;
	  }	     	

	ll_pt.erase(ll_pt.begin() + id_remove.at(j));
	
	if( is_electron )
	  ntElectron->erase(ntElectron->begin() + index);
	else
	  ntMuon->erase(ntMuon->begin() + index);
     }   

   bool l1_e = false;
   bool l1_m = false;

   bool l2_e = false;
   bool l2_m = false;

   bool l3_e = false;
   bool l3_m = false;

   bool l4_e = false;
   bool l4_m = false;
   
   int idx1 = -1;
   int idx2 = -1;
   int idx3 = -1;
   int idx4 = -1;

   // descending
   std::sort( ll_pt.begin(), ll_pt.end(),
	      boost::bind(&std::pair<int, double>::second, _1) >
	      boost::bind(&std::pair<int, double>::second, _2));

   if( ll_pt.size() > 1 )
     {
	if( ll_pt.at(0).first < 1000 ) l1_e = true;
	if( ll_pt.at(0).first >= 1000 ) l1_m = true;
	
	if( ll_pt.at(1).first < 1000 ) l2_e = true;
	if( ll_pt.at(1).first >= 1000 ) l2_m = true;

	if( ll_pt.size() > 2 )
	  {
	     if( ll_pt.at(2).first < 1000 ) l3_e = true;
	     if( ll_pt.at(2).first >= 1000 ) l3_m = true;	     
	  }	

	if( ll_pt.size() > 3 )
	  {
	     if( ll_pt.at(3).first < 1000 ) l4_e = true;
	     if( ll_pt.at(3).first >= 1000 ) l4_m = true;	     
	  }	
	
	// -1 - < 2 leptons
	//  0 - ee
	//  1 - mm
	//  2 - em
   	 
	// -1 - < 2 leptons
	//  0 - SS
	//  1 - OS

	if( l1_e && l2_e ) llChan = 0;
	if( l1_m && l2_m ) llChan = 1;
	if( (l1_e && l2_m) || 
	    (l1_m && l2_e) ) llChan = 2;
	
	idx1 = (ll_pt.at(0).first < 1000) ? ll_pt.at(0).first : ll_pt.at(0).first - 1000;
	idx2 = (ll_pt.at(1).first < 1000) ? ll_pt.at(1).first : ll_pt.at(1).first - 1000;
	if( ll_pt.size() > 2 ) idx3 = (ll_pt.at(2).first < 1000) ? ll_pt.at(2).first : ll_pt.at(2).first - 1000;
	if( ll_pt.size() > 3 ) idx4 = (ll_pt.at(3).first < 1000) ? ll_pt.at(3).first : ll_pt.at(3).first - 1000;

	int cha1 = -1;
	int cha2 = -1;
	
	if( idx1 >= 0 )
	  {	     
	     if( l1_e ) cha1 = ntElectron->at(idx1).charge();
	     if( l1_m ) cha1 = ntMuon->at(idx1).charge();
	  }
	if( idx2 >= 0 )
	  {	     
	     if( l2_e ) cha2 = ntElectron->at(idx2).charge();
	     if( l2_m ) cha2 = ntMuon->at(idx2).charge();
	  }	
	
	if( idx1 >= 0 && idx2 >= 0 )
	  {
	     llChar = (cha1 * cha2 == 1) ? 0 : 1;
	  }
     }

   std::vector<int> res;
   res.clear();
   res.push_back(llChan);
   res.push_back(llChar);
   res.push_back(idx1);
   res.push_back(idx2);
   res.push_back(l1_e);
   res.push_back(l2_e);

   res.push_back(idx3);
   res.push_back(idx4);
   res.push_back(l3_e);
   res.push_back(l4_e);
   
   return res;
}

/*
float SYSANA::Hist::getVar(std::string sys,int ijet,std::string varName,int ibin)
{
   float var = -666;

//   if( ibin-1 < 0 ) return var;
   
   if( ijet <= ntP->nJet )
     {	   
	if( strcmp(varName.c_str(),"h_j1_JP_") == 0 )
	  {
	     var = Jet_Proba_New;
	     if( var == 0. ) var = -666.;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_BJP_") == 0 )
	  {		  
	     var = ntP->Jet_Bprob[ijet];
	     if( var == 0. ) var = -666.;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_pthat_") == 0 )
	  {		  
	     var = ntP->pthat/1000.;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_ntrk_") == 0 )
	  {		  
	     var = ntP->Jet_ntracks[ijet];
	  }	     
	if( strcmp(varName.c_str(),"h_j1_nseltrk_") == 0 )
	  {		  
	     var = ntP->Jet_nseltracks[ijet];
	  }	     
	if( strcmp(varName.c_str(),"h_j1_eta_") == 0 )
	  {		  
	     var = ntP->Jet_eta[ijet];
	  }	     
	if( strcmp(varName.c_str(),"h_j1_phi_") == 0 )
	  {		  
	     var = ntP->Jet_phi[ijet];
	  }	     
	if( strcmp(varName.c_str(),"h_j1_mass_") == 0 )
	  {		  
	     var = ntP->Jet_mass[ijet];
	  }	     
	if( strcmp(varName.c_str(),"h_j1_mupt_") == 0 )
	  {
	     if( muidx.size() > 0 )
	       var = ntP->Muon_pt[muidx[0]];
	     else fillThis = false;
//	     for(int im=0;im<muidx.size();im++)
//	       {		  
//		  var = ntP->Muon_pt[muidx[im]];
//	       }	     
	  }	     
	if( strcmp(varName.c_str(),"h_j1_mueta_") == 0 )
	  {
	     if( muidx.size() > 0 )
	       var = ntP->Muon_eta[muidx[0]];
	     else fillThis = false;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_muphi_") == 0 )
	  {
	     if( muidx.size() > 0 )
	       var = ntP->Muon_phi[muidx[0]];
	     else fillThis = false;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_muptrel_") == 0 )
	  {
	     if( muidx.size() > 0 )
	       var = ntP->Muon_ptrel[muidx[0]];
	     else fillThis = false;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_npv_") == 0 )
	  {
	     if( ijet == 0 )
	       var = ntP->nPV;
	     else fillThis = false;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_npu_") == 0 )
	  {
	     if( ijet == 0 )
	       var = ntP->nPU;
	     else fillThis = false;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_njet_") == 0 )
	  {
	     if( ijet == 0 )
	       var = ntP->nJet;
	     else fillThis = false;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_nmuon_") == 0 )
	  {
	     if( ijet == 0 )
	       var = ntP->nMuon;
	     else fillThis = false;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_nsv_") == 0 )
	  {
	     if( ijet == 0 )
	       var = ntP->nSV;
	     else fillThis = false;
	  }	     
	if( strcmp(varName.c_str(),"h_j1_nsvj_") == 0 )
	  {
	     var = ntP->Jet_nLastSV[ijet];
	  }	     

	     if( strcmp(varName.c_str(),"h_j1_pt_") == 0 )
	       {	
		  double jptmax = 1250.;
		  double jptmin = 0.;
		  if( ibin > 0 ) 
		    {
		       jptmax = jPtMax[ibin-1];
		       jptmin = jPtMin[ibin-1];
		    }	     
		  
		  var = (getPt(sys)-jptmin)/(jptmax-jptmin);
	       }	     
     }   

   return var;
}

std::vector<float> SYSANA::Hist::getVarVec(std::string sys,int ijet,std::string varName,int ibin)
{
   std::vector<float> var;
   var.clear();

   if( ibin-1 < 0 ) return var;
   
   if( ijet <= ntP->nJet )
     {	   
	if( strcmp(varName.c_str(),"h_j1_svntrk_") == 0 )
	  {
	     if( ijet == 0 )
	       {
		  for(int isv=0;isv<ntP->nSV;isv++)
		    var.push_back(ntP->SV_nTrk[isv]);
	       }	     
	     else fillThis = false;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_svmass_") == 0 )
	  {
	     if( ijet == 0 )
	       {
		  for(int isv=0;isv<ntP->nSV;isv++)
		    var.push_back(ntP->SV_mass[isv]);
	       }	     
	     else fillThis = false;
	  }	     
     }   

   return var;
}

std::pair<float,float> SYSANA::Hist::getVar2d(std::string sys,int ijet,std::string varName,int ibin)
{
   float varX = -666;
   float varY = -666;
   
   if( ijet <= ntP->nJet )
     {	   
	if( strcmp(varName.c_str(),"h_j1_JP_vs_nseltrk_") == 0 )
	  {		  
	     varX = ntP->Jet_nseltracks[ijet];
	     varY = Jet_Proba_New;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_nseltrk_vs_ntrkgen_") == 0 )
	  {		  
	     varY = ntP->Jet_nseltracks[ijet];
	     varX = ntrkgen;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_JP_vs_njet_") == 0 )
	  {
	     varX = ntP->nJet;
	     varY = Jet_Proba_New;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_JP_vs_jeta_") == 0 )
	  {
	     varX = ntP->Jet_eta[ijet];
	     varY = Jet_Proba_New;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_JP_vs_nsv_") == 0 )
	  {
	     varX = ntP->nSV;
	     varY = Jet_Proba_New;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_JP_vs_pt_") == 0 )
	  {	     
	     double jptmax = 1250.;
	     double jptmin = 0.;
	     if( ibin > 0 ) 
	       {
		  jptmax = jPtMax[ibin-1];
		  jptmin = jPtMin[ibin-1];
	       }	     

	     varX = (getPt(sys)-jptmin)/(jptmax-jptmin);
	     varY = Jet_Proba_New;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_JP_vs_CSV_") == 0 )
	  {
	     varX = ntP->Jet_CombSvx[ijet];
	     varY = Jet_Proba_New;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_BJP_vs_CSV_") == 0 )
	  {
	     varX = ntP->Jet_CombSvx[ijet];
	     varY = ntP->Jet_Bprob[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_pt_vs_nseltrk_") == 0 )
	  {	    
	     double jptmax = 1250.;
	     double jptmin = 0.;
	     if( ibin > 0 ) 
	       {
		  jptmax = jPtMax[ibin-1];
		  jptmin = jPtMin[ibin-1];
	       }	     

	     varY = (getPt(sys)-jptmin)/(jptmax-jptmin);
	     varX = ntP->Jet_nseltracks[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_pt_vs_njet_") == 0 )
	  {	    
	     double jptmax = 1250.;
	     double jptmin = 0.;
	     if( ibin > 0 ) 
	       {
		  jptmax = jPtMax[ibin-1];
		  jptmin = jPtMin[ibin-1];
	       }	     

	     varY = (getPt(sys)-jptmin)/(jptmax-jptmin);
	     varX = ntP->nJet;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_pt_vs_jeta_") == 0 )
	  {	    
	     double jptmax = 1250.;
	     double jptmin = 0.;
	     if( ibin > 0 ) 
	       {
		  jptmax = jPtMax[ibin-1];
		  jptmin = jPtMin[ibin-1];
	       }	     

	     varY = (getPt(sys)-jptmin)/(jptmax-jptmin);
	     varX = ntP->Jet_eta[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_JP_vs_mupt_") == 0 )
	  {		  
	     // highest-Pt muon in a jet
	     if( muidx.size() > 0 )
	       varX = ntP->Muon_pt[muidx[0]];
	     else fillThis = false;
	     varY = Jet_Proba_New;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_njet_vs_mupt_") == 0 )
	  {		  
	     // highest-Pt muon in a jet
	     if( muidx.size() > 0 )
	       varX = ntP->Muon_pt[muidx[0]];
	     else fillThis = false;
	     varY = ntP->nJet;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_mupt_vs_jeta_") == 0 )
	  {		  
	     // highest-Pt muon in a jet
	     if( muidx.size() > 0 )
	       varY = ntP->Muon_pt[muidx[0]];
	     else fillThis = false;
	     varX = ntP->Jet_eta[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_njet_vs_jeta_") == 0 )
	  {		  
	     varX = ntP->Jet_eta[ijet];
	     varY = ntP->nJet;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_nseltrk_vs_mupt_") == 0 )
	  {		  
	     if( muidx.size() > 0 )
	       varX = ntP->Muon_pt[muidx[0]];
	     else fillThis = false;
	     varY = ntP->Jet_nseltracks[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_nseltrk_vs_njet_") == 0 )
	  {		  
	     varX = ntP->nJet;
	     varY = ntP->Jet_nseltracks[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_nseltrk_vs_jeta_") == 0 )
	  {		  
	     varX = ntP->Jet_eta[ijet];
	     varY = ntP->Jet_nseltracks[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_pt_vs_mupt_") == 0 )
	  {	     	     
	     // highest-Pt muon in a jet
	     if( muidx.size() > 0 )
	       varX = ntP->Muon_pt[muidx[0]];
	     else fillThis = false;

	     double jptmax = 1250.;
	     double jptmin = 0.;
	     if( ibin > 0 )
	       {
		  jptmax = jPtMax[ibin-1];
		  jptmin = jPtMin[ibin-1];
	       }	     
	     
	     varY = (getPt(sys)-jptmin)/(jptmax-jptmin);
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_nseltrk_vs_nsv_") == 0 )
	  {		  
	     varX = ntP->nSV;
	     varY = ntP->Jet_nseltracks[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_njet_vs_nsv_") == 0 )
	  {
	     if( ijet == 0 )
	       {		  
		  varX = ntP->nSV;
		  varY = ntP->nJet;
	       }
	     else fillThis = false;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_mupt_vs_nsv_") == 0 )
	  {		  
	     // highest-Pt muon in a jet
	     if( muidx.size() > 0 )
	       varY = ntP->Muon_pt[muidx[0]];
	     else fillThis = false;
	     varX = ntP->nSV;
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_jeta_vs_nsv_") == 0 )
	  {		  
	     varX = ntP->nSV;
	     varY = ntP->Jet_eta[ijet];
	  }	     
	else if( strcmp(varName.c_str(),"h_j1_pt_vs_nsv_") == 0 )
	  {	    
	     double jptmax = 1250.;
	     double jptmin = 0.;
	     if( ibin > 0 ) 
	       {
		  jptmax = jPtMax[ibin-1];
		  jptmin = jPtMin[ibin-1];
	       }	     

	     varY = (getPt(sys)-jptmin)/(jptmax-jptmin);
	     varX = ntP->nSV;
	  }	     
	else
	  {
	     std::cout << "Not found:" << varName << std::endl;
	     exit(1);
	  }	
     }   
   
   return std::make_pair(varX,varY);
}

std::vector<float> SYSANA::Hist::getVar3d(std::string sys,int ijet,std::string varName,int ibin)
{
   float varX = -666;
   float varY = -666;
   float varZ = -666;
   
   if( ijet <= ntP->nJet )
     {	   
	if( strcmp(varName.c_str(),"h_j1_ntrk_vs_pt_vs_mupt_") == 0 )
	  {		  
	     varZ = ntP->Jet_ntracks[ijet];
	     
	     double jptmax = 1250.;
	     double jptmin = 0.;
	     if( ibin > 0 ) 
	       {
		  jptmax = jPtMax[ibin-1];
		  jptmin = jPtMin[ibin-1];
	       }	     
	     varY = (getPt(sys)-jptmin)/(jptmax-jptmin);

	     // highest-Pt muon in a jet
	     if( muidx.size() > 0 )
	       varX = ntP->Muon_pt[muidx[0]];
	  }	     
     }   

   std::vector<float> res;
   res.push_back(varX);
   res.push_back(varY);
   res.push_back(varZ);
   return res;
}

float SYSANA::Hist::getPt(std::string sys)
{
   float var = -666;
   
   if( (sys == "" ||
	sys == "_sys_up" || sys == "_sys_low" ||
	sys == "_pu_up" || sys == "_pu_low" ||
	sys == "_gsplit_up" || sys == "_gsplit_low" ||
	sys == "_bfrag_up" || sys == "_bfrag_low" ||
	sys == "_cdfrag_up" || sys == "_cdfrag_low" ||
	sys == "_cfrag_up" || sys == "_cfrag_low" ||
	sys == "_ksl_up" || sys == "_ksl_low" ||
	sys == "_ntrkgen_up" || sys == "_ntrkgen_low"
       )
       ||
       isdata )
     {
	var = v_jet->Pt();
     }
   else if( sys == "_jes_up" )
     {
	var = v_jet_sys_jesTotalUp->Pt();
     }	
   else if( sys == "_jes_low" )
     {
	var = v_jet_sys_jesTotalLow->Pt();
     }
   else if( sys == "_jer_up" )
     {
	var = v_jet_sys_jerTotalUp->Pt();
     }	
   else if( sys == "_jer_low" )
     {
	var = v_jet_sys_jerTotalLow->Pt();
     }	
   else
     {
	std::cout << "No jet Pt for: " << sys << std::endl;
	exit(1);
     }	

   return var;
}
*/

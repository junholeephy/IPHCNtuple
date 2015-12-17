void draw()
{
   gROOT->SetBatch(1);
   
   gROOT->ProcessLine(".L PlotStyle.C");
   
   SetPlotStyle();

   gStyle->SetOptFit();
   
   TFile *f = TFile::Open("../hist/output.root");
   TTree *tr = (TTree*)f->Get("tr");

   int lepTruth_n, lepTruth_label[1000];
   float lepTruth_pt[1000], lepTruth_eta[1000], lepTruth_phi[1000], lepTruth_E[1000];

   int lepRec_n, lepRec_label[1000];
   float lepRec_pt[1000], lepRec_eta[1000], lepRec_phi[1000], lepRec_E[1000];
   
   tr->SetBranchAddress("lepTruth_n",&lepTruth_n);
   tr->SetBranchAddress("lepTruth_pt",&lepTruth_pt);
   tr->SetBranchAddress("lepTruth_eta",&lepTruth_eta);
   tr->SetBranchAddress("lepTruth_phi",&lepTruth_phi);
   tr->SetBranchAddress("lepTruth_E",&lepTruth_E);
   tr->SetBranchAddress("lepTruth_label",&lepTruth_label);

   tr->SetBranchAddress("lepRec_n",&lepRec_n);
   tr->SetBranchAddress("lepRec_pt",&lepRec_pt);
   tr->SetBranchAddress("lepRec_eta",&lepRec_eta);
   tr->SetBranchAddress("lepRec_phi",&lepRec_phi);
   tr->SetBranchAddress("lepRec_E",&lepRec_E);
   tr->SetBranchAddress("lepRec_label",&lepRec_label);
   
   TCanvas *c1 = new TCanvas();

   TH1F *h_lep_eta_pos = new TH1F("h_lep_eta_pos","h_lep_eta_pos",20,-2.5,2.5);
   h_lep_eta_pos->Sumw2();
   TH1F *h_lep_eta_neg = new TH1F("h_lep_eta_neg","h_lep_eta_neg",20,-2.5,2.5);
   h_lep_eta_neg->Sumw2();

   TH1F *h_lepRec_eta_pos = new TH1F("h_lepRec_eta_pos","h_lepRec_eta_pos",20,-2.5,2.5);
   h_lepRec_eta_pos->Sumw2();
   TH1F *h_lepRec_eta_neg = new TH1F("h_lepRec_eta_neg","h_lepRec_eta_neg",20,-2.5,2.5);
   h_lepRec_eta_neg->Sumw2();

   TH1F *h_lep_dr = new TH1F("h_lep_dr","h_lep_dr",30,0.0,1.0);
   h_lep_dr->Sumw2();
   
   int nEvents = tr->GetEntries();
   for(int i=0;i<nEvents;i++)
     {
	tr->GetEntry(i);

	for(int il=0;il<lepTruth_n;il++)
	  {
	     if( lepTruth_pt[il] > 40. && fabs(lepTruth_eta[il]) < 2.5 )
	       {
		  if( lepTruth_label[il] == -11 || lepTruth_label[il] == -13 )
		    h_lep_eta_pos->Fill(lepTruth_eta[il]);
		  else if( lepTruth_label[il] == 11 || lepTruth_label[il] == 13 )
		    h_lep_eta_neg->Fill(lepTruth_eta[il]);
	       }
	  }

	for(int il=0;il<lepRec_n;il++)
	  {
	     if( lepRec_pt[il] > 40. && fabs(lepRec_eta[il]) < 2.5 )
	       {
		  float dr = 666;
		  if( lepRec_n > 0 && lepTruth_n > 0 )
		    {
		       TLorentzVector v1; v1.SetPtEtaPhiE(lepRec_pt[0],lepRec_eta[0],lepRec_phi[0],lepRec_E[0]);
		       TLorentzVector v2; v2.SetPtEtaPhiE(lepTruth_pt[0],lepTruth_eta[0],lepTruth_phi[0],lepTruth_E[0]);
		       dr = v1.DeltaR(v2);
		    }	

		  if( dr < 0.1 )
		    {		       
		       if( lepRec_label[il] == -11 || lepRec_label[il] == -13 )
			 h_lepRec_eta_pos->Fill(lepRec_eta[il]);
		       else if( lepRec_label[il] == 11 || lepRec_label[il] == 13 )
			 h_lepRec_eta_neg->Fill(lepRec_eta[il]);
		    }		  
	       }
	  }
	
//	if( lepRec_n > 0 && lepTruth_n > 0 )
//	  {
//	     TLorentzVector v1; v1.SetPtEtaPhiE(lepRec_pt[0],lepRec_eta[0],lepRec_phi[0],lepRec_E[0]);
//	     TLorentzVector v2; v2.SetPtEtaPhiE(lepTruth_pt[0],lepTruth_eta[0],lepTruth_phi[0],lepTruth_E[0]);
//	     float dr = v1.DeltaR(v2);
//	     h_lep_dr->Fill(dr);
//	  }	
     }   
	
   addbin(h_lep_eta_pos);
   addbin(h_lep_eta_neg);

   addbin(h_lepRec_eta_pos);
   addbin(h_lepRec_eta_neg);
   
   addbin(h_lep_dr);
   
   h_lep_eta_pos->Divide(h_lep_eta_pos,h_lep_eta_neg);
   
   h_lepRec_eta_pos->Divide(h_lepRec_eta_pos,h_lepRec_eta_neg);
   
   std::cout << h_lepRec_eta_pos->GetEntries() << std::endl;
   
   h_lep_eta_pos->Draw("hist e1p");
   h_lep_eta_pos->GetYaxis()->SetTitle("ratio");
   h_lep_eta_pos->GetXaxis()->SetTitle("#eta");
   c1->Print("pics/lep_eta_ratio.eps");
   c1->Clear();

   h_lepRec_eta_pos->Draw("hist e1p");
   h_lepRec_eta_pos->GetYaxis()->SetTitle("ratio");
   h_lepRec_eta_pos->GetXaxis()->SetTitle("#eta");
   c1->Print("pics/lepRec_eta_ratio.eps");
   c1->Clear();

   h_lep_dr->Draw("hist e1p");
   h_lep_dr->GetYaxis()->SetTitle("#DeltaR");
   h_lep_dr->GetXaxis()->SetTitle("Number of events");
   c1->Print("pics/lep_dr.eps");
   c1->Clear();
   
   gApplication->Terminate();
}

void addbin(TH1F *h)
{   
   // Add overflow and underflow bins
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
   h->SetBinError(1,TMath::Sqrt(pow(h->GetBinError(0),2)+pow(h->GetBinError(1),2)));
   h->SetBinContent(x_nbins,h->GetBinContent(x_nbins)+h->GetBinContent(x_nbins+1));
   h->SetBinError(x_nbins,TMath::Sqrt(pow(h->GetBinError(x_nbins),2)+
				      pow(h->GetBinError(x_nbins+1),2)));
   // Set overflow and underflow bins to 0
   h->SetBinContent(0,0.);
   h->SetBinError(0,0.);
   h->SetBinContent(x_nbins+1,0.);
   h->SetBinError(x_nbins+1,0.);
}

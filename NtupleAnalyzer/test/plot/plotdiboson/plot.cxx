{

    // Setup a style
    using namespace RooFit;

    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleSize(0.07, "XYZ");
    gStyle->SetTitleFont(22,"X");
    gStyle->SetTitleFont(22,"Y");
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetHistLineWidth(2);

    // Grab PDF's from histos

    THStack hs("hs","test stacked histograms");

    // WZ
    TFile *_file ;
    _file = TFile::Open("./Input/output_WZ_MC.root");
    _file->cd();

    TH1F*  hist0 = (TH1F*)_file->Get("InvMassRemainingLepton_FinalCut_WZ_CR__TTbarHiggs");
    hist0->Sumw2();

    // Other contributions
    // Sorted by number of events: WZ, ZZ, VVV, DY, TT, TTW, TTH

    TFile *_file01 ;
    _file01 = TFile::Open("./Input/output_ZZ_MC.root");
    _file01->cd();

    TH1F*  hist1 = (TH1F*)_file01->Get("InvMassRemainingLepton_FinalCut_WZ_CR__TTbarHiggs");
    hist1->Sumw2();

    TFile *_file02 ;
    _file02 = TFile::Open("./Input/output_DY_MC.root");
    _file02->cd();

    TH1F*  hist2 = (TH1F*)_file02->Get("InvMassRemainingLepton_FinalCut_WZ_CR__TTbarHiggs");
    hist2->Sumw2();

    //TH1F*  hist1 = (TH1F*)_file->Get("mtw_TTZ");
    //hist2->Sumw2();hist2->Rebin(2);    
    //TH1F*  hist2 = (TH1F*)_file->Get("mtw_DY");
    //hist2->Sumw2();hist2->Rebin(2);
    //TH1F*  hist3 = (TH1F*)_file->Get("mtw_ZZ");
    //hist3->Sumw2();hist3->Rebin(2);
    //TH1F*  hist4 = (TH1F*)_file->Get("mtw_VVV");
    //hist4->Sumw2();hist4->Rebin(2);
    //TH1F*  hist5 = (TH1F*)_file->Get("mtw_TT");
    //hist5->Sumw2();hist5->Rebin(2);
    //TH1F*  hist6 = (TH1F*)_file->Get("mtw_TT");
    //hist6->Sumw2();hist6->Rebin(2);

    TH1F*  histBkg = hist1->Clone(); //histBkg->Add(hist1,1);
    //histBkg->Add(hist2);             //histBkg->Add(hist2,1);
    //histBkg->Add(hist3);         //histBkg->Add(hist3,1);
    //histBkg->Add(hist4);         //histBkg->Add(hist4,1);
    //histBkg->Add(hist5);         //histBkg->Add(hist5,1);
    //histBkg->Add(hist6);         //histBkg->Add(hist6,1);

    // Get Data histo

    TFile *_file1 ;
    _file1 = TFile::Open("./Input/output_DATA.root");
    _file1->cd();

    TH1F*  histData = (TH1F*)_file1->Get("InvMassRemainingLepton_FinalCut_WZ_CR__TTbarHiggs");
    histData->Sumw2();
    histData->SetMarkerStyle(20);
    RooRealVar
        x("x","x",histData->GetXaxis()->GetXmin(),histData->GetXaxis()->GetXmax());
    RooDataHist* data = new RooDataHist("data","data",x, histData);

    // Define MC PDFs
    RooDataHist * rooHist0   = new RooDataHist("WZ",  "WZ",  x, hist0  );
    RooDataHist * rooHistBkg = new RooDataHist("BKG", "BKG", x, histBkg);

    RooHistPdf histpdf0("histpdf0",   "histpdf0",   x, *rooHist0  );
    RooHistPdf histpdf5("histpdfBkg", "histpdfBkg", x, *rooHistBkg);


    // Variable to be fitted nWZ, nBkg set to constant
    RooRealVar nWZ("nWZ", "nWZ", 100., 0., 200.);

    double yieldBkg = histBkg->Integral();
    RooRealVar nBkg("nBkg", "nBkg", yieldBkg);
    nBkg.setConstant(kTRUE);
    cout << "yieldBkg " << yieldBkg << endl;

    // Combine PDFs
    RooArgList shapes;
    RooArgList yields;

    shapes.add(histpdf0); yields.add(nWZ);
    shapes.add(histpdf5); yields.add(nBkg);

    RooAddPdf *pdf = new RooAddPdf("pdf"," pdf", shapes, yields);

    // Do the fit.
    RooFitResult* fitResult = pdf->fitTo(*data, Save()) ;

    //std::cout << nWZ.getVal() << std::endl;
    std::cout << nWZ.Print()  << std::endl;

    // Plot pdf's

    RooPlot* mframe =
        x.frame(histData->GetXaxis()->GetXmin(),histData->GetXaxis()->GetXmax()) ;
    data->plotOn(mframe,Name("data")) ;
    pdf->plotOn(mframe,LineColor(kBlack),LineStyle(kDashed));
    pdf->plotOn(mframe,Name("pdf")) ;
    mframe->Draw() ;

    // Calculate the various contributions

    /*cout << " the WZ " << c0.getVal()*histData->Integral() << " ± " <<
      c0.getError()*histData->Integral() << endl;
      cout << " the BKG        " << (1.-c0.getVal())*histData->Integral()
      << " ± " << c0.getError()*histData->Integral() << endl;
      cout << " the Data "            << histData->Integral()    << endl;
      cout << " the MC WZ " << hist0->Integral() << endl;

      cout << "scale factor WZ " <<
      c0.getVal()*histData->Integral()/hist0->Integral() << " ± " <<
      c0.getError()*histData->Integral()/hist0->Integral()<< endl;
      cout << "scale factor Bkg "        <<
      (1-c0.getVal())*histData->Integral()/histBkg->Integral() << " ± " <<
      c0.getError()*histData->Integral()/histBkg->Integral()<< endl;
      cout <<" " <<endl;
      cout << "chisquare " <<mframe->chiSquare("pdf","data",1)<<endl;
      cout << "pull "
      <<(c0.getVal()*histData->Integral()-hist0->Integral())/hist0->Integral()
      << std::endl;
      */
    // Plot the stuff

    //

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! A DECOMMENTER
    //histBkg->Scale((1-c0.getVal())*histData->Integral()/histBkg->Integral());


    //TEST
    //histBkg->Scale((1-nWZ.getVal())*histData->Integral()/histBkg->Integral());
    //histBkg->Scale((nWZ.getVal())*histData->Integral()/histBkg->Integral());
    //histBkg->Scale(histData->Integral()/histBkg->Integral());

    std::cout << "scale factor : "<< nWZ.getVal()/hist0->Integral() << std::endl;
    std::cout << "nWZ MC : "      << hist0->Integral()              << std::endl;
    std::cout << "yieldBkg "      << yieldBkg                       << std::endl;

    //hist0->Scale(nWZ.getVal()/hist0->Integral());

    histBkg->SetFillColor(kRed+2);
    hist0->SetFillColor(kBlue+2);

    histData->GetXaxis()->SetTitle("m_TW(l) [GEV]");
    histData->GetYaxis()->SetTitle("Events");
    histData->GetYaxis()->SetRangeUser(0,20);

    histData->Draw("e");
    hs.Add(histBkg);
    hs.Add(hist0);
    hs.Draw("histosame");
    histData->SetLineColor(kRed+2);
    histData->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData, "Data",                  "p");
    qw->AddEntry(hist0,    "WZ",                    "f");
    qw->AddEntry(histBkg,  "residual backgrounds",  "f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"CMS, #sqrt{s}=13 TeV, 2.2 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.18);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_mTW.pdf");

    fitResult->Print("") ;


    // #################################
    // # the selected jet multiplicity #
    // #################################

     _file->cd();
    TH1F*  hist0_nJet25 = (TH1F*)_file->Get("JetMultiplicity_FinalCut_WZ_CR__TTbarHiggs");
    hist0_nJet25->Sumw2();
     _file01->cd();
    TH1F*  hist1_nJet25 = (TH1F*)_file01->Get("JetMultiplicity_FinalCut_WZ_CR__TTbarHiggs");
    hist1_nJet25->Sumw2();
    _file02->cd();
    TH1F*  hist2_nJet25 = (TH1F*)_file02->Get("JetMultiplicity_FinalCut_WZ_CR__TTbarHiggs");
    hist2_nJet25->Sumw2();
    //TH1F*  hist3_nJet25 = (TH1F*)_file->Get("nJet25_ZZ");
    //hist3_nJet25->Sumw2();
    //TH1F*  hist4_nJet25 = (TH1F*)_file->Get("nJet25_VVV");
    //hist4_nJet25->Sumw2();
    //TH1F*  hist5_nJet25 = (TH1F*)_file->Get("nJet25_TT");
    //hist5_nJet25->Sumw2();
    //TH1F*  hist6_nJet25 = (TH1F*)_file->Get("nJet25_TTW");
    //hist6_nJet25->Sumw2();


    TH1F*  histBkg_nJet25 = hist1_nJet25->Clone();
    //histBkg_nJet25->Add(hist2_nJet25);
    //histBkg_nJet25->Add(hist3_nJet25);
    //histBkg_nJet25->Add(hist4_nJet25);
    //histBkg_nJet25->Add(hist5_nJet25);
    //histBkg_nJet25->Add(hist6_nJet25);

     _file1->cd();
    TH1F*  histData_nJet25 = (TH1F*)_file1->Get("JetMultiplicity_FinalCut_WZ_CR__TTbarHiggs");
    histData_nJet25->Sumw2();
    histData_nJet25->SetMarkerStyle(20);

    //hist0_nJet25->Scale(nWZ.getVal()/hist0_nJet25->Integral());

    THStack hs_nJet25("hs_nJet25","test stacked histograms");

    histBkg_nJet25->SetFillColor(kRed+2);
    hist0_nJet25->SetFillColor(kBlue+2);

    histData_nJet25->GetXaxis()->SetTitle("N(jet, p_{T} > 25)");
    histData_nJet25->GetYaxis()->SetTitle("Events");
    histData_nJet25->GetYaxis()->SetRangeUser(0,40);

    histData_nJet25->Draw("e");
    hs_nJet25.Add(histBkg_nJet25);
    hs_nJet25.Add(hist0_nJet25);
    hs_nJet25.Draw("histosame");
    histData_nJet25->SetLineColor(kRed+2);
    histData_nJet25->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_nJet25, "Data" ,                "p");
    qw->AddEntry(hist0_nJet25,    "WZ" ,                  "f");
    qw->AddEntry(histBkg_nJet25,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"CMS, #sqrt{s}=13 TeV, 2.2 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.18);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plotsWZ/WZ_3l_nJet25.pdf");

    // #########################
    // # Sum of lepton charges #
    // #########################

     _file->cd();
    TH1F*  hist0_SumLepCharges = (TH1F*)_file->Get("SumOfLeptonCharges_FinalCut_WZ_CR__TTbarHiggs");
    hist0_SumLepCharges->Sumw2();
     _file01->cd();
    TH1F*  hist1_SumLepCharges = (TH1F*)_file01->Get("SumOfLeptonCharges_FinalCut_WZ_CR__TTbarHiggs");
    hist1_SumLepCharges->Sumw2();
    _file02->cd();
    TH1F*  hist2_SumLepCharges = (TH1F*)_file02->Get("SumOfLeptonCharges_FinalCut_WZ_CR__TTbarHiggs");
    hist2_SumLepCharges->Sumw2();
    //TH1F*  hist3_SumLepCharges = (TH1F*)_file->Get("SumLepCharges_ZZ");
    //hist3_SumLepCharges->Sumw2();
    //TH1F*  hist4_SumLepCharges = (TH1F*)_file->Get("SumLepCharges_VVV");
    //hist4_SumLepCharges->Sumw2();
    //TH1F*  hist5_SumLepCharges = (TH1F*)_file->Get("SumLepCharges_TT");
    //hist5_SumLepCharges->Sumw2();
    //TH1F*  hist6_SumLepCharges = (TH1F*)_file->Get("SumLepCharges_TTW");
    //hist6_SumLepCharges->Sumw2();

    TH1F*  histBkg_SumLepCharges = hist1_SumLepCharges->Clone();
    //histBkg_SumLepCharges->Add(hist2_SumLepCharges);
    ///histBkg_SumLepCharges->Add(hist3_SumLepCharges);
    //histBkg_SumLepCharges->Add(hist4_SumLepCharges);
    //histBkg_SumLepCharges->Add(hist5_SumLepCharges);
    //histBkg_SumLepCharges->Add(hist6_SumLepCharges);

    _file1->cd();
    TH1F*  histData_SumLepCharges = (TH1F*)_file1->Get("SumOfLeptonCharges_FinalCut_WZ_CR__TTbarHiggs");
    histData_SumLepCharges->Sumw2();
    histData_SumLepCharges->SetMarkerStyle(20);

    //hist0_SumLepCharges->Scale(nWZ.getVal()/hist0_SumLepCharges->Integral());
    //std::cout << nWZ.getVal()/hist0_SumLepCharges->Integral()<< std::endl;

    THStack hs_SumLepCharges("hs_SumLepCharges","test stacked histograms");

    histBkg_SumLepCharges->SetFillColor(kRed+2);
    hist0_SumLepCharges->SetFillColor(kBlue+2);

    histData_SumLepCharges->GetXaxis()->SetTitle("sum lepton charges");
    histData_SumLepCharges->GetYaxis()->SetTitle("Events");
    histData_SumLepCharges->GetYaxis()->SetRangeUser(0,40);

    histData_SumLepCharges->Draw("e");
    hs_SumLepCharges.Add(histBkg_SumLepCharges);
    hs_SumLepCharges.Add(hist0_SumLepCharges);
    hs_SumLepCharges.Draw("histosame");
    histData_SumLepCharges->SetLineColor(kRed+2);
    histData_SumLepCharges->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.2,0.657163,0.45,0.85,NULL,"brNDC");

    qw->AddEntry(histData_SumLepCharges, "Data",                  "p");
    qw->AddEntry(hist0_SumLepCharges,    "WZ",                    "f");
    qw->AddEntry(histBkg_SumLepCharges,  "residual backgrounds",  "f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"CMS, #sqrt{s}=13 TeV, 2.26 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.18);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_sumJetCharges.pdf");

    // ###################################################################
    // # vectorial sum of the transverse momenta of the selected leptons #
    // ###################################################################

    _file->cd();
    TH1F*  hist0_SumptVec3l = (TH1F*)_file->Get("SumVecPtLeptons_FinalCut_WZ_CR__TTbarHiggs");
    hist0_SumptVec3l->Sumw2();
    _file01->cd();
    TH1F*  hist1_SumptVec3l = (TH1F*)_file01->Get("SumVecPtLeptons_FinalCut_WZ_CR__TTbarHiggs");
    hist1_SumptVec3l->Sumw2();
    _file02->cd();
    TH1F*  hist2_SumptVec3l = (TH1F*)_file02->Get("SumVecPtLeptons_FinalCut_WZ_CR__TTbarHiggs");
    hist2_SumptVec3l->Sumw2();
    //TH1F*  hist3_SumptVec3l = (TH1F*)_file->Get("SumptVec3l_ZZ");
    //hist3_SumptVec3l->Sumw2();hist3_SumptVec3l->Rebin(2);
    //TH1F*  hist4_SumptVec3l = (TH1F*)_file->Get("SumptVec3l_VVV");
    //hist4_SumptVec3l->Sumw2();hist4_SumptVec3l->Rebin(2);
    //TH1F*  hist5_SumptVec3l = (TH1F*)_file->Get("SumptVec3l_TT");
    //hist5_SumptVec3l->Sumw2();hist5_SumptVec3l->Rebin(2);
    //TH1F*  hist6_SumptVec3l = (TH1F*)_file->Get("SumptVec3l_TTW");
    //hist6_SumptVec3l->Sumw2();hist6_SumptVec3l->Rebin(2);

    TH1F*  histBkg_SumptVec3l = hist1_SumptVec3l->Clone();
    //histBkg_SumptVec3l->Add(hist2_SumptVec3l);
    //histBkg_SumptVec3l->Add(hist3_SumptVec3l);
    //histBkg_SumptVec3l->Add(hist4_SumptVec3l);
    //histBkg_SumptVec3l->Add(hist5_SumptVec3l);
    //histBkg_SumptVec3l->Add(hist6_SumptVec3l);

    _file1->cd();
    TH1F*  histData_SumptVec3l = (TH1F*)_file1->Get("SumVecPtLeptons_FinalCut_WZ_CR__TTbarHiggs");
    histData_SumptVec3l->Sumw2();
    histData_SumptVec3l->SetMarkerStyle(20);

    //hist0_SumptVec3l->Scale(nWZ.getVal()/hist0_SumptVec3l->Integral());

    THStack hs_SumptVec3l("hs_SumptVec3l","test stacked histograms");

    histBkg_SumptVec3l->SetFillColor(kRed+2);
    hist0_SumptVec3l->SetFillColor(kBlue+2);

    histData_SumptVec3l->GetXaxis()->SetTitle("SumVec p_{T} leptons [GEV]");
    histData_SumptVec3l->GetYaxis()->SetTitle("Events");
    histData_SumptVec3l->GetYaxis()->SetRangeUser(0,15);

    histData_SumptVec3l->Draw("e");
    hs_SumptVec3l.Add(histBkg_SumptVec3l);
    hs_SumptVec3l.Add(hist0_SumptVec3l);
    hs_SumptVec3l.Draw("histosame");
    histData_SumptVec3l->SetLineColor(kRed+2);
    histData_SumptVec3l->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_SumLepCharges, "Data",                    "p");
    qw->AddEntry(hist0_SumLepCharges,    "WZ",                      "f");
    qw->AddEntry(histBkg_SumLepCharges,  "residual backgrounds",    "f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"CMS, #sqrt{s}=13 TeV, 2.2 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.18);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_sumVec3l.pdf");

    // ##############################################
    // # the invariant mass of the selected leptons #
    // ##############################################

    _file->cd();
    TH1F*  hist0_m3l = (TH1F*)_file->Get("InvariantMassOfSelectedLeptons_FinalCut_WZ_CR__TTbarHiggs");
    hist0_m3l->Sumw2();
    _file01->cd();
    TH1F*  hist1_m3l = (TH1F*)_file01->Get("InvariantMassOfSelectedLeptons_FinalCut_WZ_CR__TTbarHiggs");
    hist1_m3l->Sumw2();
    _file02->cd();
    TH1F*  hist2_m3l = (TH1F*)_file02->Get("InvariantMassOfSelectedLeptons_FinalCut_WZ_CR__TTbarHiggs");
    hist2_m3l->Sumw2();
    //TH1F*  hist3_m3l = (TH1F*)_file->Get("m3l_ZZ");
    //hist3_m3l->Sumw2();
    //TH1F*  hist4_m3l = (TH1F*)_file->Get("m3l_VVV");
    //hist4_m3l->Sumw2();
    //TH1F*  hist5_m3l = (TH1F*)_file->Get("m3l_TT");
    //hist5_m3l->Sumw2();
    //TH1F*  hist6_m3l = (TH1F*)_file->Get("m3l_TTW");
    //hist6_m3l->Sumw2();

    TH1F*  histBkg_m3l = hist1_m3l->Clone();
    //histBkg_m3l->Add(hist2_m3l);
    //histBkg_m3l->Add(hist3_m3l);
    //histBkg_m3l->Add(hist4_m3l);
    //histBkg_m3l->Add(hist5_m3l);
    //histBkg_m3l->Add(hist6_m3l);

    _file1->cd();
    TH1F*  histData_m3l = (TH1F*)_file1->Get("InvariantMassOfSelectedLeptons_FinalCut_WZ_CR__TTbarHiggs");
    histData_m3l->Sumw2();
    histData_m3l->SetMarkerStyle(20);

    hist0_m3l->Scale(nWZ.getVal()/hist0_m3l->Integral());

    THStack hs_m3l("hs_m3l","test stacked histograms");

    histBkg_m3l->SetFillColor(kRed+2);
    hist0_m3l->SetFillColor(kBlue+2);

    histData_m3l->GetXaxis()->SetTitle("m3l [GeV]");
    histData_m3l->GetYaxis()->SetTitle("Event Y");
    histData_m3l->GetYaxis()->SetRangeUser(0,20);

    histData_m3l->Draw("e");
    hs_m3l.Add(histBkg_m3l);
    hs_m3l.Add(hist0_m3l);
    hs_m3l.Draw("histosame");
    histData_m3l->SetLineColor(kRed+2);
    histData_m3l->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_m3l, "Data",                   "p");
    qw->AddEntry(hist0_m3l,    "WZ",                     "f");
    qw->AddEntry(histBkg_m3l,  "residual backgrounds",   "f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"CMS, #sqrt{s}=13 TeV, 2.2 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.18);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_m3l.pdf");

    // ##########
    // # met LD #
    // ##########

    _file->cd();
    TH1F*  hist0_metLD = (TH1F*)_file->Get("MetLD_FinalCut_WZ_CR__TTbarHiggs");
    hist0_metLD->Sumw2();
    _file01->cd();
    TH1F*  hist1_metLD = (TH1F*)_file01->Get("MetLD_FinalCut_WZ_CR__TTbarHiggs");
    hist1_metLD->Sumw2();
    _file02->cd();
    TH1F*  hist2_metLD = (TH1F*)_file02->Get("MetLD_FinalCut_WZ_CR__TTbarHiggs");
    hist2_metLD->Sumw2();
    //TH1F*  hist3_metLD = (TH1F*)_file->Get("metLD_ZZ");
    //hist3_metLD->Sumw2();
    //TH1F*  hist4_metLD = (TH1F*)_file->Get("metLD_VVV");
    //hist4_metLD->Sumw2();
    //TH1F*  hist5_metLD = (TH1F*)_file->Get("metLD_TT");
    //hist5_metLD->Sumw2();
    //TH1F*  hist6_metLD = (TH1F*)_file->Get("metLD_TTW");
    //hist6_metLD->Sumw2();


    TH1F*  histBkg_metLD = hist1_metLD->Clone();
    //histBkg_metLD->Add(hist2_metLD);
    //histBkg_metLD->Add(hist3_metLD);
    //histBkg_metLD->Add(hist4_metLD);
    //histBkg_metLD->Add(hist5_metLD);
    //histBkg_metLD->Add(hist6_metLD);

    _file1->cd();
    TH1F*  histData_metLD = (TH1F*)_file1->Get("MetLD_FinalCut_WZ_CR__TTbarHiggs");
    histData_metLD->Sumw2();
    histData_metLD->SetMarkerStyle(20);

    hist0_metLD->Scale(nWZ.getVal()/hist0_metLD->Integral());

    THStack hs_metLD("hs_metLD","test stacked histograms");

    histBkg_metLD->SetFillColor(kRed+2);
    hist0_metLD->SetFillColor(kBlue+2);

    histData_metLD->GetXaxis()->SetTitle("MET LD");
    histData_metLD->GetYaxis()->SetTitle("Events");
    histData_metLD->GetYaxis()->SetRangeUser(0,20);

    histData_metLD->Draw("e");
    hs_metLD.Add(histBkg_metLD);
    hs_metLD.Add(hist0_metLD);
    hs_metLD.Draw("histosame");
    histData_metLD->SetLineColor(kRed+2);
    histData_metLD->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_metLD, "Data" ,          "p");
    qw->AddEntry(hist0_metLD,    "WZ" ,"f");
    qw->AddEntry(histBkg_metLD,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"CMS, #sqrt{s}=13 TeV, 2.2 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.18);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/ttH_3l_metLD.pdf");

    // ############################################
    // # the reconstructed Z boson invariant mass #
    // ############################################

    _file->cd();
    TH1F*  hist0_mZ1 = (TH1F*)_file->Get("ZCandidateInvariantMass_FinalCut_WZ_CR__TTbarHiggs");
    hist0_mZ1->Sumw2();
    _file01->cd();
    TH1F*  hist1_mZ1 = (TH1F*)_file01->Get("ZCandidateInvariantMass_FinalCut_WZ_CR__TTbarHiggs");
    hist1_mZ1->Sumw2();
    _file02->cd();
    TH1F*  hist2_mZ1 = (TH1F*)_file02->Get("ZCandidateInvariantMass_FinalCut_WZ_CR__TTbarHiggs");
    hist2_mZ1->Sumw2();
    //TH1F*  hist3_mZ1 = (TH1F*)_file->Get("mZ1_ZZ");
    //hist3_mZ1->Sumw2();
    //TH1F*  hist4_mZ1 = (TH1F*)_file->Get("mZ1_VVV");
    //hist4_mZ1->Sumw2();
    //TH1F*  hist5_mZ1 = (TH1F*)_file->Get("mZ1_TT");
    //hist5_mZ1->Sumw2();
    //TH1F*  hist6_mZ1 = (TH1F*)_file->Get("mZ1_TTW");
    //hist6_mZ1->Sumw2();

    TH1F*  histBkg_mZ1 = hist1_mZ1->Clone();
    //histBkg_mZ1->Add(hist2_mZ1);
    //histBkg_mZ1->Add(hist3_mZ1);
    //histBkg_mZ1->Add(hist4_mZ1);
    //histBkg_mZ1->Add(hist5_mZ1);
    //histBkg_mZ1->Add(hist6_mZ1);

    _file1->cd();
    TH1F*  histData_mZ1 = (TH1F*)_file1->Get("ZCandidateInvariantMass_FinalCut_WZ_CR__TTbarHiggs");
    histData_mZ1->Sumw2();
    histData_mZ1->SetMarkerStyle(20);

    hist0_mZ1->Scale(nWZ.getVal()/hist0_mZ1->Integral());
    THStack hs_mZ1("hs_mZ1","test stacked histograms");

    histBkg_mZ1->SetFillColor(kRed+2);
    hist0_mZ1->SetFillColor(kBlue+2);

    histData_mZ1->GetXaxis()->SetTitle("best m(l^{+} l^{-}) [GEV]");
    histData_mZ1->GetYaxis()->SetTitle("Events");
    histData_mZ1->GetYaxis()->SetRangeUser(0,20);

    histData_mZ1->Draw("e");
    hs_mZ1.Add(histBkg_mZ1);
    hs_mZ1.Add(hist0_mZ1);
    hs_mZ1.Draw("histosame");
    histData_mZ1->SetLineColor(kRed+2);
    histData_mZ1->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_mZ1, "Data",                 "p");
    qw->AddEntry(hist0_mZ1,    "WZ",                   "f");
    qw->AddEntry(histBkg_mZ1,  "residual backgrounds", "f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"CMS, #sqrt{s}=13 TeV, 2.2 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.18);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/ttH_3l_mZ1.pdf");


    // ########################################################
    // # the transverse momentum of the reconstructed Z boson #
    // ########################################################

    _file->cd();
    TH1F*  hist0_pTZ = (TH1F*)_file->Get("ZCandidateTransverseMomentum_FinalCut_WZ_CR__TTbarHiggs");
    hist0_pTZ->Sumw2();
    _file01->cd();
    TH1F*  hist1_pTZ = (TH1F*)_file01->Get("ZCandidateTransverseMomentum_FinalCut_WZ_CR__TTbarHiggs");
    hist1_pTZ->Sumw2();
    _file02->cd();
    TH1F*  hist2_pTZ = (TH1F*)_file02->Get("ZCandidateTransverseMomentum_FinalCut_WZ_CR__TTbarHiggs");
    hist2_pTZ->Sumw2();
    //TH1F*  hist3_pTZ = (TH1F*)_file->Get("pTZ_ZZ");
    //hist3_pTZ->Sumw2();hist3_pTZ->Rebin(3);
    //TH1F*  hist4_pTZ = (TH1F*)_file->Get("pTZ_VVV");
    //hist4_pTZ->Sumw2();hist4_pTZ->Rebin(3);
    //TH1F*  hist5_pTZ = (TH1F*)_file->Get("pTZ_TT");
    //hist5_pTZ->Sumw2();hist5_pTZ->Rebin(3);
    //TH1F*  hist6_pTZ = (TH1F*)_file->Get("pTZ_TTW");
    //hist6_pTZ->Sumw2();hist6_pTZ->Rebin(3);

    TH1F*  histBkg_pTZ = hist1_pTZ->Clone();
    //histBkg_pTZ->Add(hist2_pTZ);
    //histBkg_pTZ->Add(hist3_pTZ);
    //histBkg_pTZ->Add(hist4_pTZ);
    //histBkg_pTZ->Add(hist5_pTZ);
    //histBkg_pTZ->Add(hist6_pTZ);

    _file1->cd();
    TH1F*  histData_pTZ = (TH1F*)_file1->Get("ZCandidateTransverseMomentum_FinalCut_WZ_CR__TTbarHiggs");
    histData_pTZ->Sumw2();
    histData_pTZ->SetMarkerStyle(20);

    hist0_pTZ->Scale(nWZ.getVal()/hist0_pTZ->Integral());

    THStack hs_pTZ("hs_pTZ","test stacked histograms");

    histBkg_pTZ->SetFillColor(kRed+2);
    hist0_pTZ->SetFillColor(kBlue+2);

    histData_pTZ->GetXaxis()->SetTitle("pT(Z) [GEV]");
    histData_pTZ->GetYaxis()->SetTitle("Events");
    histData_pTZ->GetYaxis()->SetRangeUser(0,15);

    histData_pTZ->Draw("e");
    hs_pTZ.Add(histBkg_pTZ);
    hs_pTZ.Add(hist0_pTZ);
    hs_pTZ.Draw("histosame");
    histData_pTZ->SetLineColor(kRed+2);
    histData_pTZ->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_pTZ, "Data",                 "p");
    qw->AddEntry(hist0_pTZ,    "WZ" ,                  "f");
    qw->AddEntry(histBkg_pTZ,  "residual backgrounds", "f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"CMS, #sqrt{s}=13 TeV, 2.2 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.18);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_pTZ.pdf");
}

{

    // Setup a style
    using namespace RooFit;

    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleSize(0.07, "XYZ");
    gStyle->SetTitleFont(42,"X"); // was 22
    gStyle->SetTitleFont(42,"Y"); // was 22
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetHistLineWidth(2);

    // Grab PDF's from histos

    THStack hs("hs","test stacked histograms");

    // WZ
    TFile *_fileWZ ;
    _fileWZ = TFile::Open("./Input/output_WZ_MC.root");
    _fileWZ->cd();

    TH1F*  histWZ = (TH1F*)_fileWZ->Get("MTW_FinalCut_WZZZ_CR__");
    histWZ->Sumw2();

    // Residual background

    // ZZ
    TFile *_fileZZ ;
    _fileZZ = TFile::Open("./Input/output_ZZ_MC.root");
    _fileZZ->cd();

    TH1F*  histZZ = (TH1F*)_fileZZ->Get("MTW_FinalCut_WZZZ_CR__");
    histZZ->Sumw2();

    // DY
    TFile *_fileDY ;
    _fileDY = TFile::Open("./Input/output_DY_MC.root");
    _fileDY->cd();

    TH1F*  histDY = (TH1F*)_fileDY->Get("MTW_FinalCut_WZZZ_CR__");
    histDY->Sumw2();

    // TTZ
    TFile *_fileTTZ ;
    _fileTTZ = TFile::Open("./Input/output_TTZ_MC.root");
    _fileTTZ->cd();

    TH1F*  histTTZ = (TH1F*)_fileTTZ->Get("MTW_FinalCut_WZZZ_CR__");
    histTTZ->Sumw2();

    // ZZZ
    TFile *_fileZZZ ;
    _fileZZZ = TFile::Open("./Input/output_ZZZ_MC.root");
    _fileZZZ->cd();

    TH1F*  histZZZ = (TH1F*)_fileZZZ->Get("MTW_FinalCut_WZZZ_CR__");
    histZZZ->Sumw2();

    // WWZ
    TFile *_fileWWZ ;
    _fileWWZ = TFile::Open("./Input/output_TTZ_MC.root");
    _fileWWZ->cd();

    TH1F*  histWWZ = (TH1F*)_fileWWZ->Get("MTW_FinalCut_WZZZ_CR__");
    histWWZ->Sumw2();

    // Adding residual backgrounds
    TH1F*  histBkg = histZZ->Clone();
    histBkg->Add(histDY);
    histBkg->Add(histTTZ);
    histBkg->Add(histZZZ);
    histBkg->Add(histWWZ);

    // Get Data histo

    TFile *_fileDATA ;
    _fileDATA = TFile::Open("./Input/output_DATA.root");
    _fileDATA->cd();

    TH1F*  histData = (TH1F*)_fileDATA->Get("MTW_FinalCut_WZZZ_CR__");
    
    histData->Sumw2();
    histData->SetMarkerStyle(20);
    RooRealVar
        x("x","x",histData->GetXaxis()->GetXmin(),histData->GetXaxis()->GetXmax());
    RooDataHist* data = new RooDataHist("data","data",x, histData);

    // Define MC PDFs
    RooDataHist * rooHistWZ  = new RooDataHist("WZ",  "WZ",  x, histWZ );
    RooDataHist * rooHistBkg = new RooDataHist("BKG", "BKG", x, histBkg);

    RooHistPdf histpdf0("histpdf0",   "histpdf0",   x, *rooHistWZ );
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
    RooFitResult* fitResult = pdf->fitTo(*data, Save());

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
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! A DECOMMENTER

    //histBkg->Scale((1-nWZ.getVal())*histData->Integral()/histBkg->Integral());

    std::cout << "scale factor : "<< nWZ.getVal()/histWZ->Integral() << std::endl;
    std::cout << "nWZ MC : "      << histWZ->Integral()             << std::endl;
    std::cout << "yieldBkg "      << yieldBkg                       << std::endl;

    std::cout << "nWZ : " << nWZ.getVal() << " ± " << nWZ.getError() << std::endl;

    //histWZ->Scale(nWZ.getVal()/histWZ->Integral());

    histWZ->SetFillColor(kBlue+2);
    histZZ->SetFillColor(kRed+2);
    histBkg->SetFillColor(kGreen+2);

    TH1F* hmc = (TH1F*) histWZ->Clone();
    hmc->Add(histBkg);
    
    TH1F* hmcHatchedArea = (TH1F*) hmc->Clone();
    hmcHatchedArea->SetLineColor(1);
    hmcHatchedArea->SetFillStyle(3004);
    hmcHatchedArea->SetFillColor(1);

    for(int j=0; j<hmcHatchedArea->GetNbinsX(); j++)
    {
        hmcHatchedArea->SetBinError( j+1, sqrt(pow(hmc->GetBinError(j+1),2)) ); //+pow(3.,2)));
    }

    histData->GetXaxis()->SetTitle("m_{T} W(l)");
    histData->GetXaxis()->SetTitleSize(.05);
    histData->GetYaxis()->SetTitle("Events");
    histData->GetYaxis()->SetTitleSize(.05);
    histData->GetYaxis()->SetRangeUser(0,10);

    histData->Draw("e");
    hs.Add(histBkg);
    hs.Add(histWZ);
    hs.Draw("histo same");
    hmcHatchedArea->Draw("e2same");
    histData->SetLineColor(kRed+2);
    histData->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData, "Data",                  "p");
    qw->AddEntry(histWZ,   "WZ",                    "f");
    qw->AddEntry(histBkg,  "residual backgrounds",  "f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_MTW.pdf");

    fitResult->Print("") ;

    // ##########
    // # MET LD #
    // ##########

    _fileWZ->cd();
    TH1F*  histWZ_METLD = (TH1F*)_fileWZ->Get("MetLD_FinalCut_WZZZ_CR__");
    histWZ_METLD->Sumw2();
    
    _fileZZ->cd();
    TH1F*  histZZ_METLD = (TH1F*)_fileZZ->Get("MetLD_FinalCut_WZZZ_CR__");
    histZZ_METLD->Sumw2();
    
    _fileDY->cd();
    TH1F*  histDY_METLD = (TH1F*)_fileDY->Get("MetLD_FinalCut_WZZZ_CR__");
    histDY_METLD->Sumw2();

    _fileTTZ->cd();
    TH1F*  histTTZ_METLD = (TH1F*)_fileTTZ->Get("MetLD_FinalCut_WZZZ_CR__");
    histTTZ_METLD->Sumw2();

    _fileZZZ->cd();
    TH1F*  histZZZ_METLD = (TH1F*)_fileZZZ->Get("MetLD_FinalCut_WZZZ_CR__");
    histZZZ_METLD->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ_METLD = (TH1F*)_fileWWZ->Get("MetLD_FinalCut_WZZZ_CR__");
    histWWZ_METLD->Sumw2();

    TH1F*  histBkg_METLD = histZZ_METLD->Clone();
    histBkg_METLD->Add(histDY_METLD);
    histBkg_METLD->Add(histTTZ_METLD);
    histBkg_METLD->Add(histZZZ_METLD);
    histBkg_METLD->Add(histWWZ_METLD);
     
    _fileDATA->cd();
    TH1F*  histData_METLD = (TH1F*)_fileDATA->Get("MetLD_FinalCut_WZZZ_CR__");
    histData_METLD->Sumw2();
    histData_METLD->SetMarkerStyle(20);

    //histWZ_METLD->Scale(nWZ.getVal()/histWZ_METLD->Integral());

    THStack hs_METLD("hs_nJet25","test stacked histograms");

    histWZ_METLD->SetFillColor(kBlue+2);
    histZZ_METLD->SetFillColor(kRed+2);
    histBkg_METLD->SetFillColor(kGreen+2);

    histData_METLD->GetXaxis()->SetTitle("E_{T}^{miss LD}");
    histData_METLD->GetXaxis()->SetTitleSize(.05);
    histData_METLD->GetYaxis()->SetTitle("Events");
    histData_METLD->GetYaxis()->SetTitleSize(.05);
    histData_METLD->GetYaxis()->SetRangeUser(0,10);

    histData_METLD->Draw("e");
    hs_METLD.Add(histBkg_METLD);
    hs_METLD.Add(histWZ_METLD);
    hs_METLD.Draw("histosame");
    histData_METLD->SetLineColor(kRed+2);
    histData_METLD->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_METLD, "Data" ,                "p");
    qw->AddEntry(histWZ_METLD,   "WZ" ,                  "f");
    qw->AddEntry(histBkg_METLD,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_METLD.pdf");

    // #################################
    // # the selected jet multiplicity #
    // #################################

    _fileWZ->cd();
    TH1F*  histWZ_nJet25 = (TH1F*)_fileWZ->Get("JetMultiplicity_FinalCut_WZZZ_CR__");
    histWZ_nJet25->Sumw2();
    
    _fileZZ->cd();
    TH1F*  histZZ_nJet25 = (TH1F*)_fileZZ->Get("JetMultiplicity_FinalCut_WZZZ_CR__");
    histZZ_nJet25->Sumw2();
    
    _fileDY->cd();
    TH1F*  histDY_nJet25 = (TH1F*)_fileDY->Get("JetMultiplicity_FinalCut_WZZZ_CR__");
    histDY_nJet25->Sumw2();

    _fileTTZ->cd();
    TH1F*  histTTZ_nJet25 = (TH1F*)_fileTTZ->Get("JetMultiplicity_FinalCut_WZZZ_CR__");
    histTTZ_nJet25->Sumw2();

    _fileZZZ->cd();
    TH1F*  histZZZ_nJet25 = (TH1F*)_fileZZZ->Get("JetMultiplicity_FinalCut_WZZZ_CR__");
    histZZZ_nJet25->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ_nJet25 = (TH1F*)_fileWWZ->Get("JetMultiplicity_FinalCut_WZZZ_CR__");
    histWWZ_nJet25->Sumw2();

    TH1F*  histBkg_nJet25 = histZZ_nJet25->Clone();
    histBkg_nJet25->Add(histDY_nJet25);
    histBkg_nJet25->Add(histTTZ_nJet25);
    histBkg_nJet25->Add(histZZZ_nJet25);
    histBkg_nJet25->Add(histWWZ_nJet25);
     
    _fileDATA->cd();
    TH1F*  histData_nJet25 = (TH1F*)_fileDATA->Get("JetMultiplicity_FinalCut_WZZZ_CR__");
    histData_nJet25->Sumw2();
    histData_nJet25->SetMarkerStyle(20);

    //histWZ_nJet25->Scale(nWZ.getVal()/histWZ_nJet25->Integral());

    THStack hs_nJet25("hs_nJet25","test stacked histograms");

    histWZ_nJet25->SetFillColor(kBlue+2);
    histZZ_nJet25->SetFillColor(kRed+2);
    histBkg_nJet25->SetFillColor(kGreen+2);

    histData_nJet25->GetXaxis()->SetTitle("N(jet, p_{T} > 25)");
    histData_nJet25->GetXaxis()->SetTitleSize(.05);
    histData_nJet25->GetYaxis()->SetTitle("Events");
    histData_nJet25->GetYaxis()->SetTitleSize(.05);
    histData_nJet25->GetYaxis()->SetRangeUser(0,25);

    histData_nJet25->Draw("e");
    hs_nJet25.Add(histBkg_nJet25);
    hs_nJet25.Add(histWZ_nJet25);
    hs_nJet25.Draw("histosame");
    histData_nJet25->SetLineColor(kRed+2);
    histData_nJet25->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_nJet25, "Data" ,                "p");
    qw->AddEntry(histWZ_nJet25,   "WZ" ,                  "f");
    qw->AddEntry(histBkg_nJet25,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_nJet25.pdf");

    // #########################
    // # sum of lepton charges #
    // #########################

     _fileWZ->cd();
    TH1F*  histWZ_SumLepCharges = (TH1F*)_fileWZ->Get("SumOfSelectedLeptonsCharges_FinalCut_WZZZ_CR__");
    histWZ_nJet25->Sumw2();
    
    _fileZZ->cd();
    TH1F*  histZZ_SumLepCharges = (TH1F*)_fileZZ->Get("SumOfSelectedLeptonsCharges_FinalCut_WZZZ_CR__");
    histZZ_nJet25->Sumw2();
    
    _fileDY->cd();
    TH1F*  histDY_SumLepCharges = (TH1F*)_fileDY->Get("SumOfSelectedLeptonsCharges_FinalCut_WZZZ_CR__");
    histDY_nJet25->Sumw2();

    _fileTTZ->cd();
    TH1F*  histTTZ_SumLepCharges = (TH1F*)_fileTTZ->Get("SumOfSelectedLeptonsCharges_FinalCut_WZZZ_CR__");
    histTTZ_nJet25->Sumw2();

    _fileZZZ->cd();
    TH1F*  histZZZ_SumLepCharges = (TH1F*)_fileZZZ->Get("SumOfSelectedLeptonsCharges_FinalCut_WZZZ_CR__");
    histZZZ_nJet25->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ_SumLepCharges = (TH1F*)_fileWWZ->Get("SumOfSelectedLeptonsCharges_FinalCut_WZZZ_CR__");
    histWWZ_nJet25->Sumw2();

    TH1F*  histBkg_SumLepCharges = histZZ_SumLepCharges->Clone();
    histBkg_nJet25->Add(histDY_nJet25);
    histBkg_nJet25->Add(histTTZ_nJet25);
    histBkg_nJet25->Add(histZZZ_nJet25);
    histBkg_nJet25->Add(histWWZ_nJet25);

    _fileDATA->cd();
    TH1F*  histData_SumLepCharges = (TH1F*)_fileDATA->Get("SumOfSelectedLeptonsCharges_FinalCut_WZZZ_CR__");
    histData_SumLepCharges->Sumw2();
    histData_SumLepCharges->SetMarkerStyle(20);

    //histWZ_SumLepCharges->Scale(nWZ.getVal()/histWZ_SumLepCharges->Integral());

    THStack hs_SumLepCharges("hs_SumLepCharges","test stacked histograms");

    histWZ_SumLepCharges->SetFillColor(kBlue+2);
    histZZ_SumLepCharges->SetFillColor(kRed+2);
    histBkg_SumLepCharges->SetFillColor(kGreen+2);

    histData_SumLepCharges->GetXaxis()->SetTitle("sum lepton charges");
    histData_SumLepCharges->GetXaxis()->SetTitleSize(.05);
    histData_SumLepCharges->GetYaxis()->SetTitle("Events");
    histData_SumLepCharges->GetYaxis()->SetTitleSize(.05);
    histData_SumLepCharges->GetYaxis()->SetRangeUser(0,20);

    histData_SumLepCharges->Draw("e");
    hs_SumLepCharges.Add(histBkg_SumLepCharges);
    hs_SumLepCharges.Add(histWZ_SumLepCharges);
    hs_SumLepCharges.Draw("histosame");
    histData_SumLepCharges->SetLineColor(kRed+2);
    histData_SumLepCharges->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_SumLepCharges, "Data" ,                "p");
    qw->AddEntry(histWZ_SumLepCharges,   "WZ" ,                  "f");
    qw->AddEntry(histBkg_SumLepCharges,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_sumLepCharges.pdf");

    // #################
    // # AND THEN SOME #
    // #################

    // ##############################################
    // # vectorial sum of lepton transverse momenta #
    // ##############################################

    _fileWZ->cd();
    TH1F*  histWZ_SVLpt = (TH1F*)_fileWZ->Get("SumVecPtSelectedLeptons_FinalCut_WZZZ_CR__");
    histWZ_SVLpt->Sumw2();
    
    _fileZZ->cd();
    TH1F*  histZZ_SVLpt = (TH1F*)_fileZZ->Get("SumVecPtSelectedLeptons_FinalCut_WZZZ_CR__");
    histZZ_SVLpt->Sumw2();
    
    _fileDY->cd();
    TH1F*  histDY_SVLpt = (TH1F*)_fileDY->Get("SumVecPtSelectedLeptons_FinalCut_WZZZ_CR__");
    histDY_SVLpt->Sumw2();

    _fileTTZ->cd();
    TH1F*  histTTZ_SVLpt = (TH1F*)_fileTTZ->Get("SumVecPtSelectedLeptons_FinalCut_WZZZ_CR__");
    histTTZ_SVLpt->Sumw2();

    _fileZZZ->cd();
    TH1F*  histZZZ_SVLpt = (TH1F*)_fileZZZ->Get("SumVecPtSelectedLeptons_FinalCut_WZZZ_CR__");
    histZZZ_SVLpt->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ_SVLpt = (TH1F*)_fileWWZ->Get("SumVecPtSelectedLeptons_FinalCut_WZZZ_CR__");
    histWWZ_SVLpt->Sumw2();

    TH1F*  histBkg_SVLpt = histZZ_SVLpt->Clone();
    histBkg_SVLpt->Add(histDY_SVLpt);
    histBkg_SVLpt->Add(histTTZ_SVLpt);
    histBkg_SVLpt->Add(histZZZ_SVLpt);
    histBkg_SVLpt->Add(histWWZ_SVLpt);
     
    _fileDATA->cd();
    TH1F*  histData_SVLpt = (TH1F*)_fileDATA->Get("SumVecPtSelectedLeptons_FinalCut_WZZZ_CR__");
    histData_SVLpt->Sumw2();
    histData_SVLpt->SetMarkerStyle(20);

    //histWZ_SVLpt->Scale(nWZ.getVal()/histWZ_SVLpt->Integral());

    THStack hs_SVLpt("hs_SVLpt","test stacked histograms");

    histWZ_SVLpt->SetFillColor(kBlue+2);
    histZZ_SVLpt->SetFillColor(kRed+2);
    histBkg_SVLpt->SetFillColor(kGreen+2);

    histData_SVLpt->GetXaxis()->SetTitle("SumVec p_{T} leptons [GeV]");
    histData_SVLpt->GetXaxis()->SetTitleSize(.05);
    histData_SVLpt->GetYaxis()->SetTitle("Events");
    histData_SVLpt->GetYaxis()->SetTitleSize(.05);
    histData_SVLpt->GetYaxis()->SetRangeUser(0,5);

    histData_SVLpt->Draw("e");
    hs_SVLpt.Add(histBkg_SVLpt);
    hs_SVLpt.Add(histWZ_SVLpt);
    hs_SVLpt.Draw("histosame");
    histData_SVLpt->SetLineColor(kRed+2);
    histData_SVLpt->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_SVLpt, "Data" ,                "p");
    qw->AddEntry(histWZ_SVLpt,   "WZ" ,                  "f");
    qw->AddEntry(histBkg_SVLpt,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_SVLpt.pdf");

    // ######################################
    // # invariant mass of selected leptons #
    // ######################################

    _fileWZ->cd();
    TH1F*  histWZ_Mnl = (TH1F*)_fileWZ->Get("InvariantMassOfSelectedLeptons_FinalCut_WZZZ_CR__");
    histWZ_Mnl->Sumw2();
    
    _fileZZ->cd();
    TH1F*  histZZ_Mnl = (TH1F*)_fileZZ->Get("InvariantMassOfSelectedLeptons_FinalCut_WZZZ_CR__");
    histZZ_Mnl->Sumw2();
    
    _fileDY->cd();
    TH1F*  histDY_Mnl = (TH1F*)_fileDY->Get("InvariantMassOfSelectedLeptons_FinalCut_WZZZ_CR__");
    histDY_Mnl->Sumw2();

    _fileTTZ->cd();
    TH1F*  histTTZ_Mnl = (TH1F*)_fileTTZ->Get("InvariantMassOfSelectedLeptons_FinalCut_WZZZ_CR__");
    histTTZ_Mnl->Sumw2();

    _fileZZZ->cd();
    TH1F*  histZZZ_Mnl = (TH1F*)_fileZZZ->Get("InvariantMassOfSelectedLeptons_FinalCut_WZZZ_CR__");
    histZZZ_Mnl->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ_Mnl = (TH1F*)_fileWWZ->Get("InvariantMassOfSelectedLeptons_FinalCut_WZZZ_CR__");
    histWWZ_Mnl->Sumw2();

    TH1F*  histBkg_Mnl = histZZ_Mnl->Clone();
    histBkg_Mnl->Add(histDY_Mnl);
    histBkg_Mnl->Add(histTTZ_Mnl);
    histBkg_Mnl->Add(histZZZ_Mnl);
    histBkg_Mnl->Add(histWWZ_Mnl);
     
    _fileDATA->cd();
    TH1F*  histData_Mnl = (TH1F*)_fileDATA->Get("InvariantMassOfSelectedLeptons_FinalCut_WZZZ_CR__");
    histData_Mnl->Sumw2();
    histData_Mnl->SetMarkerStyle(20);

    //histWZ_Mnl->Scale(nWZ.getVal()/histWZ_Mnl->Integral());

    THStack hs_Mnl("hs_Mnl","test stacked histograms");

    histWZ_Mnl->SetFillColor(kBlue+2);
    histZZ_Mnl->SetFillColor(kRed+2);
    histBkg_Mnl->SetFillColor(kGreen+2);

    histData_Mnl->GetXaxis()->SetTitle("ml [GeV]");
    histData_Mnl->GetXaxis()->SetTitleSize(.05);
    histData_Mnl->GetYaxis()->SetTitle("Events");
    histData_Mnl->GetYaxis()->SetTitleSize(.05);
    histData_Mnl->GetYaxis()->SetRangeUser(0,10);

    histData_Mnl->Draw("e");
    hs_Mnl.Add(histBkg_Mnl);
    hs_Mnl.Add(histWZ_Mnl);
    hs_Mnl.Draw("histosame");
    histData_Mnl->SetLineColor(kRed+2);
    histData_Mnl->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_Mnl, "Data" ,                "p");
    qw->AddEntry(histWZ_Mnl,   "WZ" ,                  "f");
    qw->AddEntry(histBkg_Mnl,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_Mnl.pdf");

    // ########################################
    // # reconstructed Z boson invariant mass #
    // ########################################

    _fileWZ->cd();
    TH1F*  histWZ_MZpeak = (TH1F*)_fileWZ->Get("ZCandidateInvariantMass_FinalCut_WZZZ_CR__");
    histWZ_MZpeak->Sumw2();
    
    _fileZZ->cd();
    TH1F*  histZZ_MZpeak = (TH1F*)_fileZZ->Get("ZCandidateInvariantMass_FinalCut_WZZZ_CR__");
    histZZ_MZpeak->Sumw2();
    
    _fileDY->cd();
    TH1F*  histDY_MZpeak = (TH1F*)_fileDY->Get("ZCandidateInvariantMass_FinalCut_WZZZ_CR__");
    histDY_MZpeak->Sumw2();

    _fileTTZ->cd();
    TH1F*  histTTZ_MZpeak = (TH1F*)_fileTTZ->Get("ZCandidateInvariantMass_FinalCut_WZZZ_CR__");
    histTTZ_MZpeak->Sumw2();

    _fileZZZ->cd();
    TH1F*  histZZZ_MZpeak = (TH1F*)_fileZZZ->Get("ZCandidateInvariantMass_FinalCut_WZZZ_CR__");
    histZZZ_MZpeak->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ_MZpeak = (TH1F*)_fileWWZ->Get("ZCandidateInvariantMass_FinalCut_WZZZ_CR__");
    histWWZ_MZpeak->Sumw2();

    TH1F*  histBkg_MZpeak = histZZ_MZpeak->Clone();
    histBkg_MZpeak->Add(histDY_MZpeak);
    histBkg_MZpeak->Add(histTTZ_MZpeak);
    histBkg_MZpeak->Add(histZZZ_MZpeak);
    histBkg_MZpeak->Add(histWWZ_MZpeak);
     
    _fileDATA->cd();
    TH1F*  histData_MZpeak = (TH1F*)_fileDATA->Get("ZCandidateInvariantMass_FinalCut_WZZZ_CR__");
    histData_MZpeak->Sumw2();
    histData_MZpeak->SetMarkerStyle(20);

    //histWZ_MZpeak->Scale(nWZ.getVal()/histWZ_MZpeak->Integral());

    THStack hs_MZpeak("hs_MZpeak","test stacked histograms");

    histWZ_MZpeak->SetFillColor(kBlue+2);
    histZZ_MZpeak->SetFillColor(kRed+2);
    histBkg_MZpeak->SetFillColor(kGreen+2);

    histData_MZpeak->GetXaxis()->SetTitle("best Z candidate mass [GeV]");
    histData_MZpeak->GetXaxis()->SetTitleSize(.05);
    histData_MZpeak->GetYaxis()->SetTitle("Events");
    histData_MZpeak->GetYaxis()->SetTitleSize(.05);
    histData_MZpeak->GetYaxis()->SetRangeUser(0,25);

    histData_MZpeak->Draw("e");
    hs_MZpeak.Add(histBkg_MZpeak);
    hs_MZpeak.Add(histWZ_MZpeak);
    hs_MZpeak.Draw("histosame");
    histData_MZpeak->SetLineColor(kRed+2);
    histData_MZpeak->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_MZpeak, "Data" ,                "p");
    qw->AddEntry(histWZ_MZpeak,   "WZ" ,                  "f");
    qw->AddEntry(histBkg_MZpeak,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_MZpeak.pdf");

    // ##########################################
    // # transverse momentum of the Z candidate #
    // ##########################################

    _fileWZ->cd();
    TH1F*  histWZ_ptZ = (TH1F*)_fileWZ->Get("ZCandidateTransverseMomentum_FinalCut_WZZZ_CR__");
    histWZ_ptZ->Sumw2();
    
    _fileZZ->cd();
    TH1F*  histZZ_ptZ = (TH1F*)_fileZZ->Get("ZCandidateTransverseMomentum_FinalCut_WZZZ_CR__");
    histZZ_ptZ->Sumw2();
    
    _fileDY->cd();
    TH1F*  histDY_ptZ = (TH1F*)_fileDY->Get("ZCandidateTransverseMomentum_FinalCut_WZZZ_CR__");
    histDY_ptZ->Sumw2();

    _fileTTZ->cd();
    TH1F*  histTTZ_ptZ = (TH1F*)_fileTTZ->Get("ZCandidateTransverseMomentum_FinalCut_WZZZ_CR__");
    histTTZ_ptZ->Sumw2();

    _fileZZZ->cd();
    TH1F*  histZZZ_ptZ = (TH1F*)_fileZZZ->Get("ZCandidateTransverseMomentum_FinalCut_WZZZ_CR__");
    histZZZ_ptZ->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ_ptZ = (TH1F*)_fileWWZ->Get("ZCandidateTransverseMomentum_FinalCut_WZZZ_CR__");
    histWWZ_ptZ->Sumw2();

    TH1F*  histBkg_ptZ = histZZ_ptZ->Clone();
    histBkg_ptZ->Add(histDY_ptZ);
    histBkg_ptZ->Add(histTTZ_ptZ);
    histBkg_ptZ->Add(histZZZ_ptZ);
    histBkg_ptZ->Add(histWWZ_ptZ);
     
    _fileDATA->cd();
    TH1F*  histData_ptZ = (TH1F*)_fileDATA->Get("ZCandidateTransverseMomentum_FinalCut_WZZZ_CR__");
    histData_ptZ->Sumw2();
    histData_ptZ->SetMarkerStyle(20);

    //histWZ_ptZ->Scale(nWZ.getVal()/histWZ_ptZ->Integral());

    THStack hs_ptZ("hs_ptZ","test stacked histograms");

    histWZ_ptZ->SetFillColor(kBlue+2);
    histZZ_ptZ->SetFillColor(kRed+2);
    histBkg_ptZ->SetFillColor(kGreen+2);

    histData_ptZ->GetXaxis()->SetTitle("p_{T}( Z candidate )");
    histData_ptZ->GetXaxis()->SetTitleSize(.05);
    histData_ptZ->GetYaxis()->SetTitle("Events");
    histData_ptZ->GetYaxis()->SetTitleSize(.05);
    histData_ptZ->GetYaxis()->SetRangeUser(0,15);

    histData_ptZ->Draw("e");
    hs_ptZ.Add(histBkg_ptZ);
    hs_ptZ.Add(histWZ_ptZ);
    hs_ptZ.Draw("histosame");
    histData_ptZ->SetLineColor(kRed+2);
    histData_ptZ->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_ptZ, "Data" ,                "p");
    qw->AddEntry(histWZ_ptZ,   "WZ" ,                  "f");
    qw->AddEntry(histBkg_ptZ,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_ptZ.pdf");

    // ######################
    // # AND THEN EVEN MORE #
    // ######################

    // ############
    // # cut flow #
    // ############

    _fileWZ->cd();
    TH1F*  histWZ_CutFlow = (TH1F*)_fileWZ->Get("CutFlow_noSel_WZZZ_CR__");
    histWZ_CutFlow->Sumw2();
    
    _fileZZ->cd();
    TH1F*  histZZ_CutFlow = (TH1F*)_fileZZ->Get("CutFlow_noSel_WZZZ_CR__");
    histZZ_CutFlow->Sumw2();
    
    _fileDY->cd();
    TH1F*  histDY_CutFlow = (TH1F*)_fileDY->Get("CutFlow_noSel_WZZZ_CR__");
    histDY_CutFlow->Sumw2();

    _fileTTZ->cd();
    TH1F*  histTTZ_CutFlow = (TH1F*)_fileTTZ->Get("CutFlow_noSel_WZZZ_CR__");
    histTTZ_CutFlow->Sumw2();

    _fileZZZ->cd();
    TH1F*  histZZZ_CutFlow = (TH1F*)_fileZZZ->Get("CutFlow_noSel_WZZZ_CR__");
    histZZZ_CutFlow->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ_CutFlow = (TH1F*)_fileWWZ->Get("CutFlow_noSel_WZZZ_CR__");
    histWWZ_CutFlow->Sumw2();

    TH1F*  histBkg_CutFlow = histZZ_CutFlow->Clone();
    histBkg_CutFlow->Add(histDY_CutFlow);
    histBkg_CutFlow->Add(histTTZ_CutFlow);
    histBkg_CutFlow->Add(histZZZ_CutFlow);
    histBkg_CutFlow->Add(histWWZ_CutFlow);
     
    _fileDATA->cd();
    TH1F*  histData_CutFlow = (TH1F*)_fileDATA->Get("CutFlow_noSel_WZZZ_CR__");
    histData_CutFlow->Sumw2();
    histData_CutFlow->SetMarkerStyle(20);

    //histWZ_CutFlow->Scale(nWZ.getVal()/histWZ_CutFlow->Integral());

    THStack hs_CutFlow("hs_CutFlow","test stacked histograms");

    histWZ_CutFlow->SetFillColor(kBlue+2);
    histZZ_CutFlow->SetFillColor(kRed+2);
    histDY_CutFlow->SetFillColor(kGray+1);
    histTTZ_CutFlow->SetFillColor(kGray+2);
    histZZZ_CutFlow->SetFillColor(kGray+3);
    histWWZ_CutFlow->SetFillColor(kGray+4);
    histBkg_CutFlow->SetFillColor(kGreen+2);

    histData_CutFlow->GetXaxis()->SetTitle("Cut #");
    histData_CutFlow->GetXaxis()->SetTitleSize(.05);
    histData_CutFlow->GetYaxis()->SetTitle("Events");
    histData_CutFlow->GetYaxis()->SetTitleSize(.05);
    histData_CutFlow->GetYaxis()->SetRangeUser(0,60);

    histData_CutFlow->Draw("e");
    //hs_CutFlow.Add(histBkg_CutFlow);
    hs_CutFlow.Add(histWZ_CutFlow);
    hs_CutFlow.Add(histZZ_CutFlow);
    hs_CutFlow.Add(histDY_CutFlow);
    hs_CutFlow.Add(histTTZ_CutFlow);
    hs_CutFlow.Add(histZZZ_CutFlow);
    hs_CutFlow.Add(histWWZ_CutFlow);
    hs_CutFlow.Draw("histosame");
    histData_CutFlow->SetLineColor(kRed+2);
    histData_CutFlow->Draw("esame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData_CutFlow, "Data" ,                "p");
    qw->AddEntry(histWZ_CutFlow,   "WZ" ,                  "f");
    qw->AddEntry(histZZ_CutFlow,   "ZZ" ,                  "f");
    qw->AddEntry(histDY_CutFlow,   "DY" ,                  "f");
    qw->AddEntry(histTTZ_CutFlow,  "TTZ" ,                 "f");
    qw->AddEntry(histZZZ_CutFlow,  "ZZZ" ,                 "f");
    qw->AddEntry(histWWZ_CutFlow,  "WWZ" ,                 "f");
    //qw->AddEntry(histBkg_CutFlow,  "residual backgrounds" ,"f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    //text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                  2.26 fb^{-1} (13TeV)");
    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       2.26 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    //text1->SetX(0.18);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("plots/WZ_3l_CutFlow.pdf");


}

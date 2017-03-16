{
    // ##########
    // # MET LD #
    // ##########
    
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

 
     TFile *_fileWZ ;
    _fileWZ = TFile::Open("./Input/output_WZ.root");
     TFile *_fileZZ ;
    _fileZZ = TFile::Open("./Input/output_ZZ.root");
     TFile *_fileDY ;
    _fileDY = TFile::Open("./Input/output_DY.root");
     TFile *_fileWWZ ;
    _fileWWZ = TFile::Open("./Input/output_WWZ.root");
     TFile *_fileTTZ ;
    _fileTTZ = TFile::Open("./Input/output_ttZ.root");
     TFile *_fileDATA ;
    _fileDATA = TFile::Open("./Input/output_data.root");


    _fileWZ->cd();
    TH1F*  histWZ = (TH1F*)_fileWZ->Get("MetLD_AN_finalSel_WZ_CR__");
    histWZ->Sumw2();
    histWZ->SetFillColor(kMagenta-6);
    histWZ->SetLineColor(kMagenta-6);
 
    _fileZZ->cd();
    TH1F*  histZZ = (TH1F*)_fileZZ->Get("MetLD_AN_finalSel_WZ_CR__");
    histZZ->Sumw2();
    histZZ->SetFillColor(kAzure+6);
    histZZ->SetLineColor(kAzure+6);
 
    _fileDY->cd();
    TH1F*  histDY = (TH1F*)_fileDY->Get("MetLD_AN_finalSel_WZ_CR__");
    histDY->Sumw2();
    histDY->SetFillColor(kAzure+6);
    histDY->SetLineColor(kAzure+6);
 
    
    _fileTTZ->cd();
    TH1F*  histTTZ = (TH1F*)_fileTTZ->Get("MetLD_AN_finalSel_WZ_CR__");
    histTTZ->Sumw2();
    histTTZ->SetFillColor(kAzure+6);
    histTTZ->SetLineColor(kAzure+6);
 
   // _fileZZZ->cd();
   // TH1F*  histZZZ = (TH1F*)_fileZZZ->Get("MetLD_AN_finalSel_WZ_CR__");
    //histZZZ->Sumw2();

    _fileWWZ->cd();
    TH1F*  histWWZ = (TH1F*)_fileWWZ->Get("MetLD_AN_finalSel_WZ_CR__");
    histWWZ->Sumw2();
    histWWZ->SetFillColor(kAzure+6);
    histWWZ->SetLineColor(kAzure+6);
 
    //TH1F*  histBkg = histZZ->Clone();
    //histBkg->Add(histDY);
    //histBkg->Add(histTTZ);
    //histBkg->Add(histZZZ);
    //histBkg->Add(histWWZ);
     
    _fileDATA->cd();
    TH1F*  histData = (TH1F*)_fileDATA->Get("MetLD_AN_finalSel_WZ_CR__");
    histData->Sumw2();
    histData->SetMarkerStyle(20);

    histWZ->Scale(0.96);
    
    TH1F*  histDataFake = (TH1F*)_fileDATA->Get("MetLD_AN_Fakes_finalSel_WZ_CR__");
    histDataFake->Sumw2();
    histDataFake->SetMarkerStyle(20);
    histDataFake->SetFillColor(1); 
    histDataFake->SetLineColor(1);
    histDataFake->SetFillStyle(3325);
    
    
    TH1F* hmc = (TH1F*) histWZ->Clone();
    hmc->Add(histDataFake);
    hmc->Add(histWWZ);
    hmc->Add(histZZ); 
    hmc->Add(histDY);  
    hmc->Add(histTTZ);
    
    TH1F* hmcHatchedArea = (TH1F*) hmc->Clone();
    hmcHatchedArea->SetLineColor(1);
    hmcHatchedArea->SetFillStyle(3013);
    hmcHatchedArea->SetFillColor(1);

    for(int j=0; j<hmcHatchedArea->GetNbinsX(); j++)
    {
        hmcHatchedArea->SetBinError( j+1, sqrt(pow(hmc->GetBinError(j+1),2)) ); //+pow(3.,2)));
    }

    histData->GetXaxis()->SetTitle("E_{T}^{miss LD} in 3l + 2 non-b jets");
    histData->GetXaxis()->SetTitleSize(.05);
    histData->GetYaxis()->SetTitle("Events");
    histData->GetYaxis()->SetTitleSize(.05);
    histData->GetYaxis()->SetRangeUser(0,100);
    
    histData->SetMaximum(150);
    
    histData->Draw("e"); 
    THStack hs("hs","test stacked histograms");
    hs.Add(histDataFake);
    hs.Add(histWWZ);
    hs.Add(histZZ); 
    hs.Add(histDY);  
    hs.Add(histTTZ);
    hs.Add(histWZ); 
    hs.Draw("histo same");
    hmcHatchedArea->Draw("e2same");
    histData->SetLineColor(kBlack);
    histData->Draw("esame");

  
    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(histData, "Data",                  "p");
    qw->AddEntry(histWZ,   "WZ",                    "f");
    qw->AddEntry(histWWZ,  "rares",  "f");
    qw->AddEntry(histDataFake,  "fakes",  "f");

    qw->SetFillColor(0);
    qw->SetTextFont(42);
    qw->SetLineWidth(0);
    qw->SetBorderSize(0);

    qw->Draw();

    text1 = new TLatex(0.15,0.93,"#bf{CMS} #it{Preliminary},                       35.9 fb^{-1} (13TeV)");
    text1->SetNDC();
    text1->SetTextAlign(12);
    text1->SetX(0.16);
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->Draw();

    c1->SaveAs("WZ_3l_METLD.pdf");
    

}

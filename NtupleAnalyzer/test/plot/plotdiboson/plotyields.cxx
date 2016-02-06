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

    TH1F*  hist0 = (TH1F*)_file->Get("CutFlow_noSel_WZ_CR__TTbarHiggs");
    hist0->Sumw2();

    // Other contributions
    // Sorted by number of events: WZ, ZZ, VVV, DY, TT, TTW, TTH

    TFile *_file01 ;
    _file01 = TFile::Open("./Input/output_ZZ_MC.root");
    _file01->cd();

    TH1F*  hist1 = (TH1F*)_file01->Get("CutFlow_noSel_WZ_CR__TTbarHiggs");
    hist1->Sumw2();

    TFile *_file02 ;
    _file02 = TFile::Open("./Input/output_ZZ_MC.root");
    _file02->cd();

    TH1F*  hist2 = (TH1F*)_file02->Get("CutFlow_noSel_WZ_CR__TTbarHiggs");
    hist2->Sumw2();

    TH1F*  histBkg = hist1->Clone();
    //histBkg->Add(hist2);

    histBkg->SetFillColor(kRed+2);
    hist0->SetFillColor(kBlue+2);
    hist1->SetFillColor(kGray+2);
    hist2->SetFillColor(kGreen+2);

    histBkg->GetXaxis()->SetTitle("Cuts");
    histBkg->GetXaxis()->SetRangeUser(0,10);
    histBkg->GetYaxis()->SetTitle("Yields (#events)");
    histBkg->GetYaxis()->SetRangeUser(0,1.1);

    TH1F*  histNorm = hist0->Clone();
    histNorm->Add(hist1);
    histNorm->Add(hist2);

    //histBkg->Draw("e");
    //hs.Add(hist1);
    //hs.Draw("histosame");
    //hs.Add(hist2);
    //hs.Draw("histosame");
    //hs.Add(hist0);
    //hs.Draw("histosame");
    //histBkg->Draw("esame");

    histBkg.Divide(histNorm);
    hist0.Divide(histNorm);
    hist1.Divide(histNorm);
    hist2.Divide(histNorm);

    histBkg->Draw("hbar");
    hs.Add(hist2);
    hs.Add(hist1);
    hs.Add(hist0);
    hs.Draw("hbarsame");

    TLegend* qw = 0;
    qw = new TLegend(0.6,0.657163,0.85,0.85,NULL,"brNDC");

    qw->AddEntry(hist0,    "WZ",                    "f");
    qw->AddEntry(hist1,    "ZZ",                    "f");
    qw->AddEntry(hist2,    "DY",                    "f");

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

    c1->SaveAs("plots/WZ_yields.pdf");
}

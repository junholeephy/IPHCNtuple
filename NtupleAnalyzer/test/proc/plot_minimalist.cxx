#include "PlotStyle.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TMath.h"
#include "TFile.h"
#include "TLegend.h"

#include <vector>
#include <iostream>

void addbin(TH1F*);
void PUreweighting(TH1F *h_data, TH1F *h_MC);

int main()
{
    gROOT->SetBatch(1);

    SetPlotStyle();

    gStyle->SetOptFit();

    TFile *fout = new TFile("hist.root","RECREATE");

    // MC
        // files

    TFile *f_WZJets;
    f_WZJets          = TFile::Open("./InputFile/output_WZJets_MC.root");
        std::cout << "File MC opened" << std::endl;

        // histos
    
    TH1F *h_WZJets_nPV;
    TH1F *h_WZJets;
    
    h_WZJets_nPV = (TH1F*)f_WZJets->Get("NomberOfPrimaryVertex_noSel_emu__TTbarHiggs");
    h_WZJets     = (TH1F*)f_WZJets->Get("ZCandidateInvariantMass_noSel_emu__TTbarHiggs");


    // Data
        // files
    
    TFile *f_sig_doubleEG;
    //f_sig_doubleEG    = TFile::Open("./InputFile/output_DoubleEG_DATA.root");
    f_sig_doubleEG    = TFile::Open("./InputFile/output_DoubleMuon_DATA.root");
        std::cout << "File DATA opened" << std::endl;

        // histos

    TH1F *h_sig_doubleEG_nPV;
    TH1F *h_sig_doubleEG;

    h_sig_doubleEG_nPV = (TH1F*)f_sig_doubleEG->Get("NomberOfPrimaryVertex_noSel_emu__TTbarHiggs");    
    h_sig_doubleEG     = (TH1F*)f_sig_doubleEG->Get("ZCandidateInvariantMass_noSel_emu__TTbarHiggs");

    fout->cd();

    TCanvas *c1 = new TCanvas();

    // distributions WZ

    h_WZJets_nPV->SetMarkerSize(0);
    h_WZJets_nPV->SetMarkerColor(kGreen-1);
    h_WZJets_nPV->SetLineColor(kGreen-1);
    h_WZJets_nPV->SetFillColor(kGreen-1);
    h_WZJets_nPV->SetLineStyle(1);

    h_WZJets->SetMarkerSize(0);
    h_WZJets->SetMarkerColor(kGreen-1);
    h_WZJets->SetLineColor(kGreen-1);
    h_WZJets->SetFillColor(kGreen-1);
    h_WZJets->SetLineStyle(1);

    std::cout << "Fin de la cosmétique h_WZJets" << std::endl;

    // distributions double EG

    h_sig_doubleEG_nPV->SetMarkerSize(0);
    h_sig_doubleEG_nPV->SetMarkerColor(kGray);
    h_sig_doubleEG_nPV->SetLineColor(kGray);
    h_sig_doubleEG_nPV->SetFillColor(kGray);
    h_sig_doubleEG_nPV->SetLineStyle(1);

    h_sig_doubleEG->SetMarkerSize(0);
    h_sig_doubleEG->SetMarkerColor(kGray);
    h_sig_doubleEG->SetLineColor(kGray);
    h_sig_doubleEG->SetFillColor(kGray);
    h_sig_doubleEG->SetLineStyle(1);

    std::cout << "Fin de la cosmétique h_sig_doubleEG" << std::endl;
   
    // PU reweighting

    TCanvas c2("c2", "PU", 300, 200);
    c2.Divide(2,2);
    c2.cd(1);
    h_sig_doubleEG_nPV->Draw();
    c2.cd(2);
    h_WZJets_nPV->Draw();
    c2.cd(3);
    PUreweighting(h_sig_doubleEG_nPV, h_WZJets_nPV);
    h_sig_doubleEG_nPV->Draw();
    c2.Print("PUreweighting.pdf");

    float WZ_int = h_WZJets->Integral(0,h_WZJets->GetXaxis()->GetNbins()+1);

    // Normalisation background

    float bg_int = WZ_int;

    std::string h_WZJets_name = "h_WZJets_whatever";
    TH1F * h_WZJets_orig = (TH1F*)h_WZJets->Clone(h_WZJets_name.c_str());
    std::string h_WZJets_name_norm = "h_WZJets_whatever_norm";
    TH1F * h_WZJets_norm = (TH1F*)h_WZJets->Clone(h_WZJets_name_norm.c_str());
    if(bg_int > 0) h_WZJets_norm->Scale(1./bg_int);

    std::cout << "End of background normalisation" << std::endl;

    // Normalisation signal

    float sig_int = h_sig_doubleEG->Integral(0,h_sig_doubleEG->GetXaxis()->GetNbins()+1);

    std::string h_sig_doubleEG_name = "h_sig_doubleEG_whatever";
    TH1F * h_sig_doubleEG_orig = (TH1F*)h_sig_doubleEG->Clone(h_sig_doubleEG_name.c_str());
    std::string h_sig_doubleEG_name_norm = "h_sig_doubleEG_whatever_norm";
    TH1F * h_sig_doubleEG_norm = (TH1F*)h_sig_doubleEG->Clone(h_sig_doubleEG_name_norm.c_str());
    if(sig_int > 0) h_sig_doubleEG_norm->Scale(1./sig_int);

    std::string h_bg_name_merged = "bg_whatever_merged";
    TH1F *h_bg_merged = (TH1F*)h_WZJets->Clone(h_bg_name_merged.c_str());
    //h_bg_merged->Add(h_WZJets);

    std::cout << "End of signal normalisation" << std::endl;

    THStack *h_bg_norm = new THStack();
    h_bg_norm->Add(h_WZJets);

    THStack *h_bg = new THStack();
    h_bg->Add(h_WZJets);

    // IL normalised
    h_bg->Draw("hist e1");
    h_bg->GetXaxis()->SetTitle("Invariant Mass");
    h_bg->GetYaxis()->SetTitle("Number of events");

    h_sig_doubleEG->Scale(1);
    h_sig_doubleEG->Draw("hist e1 same");

    h_sig_doubleEG->Scale(1);
    h_sig_doubleEG->Draw("hist e1 same");

    std::string h_sig_name_merged = "sig_wathever_merged";
    TH1F *h_sig_merged = (TH1F*)h_sig_doubleEG->Clone(h_sig_name_merged.c_str());
    //h_sig_merged->Add(h_sig_hut_ttbar);

    TLegend *leg = new TLegend(0.70,0.92,0.92,0.60);
    leg->SetFillColor(253);
    leg->SetBorderSize(0);
    leg->AddEntry(h_WZJets,        "WZ",                      "f");
    leg->AddEntry(h_sig_doubleEG,  "Double EG",               "lp");
    leg->Draw();

    c1->SetLogy(1); 
    h_bg->SetMinimum(0.001);

    std::string pname = "pics/whatever.eps";
    c1->Print(pname.c_str());
    c1->Clear();

    // unity normalised
    h_bg_norm->Draw("hist e1");
    h_bg_norm->GetXaxis()->SetTitle("Invariant mass?");
    h_bg_norm->GetYaxis()->SetTitle("");

    h_sig_doubleEG_norm->Draw("hist e1 same");
    h_sig_doubleEG_norm->Draw("hist e1 same");

    float bg_norm_max = h_bg_norm->GetMaximum();
    float h_sig_doubleEG_max = h_sig_doubleEG_norm->GetMaximum();
    //float h_sig_doubleEG_max = h_sig_doubleEG_norm->GetMaximum();

    float max_sig = std::max(h_sig_doubleEG_max,h_sig_doubleEG_max);
    float max = std::max(bg_norm_max,max_sig);
    h_bg_norm->SetMaximum(1.2*max);

    leg->Draw();

    c1->SetLogy(0);

    std::string pname_norm = "pics/whatever_norm.eps";
    c1->Print(pname_norm.c_str());
    c1->Clear();

    delete leg;
    delete c1;

    fout->Write();   
    fout->Close();
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

void PUreweighting(TH1F *h_data, TH1F *h_MC)
{
    h_data->Divide(h_MC);
}

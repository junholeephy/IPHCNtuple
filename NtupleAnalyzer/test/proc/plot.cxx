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

int main()
{
    gROOT->SetBatch(1);

    SetPlotStyle();

    gStyle->SetOptFit();

    TFile *fout = new TFile("hist.root","RECREATE");

    std::vector<std::string> ttbar;
    std::vector<std::string> stop;
    std::vector<std::string> higgs;
    std::vector<std::string> ttv;
    std::vector<std::string> thq;
    std::vector<std::string> tzq;
    std::vector<std::string> wjets;
    std::vector<std::string> zjets;

    std::vector<std::string> sig_hut_stop;
    std::vector<std::string> sig_hct_stop;
    std::vector<std::string> sig_gut;
    std::vector<std::string> sig_gct;

    std::vector<std::string> sig_hut_ttbar;
    std::vector<std::string> sig_hct_ttbar;

    std::vector<std::string> sig_hut;
    std::vector<std::string> sig_hct;

    ttbar.push_back("TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola");

    stop.push_back("TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola");
    stop.push_back("TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola");
    stop.push_back("T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola");
    stop.push_back("Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola");
    stop.push_back("TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola");
    stop.push_back("TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola");

    higgs.push_back("VBF_HToBB_M-125_13TeV-powheg-pythia6");
    higgs.push_back("WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp");
    ////   higgs.push_back("ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp");

    ttv.push_back("TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola");
    ttv.push_back("TTWJets_Tune4C_13TeV-madgraph-tauola");
    ttv.push_back("TTZJets_Tune4C_13TeV-madgraph-tauola");

    thq.push_back("THQtoTBBQ_1L");

    tzq.push_back("TZQtoTQQQ_1L");

    wjets.push_back("WJetsToLNu_13TeV-madgraph-pythia8-tauola");

    zjets.push_back("DYJetsToLL_M-50_13TeV-madgraph-pythia8");

    sig_hut_stop.push_back("tHToBB_1L_Kappa_hut");
    sig_hct_stop.push_back("tHToBB_1L_Kappa_hct");

    sig_gut.push_back("tHToBB_1L_Kappa_gut");
    sig_gct.push_back("tHToBB_1L_Kappa_gct");

    sig_hut_ttbar.push_back("TTtoTHToBB_1L_Kappa_hut");
    sig_hct_ttbar.push_back("TTtoTHToBB_1L_Kappa_hct");

    sig_hut.push_back("tHToBB_1L_Kappa_hut");
    sig_hut.push_back("TTtoTHToBB_1L_Kappa_hut");

    sig_hct.push_back("tHToBB_1L_Kappa_hct");
    sig_hct.push_back("TTtoTHToBB_1L_Kappa_hct");   

    float ksig = 1;

    const int nvar = 34;
    std::string var[nvar] =
    {
        "H_m",
        "H_pt",
        "H_eta",
        "top_m",
        "top_pt",
        "top_eta",
        "W_m",
        "W_pt",
        "W_eta",
        "l_pt",
        "l_eta",
        "nu_pt",
        "nu_eta",
        "Hb_pt",
        "Hb_eta",
        "topb_pt",
        "topb_eta",
        "Hb1_Hb2_dr",
        "H_top_dr",
        "l_nu_dr",
        "W_topb_dr",
        "Hb1_Hb2_costheta",
        "H_top_costheta",
        "l_nu_costheta",
        "W_topb_costheta",
        "HT",
        "MET",
        "njet",
        "nu_phi",
        "MET_phi",
        "H_nu_dr",
        "H_l_dr",
        "chi2",
        "l_charge"
    };   

    std::string xtit[nvar] =
    {
        "m(H) [GeV]",
        "p_{T} (H) [GeV]",
        "#eta (H)",
        "m(top) [GeV]",
        "p_{T} (top) [GeV]",
        "#eta (top)",
        "m(W) [GeV]",
        "p_{T} (W) [GeV]",
        "#eta (W)",
        "p_{T} (l) [GeV]",
        "#eta (l)",
        "p_{T} (#nu) [GeV]",
        "#eta (#nu)",
        "p_{T} (b_{H}) [GeV]",
        "#eta (b_{H})",
        "p_{T} (b_{top}) [GeV]",
        "#eta (b_{top})",
        "#Delta R (b_{H1},b_{H2})",
        "#Delta R (H,top)",
        "#Delta R (l,#nu)",
        "#Delta R (W,b_{top})",
        "cos #theta (b_{H1},b_{H2})",
        "cos #theta (H,top)",
        "cos #theta (l,#nu)",
        "cos #theta (W,b_{top})",
        "H_{T} [GeV]",
        "E_{T}^{miss} [GeV]",
        "Number of jets",
        "#phi (#nu)",
        "#phi (MET)",
        "#Delta R (H,#nu)",
        "#Delta R (H,l)",
        "#chi^{2}",
        "lepton charge"
    };   

    const int nchan = 2;
    std::string chan[nchan] = 
    {
        "e",
        "m"
    };   

    const int njet = 1;
    std::string jet[njet] =
    {
        "bjge3"
    };   

    const int nsel = 1;
    std::string sel[nsel] =
    {
        "nosel"
    };   

    std::string hname[1000];
    int varid[1000];

    int idx = 0;
    for(int iv=0;iv<nvar;iv++)
    {
        for(int ic=0;ic<nchan;ic++)
        {
            for(int ij=0;ij<njet;ij++)
            {
                for(int is=0;is<nsel;is++)
                {
                    hname[idx] = "h_"+var[iv]+"_"+chan[ic]+"_"+jet[ij]+"_"+sel[is]+"_ALL";
                    varid[idx] = iv;
                    idx++;
                }
            }
        } 
    }

    TH1F *h_ttbar;   
    TFile *f_ttbar;   

    TH1F *h_stop;
    TFile *f_stop;   

    TH1F *h_ttv;
    TFile *f_ttv;   

    TH1F *h_zjets;
    TFile *f_zjets;   

    TH1F *h_wjets;
    TFile *f_wjets;   

    TH1F *h_higgs;
    TFile *f_higgs;   

    TH1F *h_thq;
    TFile *f_thq;   

    TH1F *h_tzq;
    TFile *f_tzq;   

    TH1F *h_sig_hut;
    TFile *f_sig_hut;

    TH1F *h_sig_hut_stop;
    TFile *f_sig_hut_stop;

    TH1F *h_sig_hut_ttbar;
    TFile *f_sig_hut_ttbar;

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_ttbar = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_ttbar->Get(hname.c_str());
    h_ttbar = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_zjets = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_zjets->Get(hname.c_str());
    h_zjets = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_wjets = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_wjets->Get(hname.c_str());
    h_wjets = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_stop = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_stop->Get(hname.c_str());
    h_stop = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_ttv = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_ttv->Get(hname.c_str());
    h_ttv = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_higgs = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_higgs->Get(hname.c_str());
    h_higgs = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_thq = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_thq->Get(hname.c_str());
    h_thq = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_tzq = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_tzq->Get(hname.c_str());
    h_tzq = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_sig_hut = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_sig_hut->Get(hname.c_str());
    h_sig_hut = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_sig_hut_stop = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_sig_hut_stop->Get(hname.c_str());
    h_sig_hut_stop = (TH1F*)h->Clone(hname_clone.c_str());

    std::string fname = "./InputFile/output_WZJets_MC.root";
    f_sig_hut_ttbar = TFile::Open(fname.c_str());
    std::string hname_clone = hname+"_clone";
    TH1F *h = (TH1F*)f_sig_hut_ttbar->Get(hname.c_str());
    h_sig_hut_ttbar = (TH1F*)h->Clone(hname_clone.c_str());

    fout->cd();

    TCanvas *c1 = new TCanvas();

    addbin(h_ttbar);
    addbin(h_stop);
    addbin(h_ttv);
    addbin(h_higgs);
    addbin(h_zjets);
    addbin(h_wjets);
    addbin(h_thq);
    addbin(h_tzq);
    addbin(h_sig_hut);
    addbin(h_sig_hut_stop);
    addbin(h_sig_hut_ttbar);

    h_zjets->SetMarkerSize(0);
    h_zjets->SetMarkerColor(94);
    h_zjets->SetLineColor(94);
    h_zjets->SetFillColor(94);
    h_zjets->SetLineStyle(1);

    h_ttbar->SetMarkerSize(0);
    h_ttbar->SetMarkerColor(kAzure+1);
    h_ttbar->SetLineColor(kAzure+1);
    h_ttbar->SetFillColor(kAzure+1);
    h_ttbar->SetLineStyle(1);

    h_stop->SetMarkerSize(0);
    h_stop->SetMarkerColor(kGreen);
    h_stop->SetLineColor(kGreen);
    h_stop->SetFillColor(kGreen);
    h_stop->SetLineStyle(1);

    h_ttv->SetMarkerSize(0);
    h_ttv->SetMarkerColor(kOrange);
    h_ttv->SetLineColor(kOrange);
    h_ttv->SetFillColor(kOrange);
    h_ttv->SetLineStyle(1);

    h_higgs->SetMarkerSize(0);
    h_higgs->SetMarkerColor(kYellow);
    h_higgs->SetLineColor(kYellow);
    h_higgs->SetFillColor(kYellow);
    h_higgs->SetLineStyle(1);

    h_wjets->SetMarkerSize(0);
    h_wjets->SetMarkerColor(kMagenta-7);
    h_wjets->SetLineColor(kMagenta-7);
    h_wjets->SetFillColor(kMagenta-7);
    h_wjets->SetLineStyle(1);

    h_thq->SetMarkerSize(0);
    h_thq->SetMarkerColor(kGreen+1);
    h_thq->SetLineColor(kGreen+1);
    h_thq->SetFillColor(kGreen+1);
    h_thq->SetLineStyle(1);

    h_tzq->SetMarkerSize(0);
    h_tzq->SetMarkerColor(kGray);
    h_tzq->SetLineColor(kGray);
    h_tzq->SetFillColor(kGray);
    h_tzq->SetLineStyle(1);

    h_sig_hut->SetMarkerSize(0);
    h_sig_hut->SetMarkerColor(kRed);
    h_sig_hut->SetLineColor(kRed);
    h_sig_hut->SetFillColor(0);
    h_sig_hut->SetLineWidth(2);

    h_sig_hut_stop->SetMarkerSize(0);
    h_sig_hut_stop->SetMarkerColor(kRed+1);
    h_sig_hut_stop->SetLineColor(kRed+1);
    h_sig_hut_stop->SetFillColor(0);
    h_sig_hut_stop->SetLineWidth(2);

    h_sig_hut_ttbar->SetMarkerSize(0);
    h_sig_hut_ttbar->SetMarkerColor(kRed-7);
    h_sig_hut_ttbar->SetLineColor(kRed-7);
    h_sig_hut_ttbar->SetFillColor(0);
    h_sig_hut_ttbar->SetLineWidth(2);

    float thq_int = h_thq->Integral(0,h_thq->GetXaxis()->GetNbins()+1);
    float tzq_int = h_tzq->Integral(0,h_tzq->GetXaxis()->GetNbins()+1);
    float ttbar_int = h_ttbar->Integral(0,h_ttbar->GetXaxis()->GetNbins()+1);
    float stop_int = h_stop->Integral(0,h_stop->GetXaxis()->GetNbins()+1);
    float ttv_int = h_ttv->Integral(0,h_ttv->GetXaxis()->GetNbins()+1);
    float zjets_int = h_zjets->Integral(0,h_zjets->GetXaxis()->GetNbins()+1);
    float wjets_int = h_wjets->Integral(0,h_wjets->GetXaxis()->GetNbins()+1);
    float higgs_int = h_higgs->Integral(0,h_higgs->GetXaxis()->GetNbins()+1);

    float bg_int = thq_int+tzq_int+ttbar_int+stop_int+ttv_int+zjets_int+wjets_int+higgs_int;

    std::string h_thq_name = "thq_"+hname;
    TH1F * h_thq_orig = (TH1F*)h_thq->Clone(h_thq_name.c_str());
    std::string h_thq_name_norm = "thq_"+hname+"_norm";
    TH1F * h_thq_norm = (TH1F*)h_thq->Clone(h_thq_name_norm.c_str());
    if(bg_int > 0) h_thq_norm->Scale(1./bg_int);

    std::string h_tzq_name = "tzq_"+hname;
    TH1F * h_tzq_orig = (TH1F*)h_tzq->Clone(h_tzq_name.c_str());
    std::string h_tzq_name_norm = "tzq_"+hname+"_norm";
    TH1F * h_tzq_norm = (TH1F*)h_tzq->Clone(h_tzq_name_norm.c_str());
    if(bg_int > 0) h_tzq_norm->Scale(1./bg_int);

    std::string h_ttbar_name = "ttbar_"+hname;
    TH1F * h_ttbar_orig = (TH1F*)h_ttbar->Clone(h_ttbar_name.c_str());
    std::string h_ttbar_name_norm = "ttbar_"+hname+"_norm";
    TH1F * h_ttbar_norm = (TH1F*)h_ttbar->Clone(h_ttbar_name_norm.c_str());
    if(bg_int > 0) h_ttbar_norm->Scale(1./bg_int);

    std::string h_stop_name = "stop_"+hname;
    TH1F * h_stop_orig = (TH1F*)h_stop->Clone(h_stop_name.c_str());
    std::string h_stop_name_norm = "stop_"+hname+"_norm";
    TH1F * h_stop_norm = (TH1F*)h_stop->Clone(h_stop_name_norm.c_str());
    if(bg_int > 0) h_stop_norm->Scale(1./bg_int);

    std::string h_ttv_name = "ttv_"+hname;
    TH1F * h_ttv_orig = (TH1F*)h_ttv->Clone(h_ttv_name.c_str());
    std::string h_ttv_name_norm = "ttv_"+hname+"_norm";
    TH1F * h_ttv_norm = (TH1F*)h_ttv->Clone(h_ttv_name_norm.c_str());
    if(bg_int > 0) h_ttv_norm->Scale(1./bg_int);

    std::string h_zjets_name = "zjets_"+hname;
    TH1F * h_zjets_orig = (TH1F*)h_zjets->Clone(h_zjets_name.c_str());
    std::string h_zjets_name_norm = "zjets_"+hname+"_norm";
    TH1F * h_zjets_norm = (TH1F*)h_zjets->Clone(h_zjets_name_norm.c_str());
    if(bg_int > 0) h_zjets_norm->Scale(1./bg_int);

    std::string h_wjets_name = "wjets_"+hname;
    TH1F * h_wjets_orig = (TH1F*)h_wjets->Clone(h_wjets_name.c_str());
    std::string h_wjets_name_norm = "wjets_"+hname+"_norm";
    TH1F * h_wjets_norm = (TH1F*)h_wjets->Clone(h_wjets_name_norm.c_str());
    if(bg_int > 0) h_wjets_norm->Scale(1./bg_int);

    std::string h_higgs_name = "higgs_"+hname;
    TH1F * h_higgs_orig = (TH1F*)h_higgs->Clone(h_higgs_name.c_str());
    std::string h_higgs_name_norm = "higgs_"+hname+"_norm";
    TH1F * h_higgs_norm = (TH1F*)h_higgs->Clone(h_higgs_name_norm.c_str());
    if(bg_int > 0) h_higgs_norm->Scale(1./bg_int);

    float sig_hut_stop_int = h_sig_hut_stop->Integral(0,h_sig_hut_stop->GetXaxis()->GetNbins()+1);
    float sig_hut_ttbar_int = h_sig_hut_ttbar->Integral(0,h_sig_hut_ttbar->GetXaxis()->GetNbins()+1);

    std::string h_sig_hut_stop_name = "sig_hut_stop_"+hname;
    TH1F * h_sig_hut_stop_orig = (TH1F*)h_sig_hut_stop->Clone(h_sig_hut_stop_name.c_str());
    std::string h_sig_hut_stop_name_norm = "sig_hut_stop_"+hname+"_norm";
    TH1F * h_sig_hut_stop_norm = (TH1F*)h_sig_hut_stop->Clone(h_sig_hut_stop_name_norm.c_str());
    if(sig_hut_stop_int > 0) h_sig_hut_stop_norm->Scale(1./sig_hut_stop_int);

    std::string h_sig_hut_ttbar_name = "sig_hut_ttbar_"+hname;
    TH1F * h_sig_hut_ttbar_orig = (TH1F*)h_sig_hut_ttbar->Clone(h_sig_hut_ttbar_name.c_str());
    std::string h_sig_hut_ttbar_name_norm = "sig_hut_ttbar_"+hname+"_norm";
    TH1F * h_sig_hut_ttbar_norm = (TH1F*)h_sig_hut_ttbar->Clone(h_sig_hut_ttbar_name_norm.c_str());
    if(sig_hut_ttbar_int > 0) h_sig_hut_ttbar_norm->Scale(1./sig_hut_ttbar_int);

    std::string h_bg_name_merged = "bg_"+hname+"_merged";
    TH1F *h_bg_merged = (TH1F*)h_thq->Clone(h_bg_name_merged.c_str());
    h_bg_merged->Add(h_tzq);
    h_bg_merged->Add(h_ttbar);
    h_bg_merged->Add(h_stop);
    h_bg_merged->Add(h_ttv);
    h_bg_merged->Add(h_zjets);
    h_bg_merged->Add(h_wjets);
    h_bg_merged->Add(h_higgs);

    THStack *h_bg_norm = new THStack();
    h_bg_norm->Add(h_thq_norm);
    h_bg_norm->Add(h_tzq_norm);
    h_bg_norm->Add(h_ttbar_norm);
    h_bg_norm->Add(h_stop_norm);
    h_bg_norm->Add(h_ttv_norm);
    h_bg_norm->Add(h_zjets_norm);
    h_bg_norm->Add(h_wjets_norm);
    h_bg_norm->Add(h_higgs_norm);

    THStack *h_bg = new THStack();
    h_bg->Add(h_thq);
    h_bg->Add(h_tzq);
    h_bg->Add(h_ttbar);
    h_bg->Add(h_stop);
    h_bg->Add(h_ttv);
    h_bg->Add(h_zjets);
    h_bg->Add(h_wjets);
    h_bg->Add(h_higgs);

    //int itit = varid;
    //std::string xtitle = xtit[itit];

    // IL normalised
    h_bg->Draw("hist e1");
    h_bg->GetXaxis()->SetTitle(xtitle.c_str());
    h_bg->GetYaxis()->SetTitle("Number of events");

    h_sig_hut_stop->Scale(ksig);
    h_sig_hut_stop->Draw("hist e1 same");

    h_sig_hut_ttbar->Scale(ksig);
    h_sig_hut_ttbar->Draw("hist e1 same");

    std::string h_sig_name_merged = "sig_"+hname+"_merged";
    TH1F *h_sig_merged = (TH1F*)h_sig_hut_stop->Clone(h_sig_name_merged.c_str());
    h_sig_merged->Add(h_sig_hut_ttbar);

    TLegend *leg = new TLegend(0.70,0.92,0.92,0.60);
    leg->SetFillColor(253);
    leg->SetBorderSize(0);
    leg->AddEntry(h_ttbar,         "ttbar",                      "f");
    leg->AddEntry(h_stop,          "Single top",                 "f");
    leg->AddEntry(h_zjets,         "Z+jets",                     "f");
    leg->AddEntry(h_ttv,           "ttbarV",                     "f");
    leg->AddEntry(h_wjets,         "W+jets",                     "f");
    leg->AddEntry(h_higgs,         "Higgs",                      "f");
    leg->AddEntry(h_thq,           "tHq",                        "f");
    leg->AddEntry(h_tzq,           "tZq",                        "f");
    leg->AddEntry(h_sig_hut_stop,  "FCNC #kappa_{Hut}",         "lp");
    leg->AddEntry(h_sig_hut_ttbar, "FCNC #kappa_{Hut} (ttbar)", "lp");
    leg->Draw();

    c1->SetLogy(1); 
    h_bg->SetMinimum(0.001);

    std::string pname = "pics/"+hname+".eps";
    c1->Print(pname.c_str());
    c1->Clear();

    // unity normalised
    h_bg_norm->Draw("hist e1");
    h_bg_norm->GetXaxis()->SetTitle(xtitle.c_str());
    h_bg_norm->GetYaxis()->SetTitle("");

    h_sig_hut_stop_norm->Draw("hist e1 same");
    h_sig_hut_ttbar_norm->Draw("hist e1 same");

    float bg_norm_max = h_bg_norm->GetMaximum();
    float sig_hut_stop_max = h_sig_hut_stop_norm->GetMaximum();
    float sig_hut_ttbar_max = h_sig_hut_ttbar_norm->GetMaximum();

    float max_sig = std::max(sig_hut_stop_max,sig_hut_ttbar_max);
    float max = std::max(bg_norm_max,max_sig);
    h_bg_norm->SetMaximum(1.2*max);

    leg->Draw();

    c1->SetLogy(0);

    std::string pname_norm = "pics/"+hname+"_norm.eps";
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

#include "include/NtupleProducer.h"

ClassImp(Truth)

Truth::Truth()
{
}

Truth::~Truth()
{
}

void Truth::read()
{
    int UNINT = -666;

    int gen_n = ntP->gen_n;
    int gen_n_sel = 0;
    for(int i=0;i<gen_n;i++)
    {
        int status = ntP->gen_status->at(i);
        if( status != 1 && status != 3 ) continue;
        _gen_pt.push_back(ntP->gen_pt->at(i));
        _gen_eta.push_back(ntP->gen_eta->at(i));
        _gen_phi.push_back(ntP->gen_phi->at(i));
        _gen_m.push_back(ntP->gen_m->at(i));
        _gen_id.push_back(ntP->gen_id->at(i));
        _gen_status.push_back(status);

        int mom = ntP->gen_mother_index->at(i);

        if( mom >= ntP->gen_index->size() )
        {

            std::cout << "Problem with gen length" << std::endl;
            exit(1);
        }

        int mother_id = ntP->gen_id->at(mom);

        _gen_mother_id.push_back(mother_id);

        gen_n_sel++;
    }

    _mc_truth_n = 0;

    int mc_truth_h0W1_id = ntP->mc_truth_h0W1_id;
    if( mc_truth_h0W1_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0W1_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0W1_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0W1_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0W1_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0W1_E);
        _mc_truth_label.push_back(12);
        _mc_truth_n++;
    }
    int mc_truth_h0Wl1_id = ntP->mc_truth_h0Wl1_id;
    if( mc_truth_h0Wl1_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wl1_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wl1_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wl1_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wl1_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wl1_E);
        _mc_truth_label.push_back(120);
        _mc_truth_n++;
    }
    int mc_truth_h0Wq11_id = ntP->mc_truth_h0Wq11_id;
    if( mc_truth_h0Wq11_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wq11_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wq11_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wq11_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wq11_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wq11_E);
        _mc_truth_label.push_back(122);
        _mc_truth_n++;
    }
    int mc_truth_h0Wq21_id = ntP->mc_truth_h0Wq21_id;
    if( mc_truth_h0Wq21_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wq21_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wq21_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wq21_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wq21_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wq21_E);
        _mc_truth_label.push_back(123);
        _mc_truth_n++;
    }
    int mc_truth_h0W2_id = ntP->mc_truth_h0W2_id;
    if( mc_truth_h0W2_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0W2_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0W2_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0W2_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0W2_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0W2_E);
        _mc_truth_label.push_back(13);
        _mc_truth_n++;
    }
    int mc_truth_h0Wl2_id = ntP->mc_truth_h0Wl2_id;
    if( mc_truth_h0Wl2_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wl2_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wl2_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wl2_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wl2_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wl2_E);
        _mc_truth_label.push_back(130);
        _mc_truth_n++;
    }
    int mc_truth_h0Wq12_id = ntP->mc_truth_h0Wq12_id;
    if( mc_truth_h0Wq12_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wq12_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wq12_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wq12_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wq12_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wq12_E);
        _mc_truth_label.push_back(132);
        _mc_truth_n++;
    }
    int mc_truth_h0Wq22_id = ntP->mc_truth_h0Wq22_id;
    if( mc_truth_h0Wq22_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wq22_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wq22_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wq22_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wq22_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wq22_E);
        _mc_truth_label.push_back(133);
        _mc_truth_n++;
    }
    int mc_truth_h0Z1_id = ntP->mc_truth_h0Z1_id;
    if( mc_truth_h0Z1_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Z1_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Z1_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Z1_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Z1_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Z1_E);
        _mc_truth_label.push_back(14);
        _mc_truth_n++;
    }
    int mc_truth_h0Zl11_id = ntP->mc_truth_h0Zl11_id;
    if( mc_truth_h0Zl11_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zl11_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zl11_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zl11_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zl11_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zl11_E);
        _mc_truth_label.push_back(140);
        _mc_truth_n++;
    }
    int mc_truth_h0Zl21_id = ntP->mc_truth_h0Zl21_id;
    if( mc_truth_h0Zl21_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zl21_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zl21_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zl21_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zl21_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zl21_E);
        _mc_truth_label.push_back(141);
        _mc_truth_n++;
    }
    int mc_truth_h0Zq11_id = ntP->mc_truth_h0Zq11_id;
    if( mc_truth_h0Zq11_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zq11_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zq11_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zq11_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zq11_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zq11_E);
        _mc_truth_label.push_back(142);
        _mc_truth_n++;
    }
    int mc_truth_h0Zq21_id = ntP->mc_truth_h0Zq21_id;
    if( mc_truth_h0Zq21_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zq21_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zq21_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zq21_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zq21_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zq21_E);
        _mc_truth_label.push_back(143);
        _mc_truth_n++;
    }
    int mc_truth_h0Z2_id = ntP->mc_truth_h0Z2_id;
    if( mc_truth_h0Z2_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Z2_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Z2_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Z2_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Z2_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Z2_E);
        _mc_truth_label.push_back(15);
        _mc_truth_n++;
    }
    int mc_truth_h0Zl12_id = ntP->mc_truth_h0Zl12_id;
    if( mc_truth_h0Zl12_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zl12_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zl12_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zl12_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zl12_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zl12_E);
        _mc_truth_label.push_back(150);
        _mc_truth_n++;
    }
    int mc_truth_h0Zl22_id = ntP->mc_truth_h0Zl22_id;
    if( mc_truth_h0Zl22_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zl22_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zl22_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zl22_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zl22_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zl22_E);
        _mc_truth_label.push_back(151);
        _mc_truth_n++;
    }
    int mc_truth_h0Zq12_id = ntP->mc_truth_h0Zq12_id;
    if( mc_truth_h0Zq12_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zq12_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zq12_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zq12_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zq12_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zq12_E);
        _mc_truth_label.push_back(152);
        _mc_truth_n++;
    }
    int mc_truth_h0Zq22_id = ntP->mc_truth_h0Zq22_id;
    if( mc_truth_h0Zq22_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zq22_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zq22_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zq22_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zq22_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zq22_E);
        _mc_truth_label.push_back(153);
        _mc_truth_n++;
    }
    int mc_truth_h0tau1_id = ntP->mc_truth_h0tau1_id;
    if( mc_truth_h0tau1_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0tau1_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0tau1_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0tau1_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0tau1_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0tau1_E);
        _mc_truth_label.push_back(16);
        _mc_truth_n++;
    }
    int mc_truth_h0tau2_id = ntP->mc_truth_h0tau2_id;
    if( mc_truth_h0tau2_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0tau2_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0tau2_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0tau2_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0tau2_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0tau2_E);
        _mc_truth_label.push_back(17);
        _mc_truth_n++;	
    }
}

void Truth::init()
{
    _mc_truth_n = 0;

    _mc_truth_id.clear();
    _mc_truth_pt.clear();
    _mc_truth_eta.clear();
    _mc_truth_phi.clear();
    _mc_truth_E.clear();
    _mc_truth_label.clear();

    _gen_pt.clear();
    _gen_eta.clear();
    _gen_phi.clear();
    _gen_m.clear();
    _gen_id.clear();
    _gen_status.clear();
    _gen_mother_id.clear();
}


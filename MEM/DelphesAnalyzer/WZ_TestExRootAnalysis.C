#include <iostream>

#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"
#include "externel/ExRootAnalysis/ExRootTreeReader.h"
#include "WZ_TestExRootAnalysis.h"
#include "MultiLeptonTree.h"

using namespace std;

void WZ_TestExRootAnalysis()
{
	TChain chain("Delphyes");
	
	chain.add("");
	chain.add("");

	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64t numberOfEntries = treeReader->GetEntries();
	cout<<"Total Entry number of Events : "<<numberOfEntries<<endl;

	TClonesArray *branchEvent = treeReader->UseBranch("Event");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");        // what is the difference between this and GenMissingET????
	TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
	TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");

//	std::vector<Jet*>	vSelectedJets;
	std::vector<Lepton*>	vSelectedLeptons;
//	std::vector<GenParticle*>	vSelectedPartons;
	std::vector<GenParticle*>	vSelectedNeutrinos;
	std::vector<GenParticle*>	vSelectedGenLeptons;
	std::vector<GenParticle*>	vSelectedW;
	std::vector<GenParticle*>	vSelectedNeutrinoFromW;
	std::vector<GenParticle*>	vSelectedZ;
	
//	std::vector<std::pair<Jet*, Genparticle*>> MatchedJets;

	MissingET *mET = NULL;
	ScalarHT *msumET = NULL;

	double Ptot = 0;
	double Ptot_Px=0, Ptot_Py=0, Ptot_Pz=0, Ptot_Eta=0, Ptot_E=0, Ptot_M=0;
    TLorentzVector Ptot_P4(0,0,0,0);

	unsigned int nSelectedEventsWZwithGen = 0;
	unsigned int nSelectedEventsWZ = 0;
	unsigned int nSelectedEvent = 0;

	TFile *fOutput = new TFile("output.root","RECREATE");
	MultiLeptonTree tree;
	tree.initializeOutputTree();

	bool doWZselection = true;
	bool doTFselection = false;
	bool storeGenLevel = true;

	
	for(Int_ entry = 0; entry<numberOfEntries; ++entry)
	{
		treeReader->ReadEntry(entry);
		if(entry % 1000 == 0) cout<<"Now looping "<<entry<<"th event"<<endl;

//		vSelectedJets.clear();
		vSelectedGenLeptons.clear();
		vSelectedLeptons.clear();
//		vSelectedPartons.clear();
		vSelectedW.clear();
		vSelectedNeutrinosFromW.clear();
		vSelectedNeutrinos.clear();
		vSelectedZ.clear();

		GenParticle *genparticle;
		Ptot_Px=0, Ptot_Py=0, Ptot_Pz=0, Ptot_Eta=0, Ptot_E=0, Ptot_M=0;

		if(branchGenParticle->GetEntries()>0)
		{
			for(int ig=0; ig<branchGenParticle->GetEntries(); ig++)
			{
				genparticle = (GenParticle*)vranchGenParticle->At(ig);

				//Picking Neutrinos, and gathering its mother W boson and its by product a lepton
				if(abs(genparticle->PID)==12 || abs(genparticle->PID)==14 || abs(genparticle->PID)==16)
				{
					vSelectedNeutrinos.push_back(genparticle);
					Ptot_Px += genparticle->Px;
					Ptot_Py += genparticle->Py;
					Ptot_Pz += genparticle->Pz;
					Ptot_E  += genparticle->E;

					int m = genparticle->M1;
					GenParticle *genmother = (GenParticle*) branchGenParticle->At(m);
					if(abs(genmother->PID)==24)
					{
						vSelectedW.push_back(genmother);
						vSelectedNeutrinoFromW.push_back(genparticle);
						int md1 = genmother->D1;
						int md2 = genmother->D2;
						GenParticle *gendaughter1 = (GenParticle*)branchGenParticle->At(md1);
						GenParticle *gendaughter1 = (GenParticle*)branchGenParticle->At(md2);
						if(abs(gendaughter1->PID)==11 ||abs(gendaughter1->PID)==13 || abs(gendaughter1->PID)==15)
						{
							vSelectedGenLeptons.push_bach(gendaughter1);
							Ptot_Px += genparticle->Px;
							Ptot_Py += genparticle->Py;
							Ptot_Pz += genparticle->Pz;
							Ptot_E  += genparticle->E;
						}
						else if(abs(gendaughter2->PID)==11 || abs(gendaughter2->PID)==13 || abs(gendaughter2->PID)==15)
						{
							vSelectedGenLeptons.pusch_back(gendaughter2);
							Ptot_Px += genparticle->Px;
							Ptot_Py += genparticle->Py;
							Ptot_Pz += genparticle->Pz;
							Ptot_E  += genparticle->E;
						}
					}
				}


				//For picked Leptons, whose mother is Z boson. 
				if(abs(genparticle->PID)==11 || abs(genparticle->PID)==13 || abs(genparticle->PID)==15)
				{
					int ml = genparticle->M1;
					GenParticle *genmother = (GenParticle*)branchGenParticle->At(ml);
					if(abs(genmother->PID)==22)    // Need Z boson, ????? need to check PID of Z boson
					{
						vSelectedZ.push_back(genmother);
						int md1 = genmother->D1;
						int md2 = genmother->D2;
						GenParticle *gendaughter1 = (GenParticle*)branchGenParticle->At(md1);
						GenParticle *gendaughter2 = (GenParticle*)branchGenParticle->At(md2);
						vSelectedGenLeptons.push_back(gendaughter1);
						vSelectedGenLeptons.push_back(gendaughter2);

						Ptot_Px += gendaughter1->Px;
						Ptot_Px += gendaughter2->Px;
						Ptot_Py += gendaughter1->Py;
						Ptot_Py += gendaughter2->Py;
						Ptot_Pz += gendaughter1->Pz;
						Ptot_Pz += gendaughter2->Pz;
						Ptot_E  += gendaughter1->E;
						Ptot_E  += gendaughter2->E;
					}
				}
			}
		} //End of Looping GenParticles

		Ptot_P4.SetPxPyPzE(Ptot_Px, Ptot_Py, Ptot_Pz, Ptot_E);
		Ptot = sqrt(Ptot_Px*Ptot_Px + Ptot_Py*Ptot_Py);  // transverse momemtum


		//selectiong reco Leptons
		if(branchElectron->GetEntries() >0)
		{
			for(int ie=0; ie<branchElectron->GetEntries(); ie++)
			{
				Lepton *elec = new Lepton((Electron*)branchElectron->At(ie));
				if(elec->PT > 15 && fabs(elec->Eta)<2.4)
				{
					vSelectedLeptons.push_back(elec);
				}
			}
		}
		if(branchMuon->GetEntries() > 0)
		{
			for(int im=0; im<branchMuon->GetEntries(); im++)
			{
				Lepton *muon = new Lepton((Muon*)branchMuon->At(im));
				if(muon->PT >15 && fabs(muon->Eta)<2.4)
				{
					vSelectedLeptons.pusch_back(muon);
				}
			}
		}
		std::sort(vSelectedLeptons.begin(), vSelectedLeptons.end(), SortingLeptonPt);

		if(branchMissingET->GetEntries() == 1)
		{
			mET = (MissingET*)branchMissingET->At(0);
		}
		if(branchScalarHT->GetEntries() == 1)
		{
			msumET = (ScalarHT*)branchScalarHT->At(0);
		}

		
		int tot_charge = 0;
		int tot_id = 0;
		if(vSelectedLeptons.size()>=4)
		{
			for(unsigned int i = 0; i<vSelectedLeptons.size(); i++)
			{
				tot_charge += vSelectedLeptons.at(i)->charge;
				tot_id += vSelectedLeptons.at(i)->Id;
			}
		}

		bool is3l = false;
		if(vSelectedLeptons.size()==3) is3l = true;

		// is_Z boson
		bool is_Z = false;
		if(vSelectedLeptons.size()>=2)
		{
			for(unsigned int i=0; i<vSelectedLeptons.szie()-1; i++)
			{
				TLorentzVector p1; 
				p1.SetPtEtaPhiM(vSelectedLeptons.at(i)->PT, vSelectedLeptons.at(i)->Eta, vSelectedLeptons.at(i)->Phi, 0);
				for(unsigned int j=i+1; i<vSelectedLeptons.size(); j++)
				{
					TLorentzVector P2;
					P2.SetPtEtaPhiM(vSelectedLetons.at(j)->PT, vSelectedLeptons.at(j)->Eta, vSelectedLeptons.at(j)->Phi, 0);
					if(fabs((p1+p2).M())<12 &&(vSelectedLeptons.at(i)->Charge) == -(vSelectedLeptons.at(j)->Charge) && (vSelectedLeptons.at(i)->Id) == -(vSelectedLeptons.at(j)->Id) && (fabs(p1+p2).M() - 91.188) < 15 )
					{
						is_Z = true;
					}
				}
			}
		}


		if(doWZselection)
		{
			if(vSelectedLeptions.size()<3) continue;
			if(!is_Z) continue;

			nSelectedEventsWZ++;

			if(storeGenLevel)
			{
				if(vSelectedW.size()!=1) continue;
				if(vSelectedNeutrinosFromW.size()!=1) continue;

				if(vSelectedNeutrinosFromW.size()>0)
				{
					TLorentzVector Pneutrino;
					Pneutrino.SetPtEtaPhiM(vSelectedNeutrinosFromW.at(0)->PT, vSelectedNeutrinosFromW.at(0)->Eta, vSelectedNeutrinosFromW.at(0)->Phi, 0);
					tree.multilepton_mET_Matched = Pneutrino;
			
					if(vSelectedW.size()>=1)
					{
						TLorentzVector PW1;
						PW1.SetPtEtaPhiM(vSelectedW.at(0)->PT, vSelectedW.at(0)->Eta, vSelectedW.at(0)->Phi, vSelectedW.at(0)->Mass);
						tree.multilepton_W1_P4_Matched = PW1;
					}
					tree_multilepton_GenLepton1_P4.SetPtEtaPhiM(vSelectedGenLeptons.at(0)->PT, vSelectedGenLeptons.at(0)->Eta, vSelectedGenLeptons.at(0)->Phi, 0);
					tree_multilepton_GenLepton2_P4.SetPtEtaPhiM(vSelectedGenLeptons.at(1)->PT, vSelectedGenLeptons.at(1)->Eta, vSelectedGenLeptons.at(1)->Phi, 0);
					tree_multilepton_GenLepton3_P4.SetPtEtaPhiM(vSelectedGenLeptons.at(2)->PT, vSelectedGenLeptons.at(2)->Eta, vSelectedGenLeptons.at(2)->Phi, 0);
				}
				
				nSelectedEventsWZwithGen++;
			}

			tree.multilepton_Lepton1_Id = -999;
    		tree.multilepton_Lepton2_Id = -999;
   		 	tree.multilepton_Lepton3_Id = -999;
			if(vSelectedLeptons.size()>=3)
			{
				tree.multilepton_Lepton1_P4.SetPtEtaPhiM(vSelectedLeptons.at(0)->PT, vSelectedLeptons.at(0)->Eta, vSelectedLeptons.at(0)->Phi,0);
				tree.multilepton_Lepton1_Id = vSelectedLeptons.at(0)->Id;
				tree.multilepton_Lepton2_P4.SetPtEtaPhiM(vSelectedLeptons.at(1)->PT, vSelectedLeptons.at(1)->Eta, vSelectedLeptons.at(1)->Phi,0);
				tree.multilepton_Lepton2_Id = vSelectedLeptons.at(1)->Id;
				tree.multilepton_Lepton3_P4.SetPtEtaPhiM(vSelectedLeptons.at(2)->PT, vSelectedLeptons.at(2)->Eta, vSelectedLeptons.at(2)->Phi,0);
				tree.multilepton_Lepton3_Id = vSelectedLeptons.at(2)->Id;
			}

			tree.multilepton_mET.SetPtEtaPhiM(mET->MET, 0, mET->Phi, mET->MET);
			tree.multilepton_mHT = msumET->HT;
			tree.tOutput->Fill();
		
		}	
		nSelectedEvent++;
	
	}

	cout << "---RESULTS---"<<endl;
	cout << "nSelectedEvent="<<nSelectedEvent<<endl;
	cout << "nSelectedEventsWZ="<<nSelectedEventsWZ<<endl;
	cout << "nSelectedEventsWZwithGen="<<nSelectedEventsWZwithGen<<endl;
		
	tree.tOutput->Write();
	fOutput->Close();
}			






















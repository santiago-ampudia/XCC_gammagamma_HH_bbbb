#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TRandom.h"
//#include "vector"
#include <vector>
#include <set>
#include <stdexcept>
#endif

//------------------------------------------------------------------------------
/////function that receives eta and returns theta
double findTheta(double eta)
{
	double theta = pow(TMath::E(), -eta);
	theta = 2*TMath::ATan(theta);
	
	return theta;
}

/////Function that gets max pt to be used for the detector radius
void findDetectorRadiusXY(TClonesArray *branchEFlowTrack, TClonesArray *branchGenJet, int nTracks, int nGenJets, double& ptMax, string level)
{
	TObject *object;
  	GenParticle *particleConst;
  	Track *trackConst;
  	Tower *towerConst;
  	
  	if(level == "gen")
  	{
  		for(int i=0;i<nGenJets;i++)
		{
			Jet *genjet = (Jet*) branchGenJet->At(i); // Take ith jet
			for(int j = 0; j < genjet->Constituents.GetEntriesFast(); ++j)
		        {
				object = genjet->Constituents.At(j);
				if(object == 0) continue;
				if(object->IsA() == GenParticle::Class())
				{
					particleConst = (GenParticle*) object;
				  	double ptSen=particleConst->PT;
				  	if(ptSen>ptMax) ptMax=ptSen;
				}
			}
		}
  	}
  	else if(level == "reco")
  	{
  		for(int i=0;i<nTracks;i++)
		{
			Track *track = (Track*) branchEFlowTrack->At(i); // Take ith track
		        if(ptMax<(track->PT)) ptMax=track->PT;
		}
  	}
}

/////Function that gets max phi and theta to be used for the detector radius
void findDetectorRadiusPhiTheta(TClonesArray *branchEFlowTrack, TClonesArray *branchGenJet, int nTracks, int nGenJets, double& phiMax, double& thetaMax, string level)
{
	TObject *object;
  	GenParticle *particleConst;
  	Track *trackConst;
  	Tower *towerConst;
  	
  	if(level == "gen")
  	{
  		for(int i=0;i<nGenJets;i++)
		{
			Jet *genjet = (Jet*) branchGenJet->At(i); // Take ith jet
			for(int j = 0; j < genjet->Constituents.GetEntriesFast(); ++j)
		        {
				object = genjet->Constituents.At(j);
				if(object == 0) continue;
				if(object->IsA() == GenParticle::Class())
				{
					particleConst = (GenParticle*) object;
				  	double phiSen=particleConst->Phi;
				  	if(phiSen>phiMax) phiMax=phiSen;
				  	double thetaSen=particleConst->Eta;
				  	thetaSen = findTheta(thetaSen);
				  	if(thetaSen>thetaMax) thetaMax=thetaSen;
				}
			}
		}
  	}
  	else if(level == "reco")
  	{
  		for(int i=0;i<nTracks;i++)
		{
			Track *track = (Track*) branchEFlowTrack->At(i); // Take ith track
		        if(phiMax<(track->Phi)) phiMax=track->Phi;
		        double thetaSen = track->Eta;
		        thetaSen = findTheta(thetaSen);
		        if(thetaMax<thetaSen) thetaMax=thetaSen;
		}
  	}
}

//Function that receives the true particles branch, a particle, and a number 1 or 2 saying if it should return the PID of the daughter 1 or the daughter 2. If it has no daughter for the number requested, it reutrns -999.
int findDaughter(TClonesArray *branchParticle, GenParticle *particle, int d) 
{	
  	int dIndex;
  	
  	if(d==1) 
  	{
  		if(particle->D1 == -1) return -999;
  		else dIndex=particle->D1;
  	}
  	if(d==2) 
  	{
  		if(particle->D2 == -1) return -999;
  		else dIndex=particle->D2;
  	}
  	
  	return abs(static_cast<GenParticle*>(branchParticle->At(dIndex))->PID);
}

//Function that receives the true particles branch, a particle, and a number 1 or 2 saying if it should return the index of the daughter 1 or the daughter 2. If it has no daughter for the number requested, it reutrns -999.
int findDaughterIndex(TClonesArray *branchParticle, GenParticle *particle, int d) 
{	
  	int dIndex;
  	
  	if(d==1) 
  	{
  		if(particle->D1 == -1) return -999;
  		else dIndex=particle->D1;
  	}
  	if(d==2) 
  	{
  		if(particle->D2 == -1) return -999;
  		else dIndex=particle->D2;
  	}
  	
  	return dIndex;
}

//Function that receives the true particles branch, a particle, and a number 1 or 2 saying if it should return the PID of the mother 1 or the mother 2. If it has no mother for the number requested, it reutrns -999.
int findMother(TClonesArray *branchParticle, GenParticle *particle, int m) 
{	
  	int mIndex;
  	
  	if(m==1) 
  	{
  		if(particle->M1 == -1) return -999;
  		else mIndex=particle->M1;
  	}
  	if(m==2) 
  	{
  		if(particle->M2 == -1) return -999;
  		else mIndex=particle->M2;
  	}
  	
  	return abs(static_cast<GenParticle*>(branchParticle->At(mIndex))->PID);
}

//Function that receives the true particles branch, a particle, and a number 1 or 2 saying if it should return the index of the mother 1 or the mother 2. If it has no mother for the number requested, it reutrns -999.
int findMotherIndex(TClonesArray *branchParticle, GenParticle *particle, int m) 
{	
  	int mIndex;
  	
  	if(m==1) 
  	{
  		if(particle->M1 == -1) return -999;
  		else mIndex=particle->M1;
  	}
  	else if(m==2) 
  	{
  		if(particle->M2 == -1) return -999;
  		else mIndex=particle->M2;
  	}
  	
  	return mIndex;
}

///Function that checks if a particle is coming from a photon.
void isGammaDaughter(TClonesArray *branchParticle, GenParticle *particle, bool& isGammaDaughterN)
{
	int pid = abs(particle->PID);
  	int m1 = abs(findMother(branchParticle, particle, 1));
  	
  	int m1Index = findMotherIndex(branchParticle, particle, 1);
  	//cout<<endl<<"isGammaDaughter func: "<<endl<<"pid: "<<pid<<"    m1: "<<m1;
  	if(pid != m1 && m1 != 22)
  	{
  		isGammaDaughterN = false;
  		return;
  	}
  	else if(m1 == 22) return;
  	else if(pid == m1)
  	{
  		if(m1Index != -999)
	  	{
	  		GenParticle *m1 = (GenParticle*) branchParticle->At(m1Index);
	  		//cout<<"Entra m1"<<endl;
	  		isGammaDaughter(branchParticle, m1, isGammaDaughterN);
	  	}
  	}
  	
  	return;
}

///Function that fills out a vector with the indices in branchParticle of the immediate decays of the gamma gamma decay. 
void findInitialMothersIndeces(TClonesArray *branchParticle, vector<int>& initialMothersIndices, string topology, string mothers)
{
	set<int> gammaGammaDecaysPIDs;
	vector<int> particlePIDs;
	if(topology == "HH")
	{ 
		gammaGammaDecaysPIDs.insert(25);
		gammaGammaDecaysPIDs.insert(36);	
	}
	else if(topology == "qq")
	{ 
		gammaGammaDecaysPIDs.insert(4);
		gammaGammaDecaysPIDs.insert(5);	
		gammaGammaDecaysPIDs.insert(6);	
		//////ADD ALL QUARKS
	}
	else if(topology == "ttbar") gammaGammaDecaysPIDs.insert(6);
	else if(topology == "ZZ") gammaGammaDecaysPIDs.insert(23);
	else if(topology == "WW") gammaGammaDecaysPIDs.insert(24);
	
	int pid=11, m1Pid=-999, particleIndex=0;
	////to find from which b the constituents are coming from
	if(mothers == "b")
	{
		while(pid == 11 || m1Pid == 11 || gammaGammaDecaysPIDs.find(m1Pid) != gammaGammaDecaysPIDs.end())
		{
			GenParticle *particle = (GenParticle*) branchParticle->At(particleIndex); 
			pid = abs(particle->PID);
			particlePIDs.push_back(pid);
			m1Pid = findMother(branchParticle, particle, 1);
			if(gammaGammaDecaysPIDs.find(m1Pid) != gammaGammaDecaysPIDs.end()) initialMothersIndices.push_back(particleIndex);
			particleIndex++;
		}
	}
	
	//////to find from which gammaGamma decay the consts. are coming from
	if(mothers == "gammaGammaDecay")
	{
		while(pid == 11 || m1Pid == 11)
		{
			GenParticle *particle = (GenParticle*) branchParticle->At(particleIndex); 
			pid = abs(particle->PID);
			particlePIDs.push_back(pid);
			m1Pid = findMother(branchParticle, particle, 1);
			if(m1Pid == 11 && pid != 11) initialMothersIndices.push_back(particleIndex);
			if(initialMothersIndices.size() == 1) initialMothersIndices.push_back(0); ////to make second gamma gamma decay be in pos=2 and align with color
			particleIndex++;
		}
	}
	
	cout<<"initialMothersIndices: "; 
	for(int i=0; i<initialMothersIndices.size(); i++) cout<<initialMothersIndices[i]<<" ";
	cout<<endl;
	cout<<"particlePIDs: "; 
	for(int i=0; i<10; i++) cout<<particlePIDs[i]<<" ";
	cout<<endl;
	
	return;
}

///Function that tells you the index in the vector of mothers of the initial mother of a particle. 
void findInitialMotherIndex(TClonesArray *branchParticle, GenParticle *particle, vector<int>& initialMothersIndices, int& initialMotherIndex)
{
	set<int> initialMothersIndicesSet;
	for(int i=0; i<initialMothersIndices.size(); i++) initialMothersIndicesSet.insert(initialMothersIndices[i]);
	int pid=-999, m1Index=-999, particleIndex=-999, cont=0;
	m1Index = findMotherIndex(branchParticle, particle, 1);
	while(initialMothersIndicesSet.find(m1Index) == initialMothersIndicesSet.end())
	{
		//if(cont<20) cout<<"m1Index: "<<m1Index<<endl;
		particleIndex= m1Index;
		GenParticle *motherParticle = (GenParticle*) branchParticle->At(particleIndex); 
		m1Index = findMotherIndex(branchParticle, motherParticle, 1);
		cont++;
	}
	//cout<<"m1Index: "<<m1Index<<endl;
	int index=0;
	while(m1Index != initialMothersIndices[index]) index++;
	initialMotherIndex=index;
	
	return;
}

////Function that receives two sets of angles and returns deltaR.
double findDeltaR(double eta1, double phi1, double eta2, double phi2)
{
	return sqrt(pow(eta1-eta2, 2) + pow(phi1-phi2, 2));
}

/////Function that receives a jet, the particles branch and the vector with the indices of the bs in the branchParticle and returns the position in said vector of the matched b using deltaR. 
int bMatching(Jet* jet, TClonesArray *branchParticle, vector<int>& initialMothersIndices)
{
	double minDeltaR=999999;
	double bIndex;
	for(int i=0; i<initialMothersIndices.size(); i++)
	{
		GenParticle *particle = (GenParticle*) branchParticle->At(initialMothersIndices[i]);	
		double deltaR = findDeltaR(jet->Eta, jet->Phi, particle->Eta, particle->Phi);
		if(deltaR < minDeltaR)
		{
			minDeltaR = deltaR;
			bIndex = i;
		}	
	}
	
	return bIndex;
}

////Function that prints out the decay tree of the particles
void printDecay(TClonesArray *branchParticle)
{
	int nParticles = branchParticle->GetEntries();
	//nParticles = 40;
	set<int> particlesInLevel;
	
	for(int i=0; i<nParticles; i++)
	{
		GenParticle *particle = (GenParticle*) branchParticle->At(i); 
		int pid = abs(particle->PID);
		int m1Index = findMotherIndex(branchParticle, particle, 1);
		int m2Index = findMotherIndex(branchParticle, particle, 2);
		int m1Pid = findMother(branchParticle, particle, 1);
		int m2Pid = findMother(branchParticle, particle, 2);
		int d1Index = findDaughterIndex(branchParticle, particle, 1);
		int d2Index = findDaughterIndex(branchParticle, particle, 2);
		int d1Pid = findDaughter(branchParticle, particle, 1);
		int d2Pid = findDaughter(branchParticle, particle, 2);
		
		if(particlesInLevel.find(m1Index) != particlesInLevel.end())
		{
			particlesInLevel.clear();
			cout<<endl<<endl;
		}
		else if(!(particlesInLevel.empty())) cout<<"\033[1m -/- \033[0m";
		particlesInLevel.insert(i);
		
		cout<<"\033[1mpid: \033[0m"<<pid<<"("<<i<<"); m1: "<<m1Pid<<"("<<m1Index<<"); m2: "<<m2Pid<<"("<<m2Index<<"); d1: "<<d1Pid<<"("<<d1Index<<"); d2: "<<d2Pid<<"("<<d2Index<<")";
	}
	cout<<endl<<endl;
}

////Function that receives a particle and returns its index in the branchParticle.
int findParticleIndex(TClonesArray* branchParticle, GenParticle* particleConst) {
    for (int i = 0; i < branchParticle->GetEntries(); ++i) {
        if (branchParticle->At(i) == particleConst) {
            return i;
        }
    }
    return -999; // Return -1 if the particle is not found
}

/////Function that plots genJet constituents
void plotConstituents(TClonesArray *branchGenJet, int nGenJets, string plane, TCanvas *cplotConstituents1, string color, TClonesArray *branchParticle, string topology, string mothers)
{
	cplotConstituents1->cd();
	TObject *object;
  	GenParticle *particleConst;
  	Track *trackConst;
  	Tower *towerConst;
  	
  	double phi=-999, theta=-999;
  	
  	vector<int> initialMothersIndices;
  	findInitialMothersIndeces(branchParticle, initialMothersIndices, topology, mothers);
  	
  	int cont=0, contB1=0, contB2=0, contB3=0, contB4=0, contBExtra=0;
  	if(plane == "phiTheta")
  	{
  		for(int i=0;i<nGenJets;i++)
		{
			Jet *genjet = (Jet*) branchGenJet->At(i); // Take ith jet
			vector<int> jetConstituentsIndeces;
			for(int j = 0; j < genjet->Constituents.GetEntriesFast(); ++j)
		        {
				object = genjet->Constituents.At(j);
				if(object == 0) continue;
				if(object->IsA() == GenParticle::Class())
				{
					particleConst = (GenParticle*) object;
					int constIndex = findParticleIndex(branchParticle, particleConst);
					jetConstituentsIndeces.push_back(constIndex);
				  	phi = particleConst->Phi;
				  	theta = particleConst->Eta;
				  	theta = findTheta(theta);
				  	//TLine* line = new TLine(0.0,0.0,phi,theta);
				  	TEllipse* particleDrawn = new TEllipse(phi,theta,.025,.025);
					particleDrawn->SetLineStyle(1);
					int particleColor;
					if(color == "initialParticle")
					{
						int initialMotherIndex=-999;
						findInitialMotherIndex(branchParticle, particleConst, initialMothersIndices, initialMotherIndex);
						if(initialMotherIndex == 0)
						{
							particleColor = kViolet;
							contB1++;
						}
						if(initialMotherIndex == 1)
						{
							particleColor = kMagenta;
							contB2++;
						}
						if(initialMotherIndex == 2)
						{
							particleColor = kGreen+2;
							contB3++;
						}
						if(initialMotherIndex == 3) 
						{
							particleColor = kGreen;
							contB4++;
						}
						if(initialMotherIndex > 3) 
						{
							particleColor = kBlack;
							contBExtra++;
						}
						//cout<<endl;
						findInitialMotherIndex(branchParticle, particleConst, initialMothersIndices, initialMotherIndex);
						//cout<<"initialMotherIndex: "<<initialMotherIndex<<endl;
						cont++;
					}
					else if(color == "jet")
					{
						int bMatchingIndex = bMatching(genjet, branchParticle, initialMothersIndices);
						if(bMatchingIndex == 0) particleColor = kViolet;
						if(bMatchingIndex == 1) particleColor = kMagenta;
						if(bMatchingIndex == 2) particleColor = kGreen+2;
						if(bMatchingIndex == 3) particleColor = kGreen;
						if(bMatchingIndex > 4) particleColor = kBlack;
					}
					else throw std::runtime_error("ERROR: unknown color scheme!");
					particleDrawn->SetFillColor(particleColor);  
					particleDrawn->SetLineColor(particleColor);
			   	        particleDrawn->SetLineWidth(1);
					particleDrawn->Draw(); 
	        			//line->SetLineWidth(2);
      	        			//line->Draw("same");
      	        			
				}
			}
			cout<<"Indeces of consts. in jet"<<i<<": ";
			for(int j=0; j<jetConstituentsIndeces.size(); j++) cout<<jetConstituentsIndeces[j]<<", ";
			cout<<endl;
			cout<<"jetConstituentsIndeces.size(): "<<jetConstituentsIndeces.size()<<endl<<endl;
		}
  	}
  	else if(plane == "XY")
  	{
  		
  	}
  	cout<<"jet constituents from B1: "<<contB1<<endl;
  	cout<<"jet constituents from B2: "<<contB2<<endl;
  	cout<<"jet constituents from B3: "<<contB3<<endl;
  	cout<<"jet constituents from B4: "<<contB4<<endl; 
  	cout<<"jet constituents from extra: "<<contBExtra<<endl;
}

void eventDisplayXY(const char *inputFile, int nthEvent, string jetAlgo, string genJetAlgo, string level)
{
	gSystem->Load("libDelphes");
  	TChain chain("Delphes");
  	chain.Add(inputFile);
  	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  	Long64_t numberOfEntries = treeReader->GetEntries();
  	
	TClonesArray *branchJet = treeReader->UseBranch(jetAlgo.c_str());
	TClonesArray *branchGenJet = treeReader->UseBranch(genJetAlgo.c_str());
	TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchEvent = treeReader->UseBranch("Event");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	//TClonesArray *branchNeutralHadron = treeReader->UseBranch("NeutralHadron");
	TClonesArray *branchEFlowTower = treeReader->UseBranch("Tower");
	
	treeReader->ReadEntry(nthEvent);
	
	int nJets=0;
	if(branchJet->GetEntries() > 0)  nJets =  branchJet->GetEntries();
	int nGenJets=0;
	if(branchGenJet->GetEntries() > 0)  nGenJets =  branchGenJet->GetEntries();
	int nTracks=0;
	if(branchEFlowTrack->GetEntries() > 0)  nTracks =  branchEFlowTrack->GetEntries();
	
	double ptMax=0;
	findDetectorRadiusXY(branchEFlowTrack, branchGenJet, nTracks, nGenJets, ptMax, level);
	//cout<<"ptMax: "<<ptMax<<endl<<"nJets: "<<nJets<<endl<<"nGenJets: "<<nGenJets<<endl<<endl; 
	
  	
  	TH2F *histJetPxAndPy = new TH2F("jet_px_vs_py", "jet py vs. px", 100, -ptMax*1.1, ptMax*1.1, 100, -ptMax*1.1, ptMax*1.1);
  	TCanvas *ceventDisplayXY1 = new TCanvas();
  	ceventDisplayXY1->cd();
  	histJetPxAndPy->Draw("col");
  	histJetPxAndPy->SetStats(0);
	histJetPxAndPy->GetXaxis()->SetTitle("y");
	histJetPxAndPy->GetYaxis()->SetTitle("x");
	histJetPxAndPy->SetFillColor(0);
	
	TEllipse* detector = new TEllipse(0.0,0.0,ptMax,ptMax);
        	detector->SetLineStyle(1);
                detector->SetFillStyle(4000);  
                detector->SetLineColor(kBlack);
   	        detector->SetLineWidth(1);
	        detector->Draw();
	
  	TLine* xAxis = new TLine(-ptMax,0.0,ptMax,0.0);
	        xAxis->SetLineColor(1);
	        xAxis->SetLineWidth(1);
	        xAxis->Draw("same");
	      
  	TLine* yAxis = new TLine(0.0,-ptMax,0.0,ptMax);
	        yAxis->SetLineColor(1);
	        yAxis->SetLineWidth(1);
	        yAxis->Draw("same");
	        
}

void eventDisplayPhiTheta(const char *inputFile, int nthEvent, string jetAlgo, string genJetAlgo, string level, string color, string topology)
{
	string mothers;
	if(color == "initialParticle") mothers = "gammaGammaDecay";
	else if(color == "jet") mothers = "b";
	 
	gSystem->Load("libDelphes");
  	TChain chain("Delphes");
  	chain.Add(inputFile);
  	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  	Long64_t numberOfEntries = treeReader->GetEntries();
  	
	TClonesArray *branchJet = treeReader->UseBranch(jetAlgo.c_str());
	TClonesArray *branchGenJet = treeReader->UseBranch(genJetAlgo.c_str());
	TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchEvent = treeReader->UseBranch("Event");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	//TClonesArray *branchNeutralHadron = treeReader->UseBranch("NeutralHadron");
	TClonesArray *branchEFlowTower = treeReader->UseBranch("Tower");
	
	treeReader->ReadEntry(nthEvent);
	
	int nJets=0;
	if(branchJet->GetEntries() > 0)  nJets =  branchJet->GetEntries();
	int nGenJets=0;
	if(branchGenJet->GetEntries() > 0)  nGenJets =  branchGenJet->GetEntries();
	int nTracks=0;
	if(branchEFlowTrack->GetEntries() > 0)  nTracks =  branchEFlowTrack->GetEntries();
	
	double phiMax=0, thetaMax=0;
	findDetectorRadiusPhiTheta(branchEFlowTrack, branchGenJet, nTracks, nGenJets, phiMax, thetaMax, level);
	//cout<<"ptMax: "<<ptMax<<endl<<"nJets: "<<nJets<<endl<<"nGenJets: "<<nGenJets<<endl<<endl; 
	
  	
  	TH2F *histJetPhiAndTheta = new TH2F("jet_phi_vs_Theta", "jet Phi vs. Theta", 100, -phiMax*1.1, phiMax*1.1, 100, -.5, thetaMax*1.1);
  	TCanvas *ceventDisplayPhiTheta1 = new TCanvas();
  	ceventDisplayPhiTheta1->cd();
  	histJetPhiAndTheta->Draw("col");
  	if(mothers == "gammaGammaDecay")
  	{
	  	TLegend *legendJetPhiAndThetaZ1 = new TLegend(0.15, 0.09, 0.3, 0.24);
		legendJetPhiAndThetaZ1->AddEntry(histJetPhiAndTheta, "gen. jet constituents from Z1", "");
		legendJetPhiAndThetaZ1->SetTextColor(kViolet);
		legendJetPhiAndThetaZ1->SetBorderSize(0);
		legendJetPhiAndThetaZ1->SetFillStyle(0);
		legendJetPhiAndThetaZ1->SetTextSize(0.033);
	  	legendJetPhiAndThetaZ1->Draw("same");
	  	
	  	TLegend *legendJetPhiAndThetaZ2 = new TLegend(0.50, 0.09, 0.65, 0.24);
		legendJetPhiAndThetaZ2->AddEntry(histJetPhiAndTheta, "gen. jet constituents from Z2", "");
		legendJetPhiAndThetaZ2->SetTextColor(kGreen+2);
		legendJetPhiAndThetaZ2->SetBorderSize(0);
		legendJetPhiAndThetaZ2->SetFillStyle(0);
		legendJetPhiAndThetaZ2->SetTextSize(0.033);
	  	legendJetPhiAndThetaZ2->Draw("same");
  	}
  	else if(mothers == "b")
  	{
	  	TLegend *legendJetPhiAndThetab1 = new TLegend(0.12, 0.11, 0.27, 0.26);
		legendJetPhiAndThetab1->AddEntry(histJetPhiAndTheta, "gen. jet constituents from jetb1 (from Z1)", "");
		legendJetPhiAndThetab1->SetTextColor(kViolet);
		legendJetPhiAndThetab1->SetBorderSize(0);
		legendJetPhiAndThetab1->SetFillStyle(0);
		legendJetPhiAndThetab1->SetTextSize(0.028);
	  	legendJetPhiAndThetab1->Draw("same");
	  	
	  	TLegend *legendJetPhiAndThetab2 = new TLegend(0.12, 0.07, 0.27, 0.22);
		legendJetPhiAndThetab2->AddEntry(histJetPhiAndTheta, "gen. jet constituents from jetb2 (from Z1)", "");
		legendJetPhiAndThetab2->SetTextColor(kMagenta);
		legendJetPhiAndThetab2->SetBorderSize(0);
		legendJetPhiAndThetab2->SetFillStyle(0);
		legendJetPhiAndThetab2->SetTextSize(0.028);
	  	legendJetPhiAndThetab2->Draw("same");
	  	
	  	TLegend *legendJetPhiAndThetab3 = new TLegend(0.5, 0.11, 0.65, 0.26);
		legendJetPhiAndThetab3->AddEntry(histJetPhiAndTheta, "gen. jet constituents from jetb3 (from Z2)", "");
		legendJetPhiAndThetab3->SetTextColor(kGreen+2);
		legendJetPhiAndThetab3->SetBorderSize(0);
		legendJetPhiAndThetab3->SetFillStyle(0);
		legendJetPhiAndThetab3->SetTextSize(0.028);
	  	legendJetPhiAndThetab3->Draw("same");
	  	
	  	TLegend *legendJetPhiAndThetab4 = new TLegend(0.5, 0.07, 0.65, 0.22);
		legendJetPhiAndThetab4->AddEntry(histJetPhiAndTheta, "gen. jet constituents from jetb4 (from Z2)", "");
		legendJetPhiAndThetab4->SetTextColor(kGreen);
		legendJetPhiAndThetab4->SetBorderSize(0);
		legendJetPhiAndThetab4->SetFillStyle(0);
		legendJetPhiAndThetab4->SetTextSize(0.028);
	  	legendJetPhiAndThetab4->Draw("same");
  	}
  	
  	histJetPhiAndTheta->SetStats(0);
	histJetPhiAndTheta->GetXaxis()->SetTitle("phi");
	histJetPhiAndTheta->GetYaxis()->SetTitle("theta");
	histJetPhiAndTheta->SetFillColor(0);
	
	TEllipse* detector = new TEllipse(0.0,0.0,phiMax,thetaMax);
        	detector->SetLineStyle(1);
                detector->SetFillStyle(4000);  
                detector->SetLineColor(kBlack);
   	        detector->SetLineWidth(1);
	        //detector->Draw();
	
  	TLine* xAxis = new TLine(-phiMax,0.0,phiMax,0.0);
	        xAxis->SetLineColor(1);
	        xAxis->SetLineWidth(1);
	        xAxis->Draw("same");
	      
  	TLine* yAxis = new TLine(0.0,-thetaMax,0.0,thetaMax);
	        yAxis->SetLineColor(1);
	        yAxis->SetLineWidth(1);
	        yAxis->Draw("same");
	        
	plotConstituents(branchGenJet, nGenJets, "phiTheta", ceventDisplayPhiTheta1, color, branchParticle, topology, mothers);
	/*TCanvas *ceventDisplayPhiTheta2 = new TCanvas();
  	ceventDisplayPhiTheta2->cd();
  	histJetPhiAndTheta->Draw("col");
  	histJetPhiAndTheta->SetStats(0);
	histJetPhiAndTheta->GetXaxis()->SetTitle("phi");
	histJetPhiAndTheta->GetYaxis()->SetTitle("theta");
	histJetPhiAndTheta->SetFillColor(0);
	plotJetParticles(branchGenJet, nGenJets, "phiTheta", ceventDisplayPhiTheta2, color, branchParticle, topology);*/
	
	string outputCanvasName = "eventDisplayPhiTheta"+mothers+to_string(nthEvent)+".png";
	ceventDisplayPhiTheta1->Update();
	ceventDisplayPhiTheta1->SaveAs(outputCanvasName.c_str());
	 
	printDecay(branchParticle);
}

void EventDisplay(const char *inputFile, string level, string topology, int nthEvent)
{
	string jetAlgoText = "(durham rtd_cut=10) ";
	string jetAlgo = "Jet10";
  	string genJetAlgo = "GenJet10";
  	string rtdCut = "10";
  	
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	
	//eventDisplayXY(inputFile, nthEvent, jetAlgo, genJetAlgo, level);
	//eventDisplayPhiTheta(inputFile, nthEvent, jetAlgo, genJetAlgo, level, "initialParticle", topology);
	eventDisplayPhiTheta(inputFile, nthEvent, jetAlgo, genJetAlgo, level, "jet", topology);
}








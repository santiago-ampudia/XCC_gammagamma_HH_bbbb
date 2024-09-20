
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
#endif


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

//Function that receives the true particles branch and prints the decay chain of the particles in the branch.
void printDecayHelper(TClonesArray *branchParticle)
{
	int nParticles = branchParticle->GetEntries();
	nParticles = 20;
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

void printDecay(const char *inputFile, int nthEvent)
{
	gSystem->Load("libDelphes");
  	TChain chain("Delphes");
  	chain.Add(inputFile);
  	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  	Long64_t numberOfEntries = treeReader->GetEntries();
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchJet = treeReader->UseBranch("Jet10");
	treeReader->ReadEntry(nthEvent);

	int nJets=0, nParticles=0;
	if(branchJet->GetEntries() > 0)  nJets =  branchJet->GetEntries();
	if(branchParticle->GetEntries() > 0)  nParticles =  branchParticle->GetEntries();
	cout<<"nParticles: "<<nParticles<<endl;
	cout<<"nJets: "<<nJets<<endl<<endl<<endl;

	
	printDecayHelper(branchParticle);
}

// Function to calculate the magnitude of a vector
float Magnitude(TLorentzVector vec) 
{
    return vec.P();
}
//////

//////Function that returns the maxThrust of the events
float findThrust(vector<TLorentzVector>& momenta, TLorentzVector& thrustAxis)
{
	float maxThrust = 0.0;
	
	for(float theta = 0.0; theta < TMath::Pi(); theta += 0.1) 
        {
        	for(float phi = 0.0; phi < 2.0 * TMath::Pi(); phi += 0.1) 
       		{
                	TLorentzVector axis(TMath::Sin(theta) * TMath::Cos(phi), 
                                    TMath::Sin(theta) * TMath::Sin(phi), 
                                    TMath::Cos(theta), 0.0);
                
                	float sumP = 0.0;
                	float sumP_parallel = 0.0;
		        for (auto& p : momenta) 
		        {
		            sumP += p.P();
		            sumP_parallel += TMath::Abs(p.Dot(axis) / Magnitude(axis));
		        }

		        float thrust = sumP_parallel / sumP;
		        if (thrust > maxThrust) 
		        {
		            maxThrust = thrust;
		            thrustAxis = axis;
		        }
            	}
	}
	return maxThrust;
}
//////////

void printDecayXEvents(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet10");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  vector<int> xEvents;
  

  string condition = "jetCosThetaThrust";
  double threshold = 0.995;
  int nJets=0;
  
  cout<<"Number of entries: "<<numberOfEntries<<endl<<endl;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);
	//cout<<"Event: "<<entry<<endl;

    if(condition == "jetCosThetaThrust")
    {
		if(branchJet->GetEntries() > 0)  nJets =  branchJet->GetEntries();
    	else nJets=0;
	  	vector<TLorentzVector> momenta;
    	for (int i = 0; i < branchEFlowTrack->GetEntries(); ++i) 
		{
			Track *track = (Track*) branchEFlowTrack->At(i);
			TLorentzVector track4V;
			track4V.SetPtEtaPhiM(track->PT, track->Eta, track->Phi, 0.0);
			momenta.push_back(track4V);
		}
		for (int i = 0; i < branchEFlowPhoton->GetEntries(); ++i) 
		{
			Tower *photon = (Tower*) branchEFlowPhoton->At(i);
			TLorentzVector photon4V;
			photon4V.SetPtEtaPhiM(photon->ET, photon->Eta, photon->Phi, 0.0);
			momenta.push_back(photon4V);
		}
		for (int i = 0; i < branchEFlowNeutralHadron->GetEntries(); ++i) 
		{
			Tower *neutralHadron = (Tower*) branchEFlowNeutralHadron->At(i);
			TLorentzVector neutralHadron4V;
			neutralHadron4V.SetPtEtaPhiM(neutralHadron->ET, neutralHadron->Eta, neutralHadron->Phi, 0.0);
			momenta.push_back(neutralHadron4V);
		}
		TLorentzVector thrustAxis;
		float thrust = findThrust(momenta, thrustAxis);
		double jetCosThetaThrust;
		for(int i=0;i<nJets;i++)
     	{
        	Jet *jet = (Jet*) branchJet->At(i);
			TLorentzVector jetP4 = jet->P4();
        	jetCosThetaThrust = thrustAxis.Vect().Dot(jetP4.Vect()) / (thrustAxis.P() * jetP4.P());
			//cout<<"i: "<<i<<"; jetCosThetaThrust: "<<jetCosThetaThrust<<endl;
			if(jetCosThetaThrust > threshold)
			{
				xEvents.push_back(entry);
				//cout<<"IN"<<endl;
				break;
			}
      	}
		//cout<<"OUT"<<endl;
    }
  }

  cout<<endl<<endl<<"Number of x events: "<<xEvents.size()<<endl<<endl<<endl;
  for(int i=0;i<xEvents.size();i++) 
	{
		cout<<"Event "<<i<<": "<<xEvents[i]<<endl;
		printDecay(inputFile, xEvents[i]);
		cout<<endl<<endl<<"----------------------------------------"<<endl<<endl;
	}

}


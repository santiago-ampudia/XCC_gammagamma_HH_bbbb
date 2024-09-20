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








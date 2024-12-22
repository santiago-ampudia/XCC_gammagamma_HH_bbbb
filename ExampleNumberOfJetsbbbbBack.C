/*

root -l examples/ExampleNumberOfJetsHH.C'("delphes_output.root")'
*/

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



//------------------------------------------------------------------------------

void ExampleNumberOfJetsbbbbBack(const char *inputFile)
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
  TClonesArray *branchJetAntiKt = treeReader->UseBranch("JetAntiKt");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  
  TH1 *histBJets = new TH1F("bjets", "Number of b-tagged jets in events where HH to bb and bb", 72.0, 0.0, 6.0);
  TH1 *histNonBJets = new TH1F("nonbjets", "Number of non-b-tagged jets in events where HH to bb and bb", 72.0, 0.0, 6.0);
  TH1 *histNJets = new TH1F("njets", "(antiKt R=0.5) Number of Jets in ZH (test) Events", 72.0, -1.0, 9.0);
  TH2F *histNJets2D = new TH2F("histnjetsd2", "(rtd_min=10) b-tagged vs. non-b-taged jets (bbbb) (ZH test)", 100, -1.0, 7.0, 100, -1.0, 7.0);
  TH2F *histNJets2DAntiKt = new TH2F("histnjetsd2antikt", "(antiKt R=0.5) b-tagged vs. non-b-taged jets (bbbb)", 100, -1.0, 9.0, 100, -1.0, 9.0);
  TH2F *histNJetsCompareAlgos = new TH2F("histnjjetscomparealgos", "(rtd_min=10, R=0.5) nJets in durham vs. antiKt (bbbb)", 100, -1.0, 7.0, 100, -1.0, 9.0);
  histNJetsCompareAlgos->GetXaxis()->SetTitle("Number of reco. jets with durham");
  histNJetsCompareAlgos->GetYaxis()->SetTitle("Number of reco. jets with antiKt");
  histBJets->SetLineColor(857);
  histNonBJets->SetLineColor(413);
  histNJets2D->GetXaxis()->SetTitle("b-tagged jets");
  histNJets2D->GetYaxis()->SetTitle("non-b-tagged jets");
  histNJets2DAntiKt->GetXaxis()->SetTitle("b-tagged jets");
  histNJets2DAntiKt->GetYaxis()->SetTitle("non-b-tagged jets");

  
  int njets=0, njetsAntiKt=0, nJets=0, nJetsAntiKt=0;
  int nparticles=0;
  double contEventsHH=0;
  double cont2B2NB=0;
  int contBJets=0, contBJetsAntiKt=0;
  int contNonBJets=0, contNonBJetsAntiKt=0;
  int contNon4Jets=0;
  int contEventsBBE=0;
  cout<<"Number of entries: "<<numberOfEntries<<endl<<endl;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    if(branchJet->GetEntries() > 0)  nJets =  branchJet->GetEntries();
    else nJets=0;
    if(nJets!=0 && nJets!=4) contNon4Jets++;
    
    if(branchJetAntiKt->GetEntries() > 0)  nJetsAntiKt =  branchJetAntiKt->GetEntries();
    else nJetsAntiKt=0;
    
    histNJetsCompareAlgos->Fill(nJets, nJetsAntiKt);
     	
    int contG=0;
    int contB=0;
   
    /*if(branchParticle->GetEntries() > 0) nparticles =  branchParticle->GetEntries();
    else nparticles=0;
   
    for(int i=0;i<nparticles;i++)
    {
      GenParticle *particle = (GenParticle*) branchParticle->At(i); 
      int d1Index=particle->D1;
      int d2Index=particle->D2;
      if (d1Index >= 0 && d2Index >= 0)
      {
      int d1=abs(static_cast<GenParticle*>(branchParticle->At(d1Index))->PID);
      int d2=abs(static_cast<GenParticle*>(branchParticle->At(d2Index))->PID);
      if ((particle->PID == 25 || particle->PID == 36) && (d1==5 && d2==5)) contB++; 
      }
    }*/
    
     contB=2;
     if(contB>=2) 
     {
        contEventsHH++;
        contBJets=0;
        contNonBJets=0;
        contBJetsAntiKt=0;
        contNonBJetsAntiKt=0;
        if(branchJet->GetEntries() > 0)  njets =  branchJet->GetEntries();
        else njets=0;
        histNJets->Fill(nJetsAntiKt);
     	  for(int i=0;i<njets;i++)
     	  {
        	Jet *jet = (Jet*) branchJet->At(i);
        	if(jet->BTag==1) contBJets++;
        	else contNonBJets++;
      	}
      	histNJets2D->Fill(contBJets,contNonBJets);
      	
      	if(branchJetAntiKt->GetEntries() > 0)  njetsAntiKt =  branchJetAntiKt->GetEntries();
        else njetsAntiKt=0;
        for(int i=0;i<njetsAntiKt;i++)
        {
          Jet *jet = (Jet*) branchJetAntiKt->At(i);
          if(jet->BTag==1) contBJetsAntiKt++;
          else contNonBJetsAntiKt++;
        }
      	
      	
        histNJets2DAntiKt->Fill(contBJetsAntiKt,contNonBJetsAntiKt);

        
      }
     }
     
    TCanvas *c1 = new TCanvas();
    /*TCanvas *c2 = new TCanvas();
    TCanvas *c3 = new TCanvas();*/
    TCanvas *c4 = new TCanvas();
    TCanvas *c5 = new TCanvas();
    TCanvas *c6 = new TCanvas();
    
    c1->cd();
    histNJets->Draw();
    
    /*c2->cd();
    histBJets->Draw();
    
    c3->cd();
    histNonBJets->Draw();*/
    
    gStyle->SetPaintTextFormat();
    
    c4->cd();
    histNJets2D->Draw("TEXT");
    
    c5->cd();
    histNJets2DAntiKt->Draw("TEXT");
    
    c6->cd();
    histNJetsCompareAlgos->Draw("TEXT");
    
    
    double frac=cont2B2NB/contEventsHH;
    cout<<"Events with HH to bb and bb: "<<contEventsHH<<endl<<"Events with 2 b-jets and 2 non-b-jets: "<<cont2B2NB<<endl<<"%: "<<frac<<endl;
    cout<<endl<<"Number of entries: "<<numberOfEntries<<endl<<"HH->bbbb: "<<contEventsHH<<endl<<contEventsHH/numberOfEntries<<endl;
    
    cout<<"contNon4Jets: "<<contNon4Jets<<endl;
  }


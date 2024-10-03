#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <map>
#include <algorithm>
#include <random>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include <TCanvas.h>
#include <TH2F.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TKey.h>
#include <TMath.h>
#include <TStyle.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <stdexcept>
#include <fstream>

using namespace std;

void fillOutTrees(TTree& TreeTrain, TTree& TreeTest, string fileName, string treeName)
{
    //cout<<"starting function"<<endl;

    TFile* file = new TFile(fileName.c_str());
    TTree *fullTree;
    file->GetObject(treeName.c_str(), fullTree);
    if (!fullTree) {
        cout << "Error: Tree " << treeName << " not found in file " << fileName << endl;
        return;
    }
    //else cout<<"fullTree exists"<<endl;

    Float_t aplanarity, cosThetaB1, cosThetaB2, cosThetaB3, cosThetaB4, invMassB1, invMassB2, jetB1M, jetB2M, jetB3M, jetB4M, jetB1Pt, jetB2Pt, jetB3Pt, jetB4Pt, minJetM, sphericity, sumPt, nParticles, minConstSize, jetNObjects, minJetNObjects, invMassB1AntiKt, invMassB2AntiKt, nJetsAntiKt, invMassB11Best, invMassB21Best, invMassB12Best, invMassB22Best, invMassB13Best, invMassB23Best, invMassB14Best, invMassB24Best, invMassB15Best, invMassB25Best, invMassB16Best, invMassB26Best, invMassB17Best, invMassB27Best, invMassB18Best, invMassB28Best, exclYmerge12, exclYmerge23, exclYmerge34, exclYmerge45, exclYmerge56, invMassZZ1, invMassZZ2, thrust, boostB1, boostB2, boostB3, boostB4, boostSystem, missingET, constSizeB1, constSizeB2, constSizeB3, constSizeB4, nJetsDurham15, nJetsDurham20, nJetsDurham25, nJetsDurham30, distanceZ1MinChiSquaredZZMass, distanceZ2MinChiSquaredZZMass;
    Float_t entryIndex;
    TreeTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeTrain.Branch("aplanarity",&aplanarity,"aplanarity/F");
	TreeTrain.Branch("invMassB1",&invMassB1,"invMassB1/F");
	TreeTrain.Branch("invMassB2",&invMassB2,"invMassB2/F");
    TreeTrain.Branch("minJetM",&minJetM,"minJetM/F");  
    TreeTrain.Branch("sphericity",&sphericity,"sphericity/F");
    TreeTrain.Branch("cosThetaB1",&cosThetaB1,"cosThetaB1/F");
    TreeTrain.Branch("cosThetaB2",&cosThetaB2,"cosThetaB2/F");
    TreeTrain.Branch("cosThetaB3",&cosThetaB3,"cosThetaB3/F");
    TreeTrain.Branch("cosThetaB4",&cosThetaB4,"cosThetaB4/F");
    TreeTrain.Branch("sumPt",&sumPt,"sumPt/F");
    TreeTrain.Branch("jetB1Pt",&jetB1Pt,"jetB1Pt/F");
    TreeTrain.Branch("jetB2Pt",&jetB2Pt,"jetB2Pt/F");
    TreeTrain.Branch("jetB3Pt",&jetB3Pt,"jetB3Pt/F");
    TreeTrain.Branch("jetB4Pt",&jetB4Pt,"jetB4Pt/F");
    TreeTrain.Branch("jetB1M",&jetB1M,"jetB1M/F");
    TreeTrain.Branch("jetB2M",&jetB2M,"jetB2M/F");
    TreeTrain.Branch("jetB3M",&jetB3M,"jetB3M/F");
    TreeTrain.Branch("jetB4M",&jetB4M,"jetB4M/F");
    TreeTrain.Branch("nParticles",&nParticles,"nParticles/F");
    TreeTrain.Branch("constSizeB1",&constSizeB1,"constSizeB1/F");
    TreeTrain.Branch("constSizeB2",&constSizeB2,"constSizeB2/F");
    TreeTrain.Branch("constSizeB3",&constSizeB3,"constSizeB3/F");
    TreeTrain.Branch("constSizeB4",&constSizeB4,"constSizeB4/F");
    TreeTrain.Branch("minConstSize",&minConstSize,"minConstSize/F");
    TreeTrain.Branch("jetNObjects",&jetNObjects,"jetNObjects/F");
    TreeTrain.Branch("minJetNObjects",&minJetNObjects,"minJetNObjects/F");
    TreeTrain.Branch("invMassB1AntiKt",&invMassB1AntiKt,"invMassB1AntiKt/F");
    TreeTrain.Branch("invMassB2AntiKt",&invMassB2AntiKt,"invMassB2AntiKt/F");
    TreeTrain.Branch("nJetsAntiKt",&nJetsAntiKt,"nJetsAntiKt/F");
    TreeTrain.Branch("invMassB11Best",&invMassB11Best,"invMassB11Best/F");
    TreeTrain.Branch("invMassB21Best",&invMassB21Best,"invMassB21Best/F");
    TreeTrain.Branch("invMassB12Best",&invMassB12Best,"invMassB12Best/F");
    TreeTrain.Branch("invMassB22Best",&invMassB22Best,"invMassB22Best/F");
    TreeTrain.Branch("invMassB13Best",&invMassB13Best,"invMassB13Best/F");
    TreeTrain.Branch("invMassB23Best",&invMassB23Best,"invMassB23Best/F");
    TreeTrain.Branch("invMassB14Best",&invMassB14Best,"invMassB14Best/F");
    TreeTrain.Branch("invMassB24Best",&invMassB24Best,"invMassB24Best/F");
    TreeTrain.Branch("invMassB15Best",&invMassB15Best,"invMassB15Best/F");
    TreeTrain.Branch("invMassB25Best",&invMassB25Best,"invMassB25Best/F");
    TreeTrain.Branch("invMassB16Best",&invMassB16Best,"invMassB16Best/F");
    TreeTrain.Branch("invMassB26Best",&invMassB26Best,"invMassB26Best/F");
    TreeTrain.Branch("invMassB17Best",&invMassB17Best,"invMassB17Best/F");
    TreeTrain.Branch("invMassB27Best",&invMassB27Best,"invMassB27Best/F");
    TreeTrain.Branch("invMassB18Best",&invMassB18Best,"invMassB18Best/F");
    TreeTrain.Branch("invMassB28Best",&invMassB28Best,"invMassB28Best/F");
    //TreeTrain.Branch("nJetsDurham0",&nJetsDurham0,"nJetsDurham0/F");
    //TreeTrain.Branch("nJetsDurham5",&nJetsDurham5,"nJetsDurham5/F");
    //TreeTrain.Branch("nJetsDurham10",&nJetsDurham10,"nJetsDurham10/F");
    TreeTrain.Branch("nJetsDurham15",&nJetsDurham15,"nJetsDurham15/F");
    TreeTrain.Branch("nJetsDurham20",&nJetsDurham20,"nJetsDurham20/F");
    TreeTrain.Branch("nJetsDurham25",&nJetsDurham25,"nJetsDurham25/F");
    TreeTrain.Branch("nJetsDurham30",&nJetsDurham30,"nJetsDurham30/F");
    TreeTrain.Branch("distanceZ1MinChiSquaredZZMass",&distanceZ1MinChiSquaredZZMass,"distanceZ1MinChiSquaredZZMass/F");
    TreeTrain.Branch("distanceZ2MinChiSquaredZZMass",&distanceZ2MinChiSquaredZZMass,"distanceZ2MinChiSquaredZZMass/F");
    TreeTrain.Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
    TreeTrain.Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
    TreeTrain.Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
    TreeTrain.Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
    TreeTrain.Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");
    TreeTrain.Branch("invMassZZ1",&invMassZZ1,"invMassZZ1/F");
    TreeTrain.Branch("invMassZZ2",&invMassZZ2,"invMassZZ2/F");
    TreeTrain.Branch("thrust",&thrust,"thrust/F");
    TreeTrain.Branch("boostB1",&boostB1,"boostB1/F");
    TreeTrain.Branch("boostB2",&boostB2,"boostB2/F");
    TreeTrain.Branch("boostB3",&boostB3,"boostB3/F");
    TreeTrain.Branch("boostB4",&boostB4,"boostB4/F");
    TreeTrain.Branch("boostSystem",&boostSystem,"boostSystem/F");
    TreeTrain.Branch("missingET",&missingET,"missingET/F");

    TreeTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
    TreeTest.Branch("aplanarity",&aplanarity,"aplanarity/F");
    TreeTest.Branch("invMassB1",&invMassB1,"invMassB1/F");
    TreeTest.Branch("invMassB2",&invMassB2,"invMassB2/F");
    TreeTest.Branch("minJetM",&minJetM,"minJetM/F");  
    TreeTest.Branch("sphericity",&sphericity,"sphericity/F");
    TreeTest.Branch("cosThetaB1",&cosThetaB1,"cosThetaB1/F");
    TreeTest.Branch("cosThetaB2",&cosThetaB2,"cosThetaB2/F");
    TreeTest.Branch("cosThetaB3",&cosThetaB3,"cosThetaB3/F");
    TreeTest.Branch("cosThetaB4",&cosThetaB4,"cosThetaB4/F");
    TreeTest.Branch("sumPt",&sumPt,"sumPt/F");
    TreeTest.Branch("jetB1Pt",&jetB1Pt,"jetB1Pt/F");
    TreeTest.Branch("jetB2Pt",&jetB2Pt,"jetB2Pt/F");
    TreeTest.Branch("jetB3Pt",&jetB3Pt,"jetB3Pt/F");
    TreeTest.Branch("jetB4Pt",&jetB4Pt,"jetB4Pt/F");
    TreeTest.Branch("jetB1M",&jetB1M,"jetB1M/F");
    TreeTest.Branch("jetB2M",&jetB2M,"jetB2M/F");
    TreeTest.Branch("jetB3M",&jetB3M,"jetB3M/F");
    TreeTest.Branch("jetB4M",&jetB4M,"jetB4M/F");
    TreeTest.Branch("nParticles",&nParticles,"nParticles/F");
    TreeTest.Branch("constSizeB1",&constSizeB1,"constSizeB1/F");
    TreeTest.Branch("constSizeB2",&constSizeB2,"constSizeB2/F");
    TreeTest.Branch("constSizeB3",&constSizeB3,"constSizeB3/F");
    TreeTest.Branch("constSizeB4",&constSizeB4,"constSizeB4/F");
    TreeTest.Branch("minConstSize",&minConstSize,"minConstSize/F");
    TreeTest.Branch("jetNObjects",&jetNObjects,"jetNObjects/F");
    TreeTest.Branch("minJetNObjects",&minJetNObjects,"minJetNObjects/F");
    TreeTest.Branch("invMassB1AntiKt",&invMassB1AntiKt,"invMassB1AntiKt/F");
    TreeTest.Branch("invMassB2AntiKt",&invMassB2AntiKt,"invMassB2AntiKt/F");
    TreeTest.Branch("nJetsAntiKt",&nJetsAntiKt,"nJetsAntiKt/F");
    TreeTest.Branch("invMassB11Best",&invMassB11Best,"invMassB11Best/F");
    TreeTest.Branch("invMassB21Best",&invMassB21Best,"invMassB21Best/F");
    TreeTest.Branch("invMassB12Best",&invMassB12Best,"invMassB12Best/F");
    TreeTest.Branch("invMassB22Best",&invMassB22Best,"invMassB22Best/F");
    TreeTest.Branch("invMassB13Best",&invMassB13Best,"invMassB13Best/F");
    TreeTest.Branch("invMassB23Best",&invMassB23Best,"invMassB23Best/F");
    TreeTest.Branch("invMassB14Best",&invMassB14Best,"invMassB14Best/F");
    TreeTest.Branch("invMassB24Best",&invMassB24Best,"invMassB24Best/F");
    TreeTest.Branch("invMassB15Best",&invMassB15Best,"invMassB15Best/F");
    TreeTest.Branch("invMassB25Best",&invMassB25Best,"invMassB25Best/F");
    TreeTest.Branch("invMassB16Best",&invMassB16Best,"invMassB16Best/F");
    TreeTest.Branch("invMassB26Best",&invMassB26Best,"invMassB26Best/F");
    TreeTest.Branch("invMassB17Best",&invMassB17Best,"invMassB17Best/F");
    TreeTest.Branch("invMassB27Best",&invMassB27Best,"invMassB27Best/F");
    TreeTest.Branch("invMassB18Best",&invMassB18Best,"invMassB18Best/F");
    TreeTest.Branch("invMassB28Best",&invMassB28Best,"invMassB28Best/F");
    //TreeTest.Branch("nJetsDurham0",&nJetsDurham0,"nJetsDurham0/F");
    //TreeTest.Branch("nJetsDurham5",&nJetsDurham5,"nJetsDurham5/F");
    //TreeTest.Branch("nJetsDurham10",&nJetsDurham10,"nJetsDurham10/F");
    TreeTest.Branch("nJetsDurham15",&nJetsDurham15,"nJetsDurham15/F");
    TreeTest.Branch("nJetsDurham20",&nJetsDurham20,"nJetsDurham20/F");
    TreeTest.Branch("nJetsDurham25",&nJetsDurham25,"nJetsDurham25/F");
    TreeTest.Branch("nJetsDurham30",&nJetsDurham30,"nJetsDurham30/F");
    TreeTest.Branch("distanceZ1MinChiSquaredZZMass",&distanceZ1MinChiSquaredZZMass,"distanceZ1MinChiSquaredZZMass/F");
    TreeTest.Branch("distanceZ2MinChiSquaredZZMass",&distanceZ2MinChiSquaredZZMass,"distanceZ2MinChiSquaredZZMass/F");
    TreeTest.Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
    TreeTest.Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
    TreeTest.Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
    TreeTest.Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
    TreeTest.Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");
    TreeTest.Branch("invMassZZ1",&invMassZZ1,"invMassZZ1/F");
    TreeTest.Branch("invMassZZ2",&invMassZZ2,"invMassZZ2/F");
    TreeTest.Branch("thrust",&thrust,"thrust/F");
    TreeTest.Branch("boostB1",&boostB1,"boostB1/F");
    TreeTest.Branch("boostB2",&boostB2,"boostB2/F");
    TreeTest.Branch("boostB3",&boostB3,"boostB3/F");
    TreeTest.Branch("boostB4",&boostB4,"boostB4/F");
    TreeTest.Branch("boostSystem",&boostSystem,"boostSystem/F");
    TreeTest.Branch("missingET",&missingET,"missingET/F");

    //cout<<"New trees branches assigned"<<endl;

    fullTree->SetBranchAddress("entryIndex", &entryIndex);
    fullTree->SetBranchAddress("aplanarity", &aplanarity);
    fullTree->SetBranchAddress("invMassB1", &invMassB1);
    fullTree->SetBranchAddress("invMassB2", &invMassB2);
    fullTree->SetBranchAddress("minJetM", &minJetM);
    fullTree->SetBranchAddress("sphericity", &sphericity);
    fullTree->SetBranchAddress("cosThetaB1", &cosThetaB1);
    fullTree->SetBranchAddress("cosThetaB2", &cosThetaB2);
    fullTree->SetBranchAddress("cosThetaB3", &cosThetaB3);
    fullTree->SetBranchAddress("cosThetaB4", &cosThetaB4);
    fullTree->SetBranchAddress("sumPt", &sumPt);
    fullTree->SetBranchAddress("jetB1Pt", &jetB1Pt);
    fullTree->SetBranchAddress("jetB2Pt", &jetB2Pt);
    fullTree->SetBranchAddress("jetB3Pt", &jetB3Pt);
    fullTree->SetBranchAddress("jetB4Pt", &jetB4Pt);
    fullTree->SetBranchAddress("jetB1M", &jetB1M);
    fullTree->SetBranchAddress("jetB2M", &jetB2M);
    fullTree->SetBranchAddress("jetB3M", &jetB3M);
    fullTree->SetBranchAddress("jetB4M", &jetB4M);
    fullTree->SetBranchAddress("nParticles", &nParticles);
    fullTree->SetBranchAddress("constSizeB1", &constSizeB1);
    fullTree->SetBranchAddress("constSizeB2", &constSizeB2);
    fullTree->SetBranchAddress("constSizeB3", &constSizeB3);
    fullTree->SetBranchAddress("constSizeB4", &constSizeB4);
    fullTree->SetBranchAddress("minConstSize", &minConstSize);
    fullTree->SetBranchAddress("jetNObjects", &jetNObjects);
    fullTree->SetBranchAddress("minJetNObjects", &minJetNObjects);
    fullTree->SetBranchAddress("invMassB1AntiKt", &invMassB1AntiKt);
    fullTree->SetBranchAddress("invMassB2AntiKt", &invMassB2AntiKt);
    fullTree->SetBranchAddress("nJetsAntiKt", &nJetsAntiKt);
    fullTree->SetBranchAddress("invMassB11Best", &invMassB11Best);
    fullTree->SetBranchAddress("invMassB21Best", &invMassB21Best);
    fullTree->SetBranchAddress("invMassB12Best", &invMassB12Best);
    fullTree->SetBranchAddress("invMassB22Best", &invMassB22Best);
    fullTree->SetBranchAddress("invMassB13Best", &invMassB13Best);
    fullTree->SetBranchAddress("invMassB23Best", &invMassB23Best);
    fullTree->SetBranchAddress("invMassB14Best", &invMassB14Best);
    fullTree->SetBranchAddress("invMassB24Best", &invMassB24Best);
    fullTree->SetBranchAddress("invMassB15Best", &invMassB15Best);
    fullTree->SetBranchAddress("invMassB25Best", &invMassB25Best);
    fullTree->SetBranchAddress("invMassB16Best", &invMassB16Best);
    fullTree->SetBranchAddress("invMassB26Best", &invMassB26Best);
    fullTree->SetBranchAddress("invMassB17Best", &invMassB17Best);
    fullTree->SetBranchAddress("invMassB27Best", &invMassB27Best);
    fullTree->SetBranchAddress("invMassB18Best", &invMassB18Best);
    fullTree->SetBranchAddress("invMassB28Best", &invMassB28Best);
    //fullTree->SetBranchAddress("nJetsDurham0", &nJetsDurham0);
    //fullTree->SetBranchAddress("nJetsDurham5", &nJetsDurham5);
    //fullTree->SetBranchAddress("nJetsDurham10", &nJetsDurham10);
    fullTree->SetBranchAddress("nJetsDurham15", &nJetsDurham15);
    fullTree->SetBranchAddress("nJetsDurham20", &nJetsDurham20);
    fullTree->SetBranchAddress("nJetsDurham25", &nJetsDurham25);
    fullTree->SetBranchAddress("nJetsDurham30", &nJetsDurham30);
    fullTree->SetBranchAddress("distanceZ1MinChiSquaredZZMass", &distanceZ1MinChiSquaredZZMass);
    fullTree->SetBranchAddress("distanceZ2MinChiSquaredZZMass", &distanceZ2MinChiSquaredZZMass);
    fullTree->SetBranchAddress("exclYmerge12", &exclYmerge12);
    fullTree->SetBranchAddress("exclYmerge23", &exclYmerge23);
    fullTree->SetBranchAddress("exclYmerge34", &exclYmerge34);
    fullTree->SetBranchAddress("exclYmerge45", &exclYmerge45);
    fullTree->SetBranchAddress("exclYmerge56", &exclYmerge56);
    fullTree->SetBranchAddress("invMassZZ1", &invMassZZ1);
    fullTree->SetBranchAddress("invMassZZ2", &invMassZZ2);
    fullTree->SetBranchAddress("thrust", &thrust);
    fullTree->SetBranchAddress("boostB1", &boostB1);
    fullTree->SetBranchAddress("boostB2", &boostB2);
    fullTree->SetBranchAddress("boostB3", &boostB3);
    fullTree->SetBranchAddress("boostB4", &boostB4);
    fullTree->SetBranchAddress("boostSystem", &boostSystem);
    fullTree->SetBranchAddress("missingET", &missingET);

    //cout<<"full tree branches assigned"<<endl;

    
    for(int entry=0; entry<fullTree->GetEntries(); entry++)
    {
        fullTree->GetEntry(entry);
        random_device rd;  // Will be used to obtain a seed for the random number engine
        mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        uniform_int_distribution<> distrib(1, 4); // Uniform distribution between 1 and 4
        if(distrib(gen) != 1)
        { 
            TreeTrain.Fill();
            //if(contEventsPostFilter<30) cout<<"contEventsPostFilter: "<<contEventsPostFilter<<"    NO"<<endl;
        }
        else
        { 
            TreeTest.Fill();
            //if(contEventsPostFilter<30) cout<<"contEventsPostFilter: "<<contEventsPostFilter<<"    EN CUATRO"<<endl;
        }
    }

    return;
}

void samplingTrainTest(string rtdCut = "invalid", string preselection = "", string varVersion = "invalid", string sampleName = "")
{

    // Open background files and get their trees
    std::vector<std::string> fileNames = {
        "analysis/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root",
        "analysis/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root",
        "analysis/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root",
        "analysis/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root",
        "analysis/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root",
        "analysis/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root",
        "analysis/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root",
        "analysis/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root"
    };

    std::vector<TFile*> files;
    for (const auto& fileName : fileNames) {
        files.push_back(TFile::Open(fileName.c_str()));
    }

    //cout<<"files opened"<<endl;

    vector<string> trees = {
        "TreeS",
        "TreeBqq",
        "TreeBtt",
        "TreeBZZ",
        "TreeBWW",
        "TreeBqqX",
        "TreeBqqqqX",
        "TreeBqqHX"
    };

    //cout<<"Trees names retrieved"<<endl;

    string jetAlgoOutputTreeSTrain = "analysis/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeSTest = "analysis/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	string jetAlgoOutputTreeBqqTrain = "analysis/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBqqTest = "analysis/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	string jetAlgoOutputTreeBttTrain = "analysis/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBttTest = "analysis/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	string jetAlgoOutputTreeBZZTrain = "analysis/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBZZTest = "analysis/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	string jetAlgoOutputTreeBWWTrain = "analysis/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBWWTest = "analysis/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBqqXTrain = "analysis/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBqqXTest = "analysis/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBqqqqXTrain = "analysis/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBqqqqXTest = "analysis/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBqqHXTrain = "analysis/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBqqHXTest = "analysis/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    
    TFile *outputTreeSTrain = new TFile(jetAlgoOutputTreeSTrain.c_str(), "recreate");
	TTree TreeSTrain("TreeSTrain","a simple Tree with simple variables (Train)");
	TFile *outputTreeSTest = new TFile(jetAlgoOutputTreeSTest.c_str(), "recreate");
	TTree TreeSTest("TreeSTest","a simple Tree with simple variables (Test)");
    TFile *outputTreeBqqTrain = new TFile(jetAlgoOutputTreeBqqTrain.c_str(), "recreate");
    TTree TreeBqqTrain("TreeBqqTrain","a bqqimple Tree with bqqimple variables (Train)");
    TFile *outputTreeBqqTest = new TFile(jetAlgoOutputTreeBqqTest.c_str(), "recreate");
    TTree TreeBqqTest("TreeBqqTest","a bqqimple Tree with bqqimple variables (Test)");
    TFile *outputTreeBttTrain = new TFile(jetAlgoOutputTreeBttTrain.c_str(), "recreate");
    TTree TreeBttTrain("TreeBttTrain","a bttimple Tree with bttimple variables (Train)");
    TFile *outputTreeBttTest = new TFile(jetAlgoOutputTreeBttTest.c_str(), "recreate");
    TTree TreeBttTest("TreeBttTest","a bttimple Tree with bttimple variables (Test)");
    TFile *outputTreeBZZTrain = new TFile(jetAlgoOutputTreeBZZTrain.c_str(), "recreate");
    TTree TreeBZZTrain("TreeBZZTrain","a bZZimple Tree with bZZimple variables (Train)");
    TFile *outputTreeBZZTest = new TFile(jetAlgoOutputTreeBZZTest.c_str(), "recreate");
    TTree TreeBZZTest("TreeBZZTest","a bZZimple Tree with bZZimple variables (Test)");
    TFile *outputTreeBWWTrain = new TFile(jetAlgoOutputTreeBWWTrain.c_str(), "recreate");
    TTree TreeBWWTrain("TreeBWWTrain","a bWWimple Tree with bWWimple variables (Train)");
    TFile *outputTreeBWWTest = new TFile(jetAlgoOutputTreeBWWTest.c_str(), "recreate");
    TTree TreeBWWTest("TreeBWWTest","a bWWimple Tree with bWWimple variables (Test)");
    TFile *outputTreeBqqXTrain = new TFile(jetAlgoOutputTreeBqqXTrain.c_str(), "recreate");
    TTree TreeBqqXTrain("TreeBqqXTrain","a bqqXimple Tree with bqqXimple variables (Train)");
    TFile *outputTreeBqqXTest = new TFile(jetAlgoOutputTreeBqqXTest.c_str(), "recreate");
    TTree TreeBqqXTest("TreeBqqXTest","a bqqXimple Tree with bqqXimple variables (Test)");
    TFile *outputTreeBqqqqXTrain = new TFile(jetAlgoOutputTreeBqqqqXTrain.c_str(), "recreate");
    TTree TreeBqqqqXTrain("TreeBqqqqXTrain","a bqqqqXimple Tree with bqqqqXimple variables (Train)");
    TFile *outputTreeBqqqqXTest = new TFile(jetAlgoOutputTreeBqqqqXTest.c_str(), "recreate");
    TTree TreeBqqqqXTest("TreeBqqqqXTest","a bqqqqXimple Tree with bqqqqXimple variables (Test)");
    TFile *outputTreeBqqHXTrain = new TFile(jetAlgoOutputTreeBqqHXTrain.c_str(), "recreate");
    TTree TreeBqqHXTrain("TreeBqqHXTrain","a bqqHXimple Tree with bqqHXimple variables (Train)");
    TFile *outputTreeBqqHXTest = new TFile(jetAlgoOutputTreeBqqHXTest.c_str(), "recreate");
    TTree TreeBqqHXTest("TreeBqqHXTest","a bqqHXimple Tree with bqqHXimple variables (Test)");

    TTree TreeSenTrain("TreeSenTrain","a simple Tree with simple variables (Train sen)");
    TTree TreeSenTest("TreeSenTest","a simple Tree with simple variables (Test sen)");

    //cout<<"new trees created"<<endl;

    //cout<<"calling functions"<<endl;
    fillOutTrees(TreeSTrain, TreeSTest, fileNames[0], trees[0]);
    fillOutTrees(TreeBqqTrain, TreeBqqTest, fileNames[1], trees[1]);
    fillOutTrees(TreeBttTrain, TreeBttTest, fileNames[2], trees[2]);
    fillOutTrees(TreeBZZTrain, TreeBZZTest, fileNames[3], trees[3]);
    fillOutTrees(TreeBWWTrain, TreeBWWTest, fileNames[4], trees[4]);
    fillOutTrees(TreeBqqXTrain, TreeBqqXTest, fileNames[5], trees[5]);
    fillOutTrees(TreeBqqqqXTrain, TreeBqqqqXTest, fileNames[6], trees[6]);
    fillOutTrees(TreeBqqHXTrain, TreeBqqHXTest, fileNames[7], trees[7]);
    //cout<<"functions done"<<endl;

    outputTreeSTrain->cd();
    TreeSTrain.Write();
    outputTreeSTest->cd();
    TreeSTest.Write();
    outputTreeBqqTrain->cd();
    TreeBqqTrain.Write();
    outputTreeBqqTest->cd();
    TreeBqqTest.Write();
    outputTreeBttTrain->cd();
    TreeBttTrain.Write();
    outputTreeBttTest->cd();
    TreeBttTest.Write();
    outputTreeBZZTrain->cd();
    TreeBZZTrain.Write();
    outputTreeBZZTest->cd();
    TreeBZZTest.Write();
    outputTreeBWWTrain->cd();
    TreeBWWTrain.Write();
    outputTreeBWWTest->cd();
    TreeBWWTest.Write();
    outputTreeBqqXTrain->cd();
    TreeBqqXTrain.Write();
    outputTreeBqqXTest->cd();
    TreeBqqXTest.Write();
    outputTreeBqqqqXTrain->cd();
    TreeBqqqqXTrain.Write();
    outputTreeBqqqqXTest->cd();
    TreeBqqqqXTest.Write();
    outputTreeBqqHXTrain->cd();
    TreeBqqHXTrain.Write();
    outputTreeBqqHXTest->cd();
    TreeBqqHXTest.Write();

    outputTreeSTrain->Close();
    outputTreeSTest->Close();
    outputTreeBqqTrain->Close();
    outputTreeBqqTest->Close();
    outputTreeBttTrain->Close();
    outputTreeBttTest->Close();
    outputTreeBZZTrain->Close();
    outputTreeBZZTest->Close();
    outputTreeBWWTrain->Close();
    outputTreeBWWTest->Close();
    outputTreeBqqXTrain->Close();
    outputTreeBqqXTest->Close();
    outputTreeBqqqqXTrain->Close();
    outputTreeBqqqqXTest->Close();
    outputTreeBqqHXTrain->Close();
    outputTreeBqqHXTest->Close();

    //cout<<"done"<<endl;

}

/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This executable gives an example of a very simple use of the genetic algorithm
/// of TMVA
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Executable: TMVACutsGA
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker
 
#include <iostream> // Stream declarations
#include <vector>
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include <cstdlib>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include <algorithm>
#include <random>
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
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <stdexcept>
#include <fstream>
 
using std::vector;
 
using namespace TMVA;

// Function declarations for calculating yields
Double_t CalculateSignalYield(TTree* signalTree, const vector<Double_t>& cuts, double weightHH, int& cont0);
Double_t CalculateBackgroundYield(const vector<TTree*>& backgroundTrees, const vector<Double_t>& cuts, double weightqq, double weightttbar, double weightZZ, double weightWW, double weightqqX, double weightqqqqX, double weightqqHX, double weightZH);
 
class MyFitness : public IFitterTarget {
    public:
       MyFitness(TTree* sigTree, const vector<TTree*>& bgTrees, double weightHH, double weightqq, double weightttbar, double weightZZ, double weightWW, double weightqqX, double weightqqqqX, double weightqqHX, double weightZH) 
        : signalTree(sigTree), backgroundTrees(bgTrees), weightHH(weightHH), weightqq(weightqq), weightttbar(weightttbar), weightZZ(weightZZ), weightWW(weightWW), weightqqX(weightqqX), weightqqqqX(weightqqqqX), weightqqHX(weightqqHX), weightZH(weightZH) {}
       // Estimator Function to maximize signal significance
      Double_t EstimatorFunction(std::vector<Double_t> & factors) override {
         // Calculate signal and background yields
         Double_t signalYield = CalculateSignalYield(signalTree, factors, weightHH, cont0);
         Double_t backgroundYield = CalculateBackgroundYield(backgroundTrees, factors, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, weightZH);

         // Avoid division by zero
         if (signalYield + backgroundYield == 0) return std::numeric_limits<Double_t>::max();

         // Return negative significance since we want to maximize it
         //cout<<"Significance: "<<(signalYield / sqrt(signalYield + backgroundYield))<<endl<<endl;
         return - (signalYield / sqrt(signalYield + backgroundYield));
      }
   private:
   TTree* signalTree;
   vector<TTree*> backgroundTrees;
   double weightHH;
   double weightqq;
   double weightttbar;
   double weightZZ;
   double weightWW;
   double weightqqX;
   double weightqqqqX;
   double weightqqHX;
   double weightZH;
   int cont0 = 0; 
};
 
 
class MyGA2nd : public GeneticAlgorithm {
    public:
       MyGA2nd( IFitterTarget& target, Int_t size, vector<Interval*>& ranges ) : GeneticAlgorithm(target,
       size, ranges ){
       }
};



// Function to calculate signal yield
Double_t CalculateSignalYield(TTree* signalTree, const vector<Double_t>& cuts, double weightHH, int& cont0) {
    Double_t signalYield=0;
    Float_t NNOutputs[8] = {0}; // Initialize to zero to avoid uninitialized values

    signalTree->SetBranchAddress("NN1Output", &NNOutputs[0]);
    signalTree->SetBranchAddress("NN2Output", &NNOutputs[1]);
    signalTree->SetBranchAddress("NN3Output", &NNOutputs[2]);
    signalTree->SetBranchAddress("NN4Output", &NNOutputs[3]);
    signalTree->SetBranchAddress("NN5Output", &NNOutputs[4]);
    signalTree->SetBranchAddress("NN6Output", &NNOutputs[5]);
    signalTree->SetBranchAddress("NN7Output", &NNOutputs[6]);
    signalTree->SetBranchAddress("NN8Output", &NNOutputs[7]);

    Long64_t nEntries = signalTree->GetEntries();
    int nBacks = sizeof(NNOutputs) / sizeof(NNOutputs[0]);

    bool goodCuts=true;
    /*for(int i=0; i<nBacks; i++)
    {
        if(abs(cuts[i]) == 1)
        {
            goodCuts = false;
            break;
        }
    }*/

    //cout<<"Initial --    weightHH: "<<weightHH<<"    signalYield: "<<signalYield<<"    nEntries: "<<nEntries<<"    nBacks: "<<nBacks<<"    goodCuts: "<<goodCuts;
    // Loop over all entries in the tree
    for (Long64_t i = 0; i < nEntries; ++i) {
        signalTree->GetEntry(i);
        
        //if(i<1000 && cont0 == 1) cout<<"Entry: "<<i;
        // Apply all cuts
        bool passesCuts = true;
        for (int j = 0; j < nBacks; ++j) {

            //if(i<1000 && cont0 == 1) cout<<"    NNOutputs"<<j<<": "<<NNOutputs[j]<<"    cuts"<<j<<": "<<cuts[j]<<endl;

            /*if (NNOutputs[j] == -999 || isnan(NNOutputs[j]) || isinf(NNOutputs[j])) {
                if (i < 10) cout << "Invalid NNOutputs[" << j << "] detected: " << NNOutputs[j] << endl;
                //passesCuts = false;
                //break; // Skip this entry entirely
            }*/

            /*if(NNOutputs[j] == -999) 
            {
                //if(i<10) cout<<"-999"<<endl;
                continue;
            }*/

           

            if (NNOutputs[j] < cuts[j]) {
                //if(i<1000 && cont0 == 1) cout<<"DOESNT PASS"<<endl<<endl;
                passesCuts = false;
                break;
            }
        }

        if (passesCuts) signalYield++;
    }
    /*cout<<endl<<"For cuts: ";
    for(int i=0; i<nBacks; i++) cout<<cuts[i]<<", ";
    cout<<endl<<"Surviving HH: "<<signalYield*weightHH<<" (uw: "<<signalYield<<")"<<endl;*/
    /*if(signalYield == 0 && goodCuts)
    {
        cont0++;
        if(cont0 == 1) CalculateSignalYield(signalTree, cuts, weightHH, cont0);
    }*/
    return signalYield*weightHH;
}

// Function to calculate background yield
Double_t CalculateBackgroundYield(const vector<TTree*>& backgroundTrees, const vector<Double_t>& cuts, double weightqq, double weightttbar, double weightZZ, double weightWW, double weightqqX, double weightqqqqX, double weightqqHX, double weightZH) {
    Double_t totalBackgroundYield = 0;
    int contBacks=0;
    // Loop over all background trees
    for (TTree* bgTree : backgroundTrees) {
        Float_t NNOutputs[8];
        bgTree->SetBranchAddress("NN1Output", &NNOutputs[0]);
        bgTree->SetBranchAddress("NN2Output", &NNOutputs[1]);
        bgTree->SetBranchAddress("NN3Output", &NNOutputs[2]);
        bgTree->SetBranchAddress("NN4Output", &NNOutputs[3]);
        bgTree->SetBranchAddress("NN5Output", &NNOutputs[4]);
        bgTree->SetBranchAddress("NN6Output", &NNOutputs[5]);
        bgTree->SetBranchAddress("NN7Output", &NNOutputs[6]);
        bgTree->SetBranchAddress("NN8Output", &NNOutputs[7]);

        Long64_t nEntries = bgTree->GetEntries();
        Double_t backgroundYield = 0;
        int nBacks = sizeof(NNOutputs) / sizeof(NNOutputs[0]);


        // Loop over all entries in the tree
        for (Long64_t i = 0; i < nEntries; ++i) {
            bgTree->GetEntry(i);

            // Apply all cuts
            bool passesCuts = true;
            for (int j = 0; j < nBacks; ++j) {
                if (NNOutputs[j] < cuts[j]) {
                    passesCuts = false;
                    break;
                }
            }

            if (passesCuts) backgroundYield++;
        }
        double weightBack;
        string nBack;
        if(contBacks==0)
        {
            weightBack = weightqq;
            nBack = "qq";
        } 
        else if(contBacks==1)
        {
            weightBack = weightttbar;
            nBack = "ttbar";
        }
        else if(contBacks==2)
        {
            weightBack = weightZZ;
            nBack = "ZZ";
        }
        else if(contBacks==3)
        {
            weightBack = weightWW;
            nBack = "WW";
        }
        else if(contBacks==4)
        {
            weightBack = weightqqX;
            nBack = "qqX";
        }
        else if(contBacks==5)
        {
            weightBack = weightqqqqX;
            nBack = "qqqqX";
        }
        else if(contBacks==6)
        {
            weightBack = weightqqHX;
            nBack = "qqHX";
        }
        else if(contBacks==7)
        {
            weightBack = weightZH;
            nBack = "ZH";
        }
        totalBackgroundYield += backgroundYield*weightBack;
        contBacks++;
        //cout<<"Surviving "<<nBack<<": "<<backgroundYield*weightBack<<" (uw: "<<backgroundYield<<")"<<endl;
    }

   //cout<<"total surviving back: "<<totalBackgroundYield<<endl<<endl;
    return totalBackgroundYield;
}

////Function that calculates the appropiate weights for each topology
void getWeights(double& weightHH, double& weightqq, double& weightttbar, double& weightZZ, double& weightWW, double& weightqqX, double& weightqqqqX, double& weightqqHX, double& weightZH, string rtdCut, string preselection, string varVersion, string sampleName)
{
    string inputTrainHHText="analysis/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
 	string inputTestHHText="analysis/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
 	TFile* inputTrainHH = new TFile(inputTrainHHText.c_str());
    TTree* theTreeTrainHH = (TTree*)inputTrainHH->Get("TreeSTrain");
 	TFile* inputTestHH = new TFile(inputTestHHText.c_str());
    TTree* theTreeTestHH = (TTree*)inputTestHH->Get("TreeSTest");
    double weightFactorHH = (theTreeTestHH->GetEntries()+theTreeTrainHH->GetEntries());
    weightFactorHH = weightFactorHH/(theTreeTestHH->GetEntries());
    string inputTrainqqText="analysis/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
 	string inputTestqqText="analysis/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputTrainqq = new TFile(inputTrainqqText.c_str());
    TTree* theTreeTrainqq = (TTree*)inputTrainqq->Get("TreeBqqTrain");
 	TFile* inputTestqq = new TFile(inputTestqqText.c_str());
    TTree* theTreeTestqq = (TTree*)inputTestqq->Get("TreeBqqTest");
    double weightFactorqq = (theTreeTestqq->GetEntries()+theTreeTrainqq->GetEntries());
    weightFactorqq = weightFactorqq/(theTreeTestqq->GetEntries());
    string inputTrainttbarText="analysis/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
 	string inputTestttbarText="analysis/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputTrainttbar = new TFile(inputTrainttbarText.c_str());
    TTree* theTreeTrainttbar = (TTree*)inputTrainttbar->Get("TreeBttTrain");
 	TFile* inputTestttbar = new TFile(inputTestttbarText.c_str());
    TTree* theTreeTestttbar = (TTree*)inputTestttbar->Get("TreeBttTest");
    double weightFactorttbar = (theTreeTestttbar->GetEntries()+theTreeTrainttbar->GetEntries());
    weightFactorttbar = weightFactorttbar/(theTreeTestttbar->GetEntries());
    string inputTrainZZText="analysis/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
 	string inputTestZZText="analysis/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputTrainZZ = new TFile(inputTrainZZText.c_str());
    TTree* theTreeTrainZZ = (TTree*)inputTrainZZ->Get("TreeBZZTrain");
 	TFile* inputTestZZ = new TFile(inputTestZZText.c_str());
    TTree* theTreeTestZZ = (TTree*)inputTestZZ->Get("TreeBZZTest");
    double weightFactorZZ = (theTreeTestZZ->GetEntries()+theTreeTrainZZ->GetEntries());
    weightFactorZZ = weightFactorZZ/(theTreeTestZZ->GetEntries());
    string inputTrainWWText="analysis/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
 	string inputTestWWText="analysis/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputTrainWW = new TFile(inputTrainWWText.c_str());
    TTree* theTreeTrainWW = (TTree*)inputTrainWW->Get("TreeBWWTrain");
 	TFile* inputTestWW = new TFile(inputTestWWText.c_str());
    TTree* theTreeTestWW = (TTree*)inputTestWW->Get("TreeBWWTest");
    double weightFactorWW = (theTreeTestWW->GetEntries()+theTreeTrainWW->GetEntries());
    weightFactorWW = weightFactorWW/(theTreeTestWW->GetEntries());
	string inputTrainqqXText="analysis/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	string inputTestqqXText="analysis/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile* inputTrainqqX = new TFile(inputTrainqqXText.c_str());
	TTree* theTreeTrainqqX = (TTree*)inputTrainqqX->Get("TreeBqqXTrain");
	TFile* inputTestqqX = new TFile(inputTestqqXText.c_str());
	TTree* theTreeTestqqX = (TTree*)inputTestqqX->Get("TreeBqqXTest");
	double weightFactorqqX = (theTreeTestqqX->GetEntries()+theTreeTrainqqX->GetEntries());
	weightFactorqqX = weightFactorqqX/(theTreeTestqqX->GetEntries());
	string inputTrainqqqqXText="analysis/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	string inputTestqqqqXText="analysis/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile* inputTrainqqqqX = new TFile(inputTrainqqqqXText.c_str());
	TTree* theTreeTrainqqqqX = (TTree*)inputTrainqqqqX->Get("TreeBqqqqXTrain");
	TFile* inputTestqqqqX = new TFile(inputTestqqqqXText.c_str());
	TTree* theTreeTestqqqqX = (TTree*)inputTestqqqqX->Get("TreeBqqqqXTest");
	double weightFactorqqqqX = (theTreeTestqqqqX->GetEntries()+theTreeTrainqqqqX->GetEntries());
	weightFactorqqqqX = weightFactorqqqqX/(theTreeTestqqqqX->GetEntries());
	string inputTrainqqHXText="analysis/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	string inputTestqqHXText="analysis/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile* inputTrainqqHX = new TFile(inputTrainqqHXText.c_str());
	TTree* theTreeTrainqqHX = (TTree*)inputTrainqqHX->Get("TreeBqqHXTrain");
	TFile* inputTestqqHX = new TFile(inputTestqqHXText.c_str());
	TTree* theTreeTestqqHX = (TTree*)inputTestqqHX->Get("TreeBqqHXTest");
	double weightFactorqqHX = (theTreeTestqqHX->GetEntries()+theTreeTrainqqHX->GetEntries());
	weightFactorqqHX = weightFactorqqHX/(theTreeTestqqHX->GetEntries());
    string inputTrainZHText="analysis/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	string inputTestZHText="analysis/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile* inputTrainZH = new TFile(inputTrainZHText.c_str());
	TTree* theTreeTrainZH = (TTree*)inputTrainZH->Get("TreeBZHTrain");
	TFile* inputTestZH = new TFile(inputTestZHText.c_str());
	TTree* theTreeTestZH = (TTree*)inputTestZH->Get("TreeBZHTest");
	double weightFactorZH = (theTreeTestZH->GetEntries()+theTreeTrainZH->GetEntries());
	weightFactorZH = weightFactorZH/(theTreeTestZH->GetEntries());

    cout<<"WeightFactorHH: "<<weightFactorHH<<endl;
    cout<<"WeightFactorqq: "<<weightFactorqq<<endl;
    cout<<"WeightFactorttbar: "<<weightFactorttbar<<endl;
    cout<<"WeightFactorZZ: "<<weightFactorZZ<<endl;
    cout<<"WeightFactorWW: "<<weightFactorWW<<endl;
    cout<<"WeightFactorqqX: "<<weightFactorqqX<<endl;
    cout<<"WeightFactorqqqqX: "<<weightFactorqqqqX<<endl;
    cout<<"WeightFactorqqHX: "<<weightFactorqqHX<<endl;
    cout<<"WeightFactorZH: "<<weightFactorZH<<endl;

    weightHH = 0.001225 * weightFactorHH;
    weightqq = 0.0349 * weightFactorqq;
    weightttbar = 0.503 * weightFactorttbar;
    weightZZ = 0.8167 * weightFactorZZ;
    weightWW = 0.5149 * weightFactorWW;
    weightqqX = 0.04347826 * weightFactorqqX;
    weightqqqqX = 0.04 * weightFactorqqqqX;
    weightqqHX = 0.001 * weightFactorqqHX;
    weightZH = 0.00207445 * weightFactorZH;

    cout<<"weightHH: "<<weightHH<<endl;
    cout<<"weightqq: "<<weightqq<<endl;
    cout<<"weightttbar: "<<weightttbar<<endl;
    cout<<"weightZZ: "<<weightZZ<<endl;
    cout<<"weightWW: "<<weightWW<<endl;
    cout<<"weightqqX: "<<weightqqX<<endl;
    cout<<"weightqqqqX: "<<weightqqqqX<<endl;
    cout<<"weightqqHX: "<<weightqqHX<<endl;
    cout<<"weightZH: "<<weightZH<<endl;
} 

//////function that calculates cross section and error
void findCrossSectionHHbbbb(double totalRemaining, double HHRemaining, double backRemaining, double luminosity, double& crossSection, double& errorTop, double& errorBottom, int nbin)
{ 
	double xMin = -2000.0;   // Minimum x value
   	double xMax = 10000.0;  // Maximum x value
	
	double branchingRatioHbb=0.5824; ////////////at Higgs mass of 125.00 GeV
	double totalBR = branchingRatioHbb * branchingRatioHbb; 
	double totalEffi=HHRemaining/(1462612*0.001225*totalBR);
	//cout<<"HHRemaining: "<<HHRemaining<<"    (1462612*0.001225): "<<1462612*0.001225<<"    totalEffi: "<<totalEffi<<endl;
	
	double s, b=backRemaining, n=totalRemaining; 
	cout<<"b: "<<b<<"    n: "<<n<<"    n-b: "<<n-b<<endl;
	double Lb, Lsb;
	double chiSquared, minChiSquared, minX;
	int contDebug=0;
	
	////Loop to find minChiSquared
	cout<<"Empieza minChiSquared"<<endl;
	for(double x=xMin; x<=xMax; x+=.01)
	{
		//if(contDebug<5) cout<<"Entra"<<endl<<"x: "<<x<<"    totalEffi: "<<totalEffi<<"    totalBR: "<<totalBR<<endl;
		s=x*totalEffi*totalBR;
		//if(contDebug<5) cout<<"s: "<<s<<"    b: "<<b<<"    n: "<<n<<endl;
      		//chiSquared=-2*TMath::Log((pow(TMath::E(), -(s+b)) * pow((s+b), n))/(pow(TMath::E(), -(b)) * pow((b), n)));
      		chiSquared = -2*(-s+n*TMath::Log((s+b)/b));
      		//if(contDebug<5) cout<<"Lsb: "<<Lsb<<"    Lb: "<<Lb<<"    ChiSquared: "<<chiSquared<<endl<<endl;
	      	if(chiSquared<minChiSquared)
	      	{
	      		minChiSquared=chiSquared;
	      		minX=x;
	      	}
	      	contDebug++;
	}
	TH2F *histChiSquared = new TH2F("hChiSquared", "ChiSquared as a function of cross-section * lumi", nbin, xMin, xMax, nbin, minChiSquared*1.1, 20);
	///Loop to plot chiSquared
	cout<<"Termina minChiSquared e inicia plot chiSquared"<<endl;
	for(double x=xMin; x<=xMax; x+=.01)
	{
		s=x*totalEffi*totalBR;
      		//chiSquared=-2*TMath::Log((pow(TMath::E(), -(s+b)) * pow((s+b), n))/(pow(TMath::E(), -(b)) * pow((b), n)));
      		chiSquared = -2*(-s+n*TMath::Log((s+b)/b));
      		if (chiSquared != 1000000000) histChiSquared->Fill(x, chiSquared);
	}
	cout<<"Termina plot chiSquared"<<endl;
		
	TCanvas *cfindCrossSectionHHbbbb1 = new TCanvas();
	cfindCrossSectionHHbbbb1->cd();
	histChiSquared->Draw();
	
	crossSection = minX / luminosity;
	
	cout<<"Signal efficiency: "<<totalEffi<<endl<<"Branching ratio: "<<totalBR<<endl;
	cout<<"min ChiSquared: "<<minChiSquared<<"     for x: "<<minX<<endl<<"Cross section: "<<crossSection<<endl;
	
	////To find uncertainty
	TLine* lineError = new TLine(xMin,minChiSquared+1,xMax,minChiSquared+1);
	lineError->SetLineColor(kRed);
	lineError->SetLineWidth(2);
	lineError->SetLineStyle(1);
	lineError->Draw("same");
	double xErrorLeft, errorLeft, distanceMinLeft=100000, xErrorRight, errorRight, distanceMinRight=100000; 
	for(double x=xMin; x<=xMax; x+=0.1)
	{
		s=x*totalEffi*totalBR;
      		//chiSquared=-2*TMath::Log((pow(TMath::E(), -(s+b)) * pow((s+b), n))/(pow(TMath::E(), -(b)) * pow((b), n)));
      		chiSquared = -2*(-s+n*TMath::Log((s+b)/b));
	      	if(chiSquared<minChiSquared+1.1 && chiSquared>minChiSquared+0.9)
	      	{
	      		if(x<minX)
	      		{
	      			if(abs(chiSquared-(minChiSquared+1))<distanceMinLeft)
	      			{
	      				distanceMinLeft=abs(chiSquared-(minChiSquared+1));
	      				xErrorLeft=x;
	      			}
	      		}
	      		else
	      		{
	      			if(abs(chiSquared-(minChiSquared+1))<distanceMinRight)
	      			{
	      				distanceMinRight=abs(chiSquared-(minChiSquared+1));
	      				xErrorRight=x;
	      			}
	      		}
	      	}
	      	contDebug++;
	}
	errorLeft=abs((xErrorLeft-minX)/(luminosity*crossSection));	
	errorRight=abs((xErrorRight-minX)/(luminosity*crossSection));	
	cout<<"errorLeft: "<<errorLeft<<"    for x: "<<xErrorLeft<<"    with distance: "<<distanceMinLeft<<endl;
	cout<<"errorRight: "<<errorRight<<"    for x: "<<xErrorRight<<"    with distance: "<<distanceMinRight<<endl<<endl;
	errorTop = errorRight;
	errorBottom = errorLeft;
	////////
	
}
 
void TMVACutsGA(string rtdCut = "invalid", string preselection = "", string varVersion = "invalid", string sampleName = "") 
{
 
   std::cout << "Start Test TMVACutsGA" << std::endl;
   cout<< "========================" << std::endl;


   double weightHH=0, weightqq=0, weightttbar=0, weightZZ=0, weightWW=0, weightqqX=0, weightqqqqX=0, weightqqHX=0, weightZH=0;
   getWeights(weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, weightZH, rtdCut, preselection, varVersion, sampleName);
   int cont0=0;


   // Open signal file
   string signalFileText = "analysis/outputTreeSGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
   TFile* signalFile = TFile::Open(signalFileText.c_str());
   TTree* signalTree = (TTree*)signalFile->Get("TreeSGA");

   // Open background files and get their trees
   std::vector<std::string> backgroundFileNames = {
       "analysis/outputTreeBqqGAESpreadDurham"+rtdCut+preselection+sampleName+".root",
       "analysis/outputTreeBttbarGAESpreadDurham"+rtdCut+preselection+sampleName+".root",
       "analysis/outputTreeBZZGAESpreadDurham"+rtdCut+preselection+sampleName+".root",
       "analysis/outputTreeBWWGAESpreadDurham"+rtdCut+preselection+sampleName+".root",
       "analysis/outputTreeBqqXGAESpreadDurham"+rtdCut+preselection+sampleName+".root",
       "analysis/outputTreeBqqqqXGAESpreadDurham"+rtdCut+preselection+sampleName+".root",
       "analysis/outputTreeBqqHXGAESpreadDurham"+rtdCut+preselection+sampleName+".root",
       "analysis/outputTreeBZHGAESpreadDurham"+rtdCut+preselection+sampleName+".root" 
   };

   std::vector<TFile*> backgroundFiles;
    for (const auto& fileName : backgroundFileNames) {
    backgroundFiles.push_back(TFile::Open(fileName.c_str()));
    }

    vector<TTree*> backgroundTrees = {
        (TTree*)backgroundFiles[0]->Get("TreeBqqGA"),
        (TTree*)backgroundFiles[1]->Get("TreeBttbarGA"),
        (TTree*)backgroundFiles[2]->Get("TreeBZZGA"),
        (TTree*)backgroundFiles[3]->Get("TreeBWWGA"),
        (TTree*)backgroundFiles[4]->Get("TreeBqqXGA"),
        (TTree*)backgroundFiles[5]->Get("TreeBqqqqXGA"),
        (TTree*)backgroundFiles[6]->Get("TreeBqqHXGA"), 
        (TTree*)backgroundFiles[7]->Get("TreeBZHGA")
    };


   cout<<"uw HH events: "<<signalTree->GetEntries()<<endl;
   cout<<"w HH events: "<<(signalTree->GetEntries())*weightHH<<endl;
   int nBacks=0;
   for (TTree* bgTree : backgroundTrees)
   {
        string nBack;
        double weightBack=0;
        if(nBacks==0)
        {
            nBack = "qq";
            weightBack = weightqq;
        }
        else if(nBacks==1)
        {
            nBack = "ttbar";
            weightBack = weightttbar;
        }
        else if(nBacks==2)
        {
            nBack = "ZZ";
            weightBack = weightZZ;
        }
        else if(nBacks==3)
        {
            nBack = "WW";
            weightBack = weightWW;
        }
        else if(nBacks==4)
        {
            nBack = "qqX";
            weightBack = weightqqX;
        }
        else if(nBacks==5)
        {
            nBack = "qqqqX";
            weightBack = weightqqqqX;
        }
        else if(nBacks==6)
        {
            nBack = "qqHX";
            weightBack = weightqqHX;
        }
        else if(nBacks==7)
        {
            nBack = "ZH";
            weightBack = weightZH;
        }
        nBacks++;
        cout<<"uw "<<nBack<<" events: "<<bgTree->GetEntries()<<endl;
        cout<<"w "<<nBack<<" events: "<<(bgTree->GetEntries())*weightBack<<endl;
   }

   // Define ranges for the 8 variables (NN1Output, ..., NN8Output) with min, max, bins
   vector<Interval*> ranges;
   //for (int i = 0; i < backgroundTrees.size(); ++i) {
   for (int i = 0; i < 8; ++i) {
      ranges.push_back(new Interval(-1, 1, 100)); // Using -1 to 1 as range, with 100 bins
   }
 
   // Initialize the fitness function and genetic algorithm
   IFitterTarget* myFitness = new MyFitness(signalTree, backgroundTrees, weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, weightZH);
   //MyGA2nd mg(*myFitness, 100, ranges); // Population size of 100
   MyGA2nd mg(*myFitness, 300, ranges); // Population size of 300
 
   #define CONVSTEPS 50
   #define CONVCRIT 0.0001
   #define SCSTEPS 10
   #define SCRATE 5
   #define SCFACTOR 0.95
 
   do {
      // prepares the new generation and does evolution
      mg.Init();
      // assess the quality of the individuals
      mg.CalculateFitness();
      mg.GetGeneticPopulation().Print(0);
      //std::cout << "---" << std::endl;
      // reduce the population size to the initially defined one
      mg.GetGeneticPopulation().TrimPopulation();
      // tricky thing: control the speed of how fast the "solution space" is searched through
      // this function basically influences the sigma of a gaussian around the actual value
      // of the parameter where the new value will be randomly thrown.
      // when the number of improvements within the last SCSTEPS
      // A) smaller than SCRATE: divide the preset sigma by SCFACTOR
      // B) equal to SCRATE: do nothing
      // C) greater than SCRATE: multiply the preset sigma by SCFACTOR
      // if you don't know what to do, leave it unchanged or even delete this function call
      mg.SpreadControl( SCSTEPS, SCRATE, SCFACTOR );
 
   } while (!mg.HasConverged( CONVSTEPS, CONVCRIT ));  // converged if: fitness-improvement < CONVCRIT within the last CONVSTEPS loops
 
   GeneticGenes* genes = mg.GetGeneticPopulation().GetGenes( 0 );
   std::vector<Double_t> gvec;
   gvec = genes->GetFactors();
   int n = 0;
   for( std::vector<Double_t>::iterator it = gvec.begin(); it<gvec.end(); it++ ){
      std::cout << "FACTOR " << n << " : " << (*it) << std::endl;
      n++;
   }

   cout<<endl<<endl<<endl;

   // Retrieve the best solution (genes)
   genes = mg.GetGeneticPopulation().GetGenes(0);
   std::vector<Double_t> optimalCuts = genes->GetFactors();

   // Calculate signal and background yields using the optimal cuts
   Double_t optimalSignalYield = CalculateSignalYield(signalTree, optimalCuts, weightHH, cont0);
   Double_t optimalBackgroundYield = CalculateBackgroundYield(backgroundTrees, optimalCuts, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, weightZH);
    
   // Calculate the maximum significance
   Double_t maxSignificance = optimalSignalYield / sqrt(optimalSignalYield + optimalBackgroundYield);
   
   // Output the results
   cout << "Optimal Cuts:" << endl;
   for (size_t i = 0; i < optimalCuts.size(); ++i) {
      cout << "Cut on NN" << (i+1) << "Output: " << optimalCuts[i] << endl;
   }
   cout << "Optimal Signal Yield: " << optimalSignalYield << endl;
   cout << "Optimal Background Yield: " << optimalBackgroundYield << endl;
   cout << "Maximum Significance: " << maxSignificance << endl;

   double crossSectionCombined, errorTopCombined, errorBottomCombined, luminosity=4900, nbin=10000;
   findCrossSectionHHbbbb(optimalSignalYield+optimalBackgroundYield, optimalSignalYield, optimalBackgroundYield, luminosity, crossSectionCombined, errorTopCombined, errorBottomCombined, nbin);

   /////saving significance and errors in txt
 	ofstream outFile("analysis/significanceAndErrorsFile.txt");
	if (outFile.is_open()) {
	    outFile << maxSignificance << " " << errorBottomCombined << " " << errorTopCombined << endl;
	    outFile.close();
	} else {
	    cerr << "Unable to open file for writing significance and errors" << endl;
	}
	///////saving significance and errors in txt


   // Clean up files
   signalFile->Close();
   for (TFile* bgFile : backgroundFiles) {
      bgFile->Close();
   }

}
 
 
/*int main( int argc, char** argv )
{
   TMVACutsGA();
}*/

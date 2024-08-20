/////Preselection code for the analysis

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
#include "TF1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMatrixDSymEigen.h"
#include "EventShape/Class/interface/EventShape.h"
#include "vector"
#include <vector>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <cstdlib>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <random>
#include <set>
#include <queue>
#endif

#include "epConstrainHH.h"
#include "LorentzVectorWithErrors.h"

using namespace std;

void GenarateVarsPlots()
{
	cout<<"New macro working"<<endl;
	gStyle->SetOptStat(0);
	
	TFile *fileHH = TFile::Open("analysis/outputTreeSHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	TFile *fileqq = TFile::Open("analysis/outputTreeBqqHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	TFile *filett = TFile::Open("analysis/outputTreeBttHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	TFile *fileZZ = TFile::Open("analysis/outputTreeBZZHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	TFile *fileWW = TFile::Open("analysis/outputTreeBWWHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	
	TTree *treeHH, *treeqq, *treett, *treeZZ, *treeWW;
	
	fileHH->GetObject("TreeSTrain", treeHH);
	fileqq->GetObject("TreeBqqTrain", treeqq);
	filett->GetObject("TreeBttTrain", treett);
	fileZZ->GetObject("TreeBZZTrain", treeZZ);
	fileWW->GetObject("TreeBWWTrain", treeWW);
	
	TH1 *histInvMassB1HH = new TH1F("jetb1HH", "", 142.0, -1.0, 270);
	TH1 *histInvMassB2HH = new TH1F("jetb2HH", "", 142.0, -1.0, 270);
	TH1 *histJet1CosThetaHH = new TH1F("jetcostheta1HH", "", 142.0, -1.1, 1.1);
	TH1 *histSpherHH = new TH1F("sphericity HH", "", 100.0, 0.0, 1.0);
    	histInvMassB1HH->SetFillColor(kBlue);
    	histInvMassB1HH->SetFillStyle(1001); // Solid fill style
    	histInvMassB1HH->SetXTitle("Invariant mass of leading jet pair (GeV)");
    	histInvMassB1HH->SetYTitle("Arbitrary units");
    	histInvMassB2HH->SetFillColor(kBlue);
    	histInvMassB2HH->SetFillStyle(1001); // Solid fill style
    	histInvMassB2HH->SetXTitle("Invariant mass of non-leading jet pair (GeV)");
    	histInvMassB2HH->SetYTitle("Arbitrary units");
	histJet1CosThetaHH->SetFillColor(kBlue);
    	histJet1CosThetaHH->SetFillStyle(1001); // Solid fill style
    	histJet1CosThetaHH->SetXTitle("cos(theta) of leading jet");
    	histJet1CosThetaHH->SetYTitle("Arbitrary units");
   	histSpherHH->SetFillColor(kBlue);
    	histSpherHH->SetFillStyle(1001); // Solid fill style
    	histSpherHH->SetXTitle("Sphericity");
    	histSpherHH->SetYTitle("Arbitrary units");
    	TH1 *histInvMassB1qq = new TH1F("jetb1qq", "inv. mass 1 qq", 142.0, -1.0, 270);
	TH1 *histInvMassB2qq = new TH1F("jetb2qq", "inv. mass 2 qq", 142.0, -1.0, 270);
	TH1 *histJet1CosThetaqq = new TH1F("jetcostheta1qq", "jet cosTheta 1 qq", 142.0, -1.1, 1.1);
	TH1 *histSpherqq = new TH1F("sphericity qq", "", 100.0, 0.0, 1.0);
    	histInvMassB1qq->SetFillColor(kRed+2);
    	histInvMassB1qq->SetFillStyle(3002);
    	histInvMassB1qq->SetLineColor(kRed+2);
    	histInvMassB2qq->SetFillColor(kRed+2);
    	histInvMassB2qq->SetFillStyle(3002);
    	histInvMassB2qq->SetLineColor(kRed+2);
	histJet1CosThetaqq->SetFillColor(kRed+2);
    	histJet1CosThetaqq->SetFillStyle(3002);
    	histJet1CosThetaqq->SetLineColor(kRed+2);
   	histSpherqq->SetFillColor(kRed+2);
   	histSpherqq->SetFillStyle(3002); 
   	histSpherqq->SetLineColor(kRed+2);
   	histSpherqq->SetXTitle("Sphericity");
    	histSpherqq->SetYTitle("Arbitrary units");
    	TH1 *histInvMassB1tt = new TH1F("jetb1tt", "inv. mass 1 tt", 142.0, -1.0, 270);
	TH1 *histInvMassB2tt = new TH1F("jetb2tt", "inv. mass 2 tt", 142.0, -1.0, 270);
	TH1 *histJet1CosThetatt = new TH1F("jetcostheta1tt", "jet cosTheta 1 tt", 142.0, -1.1, 1.1);
	TH1 *histSphertt = new TH1F("sphericity tt", "", 100.0, 0.0, 1.0);
    	histInvMassB1tt->SetFillColorAlpha(kPink+6, 0.9);
    	histInvMassB1tt->SetFillStyle(3002);
    	histInvMassB1tt->SetLineColor(kPink+6);
    	histInvMassB2tt->SetFillColorAlpha(kPink+6, 0.9);
    	histInvMassB2tt->SetFillStyle(3002);
    	histInvMassB2tt->SetLineColor(kPink+6);
	histJet1CosThetatt->SetFillColorAlpha(kPink+6, 0.9);
    	histJet1CosThetatt->SetFillStyle(3002);
    	histJet1CosThetatt->SetLineColor(kPink+6);
   	histSphertt->SetFillColorAlpha(kPink+6, 0.9);
   	histSphertt->SetFillStyle(3002); 
   	histSphertt->SetLineColor(kPink+6);
   	histSphertt->SetXTitle("Sphericity");
    	histSphertt->SetYTitle("Arbitrary units");
    	TH1 *histInvMassB1ZZ = new TH1F("jetb1ZZ", "inv. mass 1 ZZ", 142.0, -1.0, 270);
	TH1 *histInvMassB2ZZ = new TH1F("jetb2ZZ", "inv. mass 2 ZZ", 142.0, -1.0, 270);
	TH1 *histJet1CosThetaZZ = new TH1F("jetcostheta1ZZ", "jet cosTheta 1 ZZ", 142.0, -1.1, 1.1);
	TH1 *histSpherZZ = new TH1F("sphericity ZZ", "", 100.0, 0.0, 1.0);
    	histInvMassB1ZZ->SetFillColorAlpha(kGreen+1, 0.8);
    	histInvMassB1ZZ->SetFillStyle(3002);
    	histInvMassB1ZZ->SetLineColor(kGreen+1);
    	histInvMassB2ZZ->SetFillColorAlpha(kGreen+1, 0.8);
    	histInvMassB2ZZ->SetFillStyle(3002);
    	histInvMassB2ZZ->SetLineColor(kGreen+1);
	histJet1CosThetaZZ->SetFillColorAlpha(kGreen+1, 0.8);
    	histJet1CosThetaZZ->SetFillStyle(3002);
    	histJet1CosThetaZZ->SetLineColor(kGreen+1);
   	histSpherZZ->SetFillColorAlpha(kGreen+1, 0.8);
   	histSpherZZ->SetFillStyle(3002); 
   	histSpherZZ->SetLineColor(kGreen+1);
   	histSpherZZ->SetXTitle("Sphericity");
    	histSpherZZ->SetYTitle("Arbitrary units");
    	TH1 *histInvMassB1WW = new TH1F("jetb1WW", "inv. mass 1 WW", 142.0, -1.0, 270);
	TH1 *histInvMassB2WW = new TH1F("jetb2WW", "inv. mass 2 WW", 142.0, -1.0, 270);
	TH1 *histJet1CosThetaWW = new TH1F("jetcostheta1WW", "jet cosTheta 1 WW", 142.0, -1.1, 1.1);
	TH1 *histSpherWW = new TH1F("sphericity WW", "", 100.0, 0.0, 1.0);
    	histInvMassB1WW->SetFillColorAlpha(kYellow-2, 0.7);
    	histInvMassB1WW->SetFillStyle(3002);
    	histInvMassB1WW->SetLineColor(kYellow-2);
    	histInvMassB2WW->SetFillColorAlpha(kYellow-2, 0.7);
    	histInvMassB2WW->SetFillStyle(3002);
    	histInvMassB2WW->SetLineColor(kYellow-2);
	histJet1CosThetaWW->SetFillColorAlpha(kYellow-2, 0.7);
    	histJet1CosThetaWW->SetFillStyle(3002);
    	histJet1CosThetaWW->SetLineColor(kYellow-2);
   	histSpherWW->SetFillColorAlpha(kYellow-2, 0.7);
   	histSpherWW->SetFillStyle(3002); 
   	histSpherWW->SetLineColor(kYellow-2);
   	histSpherWW->SetXTitle("Sphericity");
    	histSpherWW->SetYTitle("Arbitrary units");
    	
    	
    	
    	////For ZZ detailed study
    	string jetAlgoText = "(durham rtd_cut=10) ";
    	string histExclYmerge12Text = jetAlgoText + "y value when merging from 1 to 2 jets"; 
  	string histExclYmerge23Text = jetAlgoText + "y value when merging from 2 to 3 jets"; 
  	string histExclYmerge34Text = jetAlgoText + "y value when merging from 3 to 4 jets"; 
  	string histExclYmerge45Text = jetAlgoText + "y value when merging from 4 to 5 jets"; 
  	string histExclYmerge56Text = jetAlgoText + "y value when merging from 5 to 6 jets"; 
  	string histMinChiSquaredZZMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZZ mass"; 
  	string histInvMassZZ1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZZ";
  	string histInvMassZZ2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZZ";
  	string histDistanceZ1MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZZ";
  	string histDistanceZ2MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZZ";
  	string histCompareAlgosB1Text = jetAlgoText + "jetPairB1 mass in durham vs. antiKt";
  	string histCompareAlgosB2Text = jetAlgoText + "jetPairB2 mass in durham vs. antiKt";
  	string histInvMassB1CorrectMinChiSquaredHHText = jetAlgoText + "jetPairB1 mass for correct pairing with min chi squared (125)";
  	string histInvMassB1IncorrectMinChiSquaredHHText = jetAlgoText + "jetPairB1 mass for incorrect pairing with min chi squared (125)";
  	string histInvMassB1CorrectMinChiSquaredZZText = jetAlgoText + "jetPairB1 mass for correct pairing with min chi squared (90)";
  	string histInvMassB1IncorrectMinChiSquaredZZText = jetAlgoText + "jetPairB1 mass for incorrect pairing with min chi squared (90)";
  	string histInvMassB1CorrectWithMinDeltaRJetB1Text = jetAlgoText + "jetPairB1 mass for correct pairing with min delta R for jetB1";
  	string histInvMassB1IncorrectWithMinDeltaRJetB1Text = jetAlgoText + "jetPairB1 mass for incorrect pairing with min delta R for jetB1";
  	string histInvMassB1CorrectWithMinDeltaRJetsText = jetAlgoText + "jetPairB1 mass for correct pairing with min delta R for any jet";
  	string histInvMassB1IncorrectWithMinDeltaRJetsText = jetAlgoText + "jetPairB1 mass for incorrect pairing with min delta R for any jet";
  	string histInvMassB2CorrectMinChiSquaredHHText = jetAlgoText + "jetPairB2 mass for correct pairing with min chi squared (125)";
  	string histInvMassB2IncorrectMinChiSquaredHHText = jetAlgoText + "jetPairB2 mass for incorrect pairing with min chi squared (125)";
  	string histInvMassB2CorrectMinChiSquaredZZText = jetAlgoText + "jetPairB2 mass for correct pairing with min chi squared (90)";
  	string histInvMassB2IncorrectMinChiSquaredZZText = jetAlgoText + "jetPairB2 mass for incorrect pairing with min chi squared (90)";
  	string histInvMassB2CorrectWithMinDeltaRJetB1Text = jetAlgoText + "jetPairB2 mass for correct pairing with min delta R for jetB1";
  	string histInvMassB2IncorrectWithMinDeltaRJetB1Text = jetAlgoText + "jetPairB2 mass for incorrect pairing with min delta R for jetB1";
  	string histInvMassB2CorrectWithMinDeltaRJetsText = jetAlgoText + "jetPairB2 mass for correct pairing with min delta R for any jet";
  	string histInvMassB2IncorrectWithMinDeltaRJetsText = jetAlgoText + "jetPairB2 mass for incorrect pairing with min delta R for any jet";
  	string histThetaStarB1Text = jetAlgoText + "theta* for b1";
  	string histThetaStarB2Text = jetAlgoText + "theta* for b2";
  	string histThetaStarB3Text = jetAlgoText + "theta* for b3";
  	string histThetaStarB4Text = jetAlgoText + "theta* for b4";
  	string histThrustText = jetAlgoText + "thrust";
    	TH1 *histExclYmerge12HH = new TH1F("ExclYmerge12HH", histExclYmerge12Text.c_str(), 142.0, -0.1, 1.1);
  	TH1 *histExclYmerge23HH = new TH1F("ExclYmerge23HH,", histExclYmerge23Text.c_str(), 142.0, -0.05, 0.35);
  	TH1 *histExclYmerge34HH = new TH1F("ExclYmerge34HH,", histExclYmerge34Text.c_str(), 142.0, -0.05, 0.20);
  	TH1 *histExclYmerge45HH = new TH1F("ExclYmerge45HH,", histExclYmerge45Text.c_str(), 142.0, -0.005, 0.04);
  	TH1 *histExclYmerge56HH = new TH1F("ExclYmerge56HH,", histExclYmerge56Text.c_str(), 142.0, -0.001, 0.015);
  	TH1 *histMinChiSquaredZZMassHH = new TH1F("MinChiSquaredZZMassHH,", histMinChiSquaredZZMassText.c_str(), 142.0, -100.0, 6000);
  	TH1 *histInvMassZZ1HH = new TH1F("InvMassZZ1HH,", histInvMassZZ1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassZZ2HH = new TH1F("InvMassZZ2HH,", histInvMassZZ2Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histDistanceZ1MinChiSquaredZZMassHH = new TH1F("distanceZ1MinChiSquaredZZMassHH,", histDistanceZ1MinChiSquaredZZMassText.c_str(), 142.0, -5.0, 150);
  	TH1 *histDistanceZ2MinChiSquaredZZMassHH = new TH1F("distanceZ2MinChiSquaredZZMassHH,", histDistanceZ2MinChiSquaredZZMassText.c_str(), 142.0, -5.0, 150);
  	TH2F *histCompareAlgosB1HH = new TH2F("histCompareAlgosB1HH", histCompareAlgosB1Text.c_str(), 142, -5.0, 270.0, 142, -5.0, 270.0);
  	TH2F *histCompareAlgosB2HH = new TH2F("histCompareAlgosB2HH", histCompareAlgosB2Text.c_str(), 142, -5.0, 270.0, 142, -5.0, 270.0);
  	TH1 *histInvMassB1CorrectMinChiSquaredHHHH = new TH1F("InvMassB1CorrectMinChiSquaredHHHH,", histInvMassB1CorrectMinChiSquaredHHText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1IncorrectMinChiSquaredHHHH = new TH1F("InvMassB1IncorrectMinChiSquaredHHHH,", histInvMassB1IncorrectMinChiSquaredHHText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1CorrectMinChiSquaredZZHH = new TH1F("InvMassB1CorrectMinChiSquaredZZHH,", histInvMassB1CorrectMinChiSquaredZZText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1IncorrectMinChiSquaredZZHH = new TH1F("InvMassB1IncorrectMinChiSquaredZZHH,", histInvMassB1IncorrectMinChiSquaredZZText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1CorrectWithMinDeltaRJetB1HH = new TH1F("InvMassB1CorrectWithMinDeltaRJetB1HH,", histInvMassB1CorrectWithMinDeltaRJetB1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1IncorrectWithMinDeltaRJetB1HH = new TH1F("InvMassB1IncorrectWithMinDeltaRJetB1HH,", histInvMassB1IncorrectWithMinDeltaRJetB1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1CorrectWithMinDeltaRJetsHH = new TH1F("InvMassB1CorrectWithMinDeltaRJetsHH,", histInvMassB1CorrectWithMinDeltaRJetsText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1IncorrectWithMinDeltaRJetsHH = new TH1F("InvMassB1IncorrectWithMinDeltaRJetsHH,", histInvMassB1IncorrectWithMinDeltaRJetsText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2CorrectMinChiSquaredHHHH = new TH1F("InvMassB2CorrectMinChiSquaredHHHH,", histInvMassB2CorrectMinChiSquaredHHText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2IncorrectMinChiSquaredHHHH = new TH1F("InvMassB2IncorrectMinChiSquaredHHHH,", histInvMassB2IncorrectMinChiSquaredHHText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2CorrectMinChiSquaredZZHH = new TH1F("InvMassB2CorrectMinChiSquaredZZHH,", histInvMassB2CorrectMinChiSquaredZZText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2IncorrectMinChiSquaredZZHH = new TH1F("InvMassB2IncorrectMinChiSquaredZZHH,", histInvMassB2IncorrectMinChiSquaredZZText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2CorrectWithMinDeltaRJetB1HH = new TH1F("InvMassB2CorrectWithMinDeltaRJetB1HH,", histInvMassB2CorrectWithMinDeltaRJetB1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2IncorrectWithMinDeltaRJetB1HH = new TH1F("InvMassB2IncorrectWithMinDeltaRJetB1HH,", histInvMassB2IncorrectWithMinDeltaRJetB1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2CorrectWithMinDeltaRJetsHH = new TH1F("InvMassB2CorrectWithMinDeltaRJetsHH,", histInvMassB2CorrectWithMinDeltaRJetsText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2IncorrectWithMinDeltaRJetsHH = new TH1F("InvMassB2IncorrectWithMinDeltaRJetsHH,", histInvMassB2IncorrectWithMinDeltaRJetsText.c_str(), 142.0, -1.0, 270);
  	TH1 *histThetaStarB1HH = new TH1F("ThetaStarB1HH,", histThetaStarB1Text.c_str(), 142.0, -0.3, 3.5);
  	TH1 *histThetaStarB2HH = new TH1F("ThetaStarB2HH,", histThetaStarB2Text.c_str(), 142.0, -0.3, 3.5);
  	TH1 *histThetaStarB3HH = new TH1F("ThetaStarB3HH,", histThetaStarB3Text.c_str(), 142.0, -0.3, 3.5);
  	TH1 *histThetaStarB4HH = new TH1F("ThetaStarB4HH,", histThetaStarB4Text.c_str(), 142.0, -0.3, 3.5);
  	TH1 *histThrustHH = new TH1F("ThrustHH,", histThrustText.c_str(), 142.0, -0.1, 1.1);
  	histExclYmerge12HH->SetFillColor(kBlue);
	histExclYmerge12HH->SetFillStyle(1001);
	histExclYmerge23HH->SetFillColor(kBlue);
	histExclYmerge23HH->SetFillStyle(1001);
	histExclYmerge34HH->SetFillColor(kBlue);
	histExclYmerge34HH->SetFillStyle(1001);
	histExclYmerge45HH->SetFillColor(kBlue);
	histExclYmerge45HH->SetFillStyle(1001);
	histExclYmerge56HH->SetFillColor(kBlue);
	histExclYmerge56HH->SetFillStyle(1001);
	histMinChiSquaredZZMassHH->SetFillColor(kBlue);
	histMinChiSquaredZZMassHH->SetFillStyle(1001);
	histInvMassZZ1HH->SetFillColor(kBlue);
	histInvMassZZ1HH->SetFillStyle(1001);
	histInvMassZZ2HH->SetFillColor(kBlue);
	histInvMassZZ2HH->SetFillStyle(1001);
	histDistanceZ1MinChiSquaredZZMassHH->SetFillColor(kBlue);
	histDistanceZ1MinChiSquaredZZMassHH->SetFillStyle(1001);
	histDistanceZ2MinChiSquaredZZMassHH->SetFillColor(kBlue);
	histDistanceZ2MinChiSquaredZZMassHH->SetFillStyle(1001);
	histCompareAlgosB1HH->SetFillColor(kBlue);
	histCompareAlgosB1HH->SetLineColor(kBlue);
	histCompareAlgosB2HH->SetFillColor(kBlue);
	histCompareAlgosB2HH->SetLineColor(kBlue);
	histCompareAlgosB1HH->GetXaxis()->SetTitle("invMassB1 with Durham");
   	histCompareAlgosB1HH->GetYaxis()->SetTitle("invMassB1 with AntiKt");
   	histCompareAlgosB2HH->GetXaxis()->SetTitle("invMassB2 with Durham");
   	histCompareAlgosB2HH->GetYaxis()->SetTitle("invMassB2 with AntiKt");
   	histInvMassB1CorrectMinChiSquaredHHHH->SetFillColor(kBlue);
	histInvMassB1CorrectMinChiSquaredHHHH->SetFillStyle(1001);
	histInvMassB1IncorrectMinChiSquaredHHHH->SetFillColor(kBlue);
	histInvMassB1IncorrectMinChiSquaredHHHH->SetFillStyle(1001);
	histInvMassB1CorrectMinChiSquaredZZHH->SetFillColor(kBlue);
	histInvMassB1CorrectMinChiSquaredZZHH->SetFillStyle(1001);
	histInvMassB1IncorrectMinChiSquaredZZHH->SetFillColor(kBlue);
	histInvMassB1IncorrectMinChiSquaredZZHH->SetFillStyle(1001);
	histInvMassB1CorrectWithMinDeltaRJetB1HH->SetFillColor(kBlue);
	histInvMassB1CorrectWithMinDeltaRJetB1HH->SetFillStyle(1001);
	histInvMassB1IncorrectWithMinDeltaRJetB1HH->SetFillColor(kBlue);
	histInvMassB1IncorrectWithMinDeltaRJetB1HH->SetFillStyle(1001);
	histInvMassB1CorrectWithMinDeltaRJetsHH->SetFillColor(kBlue);
	histInvMassB1CorrectWithMinDeltaRJetsHH->SetFillStyle(1001);
	histInvMassB1IncorrectWithMinDeltaRJetsHH->SetFillColor(kBlue);
	histInvMassB1IncorrectWithMinDeltaRJetsHH->SetFillStyle(1001);
	histInvMassB2CorrectMinChiSquaredHHHH->SetFillColor(kBlue);
	histInvMassB2CorrectMinChiSquaredHHHH->SetFillStyle(1001);
	histInvMassB2IncorrectMinChiSquaredHHHH->SetFillColor(kBlue);
	histInvMassB2IncorrectMinChiSquaredHHHH->SetFillStyle(1001);
	histInvMassB2CorrectMinChiSquaredZZHH->SetFillColor(kBlue);
	histInvMassB2CorrectMinChiSquaredZZHH->SetFillStyle(1001);
	histInvMassB2IncorrectMinChiSquaredZZHH->SetFillColor(kBlue);
	histInvMassB2IncorrectMinChiSquaredZZHH->SetFillStyle(1001);
	histInvMassB2CorrectWithMinDeltaRJetB1HH->SetFillColor(kBlue);
	histInvMassB2CorrectWithMinDeltaRJetB1HH->SetFillStyle(1001);
	histInvMassB2IncorrectWithMinDeltaRJetB1HH->SetFillColor(kBlue);
	histInvMassB2IncorrectWithMinDeltaRJetB1HH->SetFillStyle(1001);
	histInvMassB2CorrectWithMinDeltaRJetsHH->SetFillColor(kBlue);
	histInvMassB2CorrectWithMinDeltaRJetsHH->SetFillStyle(1001);
	histInvMassB2IncorrectWithMinDeltaRJetsHH->SetFillColor(kBlue);
	histInvMassB2IncorrectWithMinDeltaRJetsHH->SetFillStyle(1001);
	histThetaStarB1HH->SetFillColor(kBlue);
	histThetaStarB1HH->SetFillStyle(1001);
	histThetaStarB2HH->SetFillColor(kBlue);
	histThetaStarB2HH->SetFillStyle(1001);
	histThetaStarB3HH->SetFillColor(kBlue);
	histThetaStarB3HH->SetFillStyle(1001);
	histThetaStarB4HH->SetFillColor(kBlue);
	histThetaStarB4HH->SetFillStyle(1001);
	histThrustHH->SetFillColor(kBlue);
	histThrustHH->SetFillStyle(1001);

  	TH1 *histExclYmerge12ZZ = new TH1F("ExclYmerge12ZZ", histExclYmerge12Text.c_str(), 142.0, -0.1, 1.1);
  	TH1 *histExclYmerge23ZZ = new TH1F("ExclYmerge23ZZ,", histExclYmerge23Text.c_str(), 142.0, -0.05, 0.35);
  	TH1 *histExclYmerge34ZZ = new TH1F("ExclYmerge34ZZ,", histExclYmerge34Text.c_str(), 142.0, -0.05, 0.20);
  	TH1 *histExclYmerge45ZZ = new TH1F("ExclYmerge45ZZ,", histExclYmerge45Text.c_str(), 142.0, -0.005, 0.04);
  	TH1 *histExclYmerge56ZZ = new TH1F("ExclYmerge56ZZ,", histExclYmerge56Text.c_str(), 142.0, -0.001, 0.015);
  	TH1 *histMinChiSquaredZZMassZZ = new TH1F("MinChiSquaredZZMassZZ,", histMinChiSquaredZZMassText.c_str(), 142.0, -100.0, 6000);
  	TH1 *histInvMassZZ1ZZ = new TH1F("InvMassZZ1ZZ,", histInvMassZZ1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassZZ2ZZ = new TH1F("InvMassZZ2ZZ,", histInvMassZZ2Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histDistanceZ1MinChiSquaredZZMassZZ = new TH1F("distanceZ1MinChiSquaredZZMassZZ,", histDistanceZ1MinChiSquaredZZMassText.c_str(), 142.0, -5.0, 150);
  	TH1 *histDistanceZ2MinChiSquaredZZMassZZ = new TH1F("distanceZ2MinChiSquaredZZMassZZ,", histDistanceZ2MinChiSquaredZZMassText.c_str(), 142.0, -5.0, 150);
  	TH2F *histCompareAlgosB1ZZ = new TH2F("histCompareAlgosB1ZZ", histCompareAlgosB1Text.c_str(), 142, -5.0, 270.0, 142, -5.0, 270.0);
  	TH2F *histCompareAlgosB2ZZ = new TH2F("histCompareAlgosB2ZZ", histCompareAlgosB2Text.c_str(), 142, -5.0, 270.0, 142, -5.0, 270.0);
  	TH1 *histInvMassB1CorrectMinChiSquaredHHZZ = new TH1F("InvMassB1CorrectMinChiSquaredHHZZ,", histInvMassB1CorrectMinChiSquaredHHText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1IncorrectMinChiSquaredHHZZ = new TH1F("InvMassB1IncorrectMinChiSquaredHHZZ,", histInvMassB1IncorrectMinChiSquaredHHText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1CorrectMinChiSquaredZZZZ = new TH1F("InvMassB1CorrectMinChiSquaredZZZZ,", histInvMassB1CorrectMinChiSquaredZZText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1IncorrectMinChiSquaredZZZZ = new TH1F("InvMassB1IncorrectMinChiSquaredZZZZ,", histInvMassB1IncorrectMinChiSquaredZZText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1CorrectWithMinDeltaRJetB1ZZ = new TH1F("InvMassB1CorrectWithMinDeltaRJetB1ZZ,", histInvMassB1CorrectWithMinDeltaRJetB1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1IncorrectWithMinDeltaRJetB1ZZ = new TH1F("InvMassB1IncorrectWithMinDeltaRJetB1ZZ,", histInvMassB1IncorrectWithMinDeltaRJetB1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1CorrectWithMinDeltaRJetsZZ = new TH1F("InvMassB1CorrectWithMinDeltaRJetsZZ,", histInvMassB1CorrectWithMinDeltaRJetsText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB1IncorrectWithMinDeltaRJetsZZ = new TH1F("InvMassB1IncorrectWithMinDeltaRJetsZZ,", histInvMassB1IncorrectWithMinDeltaRJetsText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2CorrectMinChiSquaredHHZZ = new TH1F("InvMassB2CorrectMinChiSquaredHHZZ,", histInvMassB2CorrectMinChiSquaredHHText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2IncorrectMinChiSquaredHHZZ = new TH1F("InvMassB2IncorrectMinChiSquaredHHZZ,", histInvMassB2IncorrectMinChiSquaredHHText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2CorrectMinChiSquaredZZZZ = new TH1F("InvMassB2CorrectMinChiSquaredZZZZ,", histInvMassB2CorrectMinChiSquaredZZText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2IncorrectMinChiSquaredZZZZ = new TH1F("InvMassB2IncorrectMinChiSquaredZZZZ,", histInvMassB2IncorrectMinChiSquaredZZText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2CorrectWithMinDeltaRJetB1ZZ = new TH1F("InvMassB2CorrectWithMinDeltaRJetB1ZZ,", histInvMassB2CorrectWithMinDeltaRJetB1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2IncorrectWithMinDeltaRJetB1ZZ = new TH1F("InvMassB2IncorrectWithMinDeltaRJetB1ZZ,", histInvMassB2IncorrectWithMinDeltaRJetB1Text.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2CorrectWithMinDeltaRJetsZZ = new TH1F("InvMassB2CorrectWithMinDeltaRJetsZZ,", histInvMassB2CorrectWithMinDeltaRJetsText.c_str(), 142.0, -1.0, 270);
  	TH1 *histInvMassB2IncorrectWithMinDeltaRJetsZZ = new TH1F("InvMassB2IncorrectWithMinDeltaRJetsZZ,", histInvMassB2IncorrectWithMinDeltaRJetsText.c_str(), 142.0, -1.0, 270);
  	TH1 *histThetaStarB1ZZ = new TH1F("ThetaStarB1ZZ,", histThetaStarB1Text.c_str(), 142.0, -0.3, 3.5);
  	TH1 *histThetaStarB2ZZ = new TH1F("ThetaStarB2ZZ,", histThetaStarB2Text.c_str(), 142.0, -0.3, 3.5);
  	TH1 *histThetaStarB3ZZ = new TH1F("ThetaStarB3ZZ,", histThetaStarB3Text.c_str(), 142.0, -0.3, 3.5);
  	TH1 *histThetaStarB4ZZ = new TH1F("ThetaStarB4ZZ,", histThetaStarB4Text.c_str(), 142.0, -0.3, 3.5);
  	TH1 *histThrustZZ = new TH1F("ThrustZZ,", histThrustText.c_str(), 142.0, -0.1, 1.1);
  	histExclYmerge12ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histExclYmerge12ZZ->SetFillStyle(3002);
	histExclYmerge12ZZ->SetLineColor(kGreen+1);
	histExclYmerge23ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histExclYmerge23ZZ->SetFillStyle(3002);
	histExclYmerge23ZZ->SetLineColor(kGreen+1);
	histExclYmerge34ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histExclYmerge34ZZ->SetFillStyle(3002);
	histExclYmerge34ZZ->SetLineColor(kGreen+1);
	histExclYmerge45ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histExclYmerge45ZZ->SetFillStyle(3002);
	histExclYmerge45ZZ->SetLineColor(kGreen+1);
	histExclYmerge56ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histExclYmerge56ZZ->SetFillStyle(3002);
	histExclYmerge56ZZ->SetLineColor(kGreen+1);
	histMinChiSquaredZZMassZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histMinChiSquaredZZMassZZ->SetFillStyle(3002);
	histMinChiSquaredZZMassZZ->SetLineColor(kGreen+1);
	histInvMassZZ1ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassZZ1ZZ->SetFillStyle(3002);
	histInvMassZZ1ZZ->SetLineColor(kGreen+1);
	histInvMassZZ2ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassZZ2ZZ->SetFillStyle(3002);
	histInvMassZZ2ZZ->SetLineColor(kGreen+1);
	histDistanceZ1MinChiSquaredZZMassZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histDistanceZ1MinChiSquaredZZMassZZ->SetFillStyle(3002);
	histDistanceZ1MinChiSquaredZZMassZZ->SetLineColor(kGreen+1);
	histDistanceZ2MinChiSquaredZZMassZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histDistanceZ2MinChiSquaredZZMassZZ->SetFillStyle(3002);
	histDistanceZ2MinChiSquaredZZMassZZ->SetLineColor(kGreen+1);
	histCompareAlgosB1ZZ->SetFillColor(kGreen+1);
	histCompareAlgosB1ZZ->SetLineColor(kGreen+1);
	histCompareAlgosB1ZZ->SetMarkerColor(kGreen+1);
	histCompareAlgosB2ZZ->SetFillColor(kGreen+1);
	histCompareAlgosB2ZZ->SetLineColor(kGreen+1);
	histCompareAlgosB2ZZ->SetMarkerColor(kGreen+1);
	histCompareAlgosB1ZZ->GetXaxis()->SetTitle("invMassB1 with Durham");
   	histCompareAlgosB1ZZ->GetYaxis()->SetTitle("invMassB1 with AntiKt");
   	histCompareAlgosB2ZZ->GetXaxis()->SetTitle("invMassB2 with Durham");
   	histCompareAlgosB2ZZ->GetYaxis()->SetTitle("invMassB2 with AntiKt");
   	histInvMassB1CorrectMinChiSquaredHHZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB1CorrectMinChiSquaredHHZZ->SetFillStyle(3002);
	histInvMassB1CorrectMinChiSquaredHHZZ->SetLineColor(kGreen+1);
	histInvMassB1IncorrectMinChiSquaredHHZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB1IncorrectMinChiSquaredHHZZ->SetFillStyle(3002);
	histInvMassB1IncorrectMinChiSquaredHHZZ->SetLineColor(kGreen+1);
	histInvMassB1CorrectMinChiSquaredZZZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB1CorrectMinChiSquaredZZZZ->SetFillStyle(3002);
	histInvMassB1CorrectMinChiSquaredZZZZ->SetLineColor(kGreen+1);
	histInvMassB1IncorrectMinChiSquaredZZZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB1IncorrectMinChiSquaredZZZZ->SetFillStyle(3002);
	histInvMassB1IncorrectMinChiSquaredZZZZ->SetLineColor(kGreen+1);
	histInvMassB1CorrectWithMinDeltaRJetB1ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB1CorrectWithMinDeltaRJetB1ZZ->SetFillStyle(3002);
	histInvMassB1CorrectWithMinDeltaRJetB1ZZ->SetLineColor(kGreen+1);
	histInvMassB1IncorrectWithMinDeltaRJetB1ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB1IncorrectWithMinDeltaRJetB1ZZ->SetFillStyle(3002);
	histInvMassB1IncorrectWithMinDeltaRJetB1ZZ->SetLineColor(kGreen+1);
	histInvMassB1CorrectWithMinDeltaRJetsZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB1CorrectWithMinDeltaRJetsZZ->SetFillStyle(3002);
	histInvMassB1CorrectWithMinDeltaRJetsZZ->SetLineColor(kGreen+1);
	histInvMassB1IncorrectWithMinDeltaRJetsZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB1IncorrectWithMinDeltaRJetsZZ->SetFillStyle(3002);
	histInvMassB1IncorrectWithMinDeltaRJetsZZ->SetLineColor(kGreen+1);
	histInvMassB2CorrectMinChiSquaredHHZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB2CorrectMinChiSquaredHHZZ->SetFillStyle(3002);
	histInvMassB2CorrectMinChiSquaredHHZZ->SetLineColor(kGreen+1);
	histInvMassB2IncorrectMinChiSquaredHHZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB2IncorrectMinChiSquaredHHZZ->SetFillStyle(3002);
	histInvMassB2IncorrectMinChiSquaredHHZZ->SetLineColor(kGreen+1);
	histInvMassB2CorrectMinChiSquaredZZZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB2CorrectMinChiSquaredZZZZ->SetFillStyle(3002);
	histInvMassB2CorrectMinChiSquaredZZZZ->SetLineColor(kGreen+1);
	histInvMassB2IncorrectMinChiSquaredZZZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB2IncorrectMinChiSquaredZZZZ->SetFillStyle(3002);
	histInvMassB2IncorrectMinChiSquaredZZZZ->SetLineColor(kGreen+1);
	histInvMassB2CorrectWithMinDeltaRJetB1ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB2CorrectWithMinDeltaRJetB1ZZ->SetFillStyle(3002);
	histInvMassB2CorrectWithMinDeltaRJetB1ZZ->SetLineColor(kGreen+1);
	histInvMassB2IncorrectWithMinDeltaRJetB1ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB2IncorrectWithMinDeltaRJetB1ZZ->SetFillStyle(3002);
	histInvMassB2IncorrectWithMinDeltaRJetB1ZZ->SetLineColor(kGreen+1);
	histInvMassB2CorrectWithMinDeltaRJetsZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB2CorrectWithMinDeltaRJetsZZ->SetFillStyle(3002);
	histInvMassB2CorrectWithMinDeltaRJetsZZ->SetLineColor(kGreen+1);
	histInvMassB2IncorrectWithMinDeltaRJetsZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histInvMassB2IncorrectWithMinDeltaRJetsZZ->SetFillStyle(3002);
	histInvMassB2IncorrectWithMinDeltaRJetsZZ->SetLineColor(kGreen+1);
	histThetaStarB1ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histThetaStarB1ZZ->SetFillStyle(3002);
	histThetaStarB1ZZ->SetLineColor(kGreen+1);
	histThetaStarB2ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histThetaStarB2ZZ->SetFillStyle(3002);
	histThetaStarB2ZZ->SetLineColor(kGreen+1);
	histThetaStarB3ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histThetaStarB3ZZ->SetFillStyle(3002);
	histThetaStarB3ZZ->SetLineColor(kGreen+1);
	histThetaStarB4ZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histThetaStarB4ZZ->SetFillStyle(3002);
	histThetaStarB4ZZ->SetLineColor(kGreen+1);
	histThrustZZ->SetFillColorAlpha(kGreen+1, 0.8);
	histThrustZZ->SetFillStyle(3002);
	histThrustZZ->SetLineColor(kGreen+1);
    	/////For ZZ detailed study
    	
    	
	
	 // Set the branch address
	 Float_t invMassB1,invMassB2, cosThetaB1, sphericity, exclYmerge12, exclYmerge23, exclYmerge34, exclYmerge45, exclYmerge56, invMassZZ1, invMassZZ2, distanceZ1MinChiSquaredZZMass, distanceZ2MinChiSquaredZZMass, invMassB1AntiKt, invMassB2AntiKt, invMassB1CorrectMinChiSquaredHH, invMassB1IncorrectMinChiSquaredHH, invMassB1CorrectMinChiSquaredZZ, invMassB1IncorrectMinChiSquaredZZ, invMassB1CorrectWithMinDeltaRJetB1, invMassB1IncorrectWithMinDeltaRJetB1, invMassB1CorrectWithMinDeltaRJets, invMassB1IncorrectWithMinDeltaRJets, invMassB2CorrectMinChiSquaredHH, invMassB2IncorrectMinChiSquaredHH, invMassB2CorrectMinChiSquaredZZ, invMassB2IncorrectMinChiSquaredZZ, invMassB2CorrectWithMinDeltaRJetB1, invMassB2IncorrectWithMinDeltaRJetB1, invMassB2CorrectWithMinDeltaRJets, invMassB2IncorrectWithMinDeltaRJets, thetaStarB1, thetaStarB2, thetaStarB3, thetaStarB4, thrust;
	 treeHH->SetBranchAddress("invMassB1", &invMassB1);
	 treeHH->SetBranchAddress("invMassB2", &invMassB2);
	 treeHH->SetBranchAddress("invMassB1AntiKt", &invMassB1AntiKt);
	 treeHH->SetBranchAddress("invMassB2AntiKt", &invMassB2AntiKt);
	 treeHH->SetBranchAddress("cosThetaB1", &cosThetaB1);
	 treeHH->SetBranchAddress("sphericity", &sphericity);
	 treeHH->SetBranchAddress("exclYmerge12", &exclYmerge12);
	 treeHH->SetBranchAddress("exclYmerge23", &exclYmerge23);
	 treeHH->SetBranchAddress("exclYmerge34", &exclYmerge34);
	 treeHH->SetBranchAddress("exclYmerge45", &exclYmerge45);
	 treeHH->SetBranchAddress("exclYmerge56", &exclYmerge56);
	 treeHH->SetBranchAddress("invMassZZ1", &invMassZZ1);
	 treeHH->SetBranchAddress("invMassZZ2", &invMassZZ2);
	 treeHH->SetBranchAddress("distanceZ1MinChiSquaredZZMass", &distanceZ1MinChiSquaredZZMass);
	 treeHH->SetBranchAddress("distanceZ2MinChiSquaredZZMass", &distanceZ2MinChiSquaredZZMass);
	 treeHH->SetBranchAddress("invMassB1CorrectMinChiSquaredHH", &invMassB1CorrectMinChiSquaredHH);
	 treeHH->SetBranchAddress("invMassB1IncorrectMinChiSquaredHH", &invMassB1IncorrectMinChiSquaredHH);
	 treeHH->SetBranchAddress("invMassB1CorrectMinChiSquaredZZ", &invMassB1CorrectMinChiSquaredZZ);
	 treeHH->SetBranchAddress("invMassB1IncorrectMinChiSquaredZZ", &invMassB1IncorrectMinChiSquaredZZ);
	 treeHH->SetBranchAddress("invMassB1CorrectWithMinDeltaRJetB1", &invMassB1CorrectWithMinDeltaRJetB1);
	 treeHH->SetBranchAddress("invMassB1IncorrectWithMinDeltaRJetB1", &invMassB1IncorrectWithMinDeltaRJetB1);
	 treeHH->SetBranchAddress("invMassB1CorrectWithMinDeltaRJets", &invMassB1CorrectWithMinDeltaRJets);
	 treeHH->SetBranchAddress("invMassB1IncorrectWithMinDeltaRJets", &invMassB1IncorrectWithMinDeltaRJets);
	 treeHH->SetBranchAddress("invMassB2CorrectMinChiSquaredHH", &invMassB2CorrectMinChiSquaredHH);
	 treeHH->SetBranchAddress("invMassB2IncorrectMinChiSquaredHH", &invMassB2IncorrectMinChiSquaredHH);
	 treeHH->SetBranchAddress("invMassB2CorrectMinChiSquaredZZ", &invMassB2CorrectMinChiSquaredZZ);
	 treeHH->SetBranchAddress("invMassB2IncorrectMinChiSquaredZZ", &invMassB2IncorrectMinChiSquaredZZ);
	 treeHH->SetBranchAddress("invMassB2CorrectWithMinDeltaRJetB1", &invMassB2CorrectWithMinDeltaRJetB1);
	 treeHH->SetBranchAddress("invMassB2IncorrectWithMinDeltaRJetB1", &invMassB2IncorrectWithMinDeltaRJetB1);
	 treeHH->SetBranchAddress("invMassB2CorrectWithMinDeltaRJets", &invMassB2CorrectWithMinDeltaRJets);
	 treeHH->SetBranchAddress("invMassB2IncorrectWithMinDeltaRJets", &invMassB2IncorrectWithMinDeltaRJets);
	 treeHH->SetBranchAddress("thetaStarB1", &thetaStarB1);
	 treeHH->SetBranchAddress("thetaStarB2", &thetaStarB2);
	 treeHH->SetBranchAddress("thetaStarB3", &thetaStarB3);
	 treeHH->SetBranchAddress("thetaStarB4", &thetaStarB4);
	 treeHH->SetBranchAddress("thrust", &thrust);
	 
	 // Loop over the entries in the tree and fill the histogram
	 Long64_t nEntries = treeHH->GetEntries();
	 for (Long64_t i = 0; i < nEntries; ++i) {
	 	treeHH->GetEntry(i);
	 	//if(distanceZ1MinChiSquaredZZMass > 30 && distanceZ2MinChiSquaredZZMass > 30)
	 	if(1==1)
	 	{
		 	histInvMassB1HH->Fill(invMassB1);
		 	histInvMassB2HH->Fill(invMassB2);
		 	histJet1CosThetaHH->Fill(cosThetaB1);
		 	histSpherHH->Fill(sphericity);
		 	histExclYmerge12HH->Fill(exclYmerge12);
			histExclYmerge23HH->Fill(exclYmerge23);
			histExclYmerge34HH->Fill(exclYmerge34);
			histExclYmerge45HH->Fill(exclYmerge45);
			histExclYmerge56HH->Fill(exclYmerge56);
			histInvMassZZ1HH->Fill(invMassZZ1);
			histInvMassZZ2HH->Fill(invMassZZ2);
			histDistanceZ1MinChiSquaredZZMassHH->Fill(distanceZ1MinChiSquaredZZMass);
			histDistanceZ2MinChiSquaredZZMassHH->Fill(distanceZ2MinChiSquaredZZMass);
			histCompareAlgosB1HH->Fill(invMassB1, invMassB1AntiKt);
			histCompareAlgosB2HH->Fill(invMassB2, invMassB2AntiKt); 
			if(invMassB1CorrectMinChiSquaredHH != -999) histInvMassB1CorrectMinChiSquaredHHHH->Fill(invMassB1CorrectMinChiSquaredHH);
			if(invMassB1IncorrectMinChiSquaredHH != -999) histInvMassB1IncorrectMinChiSquaredHHHH->Fill(invMassB1IncorrectMinChiSquaredHH);
			if(invMassB1CorrectMinChiSquaredZZ != -999) histInvMassB1CorrectMinChiSquaredZZHH->Fill(invMassB1CorrectMinChiSquaredZZ);
			if(invMassB1IncorrectMinChiSquaredZZ != -999) histInvMassB1IncorrectMinChiSquaredZZHH->Fill(invMassB1IncorrectMinChiSquaredZZ);
			if(invMassB1CorrectWithMinDeltaRJetB1 != -999) histInvMassB1CorrectWithMinDeltaRJetB1HH->Fill(invMassB1CorrectWithMinDeltaRJetB1);
			if(invMassB1IncorrectWithMinDeltaRJetB1 != -999) histInvMassB1IncorrectWithMinDeltaRJetB1HH->Fill(invMassB1IncorrectWithMinDeltaRJetB1);
			if(invMassB1CorrectWithMinDeltaRJets != -999) histInvMassB1CorrectWithMinDeltaRJetsHH->Fill(invMassB1CorrectWithMinDeltaRJets);
			if(invMassB1IncorrectWithMinDeltaRJets != -999) histInvMassB1IncorrectWithMinDeltaRJetsHH->Fill(invMassB1IncorrectWithMinDeltaRJets);
			if(invMassB2CorrectMinChiSquaredHH != -999) histInvMassB2CorrectMinChiSquaredHHHH->Fill(invMassB2CorrectMinChiSquaredHH);
			if(invMassB2IncorrectMinChiSquaredHH != -999) histInvMassB2IncorrectMinChiSquaredHHHH->Fill(invMassB2IncorrectMinChiSquaredHH);
			if(invMassB2CorrectMinChiSquaredZZ != -999) histInvMassB2CorrectMinChiSquaredZZHH->Fill(invMassB2CorrectMinChiSquaredZZ);
			if(invMassB2IncorrectMinChiSquaredZZ != -999) histInvMassB2IncorrectMinChiSquaredZZHH->Fill(invMassB2IncorrectMinChiSquaredZZ);
			if(invMassB2CorrectWithMinDeltaRJetB1 != -999) histInvMassB2CorrectWithMinDeltaRJetB1HH->Fill(invMassB2CorrectWithMinDeltaRJetB1);
			if(invMassB2IncorrectWithMinDeltaRJetB1 != -999) histInvMassB2IncorrectWithMinDeltaRJetB1HH->Fill(invMassB2IncorrectWithMinDeltaRJetB1);
			if(invMassB2CorrectWithMinDeltaRJets != -999) histInvMassB2CorrectWithMinDeltaRJetsHH->Fill(invMassB2CorrectWithMinDeltaRJets);
			if(invMassB2IncorrectWithMinDeltaRJets != -999) histInvMassB2IncorrectWithMinDeltaRJetsHH->Fill(invMassB2IncorrectWithMinDeltaRJets);
			histThetaStarB1HH->Fill(thetaStarB1);
			histThetaStarB2HH->Fill(thetaStarB2);
			histThetaStarB3HH->Fill(thetaStarB3);
			histThetaStarB4HH->Fill(thetaStarB4);
			histThrustHH->Fill(thrust);
			
		}
	 }
	 
	 treeqq->SetBranchAddress("invMassB1", &invMassB1);
	 treeqq->SetBranchAddress("invMassB2", &invMassB2);
	 treeqq->SetBranchAddress("cosThetaB1", &cosThetaB1);
	 treeqq->SetBranchAddress("sphericity", &sphericity);
	 nEntries = treeqq->GetEntries();
	 for (Long64_t i = 0; i < nEntries; ++i) {
	 	treeqq->GetEntry(i);
	 	histInvMassB1qq->Fill(invMassB1);
	 	histInvMassB2qq->Fill(invMassB2);
	 	histJet1CosThetaqq->Fill(cosThetaB1);
	 	histSpherqq->Fill(sphericity);
	 }
	 
	 treett->SetBranchAddress("invMassB1", &invMassB1);
	 treett->SetBranchAddress("invMassB2", &invMassB2);
	 treett->SetBranchAddress("cosThetaB1", &cosThetaB1);
	 treett->SetBranchAddress("sphericity", &sphericity);
	 nEntries = treett->GetEntries();
	 for (Long64_t i = 0; i < nEntries; ++i) {
	 	treett->GetEntry(i);
	 	histInvMassB1tt->Fill(invMassB1);
	 	histInvMassB2tt->Fill(invMassB2);
	 	histJet1CosThetatt->Fill(cosThetaB1);
	 	histSphertt->Fill(sphericity);
	 }
	 
	 treeZZ->SetBranchAddress("invMassB1", &invMassB1);
	 treeZZ->SetBranchAddress("invMassB2", &invMassB2);
	 treeZZ->SetBranchAddress("invMassB1AntiKt", &invMassB1AntiKt);
	 treeZZ->SetBranchAddress("invMassB2AntiKt", &invMassB2AntiKt);
	 treeZZ->SetBranchAddress("cosThetaB1", &cosThetaB1);
	 treeZZ->SetBranchAddress("sphericity", &sphericity);
	 treeZZ->SetBranchAddress("exclYmerge12", &exclYmerge12);
	 treeZZ->SetBranchAddress("exclYmerge23", &exclYmerge23);
	 treeZZ->SetBranchAddress("exclYmerge34", &exclYmerge34);
	 treeZZ->SetBranchAddress("exclYmerge45", &exclYmerge45);
	 treeZZ->SetBranchAddress("exclYmerge56", &exclYmerge56);
	 treeZZ->SetBranchAddress("invMassZZ1", &invMassZZ1);
	 treeZZ->SetBranchAddress("invMassZZ2", &invMassZZ2);
	 treeZZ->SetBranchAddress("distanceZ1MinChiSquaredZZMass", &distanceZ1MinChiSquaredZZMass);
	 treeZZ->SetBranchAddress("distanceZ2MinChiSquaredZZMass", &distanceZ2MinChiSquaredZZMass);
	 treeZZ->SetBranchAddress("invMassB1CorrectMinChiSquaredHH", &invMassB1CorrectMinChiSquaredHH);
	 treeZZ->SetBranchAddress("invMassB1IncorrectMinChiSquaredHH", &invMassB1IncorrectMinChiSquaredHH);
	 treeZZ->SetBranchAddress("invMassB1CorrectMinChiSquaredZZ", &invMassB1CorrectMinChiSquaredZZ);
	 treeZZ->SetBranchAddress("invMassB1IncorrectMinChiSquaredZZ", &invMassB1IncorrectMinChiSquaredZZ);
	 treeZZ->SetBranchAddress("invMassB1CorrectWithMinDeltaRJetB1", &invMassB1CorrectWithMinDeltaRJetB1);
	 treeZZ->SetBranchAddress("invMassB1IncorrectWithMinDeltaRJetB1", &invMassB1IncorrectWithMinDeltaRJetB1);
	 treeZZ->SetBranchAddress("invMassB1CorrectWithMinDeltaRJets", &invMassB1CorrectWithMinDeltaRJets);
	 treeZZ->SetBranchAddress("invMassB1IncorrectWithMinDeltaRJets", &invMassB1IncorrectWithMinDeltaRJets);
	 treeZZ->SetBranchAddress("invMassB2CorrectMinChiSquaredHH", &invMassB2CorrectMinChiSquaredHH);
	 treeZZ->SetBranchAddress("invMassB2IncorrectMinChiSquaredHH", &invMassB2IncorrectMinChiSquaredHH);
	 treeZZ->SetBranchAddress("invMassB2CorrectMinChiSquaredZZ", &invMassB2CorrectMinChiSquaredZZ);
	 treeZZ->SetBranchAddress("invMassB2IncorrectMinChiSquaredZZ", &invMassB2IncorrectMinChiSquaredZZ);
	 treeZZ->SetBranchAddress("invMassB2CorrectWithMinDeltaRJetB1", &invMassB2CorrectWithMinDeltaRJetB1);
	 treeZZ->SetBranchAddress("invMassB2IncorrectWithMinDeltaRJetB1", &invMassB2IncorrectWithMinDeltaRJetB1);
	 treeZZ->SetBranchAddress("invMassB2CorrectWithMinDeltaRJets", &invMassB2CorrectWithMinDeltaRJets);
	 treeZZ->SetBranchAddress("invMassB2IncorrectWithMinDeltaRJets", &invMassB2IncorrectWithMinDeltaRJets);
	 treeZZ->SetBranchAddress("thetaStarB1", &thetaStarB1);
	 treeZZ->SetBranchAddress("thetaStarB2", &thetaStarB2);
	 treeZZ->SetBranchAddress("thetaStarB3", &thetaStarB3);
	 treeZZ->SetBranchAddress("thetaStarB4", &thetaStarB4);
	 treeZZ->SetBranchAddress("thrust", &thrust);
	 nEntries = treeZZ->GetEntries();
	 for (Long64_t i = 0; i < nEntries; ++i) {
	 	treeZZ->GetEntry(i);
	 	//if(distanceZ1MinChiSquaredZZMass > 30 && distanceZ2MinChiSquaredZZMass > 30)
	 	if(1==1)
	 	{
		 	histInvMassB1ZZ->Fill(invMassB1);
		 	histInvMassB2ZZ->Fill(invMassB2);
		 	histJet1CosThetaZZ->Fill(cosThetaB1);
		 	histSpherZZ->Fill(sphericity);
		 	histExclYmerge12ZZ->Fill(exclYmerge12);
			histExclYmerge23ZZ->Fill(exclYmerge23);
			histExclYmerge34ZZ->Fill(exclYmerge34);
			histExclYmerge45ZZ->Fill(exclYmerge45);
			histExclYmerge56ZZ->Fill(exclYmerge56);
			histInvMassZZ1ZZ->Fill(invMassZZ1);
			histInvMassZZ2ZZ->Fill(invMassZZ2);
			histDistanceZ1MinChiSquaredZZMassZZ->Fill(distanceZ1MinChiSquaredZZMass);
			histDistanceZ2MinChiSquaredZZMassZZ->Fill(distanceZ2MinChiSquaredZZMass);
			histCompareAlgosB1ZZ->Fill(invMassB1, invMassB1AntiKt);
			histCompareAlgosB2ZZ->Fill(invMassB2, invMassB2AntiKt); 
			if(invMassB1CorrectMinChiSquaredHH != -999) histInvMassB1CorrectMinChiSquaredHHZZ->Fill(invMassB1CorrectMinChiSquaredHH);
			if(invMassB1IncorrectMinChiSquaredHH != -999) histInvMassB1IncorrectMinChiSquaredHHZZ->Fill(invMassB1IncorrectMinChiSquaredHH);
			if(invMassB1CorrectMinChiSquaredZZ != -999) histInvMassB1CorrectMinChiSquaredZZZZ->Fill(invMassB1CorrectMinChiSquaredZZ);
			if(invMassB1IncorrectMinChiSquaredZZ != -999) histInvMassB1IncorrectMinChiSquaredZZZZ->Fill(invMassB1IncorrectMinChiSquaredZZ);
			if(invMassB1CorrectWithMinDeltaRJetB1 != -999) histInvMassB1CorrectWithMinDeltaRJetB1ZZ->Fill(invMassB1CorrectWithMinDeltaRJetB1);
			if(invMassB1IncorrectWithMinDeltaRJetB1 != -999) histInvMassB1IncorrectWithMinDeltaRJetB1ZZ->Fill(invMassB1IncorrectWithMinDeltaRJetB1);
			if(invMassB1CorrectWithMinDeltaRJets != -999) histInvMassB1CorrectWithMinDeltaRJetsZZ->Fill(invMassB1CorrectWithMinDeltaRJets);
			if(invMassB1IncorrectWithMinDeltaRJets != -999) histInvMassB1IncorrectWithMinDeltaRJetsZZ->Fill(invMassB1IncorrectWithMinDeltaRJets);
			if(invMassB2CorrectMinChiSquaredHH != -999) histInvMassB2CorrectMinChiSquaredHHZZ->Fill(invMassB2CorrectMinChiSquaredHH);
			if(invMassB2IncorrectMinChiSquaredHH != -999) histInvMassB2IncorrectMinChiSquaredHHZZ->Fill(invMassB2IncorrectMinChiSquaredHH);
			if(invMassB2CorrectMinChiSquaredZZ != -999) histInvMassB2CorrectMinChiSquaredZZZZ->Fill(invMassB2CorrectMinChiSquaredZZ);
			if(invMassB2IncorrectMinChiSquaredZZ != -999) histInvMassB2IncorrectMinChiSquaredZZZZ->Fill(invMassB2IncorrectMinChiSquaredZZ);
			if(invMassB2CorrectWithMinDeltaRJetB1 != -999) histInvMassB2CorrectWithMinDeltaRJetB1ZZ->Fill(invMassB2CorrectWithMinDeltaRJetB1);
			if(invMassB2IncorrectWithMinDeltaRJetB1 != -999) histInvMassB2IncorrectWithMinDeltaRJetB1ZZ->Fill(invMassB2IncorrectWithMinDeltaRJetB1);
			if(invMassB2CorrectWithMinDeltaRJets != -999) histInvMassB2CorrectWithMinDeltaRJetsZZ->Fill(invMassB2CorrectWithMinDeltaRJets);
			if(invMassB2IncorrectWithMinDeltaRJets != -999) histInvMassB2IncorrectWithMinDeltaRJetsZZ->Fill(invMassB2IncorrectWithMinDeltaRJets);
			histThetaStarB1ZZ->Fill(thetaStarB1);
			histThetaStarB2ZZ->Fill(thetaStarB2);
			histThetaStarB3ZZ->Fill(thetaStarB3);
			histThetaStarB4ZZ->Fill(thetaStarB4);
			histThrustZZ->Fill(thrust);
		}
	 }
	 
	 treeWW->SetBranchAddress("invMassB1", &invMassB1);
	 treeWW->SetBranchAddress("invMassB2", &invMassB2);
	 treeWW->SetBranchAddress("cosThetaB1", &cosThetaB1);
	 treeWW->SetBranchAddress("sphericity", &sphericity);
	 nEntries = treeWW->GetEntries();
	 for (Long64_t i = 0; i < nEntries; ++i) {
	 	treeWW->GetEntry(i);
	 	histInvMassB1WW->Fill(invMassB1);
	 	histInvMassB2WW->Fill(invMassB2);
	 	histJet1CosThetaWW->Fill(cosThetaB1);
	 	histSpherWW->Fill(sphericity);
	 }
	 
	 /*histInvMassB1HH->Scale(1.0 / histInvMassB1HH->GetMaximum());
	 histInvMassB2HH->Scale(1.0 / histInvMassB2HH->GetMaximum());
	 histJet1CosThetaHH->Scale(1.0 / histJet1CosThetaHH->GetMaximum());
	 histSpherHH->Scale(1.0 / histSpherHH->GetMaximum());
	 histInvMassB1qq->Scale(1.0 / histInvMassB1qq->GetMaximum());
	 histInvMassB2qq->Scale(1.0 / histInvMassB2qq->GetMaximum());
	 histJet1CosThetaqq->Scale(1.0 / histJet1CosThetaqq->GetMaximum());
	 histSpherqq->Scale(1.0 / histSpherqq->GetMaximum());
	 histInvMassB1tt->Scale(1.0 / histInvMassB1tt->GetMaximum());
	 histInvMassB2tt->Scale(1.0 / histInvMassB2tt->GetMaximum());
	 histJet1CosThetatt->Scale(1.0 / histJet1CosThetatt->GetMaximum());
	 histSphertt->Scale(1.0 / histSphertt->GetMaximum());
	 histInvMassB1ZZ->Scale(1.0 / histInvMassB1ZZ->GetMaximum());
	 histInvMassB2ZZ->Scale(1.0 / histInvMassB2ZZ->GetMaximum());
	 histJet1CosThetaZZ->Scale(1.0 / histJet1CosThetaZZ->GetMaximum());
	 histSpherZZ->Scale(1.0 / histSpherZZ->GetMaximum());
	 histInvMassB1WW->Scale(1.0 / histInvMassB1WW->GetMaximum());
	 histInvMassB2WW->Scale(1.0 / histInvMassB2WW->GetMaximum());
	 histJet1CosThetaWW->Scale(1.0 / histJet1CosThetaWW->GetMaximum());
	 histSpherWW->Scale(1.0 / histSpherWW->GetMaximum());*/
	 
	 
    
	/*TLegend *legendInvMassB1HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1HH->AddEntry(histInvMassB1HH, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1HH->SetTextColor(kBlue);
	legendInvMassB1HH->SetBorderSize(0);
	legendInvMassB1HH->SetFillStyle(0);
	legendInvMassB1HH->SetTextSize(0.043);
	TLegend *legendInvMassB1Backqq = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1Backqq->AddEntry(histInvMassB1HH, "#gamma#gamma#rightarrow#it{qq}", "");
	legendInvMassB1Backqq->SetTextColor(kRed+2);
	legendInvMassB1Backqq->SetBorderSize(0);
	legendInvMassB1Backqq->SetFillStyle(0);
	legendInvMassB1Backqq->SetTextSize(0.043);
	TLegend *legendInvMassB1Backtt = new TLegend(0.65, 0.55, 0.8, 0.70);
	legendInvMassB1Backtt->AddEntry(histInvMassB1HH, "#gamma#gamma#rightarrow#it{t#bar{t}}", "");
	legendInvMassB1Backtt->SetTextColor(kPink+6);
	legendInvMassB1Backtt->SetBorderSize(0);
	legendInvMassB1Backtt->SetFillStyle(0);
	legendInvMassB1Backtt->SetTextSize(0.043);
	TLegend *legendInvMassB1BackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1BackZZ->AddEntry(histInvMassB1HH, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1BackZZ->SetTextColor(kGreen+1);
	legendInvMassB1BackZZ->SetBorderSize(0);
	legendInvMassB1BackZZ->SetFillStyle(0);
	legendInvMassB1BackZZ->SetTextSize(0.043);
	TLegend *legendInvMassB1BackWW = new TLegend(0.65, 0.45, 0.8, 0.60);
	legendInvMassB1BackWW->AddEntry(histInvMassB1HH, "#gamma#gamma#rightarrow#it{WW}", "");
	legendInvMassB1BackWW->SetTextColor(kYellow-2);
	legendInvMassB1BackWW->SetBorderSize(0);
	legendInvMassB1BackWW->SetFillStyle(0);
	legendInvMassB1BackWW->SetTextSize(0.043);
	TCanvas *c1 = new TCanvas();
	c1->cd();
	histInvMassB1HH->DrawNormalized("HIST");
	//histInvMassB1qq->DrawNormalized("HIST same");
	//histInvMassB1tt->DrawNormalized("HIST same");
	histInvMassB1ZZ->DrawNormalized("HIST same");
	//histInvMassB1WW->DrawNormalized("HIST same");
	legendInvMassB1HH->Draw();
	//legendInvMassB1Backqq->Draw();
	//legendInvMassB1Backtt->Draw();
	legendInvMassB1BackZZ->Draw();
	//legendInvMassB1BackWW->Draw();
	
	c1->Update();
	c1->SaveAs("invMassB1HHZZ.png");
	
	
	TLegend *legendInvMassB2HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2HH->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2HH->SetTextColor(kBlue);
	legendInvMassB2HH->SetBorderSize(0);
	legendInvMassB2HH->SetFillStyle(0);
	legendInvMassB2HH->SetTextSize(0.043);
	TLegend *legendInvMassB2Backqq = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2Backqq->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{qq}", "");
	legendInvMassB2Backqq->SetTextColor(kRed+2);
	legendInvMassB2Backqq->SetBorderSize(0);
	legendInvMassB2Backqq->SetFillStyle(0);
	legendInvMassB2Backqq->SetTextSize(0.043);
	TLegend *legendInvMassB2Backtt = new TLegend(0.65, 0.55, 0.8, 0.70);
	legendInvMassB2Backtt->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{t#bar{t}}", "");
	legendInvMassB2Backtt->SetTextColor(kPink+6);
	legendInvMassB2Backtt->SetBorderSize(0);
	legendInvMassB2Backtt->SetFillStyle(0);
	legendInvMassB2Backtt->SetTextSize(0.043);
	TLegend *legendInvMassB2BackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2BackZZ->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2BackZZ->SetTextColor(kGreen+1);
	legendInvMassB2BackZZ->SetBorderSize(0);
	legendInvMassB2BackZZ->SetFillStyle(0);
	legendInvMassB2BackZZ->SetTextSize(0.043);
	TLegend *legendInvMassB2BackWW = new TLegend(0.65, 0.45, 0.8, 0.60);
	legendInvMassB2BackWW->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{WW}", "");
	legendInvMassB2BackWW->SetTextColor(kYellow-2);
	legendInvMassB2BackWW->SetBorderSize(0);
	legendInvMassB2BackWW->SetFillStyle(0);
	legendInvMassB2BackWW->SetTextSize(0.043);
	TCanvas *c2 = new TCanvas();
	c2->cd();
	histInvMassB2HH->DrawNormalized("HIST");
	//histInvMassB2qq->DrawNormalized("HIST same");
	//histInvMassB2tt->DrawNormalized("HIST same");
	histInvMassB2ZZ->DrawNormalized("HIST same");
	//histInvMassB2WW->DrawNormalized("HIST same");
	legendInvMassB2HH->Draw();
	//legendInvMassB2Backqq->Draw();
	//legendInvMassB2Backtt->Draw();
	//legendInvMassB2BackWW->Draw();
	legendInvMassB2BackZZ->Draw();
	
	c2->Update();
	c2->SaveAs("invMassB2HHZZ.png");

	TLegend *legendCosThetaHH = new TLegend(0.2, 0.77, 0.35, 0.92);
	legendCosThetaHH->AddEntry(legendCosThetaHH, "#gamma#gamma#rightarrow#it{HH}", "");
	legendCosThetaHH->SetTextColor(kBlue);
	legendCosThetaHH->SetBorderSize(0);
	legendCosThetaHH->SetFillStyle(0);
	legendCosThetaHH->SetTextSize(0.043);
	TLegend *legendCosThetaBackqq = new TLegend(0.2, 0.72, 0.35, 0.87);
	legendCosThetaBackqq->AddEntry(legendCosThetaHH, "#gamma#gamma#rightarrow#it{qq}", "");
	legendCosThetaBackqq->SetTextColor(kRed+2);
	legendCosThetaBackqq->SetBorderSize(0);
	legendCosThetaBackqq->SetFillStyle(0);
	legendCosThetaBackqq->SetTextSize(0.043);
	TLegend *legendCosThetaBacktt = new TLegend(0.65, 0.77, 0.8, 0.92);
	legendCosThetaBacktt->AddEntry(legendCosThetaHH, "#gamma#gamma#rightarrow#it{t#bar{t}}", "");
	legendCosThetaBacktt->SetTextColor(kPink+6);
	legendCosThetaBacktt->SetBorderSize(0);
	legendCosThetaBacktt->SetFillStyle(0);
	legendCosThetaBacktt->SetTextSize(0.043);
	TLegend *legendCosThetaBackZZ = new TLegend(0.65, 0.72, 0.8, 0.87);
	legendCosThetaBackZZ->AddEntry(legendCosThetaHH, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendCosThetaBackZZ->SetTextColor(kGreen+1);
	legendCosThetaBackZZ->SetBorderSize(0);
	legendCosThetaBackZZ->SetFillStyle(0);
	legendCosThetaBackZZ->SetTextSize(0.043);
	TLegend *legendCosThetaBackWW = new TLegend(0.65, 0.67, 0.8, 0.82);
	legendCosThetaBackWW->AddEntry(legendCosThetaHH, "#gamma#gamma#rightarrow#it{WW}", "");
	legendCosThetaBackWW->SetTextColor(kYellow-2);
	legendCosThetaBackWW->SetBorderSize(0);
	legendCosThetaBackWW->SetFillStyle(0);
	legendCosThetaBackWW->SetTextSize(0.043);
	TCanvas *c3 = new TCanvas();
	c3->cd();
	histJet1CosThetaHH->DrawNormalized("HIST");
	histJet1CosThetaqq->DrawNormalized("HIST same");
	histJet1CosThetatt->DrawNormalized("HIST same");
	histJet1CosThetaZZ->DrawNormalized("HIST same");
	histJet1CosThetaWW->DrawNormalized("HIST same");
	legendCosThetaHH->Draw();
	legendCosThetaBackqq->Draw();
	legendCosThetaBacktt->Draw();
	legendCosThetaBackZZ->Draw();
	legendCosThetaBackWW->Draw();
	
	c3->Update();
	c3->SaveAs("cosThetaB1HHqqttbarZZWW.png");
	
	TLegend *legendSpherHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendSpherHH->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{HH}", "");
	legendSpherHH->SetTextColor(kBlue);
	legendSpherHH->SetBorderSize(0);
	legendSpherHH->SetFillStyle(0);
	legendSpherHH->SetTextSize(0.043);
	TLegend *legendSpherBackqq = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendSpherBackqq->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{qq}", "");
	legendSpherBackqq->SetTextColor(kRed+2);
	legendSpherBackqq->SetBorderSize(0);
	legendSpherBackqq->SetFillStyle(0);
	legendSpherBackqq->SetTextSize(0.043);
	TLegend *legendSpherBacktt = new TLegend(0.65, 0.55, 0.8, 0.70);
	legendSpherBacktt->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{t#bar{t}}", "");
	legendSpherBacktt->SetTextColor(kPink+6);
	legendSpherBacktt->SetBorderSize(0);
	legendSpherBacktt->SetFillStyle(0);
	legendSpherBacktt->SetTextSize(0.043);
	TLegend *legendSpherBackZZ = new TLegend(0.65, 0.5, 0.8, 0.65);
	legendSpherBackZZ->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendSpherBackZZ->SetTextColor(kGreen+1);
	legendSpherBackZZ->SetBorderSize(0);
	legendSpherBackZZ->SetFillStyle(0);
	legendSpherBackZZ->SetTextSize(0.043);
	TLegend *legendSpherBackWW = new TLegend(0.65, 0.45, 0.8, 0.60);
	legendSpherBackWW->AddEntry(histInvMassB2HH, "#gamma#gamma#rightarrow#it{WW}", "");
	legendSpherBackWW->SetTextColor(kYellow-2);
	legendSpherBackWW->SetBorderSize(0);
	legendSpherBackWW->SetFillStyle(0);
	legendSpherBackWW->SetTextSize(0.043);
	TCanvas *c4 = new TCanvas();
	c4->cd();
	histSpherqq->DrawNormalized("HIST");
	histSphertt->DrawNormalized("HIST same");
	histSpherZZ->DrawNormalized("HIST same");
	histSpherWW->DrawNormalized("HIST same");
	histSpherHH->DrawNormalized("HIST same");
	histSpherqq->DrawNormalized("HIST same");
	histSphertt->DrawNormalized("HIST same");
	histSpherZZ->DrawNormalized("HIST same");
	histSpherWW->DrawNormalized("HIST same");
	legendSpherHH->Draw();
	legendSpherBackqq->Draw();
	legendSpherBacktt->Draw();
	legendSpherBackZZ->Draw();
	legendSpherBackWW->Draw();
	
	c4->Update();
	c4->SaveAs("sphericityHHqqttbarZZWW.png");
	*/
	
	
	/*TLegend *legendExclYmerge12HH = new TLegend(0.25, 0.65, 0.4, 0.8);
	legendExclYmerge12HH->AddEntry(histExclYmerge12ZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendExclYmerge12HH->SetTextColor(kBlue);
	legendExclYmerge12HH->SetBorderSize(0);
	legendExclYmerge12HH->SetFillStyle(0);
	legendExclYmerge12HH->SetTextSize(0.043);
	TLegend *legendExclYmerge12BackZZ = new TLegend(0.25, 0.6, 0.4, 0.75);
	legendExclYmerge12BackZZ->AddEntry(histExclYmerge12ZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendExclYmerge12BackZZ->SetTextColor(kGreen+1);
	legendExclYmerge12BackZZ->SetBorderSize(0);
	legendExclYmerge12BackZZ->SetFillStyle(0);
	legendExclYmerge12BackZZ->SetTextSize(0.043);
	TCanvas *c5 = new TCanvas();
	c5->cd();
	histExclYmerge12ZZ->DrawNormalized("HIST");
	histExclYmerge12HH->DrawNormalized("HIST same");
	histExclYmerge12ZZ->DrawNormalized("HIST same");
	legendExclYmerge12HH->Draw();
	legendExclYmerge12BackZZ->Draw();
	c5->Update();
	c5->SaveAs("exclYmerge12HHZZ.png");
	
	TLegend *legendExclYmerge23HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendExclYmerge23HH->AddEntry(histExclYmerge23ZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendExclYmerge23HH->SetTextColor(kBlue);
	legendExclYmerge23HH->SetBorderSize(0);
	legendExclYmerge23HH->SetFillStyle(0);
	legendExclYmerge23HH->SetTextSize(0.043);
	TLegend *legendExclYmerge23BackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendExclYmerge23BackZZ->AddEntry(histExclYmerge23ZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendExclYmerge23BackZZ->SetTextColor(kGreen+1);
	legendExclYmerge23BackZZ->SetBorderSize(0);
	legendExclYmerge23BackZZ->SetFillStyle(0);
	legendExclYmerge23BackZZ->SetTextSize(0.043);
	TCanvas *c6 = new TCanvas();
	c6->cd();
	histExclYmerge23ZZ->DrawNormalized("HIST");
	histExclYmerge23HH->DrawNormalized("HIST same");
	histExclYmerge23ZZ->DrawNormalized("HIST same");
	legendExclYmerge23HH->Draw();
	legendExclYmerge23BackZZ->Draw();
	c6->Update();
	c6->SaveAs("exclYmerge23HHZZ.png");
	
	TLegend *legendExclYmerge34HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendExclYmerge34HH->AddEntry(histExclYmerge34ZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendExclYmerge34HH->SetTextColor(kBlue);
	legendExclYmerge34HH->SetBorderSize(0);
	legendExclYmerge34HH->SetFillStyle(0);
	legendExclYmerge34HH->SetTextSize(0.043);
	TLegend *legendExclYmerge34BackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendExclYmerge34BackZZ->AddEntry(histExclYmerge34ZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendExclYmerge34BackZZ->SetTextColor(kGreen+1);
	legendExclYmerge34BackZZ->SetBorderSize(0);
	legendExclYmerge34BackZZ->SetFillStyle(0);
	legendExclYmerge34BackZZ->SetTextSize(0.043);
	TCanvas *c7 = new TCanvas();
	c7->cd();
	histExclYmerge34ZZ->DrawNormalized("HIST");
	histExclYmerge34HH->DrawNormalized("HIST same");
	histExclYmerge34ZZ->DrawNormalized("HIST same");
	legendExclYmerge34HH->Draw();
	legendExclYmerge34BackZZ->Draw();
	c7->Update();
	c7->SaveAs("exclYmerge34HHZZ.png");
	
	TLegend *legendExclYmerge45HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendExclYmerge45HH->AddEntry(histExclYmerge45ZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendExclYmerge45HH->SetTextColor(kBlue);
	legendExclYmerge45HH->SetBorderSize(0);
	legendExclYmerge45HH->SetFillStyle(0);
	legendExclYmerge45HH->SetTextSize(0.043);
	TLegend *legendExclYmerge45BackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendExclYmerge45BackZZ->AddEntry(histExclYmerge45ZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendExclYmerge45BackZZ->SetTextColor(kGreen+1);
	legendExclYmerge45BackZZ->SetBorderSize(0);
	legendExclYmerge45BackZZ->SetFillStyle(0);
	legendExclYmerge45BackZZ->SetTextSize(0.043);
	TCanvas *c8 = new TCanvas();
	c8->cd();
	histExclYmerge45ZZ->DrawNormalized("HIST");
	histExclYmerge45HH->DrawNormalized("HIST same");
	histExclYmerge45ZZ->DrawNormalized("HIST same");
	legendExclYmerge45HH->Draw();
	legendExclYmerge45BackZZ->Draw();
	c8->Update();
	c8->SaveAs("exclYmerge45HHZZ.png");
	
	TLegend *legendExclYmerge56HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendExclYmerge56HH->AddEntry(histExclYmerge56ZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendExclYmerge56HH->SetTextColor(kBlue);
	legendExclYmerge56HH->SetBorderSize(0);
	legendExclYmerge56HH->SetFillStyle(0);
	legendExclYmerge56HH->SetTextSize(0.043);
	TLegend *legendExclYmerge56BackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendExclYmerge56BackZZ->AddEntry(histExclYmerge56ZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendExclYmerge56BackZZ->SetTextColor(kGreen+1);
	legendExclYmerge56BackZZ->SetBorderSize(0);
	legendExclYmerge56BackZZ->SetFillStyle(0);
	legendExclYmerge56BackZZ->SetTextSize(0.043);
	TCanvas *c9 = new TCanvas();
	c9->cd();
	histExclYmerge56ZZ->DrawNormalized("HIST");
	histExclYmerge56HH->DrawNormalized("HIST same");
	histExclYmerge56ZZ->DrawNormalized("HIST same");
	legendExclYmerge56HH->Draw();
	legendExclYmerge56BackZZ->Draw();
	c9->Update();
	c9->SaveAs("exclYmerge56HHZZ.png");
	
	TLegend *legendInvMassZZ1HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassZZ1HH->AddEntry(histInvMassZZ1ZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassZZ1HH->SetTextColor(kBlue);
	legendInvMassZZ1HH->SetBorderSize(0);
	legendInvMassZZ1HH->SetFillStyle(0);
	legendInvMassZZ1HH->SetTextSize(0.043);
	TLegend *legendInvMassZZ1BackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassZZ1BackZZ->AddEntry(histInvMassZZ1ZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassZZ1BackZZ->SetTextColor(kGreen+1);
	legendInvMassZZ1BackZZ->SetBorderSize(0);
	legendInvMassZZ1BackZZ->SetFillStyle(0);
	legendInvMassZZ1BackZZ->SetTextSize(0.043);
	TCanvas *c10 = new TCanvas();
	c10->cd();
	//histInvMassZZ1ZZ->DrawNormalized("HIST");
	histInvMassZZ1ZZ->DrawNormalized("HIST");
	histInvMassZZ1HH->DrawNormalized("HIST same");
	histInvMassZZ1ZZ->DrawNormalized("HIST same");
	legendInvMassZZ1HH->Draw();
	legendInvMassZZ1BackZZ->Draw();
	c10->Update();
	c10->SaveAs("InvMassZZ1HHZZ.png");
	
	TLegend *legendInvMassZZ2HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassZZ2HH->AddEntry(histInvMassZZ2ZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassZZ2HH->SetTextColor(kBlue);
	legendInvMassZZ2HH->SetBorderSize(0);
	legendInvMassZZ2HH->SetFillStyle(0);
	legendInvMassZZ2HH->SetTextSize(0.043);
	TLegend *legendInvMassZZ2BackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassZZ2BackZZ->AddEntry(histInvMassZZ2ZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassZZ2BackZZ->SetTextColor(kGreen+1);
	legendInvMassZZ2BackZZ->SetBorderSize(0);
	legendInvMassZZ2BackZZ->SetFillStyle(0);
	legendInvMassZZ2BackZZ->SetTextSize(0.043);
	TCanvas *c11 = new TCanvas();
	c11->cd();
	histInvMassZZ2ZZ->DrawNormalized("HIST");
	histInvMassZZ2HH->DrawNormalized("HIST same");
	histInvMassZZ2ZZ->DrawNormalized("HIST same");
	legendInvMassZZ2HH->Draw();
	legendInvMassZZ2BackZZ->Draw();
	c11->Update();
	c11->SaveAs("InvMassZZ2HHZZ.png");
	
	TLegend *legendDistanceZ1MinChiSquaredZZMassHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendDistanceZ1MinChiSquaredZZMassHH->AddEntry(histDistanceZ1MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendDistanceZ1MinChiSquaredZZMassHH->SetTextColor(kBlue);
	legendDistanceZ1MinChiSquaredZZMassHH->SetBorderSize(0);
	legendDistanceZ1MinChiSquaredZZMassHH->SetFillStyle(0);
	legendDistanceZ1MinChiSquaredZZMassHH->SetTextSize(0.043);
	TLegend *legendDistanceZ1MinChiSquaredZZMassBackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendDistanceZ1MinChiSquaredZZMassBackZZ->AddEntry(histDistanceZ1MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendDistanceZ1MinChiSquaredZZMassBackZZ->SetTextColor(kGreen+1);
	legendDistanceZ1MinChiSquaredZZMassBackZZ->SetBorderSize(0);
	legendDistanceZ1MinChiSquaredZZMassBackZZ->SetFillStyle(0);
	legendDistanceZ1MinChiSquaredZZMassBackZZ->SetTextSize(0.043);
	TCanvas *c12 = new TCanvas();
	c12->cd();
	//histDistanceZ1MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histDistanceZ1MinChiSquaredZZMassHH->DrawNormalized("HIST");
	histDistanceZ1MinChiSquaredZZMassZZ->DrawNormalized("HIST same");
	legendDistanceZ1MinChiSquaredZZMassHH->Draw();
	legendDistanceZ1MinChiSquaredZZMassBackZZ->Draw();
	TLine* distanceZ1MassRight = new TLine(30, histDistanceZ1MinChiSquaredZZMassZZ->GetMinimum(), 30, histDistanceZ1MinChiSquaredZZMassZZ->GetMaximum());
	distanceZ1MassRight->SetLineColor(kRed);
	distanceZ1MassRight->SetLineWidth(1);
	distanceZ1MassRight->Draw("same");
	c12->Update();
	c12->SaveAs("DistanceZ1MinChiSquaredZZMassHHZZ.png");
	
	TLegend *legendDistanceZ2MinChiSquaredZZMassHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendDistanceZ2MinChiSquaredZZMassHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendDistanceZ2MinChiSquaredZZMassHH->SetTextColor(kBlue);
	legendDistanceZ2MinChiSquaredZZMassHH->SetBorderSize(0);
	legendDistanceZ2MinChiSquaredZZMassHH->SetFillStyle(0);
	legendDistanceZ2MinChiSquaredZZMassHH->SetTextSize(0.043);
	TLegend *legendDistanceZ2MinChiSquaredZZMassBackZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendDistanceZ2MinChiSquaredZZMassBackZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendDistanceZ2MinChiSquaredZZMassBackZZ->SetTextColor(kGreen+1);
	legendDistanceZ2MinChiSquaredZZMassBackZZ->SetBorderSize(0);
	legendDistanceZ2MinChiSquaredZZMassBackZZ->SetFillStyle(0);
	legendDistanceZ2MinChiSquaredZZMassBackZZ->SetTextSize(0.043);
	TCanvas *c13 = new TCanvas();
	c13->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histDistanceZ2MinChiSquaredZZMassHH->DrawNormalized("HIST");
	histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST same");
	legendDistanceZ2MinChiSquaredZZMassHH->Draw();
	legendDistanceZ2MinChiSquaredZZMassBackZZ->Draw();
	TLine* distanceZ2MassRight = new TLine(30, histDistanceZ2MinChiSquaredZZMassZZ->GetMinimum(), 30, histDistanceZ2MinChiSquaredZZMassZZ->GetMaximum());
	distanceZ2MassRight->SetLineColor(kRed);
	distanceZ2MassRight->SetLineWidth(1);
	distanceZ2MassRight->Draw("same");
	c13->Update();
	c13->SaveAs("DistanceZ2MinChiSquaredZZMassHHZZ.png");
	
	TCanvas *c14 = new TCanvas();
	c14->cd();
	histCompareAlgosB1HH->DrawNormalized("col");
	c14->Update();
	c14->SaveAs("CompareAlgosB1HH.png");
	
	TCanvas *c15 = new TCanvas();
	c15->cd();
	histCompareAlgosB1ZZ->DrawNormalized("col");
	c15->Update();
	c15->SaveAs("CompareAlgosB1ZZ.png");
	
	TCanvas *c16 = new TCanvas();
	c16->cd();
	histCompareAlgosB2HH->DrawNormalized("col");
	c16->Update();
	c16->SaveAs("CompareAlgosB2HH.png");
	
	TCanvas *c17 = new TCanvas();
	c17->cd();
	histCompareAlgosB2ZZ->DrawNormalized("col");
	c17->Update();
	c17->SaveAs("CompareAlgosB2ZZ.png");
	
	TLegend *legendInvMassB1CorrectMinChiSquaredHHHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1CorrectMinChiSquaredHHHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1CorrectMinChiSquaredHHHH->SetTextColor(kBlue);
	legendInvMassB1CorrectMinChiSquaredHHHH->SetBorderSize(0);
	legendInvMassB1CorrectMinChiSquaredHHHH->SetFillStyle(0);
	legendInvMassB1CorrectMinChiSquaredHHHH->SetTextSize(0.043);
	TLegend *legendInvMassB1CorrectMinChiSquaredHHZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1CorrectMinChiSquaredHHZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1CorrectMinChiSquaredHHZZ->SetTextColor(kGreen+1);
	legendInvMassB1CorrectMinChiSquaredHHZZ->SetBorderSize(0);
	legendInvMassB1CorrectMinChiSquaredHHZZ->SetFillStyle(0);
	legendInvMassB1CorrectMinChiSquaredHHZZ->SetTextSize(0.043);
	TCanvas *c18 = new TCanvas();
	c18->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB1CorrectMinChiSquaredHHZZ->DrawNormalized("HIST");
	histInvMassB1CorrectMinChiSquaredHHHH->DrawNormalized("HIST same");
	histInvMassB1CorrectMinChiSquaredHHZZ->DrawNormalized("HIST same");
	legendInvMassB1CorrectMinChiSquaredHHHH->Draw();
	legendInvMassB1CorrectMinChiSquaredHHZZ->Draw();
	c18->Update();
	c18->SaveAs("InvMassB1CorrectMinChiSquaredHHHHZZ.png");
	
	TLegend *legendInvMassB1IncorrectMinChiSquaredHHHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1IncorrectMinChiSquaredHHHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1IncorrectMinChiSquaredHHHH->SetTextColor(kBlue);
	legendInvMassB1IncorrectMinChiSquaredHHHH->SetBorderSize(0);
	legendInvMassB1IncorrectMinChiSquaredHHHH->SetFillStyle(0);
	legendInvMassB1IncorrectMinChiSquaredHHHH->SetTextSize(0.043);
	TLegend *legendInvMassB1IncorrectMinChiSquaredHHZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1IncorrectMinChiSquaredHHZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1IncorrectMinChiSquaredHHZZ->SetTextColor(kGreen+1);
	legendInvMassB1IncorrectMinChiSquaredHHZZ->SetBorderSize(0);
	legendInvMassB1IncorrectMinChiSquaredHHZZ->SetFillStyle(0);
	legendInvMassB1IncorrectMinChiSquaredHHZZ->SetTextSize(0.043);
	TCanvas *c19 = new TCanvas();
	c19->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB1IncorrectMinChiSquaredHHHH->DrawNormalized("HIST");
	histInvMassB1IncorrectMinChiSquaredHHZZ->DrawNormalized("HIST same");
	legendInvMassB1IncorrectMinChiSquaredHHHH->Draw();
	legendInvMassB1IncorrectMinChiSquaredHHZZ->Draw();
	c19->Update();
	c19->SaveAs("InvMassB1IncorrectMinChiSquaredHHHHZZ.png");
	
	TLegend *legendInvMassB2CorrectMinChiSquaredHHHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2CorrectMinChiSquaredHHHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2CorrectMinChiSquaredHHHH->SetTextColor(kBlue);
	legendInvMassB2CorrectMinChiSquaredHHHH->SetBorderSize(0);
	legendInvMassB2CorrectMinChiSquaredHHHH->SetFillStyle(0);
	legendInvMassB2CorrectMinChiSquaredHHHH->SetTextSize(0.043);
	TLegend *legendInvMassB2CorrectMinChiSquaredHHZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2CorrectMinChiSquaredHHZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2CorrectMinChiSquaredHHZZ->SetTextColor(kGreen+1);
	legendInvMassB2CorrectMinChiSquaredHHZZ->SetBorderSize(0);
	legendInvMassB2CorrectMinChiSquaredHHZZ->SetFillStyle(0);
	legendInvMassB2CorrectMinChiSquaredHHZZ->SetTextSize(0.043);
	TCanvas *c20 = new TCanvas();
	c20->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB2CorrectMinChiSquaredHHZZ->DrawNormalized("HIST");
	histInvMassB2CorrectMinChiSquaredHHHH->DrawNormalized("HIST same");
	histInvMassB2CorrectMinChiSquaredHHZZ->DrawNormalized("HIST same");
	legendInvMassB2CorrectMinChiSquaredHHHH->Draw();
	legendInvMassB2CorrectMinChiSquaredHHZZ->Draw();
	c20->Update();
	c20->SaveAs("InvMassB2CorrectMinChiSquaredHHHHZZ.png");
	
	TLegend *legendInvMassB2IncorrectMinChiSquaredHHHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2IncorrectMinChiSquaredHHHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2IncorrectMinChiSquaredHHHH->SetTextColor(kBlue);
	legendInvMassB2IncorrectMinChiSquaredHHHH->SetBorderSize(0);
	legendInvMassB2IncorrectMinChiSquaredHHHH->SetFillStyle(0);
	legendInvMassB2IncorrectMinChiSquaredHHHH->SetTextSize(0.043);
	TLegend *legendInvMassB2IncorrectMinChiSquaredHHZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2IncorrectMinChiSquaredHHZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2IncorrectMinChiSquaredHHZZ->SetTextColor(kGreen+1);
	legendInvMassB2IncorrectMinChiSquaredHHZZ->SetBorderSize(0);
	legendInvMassB2IncorrectMinChiSquaredHHZZ->SetFillStyle(0);
	legendInvMassB2IncorrectMinChiSquaredHHZZ->SetTextSize(0.043);
	TCanvas *c21 = new TCanvas();
	c21->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB2IncorrectMinChiSquaredHHHH->DrawNormalized("HIST");
	histInvMassB2IncorrectMinChiSquaredHHZZ->DrawNormalized("HIST same");
	legendInvMassB2IncorrectMinChiSquaredHHHH->Draw();
	legendInvMassB2IncorrectMinChiSquaredHHZZ->Draw();
	c21->Update();
	c21->SaveAs("InvMassB2IncorrectMinChiSquaredHHHHZZ.png");
	
	TLegend *legendInvMassB1CorrectMinChiSquaredZZHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1CorrectMinChiSquaredZZHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1CorrectMinChiSquaredZZHH->SetTextColor(kBlue);
	legendInvMassB1CorrectMinChiSquaredZZHH->SetBorderSize(0);
	legendInvMassB1CorrectMinChiSquaredZZHH->SetFillStyle(0);
	legendInvMassB1CorrectMinChiSquaredZZHH->SetTextSize(0.043);
	TLegend *legendInvMassB1CorrectMinChiSquaredZZZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1CorrectMinChiSquaredZZZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1CorrectMinChiSquaredZZZZ->SetTextColor(kGreen+1);
	legendInvMassB1CorrectMinChiSquaredZZZZ->SetBorderSize(0);
	legendInvMassB1CorrectMinChiSquaredZZZZ->SetFillStyle(0);
	legendInvMassB1CorrectMinChiSquaredZZZZ->SetTextSize(0.043);
	TCanvas *c22 = new TCanvas();
	c22->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB1CorrectMinChiSquaredZZZZ->DrawNormalized("HIST");
	histInvMassB1CorrectMinChiSquaredZZHH->DrawNormalized("HIST same");
	histInvMassB1CorrectMinChiSquaredZZZZ->DrawNormalized("HIST same");
	legendInvMassB1CorrectMinChiSquaredZZHH->Draw();
	legendInvMassB1CorrectMinChiSquaredZZZZ->Draw();
	c22->Update();
	c22->SaveAs("InvMassB1CorrectMinChiSquaredZZHHZZ.png");
	
	TLegend *legendInvMassB1IncorrectMinChiSquaredZZHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1IncorrectMinChiSquaredZZHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1IncorrectMinChiSquaredZZHH->SetTextColor(kBlue);
	legendInvMassB1IncorrectMinChiSquaredZZHH->SetBorderSize(0);
	legendInvMassB1IncorrectMinChiSquaredZZHH->SetFillStyle(0);
	legendInvMassB1IncorrectMinChiSquaredZZHH->SetTextSize(0.043);
	TLegend *legendInvMassB1IncorrectMinChiSquaredZZZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1IncorrectMinChiSquaredZZZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1IncorrectMinChiSquaredZZZZ->SetTextColor(kGreen+1);
	legendInvMassB1IncorrectMinChiSquaredZZZZ->SetBorderSize(0);
	legendInvMassB1IncorrectMinChiSquaredZZZZ->SetFillStyle(0);
	legendInvMassB1IncorrectMinChiSquaredZZZZ->SetTextSize(0.043);
	TCanvas *c23 = new TCanvas();
	c23->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB1IncorrectMinChiSquaredZZHH->DrawNormalized("HIST");
	histInvMassB1IncorrectMinChiSquaredZZZZ->DrawNormalized("HIST same");
	legendInvMassB1IncorrectMinChiSquaredZZHH->Draw();
	legendInvMassB1IncorrectMinChiSquaredZZZZ->Draw();
	c23->Update();
	c23->SaveAs("InvMassB1IncorrectMinChiSquaredZZHHZZ.png");
	
	TLegend *legendInvMassB2CorrectMinChiSquaredZZHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2CorrectMinChiSquaredZZHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2CorrectMinChiSquaredZZHH->SetTextColor(kBlue);
	legendInvMassB2CorrectMinChiSquaredZZHH->SetBorderSize(0);
	legendInvMassB2CorrectMinChiSquaredZZHH->SetFillStyle(0);
	legendInvMassB2CorrectMinChiSquaredZZHH->SetTextSize(0.043);
	TLegend *legendInvMassB2CorrectMinChiSquaredZZZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2CorrectMinChiSquaredZZZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2CorrectMinChiSquaredZZZZ->SetTextColor(kGreen+1);
	legendInvMassB2CorrectMinChiSquaredZZZZ->SetBorderSize(0);
	legendInvMassB2CorrectMinChiSquaredZZZZ->SetFillStyle(0);
	legendInvMassB2CorrectMinChiSquaredZZZZ->SetTextSize(0.043);
	TCanvas *c24 = new TCanvas();
	c24->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB2CorrectMinChiSquaredZZZZ->DrawNormalized("HIST");
	histInvMassB2CorrectMinChiSquaredZZHH->DrawNormalized("HIST same");
	histInvMassB2CorrectMinChiSquaredZZZZ->DrawNormalized("HIST same");
	legendInvMassB2CorrectMinChiSquaredZZHH->Draw();
	legendInvMassB2CorrectMinChiSquaredZZZZ->Draw();
	c24->Update();
	c24->SaveAs("InvMassB2CorrectMinChiSquaredZZHHZZ.png");
	
	TLegend *legendInvMassB2IncorrectMinChiSquaredZZHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2IncorrectMinChiSquaredZZHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2IncorrectMinChiSquaredZZHH->SetTextColor(kBlue);
	legendInvMassB2IncorrectMinChiSquaredZZHH->SetBorderSize(0);
	legendInvMassB2IncorrectMinChiSquaredZZHH->SetFillStyle(0);
	legendInvMassB2IncorrectMinChiSquaredZZHH->SetTextSize(0.043);
	TLegend *legendInvMassB2IncorrectMinChiSquaredZZZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2IncorrectMinChiSquaredZZZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2IncorrectMinChiSquaredZZZZ->SetTextColor(kGreen+1);
	legendInvMassB2IncorrectMinChiSquaredZZZZ->SetBorderSize(0);
	legendInvMassB2IncorrectMinChiSquaredZZZZ->SetFillStyle(0);
	legendInvMassB2IncorrectMinChiSquaredZZZZ->SetTextSize(0.043);
	TCanvas *c25 = new TCanvas();
	c25->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB2IncorrectMinChiSquaredZZHH->DrawNormalized("HIST");
	histInvMassB2IncorrectMinChiSquaredZZZZ->DrawNormalized("HIST same");
	legendInvMassB2IncorrectMinChiSquaredZZHH->Draw();
	legendInvMassB2IncorrectMinChiSquaredZZZZ->Draw();
	c25->Update();
	c25->SaveAs("InvMassB2IncorrectMinChiSquaredZZHHZZ.png");
	
	TLegend *legendInvMassB1CorrectWithMinDeltaRJetB1HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1CorrectWithMinDeltaRJetB1HH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1CorrectWithMinDeltaRJetB1HH->SetTextColor(kBlue);
	legendInvMassB1CorrectWithMinDeltaRJetB1HH->SetBorderSize(0);
	legendInvMassB1CorrectWithMinDeltaRJetB1HH->SetFillStyle(0);
	legendInvMassB1CorrectWithMinDeltaRJetB1HH->SetTextSize(0.043);
	TLegend *legendInvMassB1CorrectWithMinDeltaRJetB1ZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1CorrectWithMinDeltaRJetB1ZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1CorrectWithMinDeltaRJetB1ZZ->SetTextColor(kGreen+1);
	legendInvMassB1CorrectWithMinDeltaRJetB1ZZ->SetBorderSize(0);
	legendInvMassB1CorrectWithMinDeltaRJetB1ZZ->SetFillStyle(0);
	legendInvMassB1CorrectWithMinDeltaRJetB1ZZ->SetTextSize(0.043);
	TCanvas *c26 = new TCanvas();
	c26->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB1CorrectWithMinDeltaRJetB1ZZ->DrawNormalized("HIST");
	histInvMassB1CorrectWithMinDeltaRJetB1HH->DrawNormalized("HIST same");
	histInvMassB1CorrectWithMinDeltaRJetB1ZZ->DrawNormalized("HIST same");
	legendInvMassB1CorrectWithMinDeltaRJetB1HH->Draw();
	legendInvMassB1CorrectWithMinDeltaRJetB1ZZ->Draw();
	c26->Update();
	c26->SaveAs("InvMassB1CorrectWithMinDeltaRJetB1HHZZ.png");
	
	TLegend *legendInvMassB1IncorrectWithMinDeltaRJetB1HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1IncorrectWithMinDeltaRJetB1HH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1IncorrectWithMinDeltaRJetB1HH->SetTextColor(kBlue);
	legendInvMassB1IncorrectWithMinDeltaRJetB1HH->SetBorderSize(0);
	legendInvMassB1IncorrectWithMinDeltaRJetB1HH->SetFillStyle(0);
	legendInvMassB1IncorrectWithMinDeltaRJetB1HH->SetTextSize(0.043);
	TLegend *legendInvMassB1IncorrectWithMinDeltaRJetB1ZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1IncorrectWithMinDeltaRJetB1ZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1IncorrectWithMinDeltaRJetB1ZZ->SetTextColor(kGreen+1);
	legendInvMassB1IncorrectWithMinDeltaRJetB1ZZ->SetBorderSize(0);
	legendInvMassB1IncorrectWithMinDeltaRJetB1ZZ->SetFillStyle(0);
	legendInvMassB1IncorrectWithMinDeltaRJetB1ZZ->SetTextSize(0.043);
	TCanvas *c27 = new TCanvas();
	c27->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB1IncorrectWithMinDeltaRJetB1HH->DrawNormalized("HIST");
	histInvMassB1IncorrectWithMinDeltaRJetB1ZZ->DrawNormalized("HIST same");
	legendInvMassB1IncorrectWithMinDeltaRJetB1HH->Draw();
	legendInvMassB1IncorrectWithMinDeltaRJetB1ZZ->Draw();
	c27->Update();
	c27->SaveAs("InvMassB1IncorrectWithMinDeltaRJetB1HHZZ.png");
	
	TLegend *legendInvMassB2CorrectWithMinDeltaRJetB1HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2CorrectWithMinDeltaRJetB1HH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2CorrectWithMinDeltaRJetB1HH->SetTextColor(kBlue);
	legendInvMassB2CorrectWithMinDeltaRJetB1HH->SetBorderSize(0);
	legendInvMassB2CorrectWithMinDeltaRJetB1HH->SetFillStyle(0);
	legendInvMassB2CorrectWithMinDeltaRJetB1HH->SetTextSize(0.043);
	TLegend *legendInvMassB2CorrectWithMinDeltaRJetB1ZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2CorrectWithMinDeltaRJetB1ZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2CorrectWithMinDeltaRJetB1ZZ->SetTextColor(kGreen+1);
	legendInvMassB2CorrectWithMinDeltaRJetB1ZZ->SetBorderSize(0);
	legendInvMassB2CorrectWithMinDeltaRJetB1ZZ->SetFillStyle(0);
	legendInvMassB2CorrectWithMinDeltaRJetB1ZZ->SetTextSize(0.043);
	TCanvas *c28 = new TCanvas();
	c28->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB2CorrectWithMinDeltaRJetB1ZZ->DrawNormalized("HIST");
	histInvMassB2CorrectWithMinDeltaRJetB1HH->DrawNormalized("HIST same");
	histInvMassB2CorrectWithMinDeltaRJetB1ZZ->DrawNormalized("HIST same");
	legendInvMassB2CorrectWithMinDeltaRJetB1HH->Draw();
	legendInvMassB2CorrectWithMinDeltaRJetB1ZZ->Draw();
	c28->Update();
	c28->SaveAs("InvMassB2CorrectWithMinDeltaRJetB1HHZZ.png");
	
	TLegend *legendInvMassB2IncorrectWithMinDeltaRJetB1HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2IncorrectWithMinDeltaRJetB1HH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2IncorrectWithMinDeltaRJetB1HH->SetTextColor(kBlue);
	legendInvMassB2IncorrectWithMinDeltaRJetB1HH->SetBorderSize(0);
	legendInvMassB2IncorrectWithMinDeltaRJetB1HH->SetFillStyle(0);
	legendInvMassB2IncorrectWithMinDeltaRJetB1HH->SetTextSize(0.043);
	TLegend *legendInvMassB2IncorrectWithMinDeltaRJetB1ZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2IncorrectWithMinDeltaRJetB1ZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2IncorrectWithMinDeltaRJetB1ZZ->SetTextColor(kGreen+1);
	legendInvMassB2IncorrectWithMinDeltaRJetB1ZZ->SetBorderSize(0);
	legendInvMassB2IncorrectWithMinDeltaRJetB1ZZ->SetFillStyle(0);
	legendInvMassB2IncorrectWithMinDeltaRJetB1ZZ->SetTextSize(0.043);
	TCanvas *c29 = new TCanvas();
	c29->cd();
	histInvMassB2IncorrectWithMinDeltaRJetB1ZZ->DrawNormalized("HIST");
	histInvMassB2IncorrectWithMinDeltaRJetB1HH->DrawNormalized("HIST same");
	histInvMassB2IncorrectWithMinDeltaRJetB1ZZ->DrawNormalized("HIST same");
	legendInvMassB2IncorrectWithMinDeltaRJetB1HH->Draw();
	legendInvMassB2IncorrectWithMinDeltaRJetB1ZZ->Draw();
	c29->Update();
	c29->SaveAs("InvMassB2IncorrectWithMinDeltaRJetB1HHZZ.png");
	
	TLegend *legendInvMassB1CorrectWithMinDeltaRJetsHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1CorrectWithMinDeltaRJetsHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1CorrectWithMinDeltaRJetsHH->SetTextColor(kBlue);
	legendInvMassB1CorrectWithMinDeltaRJetsHH->SetBorderSize(0);
	legendInvMassB1CorrectWithMinDeltaRJetsHH->SetFillStyle(0);
	legendInvMassB1CorrectWithMinDeltaRJetsHH->SetTextSize(0.043);
	TLegend *legendInvMassB1CorrectWithMinDeltaRJetsZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1CorrectWithMinDeltaRJetsZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1CorrectWithMinDeltaRJetsZZ->SetTextColor(kGreen+1);
	legendInvMassB1CorrectWithMinDeltaRJetsZZ->SetBorderSize(0);
	legendInvMassB1CorrectWithMinDeltaRJetsZZ->SetFillStyle(0);
	legendInvMassB1CorrectWithMinDeltaRJetsZZ->SetTextSize(0.043);
	TCanvas *c30 = new TCanvas();
	c30->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB1CorrectWithMinDeltaRJetsZZ->DrawNormalized("HIST");
	histInvMassB1CorrectWithMinDeltaRJetsHH->DrawNormalized("HIST same");
	histInvMassB1CorrectWithMinDeltaRJetsZZ->DrawNormalized("HIST same");
	legendInvMassB1CorrectWithMinDeltaRJetsHH->Draw();
	legendInvMassB1CorrectWithMinDeltaRJetsZZ->Draw();
	c30->Update();
	c30->SaveAs("InvMassB1CorrectWithMinDeltaRJetsHHZZ.png");
	
	TLegend *legendInvMassB1IncorrectWithMinDeltaRJetsHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB1IncorrectWithMinDeltaRJetsHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB1IncorrectWithMinDeltaRJetsHH->SetTextColor(kBlue);
	legendInvMassB1IncorrectWithMinDeltaRJetsHH->SetBorderSize(0);
	legendInvMassB1IncorrectWithMinDeltaRJetsHH->SetFillStyle(0);
	legendInvMassB1IncorrectWithMinDeltaRJetsHH->SetTextSize(0.043);
	TLegend *legendInvMassB1IncorrectWithMinDeltaRJetsZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB1IncorrectWithMinDeltaRJetsZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB1IncorrectWithMinDeltaRJetsZZ->SetTextColor(kGreen+1);
	legendInvMassB1IncorrectWithMinDeltaRJetsZZ->SetBorderSize(0);
	legendInvMassB1IncorrectWithMinDeltaRJetsZZ->SetFillStyle(0);
	legendInvMassB1IncorrectWithMinDeltaRJetsZZ->SetTextSize(0.043);
	TCanvas *c31 = new TCanvas();
	c31->cd();
	histInvMassB1IncorrectWithMinDeltaRJetsZZ->DrawNormalized("HIST");
	histInvMassB1IncorrectWithMinDeltaRJetsHH->DrawNormalized("HIST same");
	histInvMassB1IncorrectWithMinDeltaRJetsZZ->DrawNormalized("HIST same");
	legendInvMassB1IncorrectWithMinDeltaRJetsHH->Draw();
	legendInvMassB1IncorrectWithMinDeltaRJetsZZ->Draw();
	c31->Update();
	c31->SaveAs("InvMassB1IncorrectWithMinDeltaRJetsHHZZ.png");
	
	TLegend *legendInvMassB2CorrectWithMinDeltaRJetsHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2CorrectWithMinDeltaRJetsHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2CorrectWithMinDeltaRJetsHH->SetTextColor(kBlue);
	legendInvMassB2CorrectWithMinDeltaRJetsHH->SetBorderSize(0);
	legendInvMassB2CorrectWithMinDeltaRJetsHH->SetFillStyle(0);
	legendInvMassB2CorrectWithMinDeltaRJetsHH->SetTextSize(0.043);
	TLegend *legendInvMassB2CorrectWithMinDeltaRJetsZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2CorrectWithMinDeltaRJetsZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2CorrectWithMinDeltaRJetsZZ->SetTextColor(kGreen+1);
	legendInvMassB2CorrectWithMinDeltaRJetsZZ->SetBorderSize(0);
	legendInvMassB2CorrectWithMinDeltaRJetsZZ->SetFillStyle(0);
	legendInvMassB2CorrectWithMinDeltaRJetsZZ->SetTextSize(0.043);
	TCanvas *c32 = new TCanvas();
	c32->cd();
	//histDistanceZ2MinChiSquaredZZMassZZ->DrawNormalized("HIST");
	histInvMassB2CorrectWithMinDeltaRJetsZZ->DrawNormalized("HIST");
	histInvMassB2CorrectWithMinDeltaRJetsHH->DrawNormalized("HIST same");
	histInvMassB2CorrectWithMinDeltaRJetsZZ->DrawNormalized("HIST same");
	legendInvMassB2CorrectWithMinDeltaRJetsHH->Draw();
	legendInvMassB2CorrectWithMinDeltaRJetsZZ->Draw();
	c32->Update();
	c32->SaveAs("InvMassB2CorrectWithMinDeltaRJetsHHZZ.png");
	
	TLegend *legendInvMassB2IncorrectWithMinDeltaRJetsHH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendInvMassB2IncorrectWithMinDeltaRJetsHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendInvMassB2IncorrectWithMinDeltaRJetsHH->SetTextColor(kBlue);
	legendInvMassB2IncorrectWithMinDeltaRJetsHH->SetBorderSize(0);
	legendInvMassB2IncorrectWithMinDeltaRJetsHH->SetFillStyle(0);
	legendInvMassB2IncorrectWithMinDeltaRJetsHH->SetTextSize(0.043);
	TLegend *legendInvMassB2IncorrectWithMinDeltaRJetsZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendInvMassB2IncorrectWithMinDeltaRJetsZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendInvMassB2IncorrectWithMinDeltaRJetsZZ->SetTextColor(kGreen+1);
	legendInvMassB2IncorrectWithMinDeltaRJetsZZ->SetBorderSize(0);
	legendInvMassB2IncorrectWithMinDeltaRJetsZZ->SetFillStyle(0);
	legendInvMassB2IncorrectWithMinDeltaRJetsZZ->SetTextSize(0.043);
	TCanvas *c33 = new TCanvas();
	c33->cd();
	histInvMassB2IncorrectWithMinDeltaRJetsZZ->DrawNormalized("HIST");
	histInvMassB2IncorrectWithMinDeltaRJetsHH->DrawNormalized("HIST same");
	histInvMassB2IncorrectWithMinDeltaRJetsZZ->DrawNormalized("HIST same");
	legendInvMassB2IncorrectWithMinDeltaRJetsHH->Draw();
	legendInvMassB2IncorrectWithMinDeltaRJetsZZ->Draw();
	c33->Update();
	c33->SaveAs("InvMassB2IncorrectWithMinDeltaRJetsHHZZ.png");*/
	
	TLegend *legendThetaStarB1HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendThetaStarB1HH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendThetaStarB1HH->SetTextColor(kBlue);
	legendThetaStarB1HH->SetBorderSize(0);
	legendThetaStarB1HH->SetFillStyle(0);
	legendThetaStarB1HH->SetTextSize(0.043);
	TLegend *legendThetaStarB1ZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendThetaStarB1ZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendThetaStarB1ZZ->SetTextColor(kGreen+1);
	legendThetaStarB1ZZ->SetBorderSize(0);
	legendThetaStarB1ZZ->SetFillStyle(0);
	legendThetaStarB1ZZ->SetTextSize(0.043);
	TCanvas *c34 = new TCanvas();
	c34->cd();
	histThetaStarB1ZZ->DrawNormalized("HIST");
	histThetaStarB1HH->DrawNormalized("HIST same");
	histThetaStarB1ZZ->DrawNormalized("HIST same");
	legendThetaStarB1HH->Draw();
	legendThetaStarB1ZZ->Draw();
	c34->Update();
	c34->SaveAs("ThetaStarB1HHZZ.png");
	
	TLegend *legendThetaStarB2HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendThetaStarB2HH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendThetaStarB2HH->SetTextColor(kBlue);
	legendThetaStarB2HH->SetBorderSize(0);
	legendThetaStarB2HH->SetFillStyle(0);
	legendThetaStarB2HH->SetTextSize(0.043);
	TLegend *legendThetaStarB2ZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendThetaStarB2ZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendThetaStarB2ZZ->SetTextColor(kGreen+1);
	legendThetaStarB2ZZ->SetBorderSize(0);
	legendThetaStarB2ZZ->SetFillStyle(0);
	legendThetaStarB2ZZ->SetTextSize(0.043);
	TCanvas *c35 = new TCanvas();
	c35->cd();
	histThetaStarB2ZZ->DrawNormalized("HIST");
	histThetaStarB2HH->DrawNormalized("HIST same");
	histThetaStarB2ZZ->DrawNormalized("HIST same");
	legendThetaStarB2HH->Draw();
	legendThetaStarB2ZZ->Draw();
	c35->Update();
	c35->SaveAs("ThetaStarB2HHZZ.png");
	
	TLegend *legendThetaStarB3HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendThetaStarB3HH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendThetaStarB3HH->SetTextColor(kBlue);
	legendThetaStarB3HH->SetBorderSize(0);
	legendThetaStarB3HH->SetFillStyle(0);
	legendThetaStarB3HH->SetTextSize(0.043);
	TLegend *legendThetaStarB3ZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendThetaStarB3ZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendThetaStarB3ZZ->SetTextColor(kGreen+1);
	legendThetaStarB3ZZ->SetBorderSize(0);
	legendThetaStarB3ZZ->SetFillStyle(0);
	legendThetaStarB3ZZ->SetTextSize(0.043);
	TCanvas *c36 = new TCanvas();
	c36->cd();
	histThetaStarB3ZZ->DrawNormalized("HIST");
	histThetaStarB3HH->DrawNormalized("HIST same");
	histThetaStarB3ZZ->DrawNormalized("HIST same");
	legendThetaStarB3HH->Draw();
	legendThetaStarB3ZZ->Draw();
	c36->Update();
	c36->SaveAs("ThetaStarB3HHZZ.png");
	
	TLegend *legendThetaStarB4HH = new TLegend(0.65, 0.65, 0.8, 0.8);
	legendThetaStarB4HH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendThetaStarB4HH->SetTextColor(kBlue);
	legendThetaStarB4HH->SetBorderSize(0);
	legendThetaStarB4HH->SetFillStyle(0);
	legendThetaStarB4HH->SetTextSize(0.043);
	TLegend *legendThetaStarB4ZZ = new TLegend(0.65, 0.6, 0.8, 0.75);
	legendThetaStarB4ZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendThetaStarB4ZZ->SetTextColor(kGreen+1);
	legendThetaStarB4ZZ->SetBorderSize(0);
	legendThetaStarB4ZZ->SetFillStyle(0);
	legendThetaStarB4ZZ->SetTextSize(0.043);
	TCanvas *c37 = new TCanvas();
	c37->cd();
	histThetaStarB4ZZ->DrawNormalized("HIST");
	histThetaStarB4HH->DrawNormalized("HIST same");
	histThetaStarB4ZZ->DrawNormalized("HIST same");
	legendThetaStarB4HH->Draw();
	legendThetaStarB4ZZ->Draw();
	c37->Update();
	c37->SaveAs("ThetaStarB4HHZZ.png");
	
	TLegend *legendThrustHH = new TLegend(0.25, 0.65, 0.4, 0.8);
	legendThrustHH->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{HH}", "");
	legendThrustHH->SetTextColor(kBlue);
	legendThrustHH->SetBorderSize(0);
	legendThrustHH->SetFillStyle(0);
	legendThrustHH->SetTextSize(0.043);
	TLegend *legendThrustZZ = new TLegend(0.25, 0.6, 0.4, 0.75);
	legendThrustZZ->AddEntry(histDistanceZ2MinChiSquaredZZMassZZ, "#gamma#gamma#rightarrow#it{ZZ}", "");
	legendThrustZZ->SetTextColor(kGreen+1);
	legendThrustZZ->SetBorderSize(0);
	legendThrustZZ->SetFillStyle(0);
	legendThrustZZ->SetTextSize(0.043);
	TCanvas *c38 = new TCanvas();
	c38->cd();
	histThrustZZ->DrawNormalized("HIST");
	histThrustHH->DrawNormalized("HIST same");
	histThrustZZ->DrawNormalized("HIST same");
	legendThrustHH->Draw();
	legendThrustZZ->Draw();
	c38->Update();
	c38->SaveAs("ThrustHHZZ.png");

}















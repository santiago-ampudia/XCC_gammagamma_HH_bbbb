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
 
 
 using namespace TMVA;
 
void FSRTMVAClassificationApplicationHHbbbbHelper( TString myMethodList, int topology, int nBack, vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, vector<double>& BDTqqXOutput, vector<double>& BDTqqqqXOutput, vector<double>& BDTqqHXOutput, int& sizeTree, int nbin, TH1F& histOutput, TH1F& histOutputW, double weight, string inputMethod, TString NNVars, string rtdCut, string sampleName, string varVersion, string preselection)
{
	if(topology == 1)
	{
		cout<<"Evaluating HH on ";
		if(nBack == 0) cout<<"NNqq"<<endl;
	 	if(nBack == 1) cout<<"NNttbar"<<endl;
	 	if(nBack == 2) cout<<"NNZZ"<<endl;
	 	if(nBack == 3) cout<<"NNWW"<<endl;
		if(nBack == 4) cout<<"NNqqX"<<endl;
		if(nBack == 5) cout<<"NNqqqqX"<<endl;
		if(nBack == 6) cout<<"NNqqHX"<<endl;
 	}
 	if(topology == 2)
	{
		cout<<"Evaluating qq on ";
		if(nBack == 0) cout<<"NNqq"<<endl;
	 	if(nBack == 1) cout<<"NNttbar"<<endl;
	 	if(nBack == 2) cout<<"NNZZ"<<endl;
	 	if(nBack == 3) cout<<"NNWW"<<endl;
		if(nBack == 4) cout<<"NNqqX"<<endl;
		if(nBack == 5) cout<<"NNqqqqX"<<endl;
		if(nBack == 6) cout<<"NNqqHX"<<endl;
 	}
 	if(topology == 3)
	{
		cout<<"Evaluating ttbar on ";
		if(nBack == 0) cout<<"NNqq"<<endl;
	 	if(nBack == 1) cout<<"NNttbar"<<endl;
	 	if(nBack == 2) cout<<"NNZZ"<<endl;
	 	if(nBack == 3) cout<<"NNWW"<<endl;
		if(nBack == 4) cout<<"NNqqX"<<endl;
		if(nBack == 5) cout<<"NNqqqqX"<<endl;
		if(nBack == 6) cout<<"NNqqHX"<<endl;
 	}
 	if(topology == 4)
	{
		cout<<"Evaluating ZZ on ";
		if(nBack == 0) cout<<"NNqq"<<endl;
	 	if(nBack == 1) cout<<"NNttbar"<<endl;
	 	if(nBack == 2) cout<<"NNZZ"<<endl;
	 	if(nBack == 3) cout<<"NNWW"<<endl;
		if(nBack == 4) cout<<"NNqqX"<<endl;
		if(nBack == 5) cout<<"NNqqqqX"<<endl;
		if(nBack == 6) cout<<"NNqqHX"<<endl;
 	}
 	if(topology == 5)
	{
		cout<<"Evaluating WW on ";
		if(nBack == 0) cout<<"NNqq"<<endl;
	 	if(nBack == 1) cout<<"NNttbar"<<endl;
	 	if(nBack == 2) cout<<"NNZZ"<<endl;
	 	if(nBack == 3) cout<<"NNWW"<<endl;
		if(nBack == 4) cout<<"NNqqX"<<endl;
		if(nBack == 5) cout<<"NNqqqqX"<<endl;
		if(nBack == 6) cout<<"NNqqHX"<<endl;
 	}
	if(topology == 6)
	{
		cout<<"Evaluating qqX on ";
		if(nBack == 0) cout<<"NNqq"<<endl;
	 	if(nBack == 1) cout<<"NNttbar"<<endl;
	 	if(nBack == 2) cout<<"NNZZ"<<endl;
	 	if(nBack == 3) cout<<"NNWW"<<endl;
		if(nBack == 4) cout<<"NNqqX"<<endl;
		if(nBack == 5) cout<<"NNqqqqX"<<endl;
		if(nBack == 6) cout<<"NNqqHX"<<endl;
 	}
	if(topology == 7)
	{
		cout<<"Evaluating qqqqX on ";
		if(nBack == 0) cout<<"NNqq"<<endl;
	 	if(nBack == 1) cout<<"NNttbar"<<endl;
	 	if(nBack == 2) cout<<"NNZZ"<<endl;
	 	if(nBack == 3) cout<<"NNWW"<<endl;
		if(nBack == 4) cout<<"NNqqX"<<endl;
		if(nBack == 5) cout<<"NNqqqqX"<<endl;
		if(nBack == 6) cout<<"NNqqHX"<<endl;
 	}
	if(topology == 8)
	{
		cout<<"Evaluating qqHX on ";
		if(nBack == 0) cout<<"NNqq"<<endl;
	 	if(nBack == 1) cout<<"NNttbar"<<endl;
	 	if(nBack == 2) cout<<"NNZZ"<<endl;
	 	if(nBack == 3) cout<<"NNWW"<<endl;
		if(nBack == 4) cout<<"NNqqX"<<endl;
		if(nBack == 5) cout<<"NNqqqqX"<<endl;
		if(nBack == 6) cout<<"NNqqHX"<<endl;
 	}
		
	
 	
 	//---------------------------------------------------------------
    // This loads the library
    TMVA::Tools::Instance();
 
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
 
    Use[inputMethod] = 1;
    /*// Cut optimisation
    Use["Cuts"]            = 0;
    Use["CutsD"]           = 0;
    Use["CutsPCA"]         = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;
    //
    // 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"]      = 0;
    Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Use["LikelihoodKDE"]   = 0;
    Use["LikelihoodMIX"]   = 0;
    //
    // Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 0;
    Use["PDERSD"]          = 0;
    Use["PDERSPCA"]        = 0;
    Use["PDEFoam"]         = 0;
    Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Use["KNN"]             = 0; // k-nearest neighbour method
    //
    // Linear Discriminant Analysis
    Use["LD"]              = 0; // Linear Discriminant identical to Fisher
    Use["Fisher"]          = 0;
    Use["FisherG"]         = 0;
    Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Use["HMatrix"]         = 0;
    //
    // Function Discriminant analysis
    Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
    Use["FDA_SA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    Use["FDA_MCMT"]        = 0;
    //
    // Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
    Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
    //
    // Support Vector Machine
    Use["SVM"]             = 0;
    //
    // Boosted Decision Trees
    Use["BDT"]             = 0; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
    //
    // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 0;
    // ---------------------------------------------------------------
    Use["Plugin"]          = 0;
    Use["Category"]        = 0;
    Use["SVM_Gauss"]       = 0;
    Use["SVM_Poly"]        = 0;
    Use["SVM_Lin"]         = 0;*/
 
   // std::cout << std::endl;
   // std::cout << "==> Start TMVAClassificationApplication" << std::endl;
 
    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
       for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
 
       std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
       for (UInt_t i=0; i<mlist.size(); i++) {
          std::string regMethod(mlist[i]);
 
          if (Use.find(regMethod) == Use.end()) {
             std::cout << "Method \"" << regMethod
                       << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
             for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
                std::cout << it->first << " ";
             }
             std::cout << std::endl;
             return;
          }
          Use[regMethod] = 1;
       }
    }
 
    // --------------------------------------------------------------------------------------------------
 
    // Create the Reader object
 
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
 
    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    Float_t aplanarity, cosThetaB1, cosThetaB2, cosThetaB3, cosThetaB4, invMassB1, invMassB2, jetB1M, jetB2M, jetB3M, jetB4M, jetB1Pt, jetB2Pt, jetB3Pt, jetB4Pt, minJetM, sphericity, sumPt, nParticles, minConstSize, jetNObjects, minJetNObjects, invMassB1AntiKt, invMassB2AntiKt, nJetsAntiKt, invMassB11Best, invMassB21Best, invMassB12Best, invMassB22Best, invMassB13Best, invMassB23Best, invMassB14Best, invMassB24Best, invMassB15Best, invMassB25Best, invMassB16Best, invMassB26Best, invMassB17Best, invMassB27Best, invMassB18Best, invMassB28Best, exclYmerge12, exclYmerge23, exclYmerge34, exclYmerge45, exclYmerge56, invMassZZ1, invMassZZ2, thrust, boostB1, boostB2, boostB3, boostB4, boostSystem, missingET;
    
    if(varVersion == "All1V")
   	{
   		reader->AddVariable( "aplanarity", &aplanarity );
		reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		reader->AddVariable( "jetB4Pt", &jetB4Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		reader->AddVariable( "sumPt", &sumPt );
	}
	else if(varVersion == "All2V")
   	{
   		reader->AddVariable( "aplanarity", &aplanarity );
		reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		reader->AddVariable( "jetB4Pt", &jetB4Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		reader->AddVariable( "sumPt", &sumPt );
		reader->AddVariable( "nParticles", &nParticles );
		reader->AddVariable( "minConstSize", &minConstSize );
    }
    else if(varVersion == "All3V")
   	{
   		reader->AddVariable( "aplanarity", &aplanarity );
		reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		reader->AddVariable( "jetB4Pt", &jetB4Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		reader->AddVariable( "sumPt", &sumPt );
		reader->AddVariable( "nParticles", &nParticles );
		reader->AddVariable( "minConstSize", &minConstSize );
		reader->AddVariable( "jetNObjects", &jetNObjects );
		reader->AddVariable( "minJetNObjects", &minJetNObjects );
   	}
    if(varVersion == "All4V")
    {
		reader->AddVariable( "aplanarity", &aplanarity );
		reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		if(nBack!=3) reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		reader->AddVariable( "jetB4Pt", &jetB4Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		reader->AddVariable( "sumPt", &sumPt );
		//reader->AddVariable( "nParticles", &nParticles );
		//reader->AddVariable( "minConstSize", &minConstSize );
		//reader->AddVariable( "jetNObjects", &jetNObjects );
		//reader->AddVariable( "minJetNObjects", &minJetNObjects );
	    reader->AddVariable("invMassB1AntiKt", &invMassB1AntiKt);
	    reader->AddVariable("invMassB2AntiKt", &invMassB2AntiKt);
		reader->AddVariable("nJetsAntiKt", &nJetsAntiKt);
		reader->AddVariable("invMassB11Best", &invMassB11Best);
		reader->AddVariable("invMassB21Best", &invMassB21Best);
		reader->AddVariable("invMassB12Best", &invMassB12Best);
		reader->AddVariable("invMassB22Best", &invMassB22Best);
		reader->AddVariable("invMassB13Best", &invMassB13Best);
		reader->AddVariable("invMassB23Best", &invMassB23Best);
		reader->AddVariable("invMassB14Best", &invMassB14Best);
		reader->AddVariable("invMassB24Best", &invMassB24Best);
		reader->AddVariable("invMassB15Best", &invMassB15Best);
		reader->AddVariable("invMassB25Best", &invMassB25Best);
		reader->AddVariable("invMassB16Best", &invMassB16Best);
		reader->AddVariable("invMassB26Best", &invMassB26Best);
		reader->AddVariable("invMassB17Best", &invMassB17Best);
		reader->AddVariable("invMassB27Best", &invMassB27Best);
		reader->AddVariable("invMassB18Best", &invMassB18Best);
		reader->AddVariable("invMassB28Best", &invMassB28Best);
	}
	if(varVersion == "All5V")
    {
		if(nBack != 0 && nBack != 3) reader->AddVariable( "aplanarity", &aplanarity );
		reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		if(nBack != 3) reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable( "jetB4Pt", &jetB4Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		reader->AddVariable( "sumPt", &sumPt );
		reader->AddVariable("invMassB1AntiKt", &invMassB1AntiKt);
		reader->AddVariable("invMassB2AntiKt", &invMassB2AntiKt);
		reader->AddVariable("nJetsAntiKt", &nJetsAntiKt);
		reader->AddVariable("invMassB11Best", &invMassB11Best);
		reader->AddVariable("invMassB21Best", &invMassB21Best);
		reader->AddVariable("invMassB12Best", &invMassB12Best);
		reader->AddVariable("invMassB22Best", &invMassB22Best);
		reader->AddVariable("invMassB13Best", &invMassB13Best);
		reader->AddVariable("invMassB23Best", &invMassB23Best);
		reader->AddVariable("invMassB14Best", &invMassB14Best);
		reader->AddVariable("invMassB24Best", &invMassB24Best);
		reader->AddVariable("invMassB15Best", &invMassB15Best);
		if(nBack != 0 && nBack != 2) reader->AddVariable("invMassB25Best", &invMassB25Best);
		reader->AddVariable("invMassB16Best", &invMassB16Best);
		if(nBack != 1) reader->AddVariable("invMassB26Best", &invMassB26Best);
		reader->AddVariable("invMassB17Best", &invMassB17Best);
		reader->AddVariable("invMassB27Best", &invMassB27Best);
		reader->AddVariable("invMassB18Best", &invMassB18Best);
		if(nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB28Best", &invMassB28Best);
	}
	if(varVersion == "All6V")
    {
		if(nBack != 0 && nBack != 2 && nBack != 3) reader->AddVariable( "aplanarity", &aplanarity );
		reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		if(nBack != 0) reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		if(nBack != 3) reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable( "jetB4Pt", &jetB4Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		if(nBack != 0 && nBack != 3) reader->AddVariable( "sumPt", &sumPt );
		reader->AddVariable("invMassB1AntiKt", &invMassB1AntiKt);
		reader->AddVariable("invMassB2AntiKt", &invMassB2AntiKt);
		if(nBack != 2) reader->AddVariable("nJetsAntiKt", &nJetsAntiKt);
		reader->AddVariable("invMassB11Best", &invMassB11Best);
		reader->AddVariable("invMassB21Best", &invMassB21Best);
		reader->AddVariable("invMassB12Best", &invMassB12Best);
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB22Best", &invMassB22Best);
		if(nBack != 1) reader->AddVariable("invMassB13Best", &invMassB13Best);
		if(nBack != 0 && nBack != 3) reader->AddVariable("invMassB23Best", &invMassB23Best);
		if(nBack != 1) reader->AddVariable("invMassB14Best", &invMassB14Best);
		if(nBack != 2) reader->AddVariable("invMassB24Best", &invMassB24Best);
		reader->AddVariable("invMassB15Best", &invMassB15Best);
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB25Best", &invMassB25Best);
		if(nBack != 1) reader->AddVariable("invMassB16Best", &invMassB16Best);
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB26Best", &invMassB26Best);
		if(nBack != 1) reader->AddVariable("invMassB17Best", &invMassB17Best);
		if(nBack != 2 && nBack != 3) reader->AddVariable("invMassB27Best", &invMassB27Best);
		reader->AddVariable("invMassB18Best", &invMassB18Best);
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB28Best", &invMassB28Best);
	}
	if(varVersion == "All7V")
    {
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable( "aplanarity", &aplanarity );
		if(nBack != 3) reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		if(nBack != 0) reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		if(nBack != 0) reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		if(nBack != 0) reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		if(nBack != 3) reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable( "jetB4Pt", &jetB4Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		if(nBack != 0 && nBack != 2 && nBack != 3) reader->AddVariable( "sumPt", &sumPt );
		reader->AddVariable("invMassB1AntiKt", &invMassB1AntiKt);
		reader->AddVariable("invMassB2AntiKt", &invMassB2AntiKt);
		if(nBack != 2 && nBack != 3) reader->AddVariable("nJetsAntiKt", &nJetsAntiKt);
		reader->AddVariable("invMassB11Best", &invMassB11Best);
		reader->AddVariable("invMassB21Best", &invMassB21Best);
		if(nBack != 1) reader->AddVariable("invMassB12Best", &invMassB12Best);
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB22Best", &invMassB22Best);
		if(nBack != 1) reader->AddVariable("invMassB13Best", &invMassB13Best);
		if(nBack != 0 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB23Best", &invMassB23Best);
		if(nBack != 1) reader->AddVariable("invMassB14Best", &invMassB14Best);
		if(nBack != 2) reader->AddVariable("invMassB24Best", &invMassB24Best);
		reader->AddVariable("invMassB15Best", &invMassB15Best);
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB25Best", &invMassB25Best);
		if(nBack != 1) reader->AddVariable("invMassB16Best", &invMassB16Best);
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB26Best", &invMassB26Best);
		if(nBack != 1) reader->AddVariable("invMassB17Best", &invMassB17Best);
		if(nBack != 2 && nBack != 3) reader->AddVariable("invMassB27Best", &invMassB27Best);
		reader->AddVariable("invMassB18Best", &invMassB18Best);
		if(nBack != 0 && nBack != 1 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB28Best", &invMassB28Best);
	}
	if(varVersion == "All8V")
    {
		if(nBack != 3) reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		if(nBack != 0) reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		if(nBack != 0) reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		if(nBack != 0) reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		if(nBack != 3) reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		if(nBack != 0 && nBack != 2 && nBack != 3) reader->AddVariable( "sumPt", &sumPt );
		reader->AddVariable("invMassB1AntiKt", &invMassB1AntiKt);
		reader->AddVariable("invMassB2AntiKt", &invMassB2AntiKt);
		if(nBack != 2 && nBack != 3) reader->AddVariable("nJetsAntiKt", &nJetsAntiKt);
		reader->AddVariable("invMassB11Best", &invMassB11Best);
		reader->AddVariable("invMassB21Best", &invMassB21Best);
		if(nBack != 1) reader->AddVariable("invMassB12Best", &invMassB12Best);
		if(nBack != 1) reader->AddVariable("invMassB13Best", &invMassB13Best);
		if(nBack != 0 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB23Best", &invMassB23Best);
		if(nBack != 1) reader->AddVariable("invMassB14Best", &invMassB14Best);
		if(nBack != 2) reader->AddVariable("invMassB24Best", &invMassB24Best);
		reader->AddVariable("invMassB15Best", &invMassB15Best);
		if(nBack != 1) reader->AddVariable("invMassB16Best", &invMassB16Best);
		if(nBack != 1) reader->AddVariable("invMassB17Best", &invMassB17Best);
		if(nBack != 2 && nBack != 3) reader->AddVariable("invMassB27Best", &invMassB27Best);
		reader->AddVariable("invMassB18Best", &invMassB18Best);
		reader->AddVariable("exclYmerge12", &exclYmerge12);
		reader->AddVariable("exclYmerge23", &exclYmerge23);
		reader->AddVariable("exclYmerge34", &exclYmerge34);
		reader->AddVariable("exclYmerge45", &exclYmerge45);
		reader->AddVariable("exclYmerge56", &exclYmerge56);
		reader->AddVariable("invMassZZ1", &invMassZZ1);
		reader->AddVariable("invMassZZ2", &invMassZZ2);
		reader->AddVariable("thrust", &thrust);
	}
	if(varVersion == "All9V")
	{
		if(nBack != 3) reader->AddVariable( "cosThetaB1", &cosThetaB1 );
		if(nBack != 0) reader->AddVariable( "cosThetaB2", &cosThetaB2 );
		if(nBack != 0) reader->AddVariable( "cosThetaB3", &cosThetaB3 );
		if(nBack != 0) reader->AddVariable( "cosThetaB4", &cosThetaB4 );
		reader->AddVariable( "invMassB1", &invMassB1 );
		reader->AddVariable( "invMassB2", &invMassB2 );
		reader->AddVariable( "jetB1M", &jetB1M );
		if(nBack != 3) reader->AddVariable( "jetB2M", &jetB2M );
		reader->AddVariable( "jetB3M", &jetB3M );
		reader->AddVariable( "jetB4M", &jetB4M );
		reader->AddVariable( "jetB1Pt", &jetB1Pt );
		reader->AddVariable( "jetB2Pt", &jetB2Pt );
		reader->AddVariable( "jetB3Pt", &jetB3Pt );
		reader->AddVariable( "minJetM", &minJetM );
		reader->AddVariable( "sphericity", &sphericity );
		if(nBack != 0 && nBack != 2 && nBack != 3) reader->AddVariable( "sumPt", &sumPt );
		reader->AddVariable("invMassB1AntiKt", &invMassB1AntiKt);
		reader->AddVariable("invMassB2AntiKt", &invMassB2AntiKt);
		if(nBack != 2 && nBack != 3) reader->AddVariable("nJetsAntiKt", &nJetsAntiKt);
		reader->AddVariable("invMassB11Best", &invMassB11Best);
		reader->AddVariable("invMassB21Best", &invMassB21Best);
		if(nBack != 1) reader->AddVariable("invMassB12Best", &invMassB12Best);
		if(nBack != 1) reader->AddVariable("invMassB13Best", &invMassB13Best);
		if(nBack != 0 && nBack != 2 && nBack != 3) reader->AddVariable("invMassB23Best", &invMassB23Best);
		if(nBack != 1) reader->AddVariable("invMassB14Best", &invMassB14Best);
		if(nBack != 2) reader->AddVariable("invMassB24Best", &invMassB24Best);
		reader->AddVariable("invMassB15Best", &invMassB15Best);
		if(nBack != 1) reader->AddVariable("invMassB16Best", &invMassB16Best);
		if(nBack != 1) reader->AddVariable("invMassB17Best", &invMassB17Best);
		if(nBack != 2 && nBack != 3) reader->AddVariable("invMassB27Best", &invMassB27Best);
		reader->AddVariable("invMassB18Best", &invMassB18Best);
		reader->AddVariable("exclYmerge12", &exclYmerge12);
		reader->AddVariable("exclYmerge23", &exclYmerge23);
		reader->AddVariable("exclYmerge34", &exclYmerge34);
		reader->AddVariable("exclYmerge45", &exclYmerge45);
		reader->AddVariable("exclYmerge56", &exclYmerge56);
		reader->AddVariable("invMassZZ1", &invMassZZ1);
		reader->AddVariable("invMassZZ2", &invMassZZ2);
		reader->AddVariable("thrust", &thrust);
		reader->AddVariable("boostB1", &boostB1);
		reader->AddVariable("boostB2", &boostB2);
		reader->AddVariable("boostB3", &boostB3);
		reader->AddVariable("boostB4", &boostB4);
		reader->AddVariable("boostSystem", &boostSystem);
		reader->AddVariable("missingET", &missingET);
	}




 
 
    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //
    string inputSText = "analysis/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputS = new TFile(inputSText.c_str());
    string inputBqqText = "analysis/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputBqq = new TFile(inputBqqText.c_str());
    string inputBttbarText = "analysis/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputBtt = new TFile(inputBttbarText.c_str());
    string inputBZZText = "analysis/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputBZZ = new TFile(inputBZZText.c_str());
    string inputBWWText = "analysis/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputBWW = new TFile(inputBWWText.c_str());
	string inputBqqXText = "analysis/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputBqqX = new TFile(inputBqqXText.c_str());
	string inputBqqqqXText = "analysis/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputBqqqqX = new TFile(inputBqqqqXText.c_str());
	string inputBqqHXText = "analysis/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputBqqHX = new TFile(inputBqqHXText.c_str());
    
    std::cout << "--- TMVAClassificationApplication   : Using input file for signal: " << inputS->GetName() << std::endl;
 
    // Event loop
 
    // Prepare the event tree
    // - Here the variable names have to corresponds to your tree
    // - You can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //
    
    TTree* theTree;
    
    if(topology == 1) theTree = (TTree*)inputS->Get("TreeSTest");
    else if(topology == 2) theTree = (TTree*)inputBqq->Get("TreeBqqTest");
    else if(topology == 3) theTree = (TTree*)inputBtt->Get("TreeBttTest");
    else if(topology == 4) theTree = (TTree*)inputBZZ->Get("TreeBZZTest");
    else if(topology == 5) theTree = (TTree*)inputBWW->Get("TreeBWWTest");
	else if(topology == 6) theTree = (TTree*)inputBqqX->Get("TreeBqqXTest");
	else if(topology == 7) theTree = (TTree*)inputBqqqqX->Get("TreeBqqqqXTest");
	else if(topology == 8) theTree = (TTree*)inputBqqHX->Get("TreeBqqHXTest");
    
	theTree->SetBranchAddress( "aplanarity", &aplanarity );
	theTree->SetBranchAddress( "cosThetaB1", &cosThetaB1 );
	theTree->SetBranchAddress( "cosThetaB2", &cosThetaB2 );
	theTree->SetBranchAddress( "cosThetaB3", &cosThetaB3 );
	theTree->SetBranchAddress( "cosThetaB4", &cosThetaB4 );
	theTree->SetBranchAddress( "invMassB1", &invMassB1 );
	theTree->SetBranchAddress( "invMassB2", &invMassB2 );
	theTree->SetBranchAddress( "jetB1M", &jetB1M );
	theTree->SetBranchAddress( "jetB2M", &jetB2M );
	theTree->SetBranchAddress( "jetB3M", &jetB3M );
	theTree->SetBranchAddress( "jetB4M", &jetB4M );
	theTree->SetBranchAddress( "jetB1Pt", &jetB1Pt );
	theTree->SetBranchAddress( "jetB2Pt", &jetB2Pt );
	theTree->SetBranchAddress( "jetB3Pt", &jetB3Pt );
	theTree->SetBranchAddress( "jetB4Pt", &jetB4Pt );
	theTree->SetBranchAddress( "minJetM", &minJetM );
	theTree->SetBranchAddress( "sphericity", &sphericity );
	theTree->SetBranchAddress( "sumPt", &sumPt );
	//theTree->SetBranchAddress( "nParticles", &nParticles );
	//theTree->SetBranchAddress( "minConstSize", &minConstSize );
	//theTree->SetBranchAddress( "jetNObjects", &jetNObjects );
	//theTree->SetBranchAddress( "minJetNObjects", &minJetNObjects );
	theTree->SetBranchAddress("invMassB1AntiKt", &invMassB1AntiKt);
	theTree->SetBranchAddress("invMassB2AntiKt", &invMassB2AntiKt);
	theTree->SetBranchAddress("nJetsAntiKt", &nJetsAntiKt);
	theTree->SetBranchAddress("invMassB11Best", &invMassB11Best);
	theTree->SetBranchAddress("invMassB21Best", &invMassB21Best);
	theTree->SetBranchAddress("invMassB12Best", &invMassB12Best);
	theTree->SetBranchAddress("invMassB22Best", &invMassB22Best);
	theTree->SetBranchAddress("invMassB13Best", &invMassB13Best);
	theTree->SetBranchAddress("invMassB23Best", &invMassB23Best);
	theTree->SetBranchAddress("invMassB14Best", &invMassB14Best);
	theTree->SetBranchAddress("invMassB24Best", &invMassB24Best);
	theTree->SetBranchAddress("invMassB15Best", &invMassB15Best);
	theTree->SetBranchAddress("invMassB25Best", &invMassB25Best);
	theTree->SetBranchAddress("invMassB16Best", &invMassB16Best);
	theTree->SetBranchAddress("invMassB26Best", &invMassB26Best);
	theTree->SetBranchAddress("invMassB17Best", &invMassB17Best);
	theTree->SetBranchAddress("invMassB27Best", &invMassB27Best);
	theTree->SetBranchAddress("invMassB18Best", &invMassB18Best);
	theTree->SetBranchAddress("invMassB28Best", &invMassB28Best);
	theTree->SetBranchAddress("exclYmerge12", &exclYmerge12);
	theTree->SetBranchAddress("exclYmerge23", &exclYmerge23);
	theTree->SetBranchAddress("exclYmerge34", &exclYmerge34);
	theTree->SetBranchAddress("exclYmerge45", &exclYmerge45);
	theTree->SetBranchAddress("exclYmerge56", &exclYmerge56);
	theTree->SetBranchAddress("invMassZZ1", &invMassZZ1);
	theTree->SetBranchAddress("invMassZZ2", &invMassZZ2);
	theTree->SetBranchAddress("thrust", &thrust);
	theTree->SetBranchAddress("boostB1", &boostB1);
	theTree->SetBranchAddress("boostB2", &boostB2);
	theTree->SetBranchAddress("boostB3", &boostB3);
	theTree->SetBranchAddress("boostB4", &boostB4);
	theTree->SetBranchAddress("boostSystem", &boostSystem);
	theTree->SetBranchAddress("missingET", &missingET);



    // Efficiency calculator for cut method
    Int_t    nSelCutsGA = 0;
    Double_t effS       = 0.7;
 
    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
 
    if(topology == 1) std::cout << "--- Processing: " << theTree->GetEntries() << " HH events" << std::endl;
    if(topology == 2) std::cout << "--- Processing: " << theTree->GetEntries() << " qq events" << std::endl;
    if(topology == 3) std::cout << "--- Processing: " << theTree->GetEntries() << " ttbar events" << std::endl;
    if(topology == 4) std::cout << "--- Processing: " << theTree->GetEntries() << " ZZ events" << std::endl;
    if(topology == 5) std::cout << "--- Processing: " << theTree->GetEntries() << " WW events" << std::endl;
	if(topology == 6) std::cout << "--- Processing: " << theTree->GetEntries() << " qqX events" << std::endl;
	if(topology == 7) std::cout << "--- Processing: " << theTree->GetEntries() << " qqqqX events" << std::endl;
	if(topology == 8) std::cout << "--- Processing: " << theTree->GetEntries() << " qqHX events" << std::endl;
    
    TStopwatch sw;
    sw.Start();
   
	   /*if(nBack == 0) cout<<"Starting application of BDTqq!"<<endl;
	   if(nBack == 1) cout<<"Starting application of BDTttbar!"<<endl;
	   if(nBack == 2) cout<<"Starting application of BDTZZ!"<<endl;
	   if(nBack == 3) cout<<"Starting application of BDTWW!"<<endl;*/
	   
	    // Book the MVA methods
	    TString dir;
	    if(nBack == 0) dir    = TString("analysis/datasetqq") + TString(rtdCut) + TString(NNVars) + TString(sampleName) + TString("/weights/");
	    else if(nBack == 1) dir    = TString("analysis/datasetttbar") + TString(rtdCut) + TString(NNVars) + TString(sampleName) + TString("/weights/");
	    else if(nBack == 2) dir    = TString("analysis/datasetZZ") + TString(rtdCut) + TString(NNVars) + TString(sampleName) + TString("/weights/");
	    else if(nBack == 3) dir    = TString("analysis/datasetWW") + TString(rtdCut) + TString(NNVars) + TString(sampleName) + TString("/weights/");
		else if(nBack == 4) dir    = TString("analysis/datasetqqX") + TString(rtdCut) + TString(NNVars) + TString(sampleName) + TString("/weights/");
		else if(nBack == 5) dir    = TString("analysis/datasetqqqqX") + TString(rtdCut) + TString(NNVars) + TString(sampleName) + TString("/weights/");
		else if(nBack == 6) dir    = TString("analysis/datasetqqHX") + TString(rtdCut) + TString(NNVars) + TString(sampleName) + TString("/weights/");

	    TString prefix = "TMVAClassification";
	    // Book method(s)
	    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
	       if (it->second) {
		  TString methodName = TString(it->first) + TString(" method");
		  TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
		  reader->BookMVA( methodName, weightfile );
	       }
	    }
    	    double sum=0;
	    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	 
	       //if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
	 
	       theTree->GetEntry(ievt);
	       
	       TString methodName = TString(inputMethod) + TString(" method");
	       
	        sum+=reader->EvaluateMVA( methodName );
	       
	       ///Printing for qualitative feeling
	       //if(topology == 1 && ievt == 0) cout<<"POS 0: "<<reader->EvaluateMVA( methodName )<<endl;
	       //if(topology == 1 && ievt == 70) cout<<"POS 70: "<<reader->EvaluateMVA( methodName )<<endl;
	       
	       ////Filling out arrays for plane of NN outputs
	       if(nBack == 0) BDTqqOutput.push_back(reader->EvaluateMVA(methodName)); 
	       if(nBack == 1) BDTttbarOutput.push_back(reader->EvaluateMVA(methodName));
	       if(nBack == 2) BDTZZOutput.push_back(reader->EvaluateMVA(methodName));
	       if(nBack == 3) BDTWWOutput.push_back(reader->EvaluateMVA(methodName));
		   if(nBack == 4) BDTqqXOutput.push_back(reader->EvaluateMVA(methodName));
		   if(nBack == 5) BDTqqqqXOutput.push_back(reader->EvaluateMVA(methodName));
		   if(nBack == 6) BDTqqHXOutput.push_back(reader->EvaluateMVA(methodName));		
	       
	       /////Filling out histograms 
	       histOutput.Fill(reader->EvaluateMVA(methodName));
	       histOutputW.Fill(reader->EvaluateMVA(methodName), weight);
	       
	       sizeTree = ievt+1;
	       
	    }
    
    
    // Get elapsed time
    sw.Stop();
    //std::cout << "--- End of event loop: "; sw.Print();
 
    // Get efficiency for cuts classifier
    if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                 << " (for a required signal efficiency of " << effS << ")" << std::endl;
 
    /*if (Use["CutsGA"]) {
 
       // test: retrieve cuts for particular signal efficiency
       // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer
       TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;
 
       if (mcuts) {
          std::vector<Double_t> cutsMin;
          std::vector<Double_t> cutsMax;
          mcuts->GetCuts( 0.7, cutsMin, cutsMax );
          std::cout << "--- -------------------------------------------------------------" << std::endl;
          std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
          for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
             std::cout << "... Cut: "
                       << cutsMin[ivar]
                       << " < \""
                       << mcuts->GetInputVar(ivar)
                       << "\" <= "
                       << cutsMax[ivar] << std::endl;
          }
          std::cout << "--- -------------------------------------------------------------" << std::endl;
       }
    }*/
 
    // Write histograms
 
    TFile *target  = new TFile( "analysis/TMVApp.root","RECREATE" );
    
    target->Close();
 
    //std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
 
    delete reader;
 
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl;
    cout<<"Average of total: "<<sum/(theTree->GetEntries())<<endl<<endl<<endl<<endl;
    
    
    /*
    if(topology == 1)
    {
	    if(nBack == 0)
	    {
	    	histBdtqqG->SetLineColor(kBlue);
	    	TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper1 = new TCanvas();
	   	cFSRTMVAClassificationApplicationHHbbbbHelper1->cd();
	    	histBdtqqG->Draw("HIST");
	    	histOutput.Add(histBdtqqG);
	    }
	    
	    if(nBack == 1)
	    {
		histBdtttbar->SetLineColor(kBlue);
		TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper2 = new TCanvas();
	    	cFSRTMVAClassificationApplicationHHbbbbHelper2->cd();
	    	histBdtttbar->Draw("HIST");
	    	histOutput.Add(histBdtttbarG);
	    }
    }	
    if(topology == 2)
    {
	    if(nBack == 0)
	    {
	    	histBdtqqG->SetLineColor(kRed+2);
	    	TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper3 = new TCanvas();
	   	cFSRTMVAClassificationApplicationHHbbbbHelper3->cd();
	    	histBdtqqG->Draw("HIST");
	    	histOutput.Add(histBdtqqG);
	    }
	    
	    if(nBack == 1)
	    {
		histBdtttbar->SetLineColor(kRed+2);
		TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper4 = new TCanvas();
	    	cFSRTMVAClassificationApplicationHHbbbbHelper4->cd();
	    	histBdtttbar->Draw("HIST");
	    	histOutput.Add(histBdtttbarG);
	    }
    }	
    if(topology == 3)
    {
	    if(nBack == 0)
	    {
	    	histBdtqq->SetLineColor(kPink+6);
	    	TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper5 = new TCanvas();
	   	cFSRTMVAClassificationApplicationHHbbbbHelper5->cd();
	    	histBdtqq->Draw("HIST");
	    	histOutput.Add(histBdtqqG);
	    }
	    
	    if(nBack == 1)
	    {
		histBdtttbar->SetLineColor(kPink+6);
		TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper6 = new TCanvas();
	    	cFSRTMVAClassificationApplicationHHbbbbHelper6->cd();
	    	histBdtttbar->Draw("HIST");
	    	histOutput.Add(histBdtttbarG);
	    }
    }
    */
   inputS->Close();
   inputBqq->Close();
   inputBtt->Close();
   inputBZZ->Close();
   inputBWW->Close();
   inputBqqX->Close();
   inputBqqqqX->Close();
   inputBqqHX->Close();
	
	return;	
 	
}
 
////Function that generates the 2D hist with the outputs of the BDTs.
void BDTOutputPlane(double bottomHistLimit, double topHistLimit, vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, int sizeHH, int sizeqq, int sizettbar, TH2F& histBDTsOutputHH)
{
 	TCanvas *cBDTOutputPlane1 = new TCanvas();
	cBDTOutputPlane1->cd();
 	histBDTsOutputHH.Draw("col");
 	for(int i=0; i<BDTqqOutput.size(); i++)
 	{
		double radius = 0.005;
		TEllipse* r = new TEllipse(BDTqqOutput[i], BDTttbarOutput[i], radius, radius); 
		if(i<sizeHH) r->SetLineColor(kBlue);	
		else if(i<(sizeHH+sizeqq)) r->SetLineColor(kRed+2);
		else if(i<(sizeHH+sizeqq+sizettbar)) r->SetLineColor(kPink+6);
		r->SetLineStyle(2);
              	r->Draw("same");
 	}
}
 
////Function that generates a histogram with the overlap of two given histograms.
void BDTOutputOverlap(double bottomHistLimit, double topHistLimit, int nbin, TH1F& histBDTHH, Int_t colorHH, TH1F& histBDTBack, Int_t colorBack, string title)
{
 	TH1F *histBdtOverlap     = new TH1F( "MVA_BDT_Overlap", title.c_str(), nbin, bottomHistLimit-.05, topHistLimit+.05);
 	histBdtOverlap->SetFillColor(kViolet);
 	
 	TH1F *histBDTHHClone = (TH1F*)histBDTHH.Clone("histBDTHHClone");
 	TH1F *histBDTBackClone = (TH1F*)histBDTBack.Clone("histBDTBackClone");
 	histBDTHHClone->SetFillColor(colorHH);
 	histBDTBackClone->SetFillColor(colorBack);

	histBDTHHClone->Scale(1.0 / histBDTHHClone->Integral());
	histBDTBackClone->Scale(1.0 / histBDTBackClone->Integral());
 	
 	
 	double contEventsHH=0, contEventsBack=0;
 	
 	for (int i=1; i<=nbin; i++) 
  	{
     		double content1 = histBDTHHClone->GetBinContent(i);
     		double content2 = histBDTBackClone->GetBinContent(i);
	     	if(content1!=0 && content2!=0)
	    	{
	     		double overlappingContent = min(content1, content2);
	     		histBdtOverlap->SetBinContent(i, overlappingContent);
	     	}
          	contEventsHH+=content1;
          	contEventsBack+=content2;
   	}
   	
   	/*TCanvas *cBDTOutputOverlap1 = new TCanvas();
   	cBDTOutputOverlap1->cd();
   	histBDTHH->Draw("HIST");
   	histBDTBack->Draw("HIST same");
   	histBdtOverlap->Draw("HIST same");*/
   	
   	
   	double maxHH = 	histBDTHHClone->GetMaximum();
   	double maxBack = histBDTBackClone->GetMaximum();
   	
   	//cout<<"maxHH: "<<maxHH<<endl<<"maxBack: "<<maxBack<<endl;
   	
   	TCanvas *cBDTOutputOverlap1 = new TCanvas();
   	cBDTOutputOverlap1->cd();
   	if(maxHH>maxBack)
   	{
   		histBDTHHClone->SetTitle(title.c_str());
   		histBDTHHClone->DrawNormalized("HIST");
   		histBDTBackClone->DrawNormalized("HIST same");
   		histBdtOverlap->Draw("HIST same");
   	}
   	else
   	{
   		histBDTBackClone->SetTitle(title.c_str());
   		histBDTBackClone->Draw("HIST");
   		histBDTHHClone->Draw("HIST same");
   		histBdtOverlap->Draw("HIST same");
   	} 
   	
   	
   	return;
}
 
 ////Function that, given two vectors, sorts the first one by increasing order and the entries of the second one to follow the first one. For example, if I have vector1 2,3,1,4,5 and vector2 7,2,9,3,1, I need vector1 to be 1,2,3,4,5 and vector2 to be 9,7,2,3,1. 
 void sortVectorPair(vector<double>& vector1, vector<double>& vector2)
 {
 	    std::vector<std::pair<double, double>> pairedVectors(vector1.size());
	    for (size_t i = 0; i < vector1.size(); ++i) {
		pairedVectors[i] = std::make_pair(vector1[i], vector2[i]);
	    }
	    // Step 2: Sort the vector of pairs based on the first element of the pair
	    std::sort(pairedVectors.begin(), pairedVectors.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
		return a.first < b.first;
	    });
	    // Step 3: Unpack the sorted pairs back into two separate vectors
	    for (size_t i = 0; i < pairedVectors.size(); ++i) {
		vector1[i] = pairedVectors[i].first;
		vector2[i] = pairedVectors[i].second;
	    }
	    return;
 }

 ////Function that, given five vectors, sorts the first one by increasing order and rearranges the entries of the second, third ones to follow the first one. For example, if I have vector1 3,1,2, vector2 9,7,8, and vector3 4,2,3, I need vector1 to be 1,2,3 and vector2 to be 7,8,9 and vector3 to be 2,3,4.
void sortVectorQuadruple(std::vector<double>& vector1, std::vector<double>& vector2, std::vector<vector<double>>& vector3, std::vector<double>& vector4)
{
    // Check that all vectors are the same size
    if (vector1.size() != vector2.size() || vector1.size() != vector3.size() || vector1.size() != vector4.size()) {
        throw std::invalid_argument("All vectors must have the same size.");
    }

    // Step 1: Pair all vectors into a vector of tuples
    std::vector<std::tuple<double, double, vector<double>, double>> pairedVectors(vector1.size());
    for (size_t i = 0; i < vector1.size(); ++i) {
        pairedVectors[i] = std::make_tuple(vector1[i], vector2[i], vector3[i], vector4[i]);
    }

    // Step 2: Sort the vector of tuples based on the first element of the tuple
    std::sort(pairedVectors.begin(), pairedVectors.end(), 
        [](const std::tuple<double, double, vector<double>, double>& a, const std::tuple<double, double, vector<double>, double>& b) {
            return std::get<0>(a) < std::get<0>(b);
        });

    // Step 3: Unpack the sorted tuples back into the four separate vectors
    for (size_t i = 0; i < pairedVectors.size(); ++i) {
        vector1[i] = std::get<0>(pairedVectors[i]);
        vector2[i] = std::get<1>(pairedVectors[i]);
        vector3[i] = std::get<2>(pairedVectors[i]);
		vector4[i] = std::get<3>(pairedVectors[i]);
    }

    return;
} 
 
 ////Function that, given five vectors, sorts the first one by increasing order and rearranges the entries of the second, third... ones to follow the first one. For example, if I have vector1 3,1,2, vector2 9,7,8, and vector3 4,2,3, I need vector1 to be 1,2,3 and vector2 to be 7,8,9 and vector3 to be 2,3,4.
void sortVectorQuintuple(std::vector<double>& vector1, std::vector<double>& vector2, std::vector<double>& vector3, std::vector<double>& vector4, std::vector<double>& vector5)
{
    // Check that all vectors are the same size
    if (vector1.size() != vector2.size() || vector1.size() != vector3.size() || vector1.size() != vector4.size() || vector1.size() != vector5.size()) {
        throw std::invalid_argument("All vectors must have the same size.");
    }

    // Step 1: Pair all vectors into a vector of tuples
    std::vector<std::tuple<double, double, double, double, double>> pairedVectors(vector1.size());
    for (size_t i = 0; i < vector1.size(); ++i) {
        pairedVectors[i] = std::make_tuple(vector1[i], vector2[i], vector3[i], vector4[i], vector5[i]);
    }

    // Step 2: Sort the vector of tuples based on the first element of the tuple
    std::sort(pairedVectors.begin(), pairedVectors.end(), 
        [](const std::tuple<double, double, double, double, double>& a, const std::tuple<double, double, double, double, double>& b) {
            return std::get<0>(a) < std::get<0>(b);
        });

    // Step 3: Unpack the sorted tuples back into the four separate vectors
    for (size_t i = 0; i < pairedVectors.size(); ++i) {
        vector1[i] = std::get<0>(pairedVectors[i]);
        vector2[i] = std::get<1>(pairedVectors[i]);
        vector3[i] = std::get<2>(pairedVectors[i]);
        vector4[i] = std::get<3>(pairedVectors[i]);
        vector5[i] = std::get<4>(pairedVectors[i]);
    }

    return;
}

////Function that, given eight vectors, sorts the first one by increasing order and rearranges the entries of the second, third... ones to follow the first one. For example, if I have vector1 3,1,2, vector2 9,7,8, and vector3 4,2,3, I need vector1 to be 1,2,3 and vector2 to be 7,8,9 and vector3 to be 2,3,4.
void sortVectorEightuple(std::vector<double>& vector1, std::vector<double>& vector2, std::vector<double>& vector3, std::vector<double>& vector4, std::vector<double>& vector5, std::vector<double>& vector6, std::vector<double>& vector7, std::vector<double>& vector8)
{
    // Check that all vectors are the same size
    if(vector1.size() != vector2.size() || vector1.size() != vector3.size() || vector1.size() != vector4.size() || vector1.size() != vector5.size() || vector1.size() != vector6.size() || vector1.size() != vector7.size() || vector1.size() != vector8.size()) {
        throw std::invalid_argument("All vectors must have the same size.");
    }

    // Step 1: Pair all vectors into a vector of tuples
    std::vector<std::tuple<double, double, double, double, double, double, double, double>> pairedVectors(vector1.size());
    for (size_t i = 0; i < vector1.size(); ++i) {
        pairedVectors[i] = std::make_tuple(vector1[i], vector2[i], vector3[i], vector4[i], vector5[i], vector6[i], vector7[i], vector8[i]);
    }

    // Step 2: Sort the vector of tuples based on the first element of the tuple
    std::sort(pairedVectors.begin(), pairedVectors.end(), 
        [](const std::tuple<double, double, double, double, double, double, double, double>& a, const std::tuple<double, double, double, double, double, double, double, double>& b) {
            return std::get<0>(a) < std::get<0>(b);
        });

    // Step 3: Unpack the sorted tuples back into the four separate vectors
    for (size_t i = 0; i < pairedVectors.size(); ++i) {
        vector1[i] = std::get<0>(pairedVectors[i]);
        vector2[i] = std::get<1>(pairedVectors[i]);
        vector3[i] = std::get<2>(pairedVectors[i]);
        vector4[i] = std::get<3>(pairedVectors[i]);
        vector5[i] = std::get<4>(pairedVectors[i]);
		vector6[i] = std::get<5>(pairedVectors[i]);
		vector7[i] = std::get<6>(pairedVectors[i]);
		vector8[i] = std::get<7>(pairedVectors[i]);
    }

    return;
}
 
 
///Function that plots the significance as a function of the cuts on BDTqq and BDTttbar; it then returns the maximum significance.
 double findMaxSignificanceAndROC(double bottomHistLimit, double topHistLimit, double nbin, TH1F& histBDTHH, TH1F& histBDTBack, double& n, double& b, double& optCut)
 {
 	 double low=bottomHistLimit, high=topHistLimit;
 	 double cut;
 	 Int_t binLow = histBDTHH.FindBin(low);
 	 Int_t binHigh = histBDTHH.FindBin(high);
 	 double totalIntegral = histBDTHH.Integral(binLow, binHigh);
 	 double totalIntegralBack = histBDTBack.Integral(binLow, binHigh);
 	 double contEvents=totalIntegral;
	 double significance=0.0, maxSignificance=0.0, signalRight=0.0, BDTEffi=0.0;
         //cout<<endl<<endl<<"Total integral for histBdt: "<<totalIntegral<<endl<<"Right integral for histBdt: "<<rightIntegral<<endl;
         
         TH2F *histROC = new TH2F("ROC", "Signal and Background Acceptance", nbin, 0.0, 1.0, nbin, 0, 1.0);
  	 TH2F *histROCRej = new TH2F("ROCRej", "Signal Acceptance and Background Rejection", nbin, 0.5, 1.0, nbin, 0, 10000.0);
  	 histROC->GetXaxis()->SetTitle("Signal accepted");
   	 histROC->GetYaxis()->SetTitle("Background rejected (1/background accepted)");
   
  	 TH2F *histSignificance     = new TH2F( "hist_significance", "Significance s/sqrt(s+b)", nbin, bottomHistLimit-.1, topHistLimit+.1, nbin, 0, 10.0);
  	 histSignificance->GetXaxis()->SetTitle("Cut value for BDT");
   	 histSignificance->GetYaxis()->SetTitle("Significance");
   
   	 for (double i=bottomHistLimit; i<=topHistLimit; i+=0.0001) 
   	 {
   		cut=i;
   		Int_t binCut = histBDTHH.FindBin(cut);
   		double rightIntegral = histBDTHH.Integral(binCut, binHigh); 
   		double rightIntegralBack = histBDTBack.Integral(binCut, binHigh);
	   	/*cout<<"Cut: "<<i<<endl;
	   	cout<<"total integral: "<<totalIntegral<<endl;
	   	cout<<"total integral back: "<<totalIntegralBack<<endl;
	   	cout<<"Right integral: "<<rightIntegral<<endl;
	   	cout<<"Right integral back: "<<rightIntegralBack<<endl<<endl;*/
	   	histROC->Fill(rightIntegral/totalIntegral, rightIntegralBack/totalIntegralBack);
	   	histROCRej->Fill(rightIntegral/totalIntegral, 1/(rightIntegralBack/totalIntegralBack));
	   	significance=rightIntegral/(sqrt(rightIntegral+rightIntegralBack));
	   	histSignificance->Fill(cut, significance);
	   	if(significance>maxSignificance) 
	   	{
	   		maxSignificance=significance;
	   		optCut=cut;
	   		n=rightIntegral+rightIntegralBack;
	   		b=rightIntegralBack;
	   		signalRight=rightIntegral;
	   	}
	   	//if((rightIntegral/totalIntegral)>.599 && (rightIntegral/totalIntegral)<.601) cout<<"For "<<(rightIntegral/totalIntegral)<<" signal acceptance, background acceptance: "<<(rightIntegralBack/totalIntegralBack)<<endl<<endl<<endl;
	  }
	  
	  TCanvas *cfindMaxSignificanceAndROC1 = new TCanvas();
	  TCanvas *cfindMaxSignificanceAndROC2 = new TCanvas();
	  TCanvas *cfindMaxSignificanceAndROC3 = new TCanvas();
	  
	  cfindMaxSignificanceAndROC1->cd();
	  histROC->Draw("BOX");
	  
	  //gStyle->SetOptStat(0);
	  cfindMaxSignificanceAndROC2->cd();
	  cfindMaxSignificanceAndROC2->SetLogy();
	  histROCRej->Draw("BOX");
	  
	  cfindMaxSignificanceAndROC3->cd();
	  histSignificance->Draw("HIST");
	  
	  return maxSignificance;
 }
 
 /////Function that generates ROC curve without using histograms and integrals, as well as calculates the max. significance
 void findSignificance(double bottomHistLimit, double topHistLimit, double nbin, vector<double> NNOutput, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, int nBack, TH1F& histBDTHH, TH1F& histBDTBack, double& defCut, double& maxSignificance, double weightHH, double weightBack, double targetFractionHH, TH2F& histROC, TH2F& histROCRej, TH2F& histSignificance, double& totalRemaining, double& HHRemaining, double& backRemaining)
 {
 	defCut=bottomHistLimit;
	double significance; 	
  	 
	double size = NNOutput.size();
	
	//cout<<endl<<endl<<endl<<endl<<endl<<"NNOutputSize: "<<size<<endl<<endl<<endl<<endl;
	
	vector<double> NNOutputHH(NNOutput.begin(), NNOutput.begin()+sizeHH);
	double bottomLimitBack, topLimitBack;
	if(nBack == 0) 
	{
		bottomLimitBack = sizeHH;
		topLimitBack = sizeHH+sizeqq;
	}
	else if(nBack == 1)
	{
		bottomLimitBack = sizeHH+sizeqq;
		topLimitBack = sizeHH+sizeqq+sizettbar;
	}
	else if(nBack == 2)
	{
		bottomLimitBack = sizeHH+sizeqq+sizettbar;
		topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ;
	}
	else if(nBack == 3)
	{
		bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ;
		topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW;
	}
	else if(nBack == 4)
	{
		bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW;
		topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX;
	}
	else if(nBack == 5)
	{
		bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX;
		topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX;
	}
	else if(nBack == 6)
	{
		bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX;
		topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX+sizeqqHX;
	}

	//cout<<"bottomLimitBack: "<<bottomLimitBack<<"      topLimitBack: "<<topLimitBack<<endl<<endl<<endl<<endl;
	vector<double> NNOutputBack(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	//cout<<"No crashea en la creacion del vector de back ouputs"<<endl<<endl<<endl;
	sort(NNOutputHH.begin(), NNOutputHH.end());
	sort(NNOutputBack.begin(), NNOutputBack.end());
	
	int sizeBack = NNOutputBack.size();
	
	//cout<<"sizeBack: "<<sizeBack<<endl<<endl<<endl;
	
	//cout<<"Size: "<<size<<endl<<"SizeHH: "<<NNOutputHH.size()<<endl<<"SizeBack: "<<NNOutputBack.size()<<endl<<endl;
	/*for(int i=0; i<100; i++) cout<<NNOutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<100; i++) cout<<NNOutputBack[i]<<", ";
	cout<<endl;*/
	
	int indexHH=0, indexBack=0;
	double elementNHH=NNOutputHH[indexHH],elementN1HH, elementNBack=NNOutputBack[indexBack], elementN1Back;
	double eventsRemainingHH=0, eventsCutHH=0, eventsRemainingBack=0, eventsCutBack=0;
	bool flagFraction=false;
	//for(double cut=bottomHistLimit; cut<=topHistLimit; cut+=0.000001)
	for(double cut=bottomHistLimit; cut<topHistLimit; cut+=0.0001)
	{
		
		
		while(cut>elementNHH && indexHH<sizeHH)
		{
			indexHH++;
			if(indexHH<sizeHH) elementNHH = NNOutputHH[indexHH];
		}
		while(cut>elementNBack && indexBack<sizeBack)
		{
			indexBack++;
			if(indexBack<sizeBack) elementNBack = NNOutputBack[indexBack];
		}
		
		eventsRemainingHH = sizeHH-indexHH;
		eventsRemainingBack = sizeBack-indexBack;
		
		double fractionHH = eventsRemainingHH/sizeHH;
		double fractionBack = eventsRemainingBack/sizeBack;
		histROCRej.Fill(fractionHH, 1/(fractionBack));
		histROC.Fill(eventsRemainingHH/sizeHH, eventsRemainingBack/sizeBack);
		
		significance=(eventsRemainingHH*weightHH)/(sqrt((eventsRemainingHH*weightHH)+(eventsRemainingBack*weightBack)));
		/*cout<<"CUT: "<<cut<<endl<<endl;
		cout<<"IndexHH: "<<indexHH<<"    elem: "<<elementNHH<<endl<<"Signal remaining: "<<eventsRemainingHH<<"    signal cut: "<<sizeHH-eventsRemainingHH<<endl<<endl;
		cout<<"IndexBack: "<<indexBack<<"    elem: "<<elementNBack<<endl<<"Back remaining: "<<eventsRemainingBack<<"    back cut: "<<sizeBack-eventsRemainingBack<<endl;*/
		//cout<<"Significance: "<<significance<<"         maxSignificnace: "<<maxSignificance<<endl<<endl;
		//cout<<"--------"<<endl;
	   	histSignificance.Fill(cut, significance);
	   	if(significance>maxSignificance) 
	   	{
	   		maxSignificance=significance;
	   		defCut=cut;
	   		backRemaining=eventsRemainingBack*weightBack;
	   		HHRemaining=eventsRemainingHH*weightHH;
	   		totalRemaining=eventsRemainingHH+eventsRemainingBack;
	   	}
		
		/*if(fractionHH < targetFractionHH+0.001 && fractionHH > targetFractionHH-0.001 && !flagFraction)
		{
			flagFraction=true;
			cout<<endl<<"(WITH ARRAYS) For "<<fractionHH*100<<"% of signal ("<<eventsRemainingHH<<" out of "<<sizeHH<<" events), background: "<<fractionBack*100<<"% ("<<eventsRemainingBack<<" out of "<<sizeBack<<" events)"<<endl;
			double low=bottomHistLimit, high=topHistLimit;
 	 		Int_t binLow = histBDTHH.FindBin(cut);
 	 		Int_t binHigh = histBDTHH.FindBin(high);
 	 		double eventsRemainingHHIntegral = histBDTHH.Integral(binLow, binHigh);
 	 		double eventsRemainingBackIntegral = histBDTBack.Integral(binLow, binHigh);
 	 		double fractionHHIntegral = eventsRemainingHHIntegral/sizeHH;
			double fractionBackIntegral = eventsRemainingBackIntegral/sizeBack;
 	 		cout<<endl<<"(WITH INTEGRALS) For "<<fractionHHIntegral*100<<"% of signal ("<<eventsRemainingHHIntegral<<" out of "<<sizeHH<<" events), background: "<<fractionBackIntegral*100<<"% ("<<eventsRemainingBackIntegral<<" out of "<<sizeBack<<" events)"<<endl;
 	 		defCut = cut;
 	 		cout<<"Cut for this %: "<<cut<<endl<<endl;
 	 		
 	 		cout<<"Significance for this: "<<significance<<endl<<endl;
		}*/
	}
	cout<<"HHRemainingUnweighted: "<<HHRemaining/weightHH<<endl<<"backRemainingUnweighted: "<<backRemaining/weightBack<<endl; 
	cout<<"HHRemaining: "<<HHRemaining<<endl<<"backRemaining: "<<backRemaining<<endl; 
	
	TCanvas *cfindSignificance1 = new TCanvas();
	TCanvas *cfindSignificance2 = new TCanvas();
	TCanvas *cfindSignificance3 = new TCanvas();
	
	cfindSignificance1->cd();
	histROC.Draw("BOX");
	  
	//gStyle->SetOptStat(0);
	cfindSignificance2->cd();
	cfindSignificance2->SetLogy();
	histROCRej.Draw("BOX");
	
	cfindSignificance3->cd();
	histSignificance.Draw("BOX");
	
 }
///Function that finds cut for NNs individually applied to all backs
 void findSignificanceAllBacks(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& NNOutput, double& defCut, double& maxSignificance, double& nBack, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, double weightHH, double weightqq, double weightttbar, double weightZZ, double weightWW, double weightqqX, double weightqqqqX, double weightqqHX, double cutStep)
 {
 	
	//cout<<"calling findSignificanceAllBacks";
	
	defCut=bottomHistLimit;
	double significance, HHRemaining, backRemaining; 	
	double size = NNOutput.size();
	double bottomLimitBack, topLimitBack;

	vector<double> NNOutputHH(NNOutput.begin(), NNOutput.begin()+sizeHH);
	bottomLimitBack = sizeHH;
	topLimitBack = sizeHH+sizeqq;
	vector<double> NNOutputBackqq(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	bottomLimitBack = sizeHH+sizeqq;
	topLimitBack = sizeHH+sizeqq+sizettbar;
	vector<double> NNOutputBackttbar(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	bottomLimitBack = sizeHH+sizeqq+sizettbar;
	topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ;
	vector<double> NNOutputBackZZ(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ;
	topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW;
	vector<double> NNOutputBackWW(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW;
	topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX;
	vector<double> NNOutputBackqqX(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX;
	topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX;
	vector<double> NNOutputBackqqqqX(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX;
	topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX+sizeqqHX;
	vector<double> NNOutputBackqqHX(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));

	sort(NNOutputHH.begin(), NNOutputHH.end());
	sort(NNOutputBackqq.begin(), NNOutputBackqq.end());
	sort(NNOutputBackttbar.begin(), NNOutputBackttbar.end());
	sort(NNOutputBackZZ.begin(), NNOutputBackZZ.end());
	sort(NNOutputBackWW.begin(), NNOutputBackWW.end());
	sort(NNOutputBackqqX.begin(), NNOutputBackqqX.end());
	sort(NNOutputBackqqqqX.begin(), NNOutputBackqqqqX.end());
	sort(NNOutputBackqqHX.begin(), NNOutputBackqqHX.end());

	int sizeBackqq = NNOutputBackqq.size();
	int sizeBackttbar = NNOutputBackttbar.size();
	int sizeBackZZ = NNOutputBackZZ.size();
	int sizeBackWW = NNOutputBackWW.size();
	int sizeBackqqX = NNOutputBackqqX.size();
	int sizeBackqqqqX = NNOutputBackqqqqX.size();
	int sizeBackqqHX = NNOutputBackqqHX.size();

	int indexHH=0, indexBackqq=0, indexBackttbar=0, indexBackZZ=0, indexBackWW=0, indexBackqqX=0, indexBackqqqqX=0, indexBackqqHX=0;
	int indexHHMaxSignificance=0, indexBackqqMaxSignificance=0, indexBackttbarMaxSignificance=0, indexBackZZMaxSignificance=0, indexBackWWMaxSignificance=0, indexBackqqXMaxSignificance=0, indexBackqqqqXMaxSignificance=0, indexBackqqHXMaxSignificance=0;
	double elementNHH=NNOutputHH[indexHH], elementN1HH, elementNBackqq=NNOutputBackqq[indexBackqq], elementN1Backqq, elementNBackttbar=NNOutputBackttbar[indexBackttbar], elementN1Backttbar, elementNBackZZ=NNOutputBackZZ[indexBackZZ], elementN1BackZZ, elementNBackWW=NNOutputBackWW[indexBackWW], elementN1BackWW, elementNBackqqX=NNOutputBackqqX[indexBackqqX], elementN1BackqqX, elementNBackqqqqX=NNOutputBackqqqqX[indexBackqqqqX], elementN1BackqqqqX, elementNBackqqHX=NNOutputBackqqHX[indexBackqqHX], elementN1BackqqHX;
	double eventsRemainingHH=0, eventsCutHH=0, eventsRemainingBack=0, eventsCutBack=0, eventsRemainingBackW=0, eventsCutBackW=0, eventsRemainingBackqq=0, eventsCutBackqq=0, eventsRemainingBackttbar=0, eventsCutBackttbar=0, eventsRemainingBackZZ=0, eventsCutBackZZ=0, eventsRemainingBackWW=0, eventsCutBackWW=0, eventsRemainingBackqqX=0, eventsCutBackqqX=0, eventsRemainingBackqqqqX=0, eventsCutBackqqqqX=0, eventsRemainingBackqqHX=0, eventsCutBackqqHX=0;
	//double cutStep=0.001;
	//cout<<endl<<"cutStep: "<<cutStep<<endl;
	for(double cut=bottomHistLimit; cut<topHistLimit; cut+=cutStep)
	{
		while(cut>elementNHH && indexHH<sizeHH)
		{
			indexHH++;
			if(indexHH<sizeHH) elementNHH = NNOutputHH[indexHH];
		}
		while(cut>elementNBackqq && indexBackqq<sizeBackqq)
		{
			indexBackqq++;
			if(indexBackqq<sizeBackqq) elementNBackqq = NNOutputBackqq[indexBackqq];
		}
		while(cut>elementNBackttbar && indexBackttbar<sizeBackttbar)
		{
			indexBackttbar++;
			if(indexBackttbar<sizeBackttbar) elementNBackttbar = NNOutputBackttbar[indexBackttbar];
		}
		while(cut>elementNBackZZ && indexBackZZ<sizeBackZZ)
		{
			indexBackZZ++;
			if(indexBackZZ<sizeBackZZ) elementNBackZZ = NNOutputBackZZ[indexBackZZ];
		}
		while(cut>elementNBackWW && indexBackWW<sizeBackWW)
		{
			indexBackWW++;
			if(indexBackWW<sizeBackWW) elementNBackWW = NNOutputBackWW[indexBackWW];
		}
		while(cut>elementNBackqqX && indexBackqqX<sizeBackqqX)
		{
			indexBackqqX++;
			if(indexBackqqX<sizeBackqqX) elementNBackqqX = NNOutputBackqqX[indexBackqqX];
		}
		while(cut>elementNBackqqqqX && indexBackqqqqX<sizeBackqqqqX)
		{
			indexBackqqqqX++;
			if(indexBackqqqqX<sizeBackqqqqX) elementNBackqqqqX = NNOutputBackqqqqX[indexBackqqqqX];
		}
		while(cut>elementNBackqqHX && indexBackqqHX<sizeBackqqHX)
		{
			indexBackqqHX++;
			if(indexBackqqHX<sizeBackqqHX) elementNBackqqHX = NNOutputBackqqHX[indexBackqqHX];
		}

		eventsRemainingHH = sizeHH-indexHH;
		eventsRemainingBackqq = sizeBackqq-indexBackqq;
		eventsRemainingBackttbar = sizeBackttbar-indexBackttbar;
		eventsRemainingBackZZ = sizeBackZZ-indexBackZZ;
		eventsRemainingBackWW = sizeBackWW-indexBackWW;
		eventsRemainingBackqqX = sizeBackqqX-indexBackqqX;
		eventsRemainingBackqqqqX = sizeBackqqqqX-indexBackqqqqX;
		eventsRemainingBackqqHX = sizeBackqqHX-indexBackqqHX;

		eventsRemainingBack = eventsRemainingBackqq+eventsRemainingBackttbar+eventsRemainingBackZZ+eventsRemainingBackWW+eventsRemainingBackqqX+eventsRemainingBackqqqqX+eventsRemainingBackqqHX;
		eventsRemainingBackW = eventsRemainingBackqq*weightqq+eventsRemainingBackttbar*weightttbar+eventsRemainingBackZZ*weightZZ+eventsRemainingBackWW*weightWW+eventsRemainingBackqqX*weightqqX+eventsRemainingBackqqqqX*weightqqqqX+eventsRemainingBackqqHX*weightqqHX;

		significance=(eventsRemainingHH*weightHH)/(sqrt((eventsRemainingHH*weightHH)+(eventsRemainingBackW)));
		/*cout<<"CUT: "<<cut<<endl<<endl;
		cout<<"IndexHH: "<<indexHH<<"    elem: "<<elementNHH<<endl<<"Signal remaining: "<<eventsRemainingHH<<"    signal cut: "<<sizeHH-eventsRemainingHH<<endl<<endl;
		cout<<"IndexBackqq: "<<indexBackqq<<"    elemqq: "<<elementNBackqq<<endl<<"Backqq remaining: "<<eventsRemainingBackqq<<"    backqq cut: "<<sizeBackqq-eventsRemainingBackqq<<endl;
		cout<<"IndexBackttbar: "<<indexBackttbar<<"    elemttbar: "<<elementNBackttbar<<endl<<"Backttbar remaining: "<<eventsRemainingBackttbar<<"    backttbar cut: "<<sizeBackttbar-eventsRemainingBackttbar<<endl;
		cout<<"IndexBackZZ: "<<indexBackZZ<<"    elemZZ: "<<elementNBackZZ<<endl<<"BackZZ remaining: "<<eventsRemainingBackZZ<<"    backZZ cut: "<<sizeBackZZ-eventsRemainingBackZZ<<endl;
		cout<<"IndexBackWW: "<<indexBackWW<<"    elemWW: "<<elementNBackWW<<endl<<"BackWW remaining: "<<eventsRemainingBackWW<<"    backWW cut: "<<sizeBackWW-eventsRemainingBackWW<<endl;
		cout<<"IndexBackqqX: "<<indexBackqqX<<"    elemqqX: "<<elementNBackqqX<<endl<<"BackqqX remaining: "<<eventsRemainingBackqqX<<"    backqqX cut: "<<sizeBackqqX-eventsRemainingBackqqX<<endl;
		cout<<"IndexBackqqqqX: "<<indexBackqqqqX<<"    elemqqqqX: "<<elementNBackqqqqX<<endl<<"BackqqqqX remaining: "<<eventsRemainingBackqqqqX<<"    backqqqqX cut: "<<sizeBackqqqqX-eventsRemainingBackqqqqX<<endl;
		cout<<"IndexBackqqHX: "<<indexBackqqHX<<"    elemqqHX: "<<elementNBackqqHX<<endl<<"BackqqHX remaining: "<<eventsRemainingBackqqHX<<"    backqqHX cut: "<<sizeBackqqHX-eventsRemainingBackqqHX<<endl;
		cout<<"Significance: "<<significance<<"         maxSignificnace: "<<maxSignificance<<endl<<endl;*/
	   	if(significance>maxSignificance) 
	   	{
	   		maxSignificance=significance;
	   		defCut=cut;
	   		backRemaining=eventsRemainingBack;
	   		HHRemaining=eventsRemainingHH*weightHH;
			indexHHMaxSignificance=indexHH;
			indexBackqqMaxSignificance=indexBackqq;
			indexBackttbarMaxSignificance=indexBackttbar;
			indexBackZZMaxSignificance=indexBackZZ;
			indexBackWWMaxSignificance=indexBackWW;
			indexBackqqXMaxSignificance=indexBackqqX;
			indexBackqqqqXMaxSignificance=indexBackqqqqX;
			indexBackqqHXMaxSignificance=indexBackqqHX;
	   	}

		//cout<<"Significance: "<<significance<<"         maxSignificnace: "<<maxSignificance<<endl<<endl;
		//cout<<"--------"<<endl;
	}
		NNOutputHH.erase(NNOutputHH.begin(), NNOutputHH.begin()+indexHHMaxSignificance);
		NNOutputBackqq.erase(NNOutputBackqq.begin(), NNOutputBackqq.begin()+indexBackqqMaxSignificance);
		NNOutputBackttbar.erase(NNOutputBackttbar.begin(), NNOutputBackttbar.begin()+indexBackttbarMaxSignificance);
		NNOutputBackZZ.erase(NNOutputBackZZ.begin(), NNOutputBackZZ.begin()+indexBackZZMaxSignificance);
		NNOutputBackWW.erase(NNOutputBackWW.begin(), NNOutputBackWW.begin()+indexBackWWMaxSignificance);
		NNOutputBackqqX.erase(NNOutputBackqqX.begin(), NNOutputBackqqX.begin()+indexBackqqXMaxSignificance);
		NNOutputBackqqqqX.erase(NNOutputBackqqqqX.begin(), NNOutputBackqqqqX.begin()+indexBackqqqqXMaxSignificance);
		NNOutputBackqqHX.erase(NNOutputBackqqHX.begin(), NNOutputBackqqHX.begin()+indexBackqqHXMaxSignificance);

	
}
 
////Function that finds the optimal cuts for all the indivivual NNs considering all backs
void findDefCuts(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, vector<double>& BDTqqXOutput, vector<double>& BDTqqqqXOutput, vector<double>& BDTqqHXOutput, double& defCutqq, double& defCutttbar, double& defCutZZ, double& defCutWW, double& defCutqqX, double& defCutqqqqX, double& defCutqqHX, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, double weightHH, double weightqq, double weightttbar, double weightZZ, double weightWW, double weightqqX, double weightqqqqX, double weightqqHX)
{
	vector<double> BDTqqOutputCopy(BDTqqOutput.begin(), BDTqqOutput.end());
	vector<double> BDTttbarOutputCopy(BDTttbarOutput.begin(), BDTttbarOutput.end());
	vector<double> BDTZZOutputCopy(BDTZZOutput.begin(), BDTZZOutput.end());
	vector<double> BDTWWOutputCopy(BDTWWOutput.begin(), BDTWWOutput.end());
	vector<double> BDTqqXOutputCopy(BDTqqXOutput.begin(), BDTqqXOutput.end());
	vector<double> BDTqqqqXOutputCopy(BDTqqqqXOutput.begin(), BDTqqqqXOutput.end());
	vector<double> BDTqqHXOutputCopy(BDTqqHXOutput.begin(), BDTqqHXOutput.end());

	vector<vector<double>> NNOutputs = {BDTqqOutputCopy, BDTttbarOutputCopy, BDTZZOutputCopy, BDTWWOutputCopy, BDTqqXOutputCopy, BDTqqqqXOutputCopy, BDTqqHXOutputCopy};
	double defCutqqSen=bottomHistLimit, defCutttbarSen=bottomHistLimit, defCutZZSen=bottomHistLimit, defCutWWSen=bottomHistLimit, defCutqqXSen=bottomHistLimit, defCutqqqqXSen=bottomHistLimit, defCutqqHXSen=bottomHistLimit;
	vector<double> NNCuts = {defCutqqSen, defCutttbarSen, defCutZZSen, defCutWWSen, defCutqqXSen, defCutqqqqXSen, defCutqqHXSen};
	vector<double> maxSignificances = {-999, -999, -999, -999, -999, -999, -999};
	vector<double> nBacks = {0, 1, 2, 3, 4, 5, 6};
	double numberOfBacks=nBacks.size();
	double cutPrev=0, cutStep=0;
	for(int i=0; i<numberOfBacks; i++)
	{
		//cout<<endl<<endl<<endl<<endl<<"Starting round "<<i<<": "<<endl;

		for(int j=0; j<NNOutputs.size(); j++)
		{
			//cout<<"In position "<<j<<": "<<endl;
			//cout<<"Back: "<<nBacks[j]<<endl<<"number of events: ";
			int contEventsRounds=0;
			for(int k=0; k<NNOutputs[j].size(); k++) if(NNOutputs[j][k] != bottomHistLimit) contEventsRounds++;
			//cout<<contEventsRounds;
			contEventsRounds=0;
			//cout<<"Before findSignificanceAllBacks: "<<endl<<"maxSignificance: "<<maxSignificances[j]<<endl<<"Cut in NN"<<NNCuts[j]<<endl<<endl;
			cutStep = 1.0 / pow(10, i + 1);
			findSignificanceAllBacks(bottomHistLimit, topHistLimit, nbin, NNOutputs[j], NNCuts[j], maxSignificances[j], nBacks[j], sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, cutStep);

			//cout<<"After findSignificanceAllBacks: "<<endl<<"maxSignificance: "<<maxSignificances[j]<<endl<<"Cut in NN: "<<NNCuts[j]<<endl<<endl;
			//cout<<"maxSignificance: "<<maxSignificances[j]<<endl<<"Cut in NN: "<<NNCuts[j]<<endl<<endl;
		}

		//cout<<endl<<endl<<"Ends round findSignificanceAllBacks and starts sortVectorQuadruple: "<<endl<<endl;
		sortVectorQuadruple(maxSignificances, NNCuts, NNOutputs, nBacks);
		//cout<<"last after sortVectorQuadruple: "<<nBacks[nBacks.size()-1]<<endl<<"first after sortVectorQuadruple: "<<nBacks[0]<<endl;
		//cout<<endl<<endl<<"Ends sortVectorQuadruple and starts event elimination: "<<endl<<endl;
		for(int j=0; j<NNOutputs[NNOutputs.size()-1].size(); j++)
		{
			
			//cout<<"Starts event elimination for back "<<j<<": "<<endl;
			//cout<<"back: "<<nBacks[j]<<endl;
			//cout<<"cut: "<<NNCuts[j]<<endl;
			//cout<<"maxSignificance: "<<maxSignificances[j]<<endl;
			//if(j==0) cout<<"range: "<<NNOutputs[NNOutputs.size()-1].size()<<endl;
			//if(j<30) cout<<"NNOutputs[NNOutputs.size()-1][j]: "<<NNOutputs[NNOutputs.size()-1][j]<<endl;
			//if(j<30) cout<<"NNCuts[NNCuts.size()-1]: "<<NNCuts[NNCuts.size()-1]<<endl;
			if(NNOutputs[NNOutputs.size()-1][j] != bottomHistLimit && NNOutputs[NNOutputs.size()-1][j] < NNCuts[(NNCuts.size()-1)])
			{
				//if(j<30) cout<<"Event eliminated"<<endl;
				for(int k=2; k<=NNOutputs.size(); k++)
				{
					//if(j<30) cout<<"Before: "<<NNOutputs[NNOutputs.size()-k][j]<<"     ";
					NNOutputs[NNOutputs.size()-k][j] = bottomHistLimit;
					//if(j<30) cout<<"After: "<<NNOutputs[NNOutputs.size()-k][j]<<endl;
				}
			} 
		}
		//cout<<endl<<endl<<"Ends event elimination and starts cut assignment: "<<endl<<endl;
		if(nBacks[nBacks.size()-1] == 0) defCutqq = NNCuts[NNCuts.size()-1];
		else if(nBacks[nBacks.size()-1] == 1) defCutttbar = NNCuts[NNCuts.size()-1];
		else if(nBacks[nBacks.size()-1] == 2) defCutZZ = NNCuts[NNCuts.size()-1];
		else if(nBacks[nBacks.size()-1] == 3) defCutWW = NNCuts[NNCuts.size()-1];
		else if(nBacks[nBacks.size()-1] == 4) defCutqqX = NNCuts[NNCuts.size()-1];
		else if(nBacks[nBacks.size()-1] == 5) defCutqqqqX = NNCuts[NNCuts.size()-1];
		else if(nBacks[nBacks.size()-1] == 6) defCutqqHX = NNCuts[NNCuts.size()-1];
		maxSignificances.pop_back();
		NNCuts.pop_back();
		NNOutputs.pop_back();
		nBacks.pop_back();
		//cout<<endl<<endl<<"Ends cut assignment and ends round: "<<endl<<endl;
	}
	/*for(int i=0; i<NNCuts.size(); i++)
	{
		cout<<"maxSignificance ("<<method<<NNCuts[i]<<"): "<<maxSignificances[i]<<endl<<"Cut in NN"<<NNCuts[i]<<" for max significance ("<<method<<NNCuts[i]<<"): "<<NNOutputs[i][0]<<endl;
	}*/
}

///Function that finds cut for NNs individually applied to eGamma backs
 void findSignificanceEGammaBacks(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& NNOutput, double& defCut, double& maxSignificance, double& nBack, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, double weightHH, double weightqq, double weightttbar, double weightZZ, double weightWW, double weightqqX, double weightqqqqX, double weightqqHX, double cutStep)
 {
 	
	//cout<<"calling findSignificanceEGammaBacks";
	
	defCut=bottomHistLimit;
	double significance, HHRemaining, backRemaining; 	
	double size = NNOutput.size();
	double bottomLimitBack, topLimitBack;

	vector<double> NNOutputHH(NNOutput.begin(), NNOutput.begin()+sizeHH);
	bottomLimitBack = sizeHH;
	bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW;
	topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX;
	vector<double> NNOutputBackqqX(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX;
	topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX;
	vector<double> NNOutputBackqqqqX(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX;
	topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX+sizeqqHX;
	vector<double> NNOutputBackqqHX(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));

	sort(NNOutputHH.begin(), NNOutputHH.end());
	sort(NNOutputBackqqX.begin(), NNOutputBackqqX.end());
	sort(NNOutputBackqqqqX.begin(), NNOutputBackqqqqX.end());
	sort(NNOutputBackqqHX.begin(), NNOutputBackqqHX.end());

	int sizeBackqqX = NNOutputBackqqX.size();
	int sizeBackqqqqX = NNOutputBackqqqqX.size();
	int sizeBackqqHX = NNOutputBackqqHX.size();

	int indexHH=0, indexBackqqX=0, indexBackqqqqX=0, indexBackqqHX=0;
	int indexHHMaxSignificance=0, indexBackqqXMaxSignificance=0, indexBackqqqqXMaxSignificance=0, indexBackqqHXMaxSignificance=0;
	double elementNHH=NNOutputHH[indexHH], elementN1HH, elementNBackqqX=NNOutputBackqqX[indexBackqqX], elementN1BackqqX, elementNBackqqqqX=NNOutputBackqqqqX[indexBackqqqqX], elementN1BackqqqqX, elementNBackqqHX=NNOutputBackqqHX[indexBackqqHX], elementN1BackqqHX;
	double eventsRemainingHH=0, eventsCutHH=0, eventsRemainingBack=0, eventsCutBack=0, eventsRemainingBackW=0, eventsCutBackW=0, eventsRemainingBackqqX=0, eventsCutBackqqX=0, eventsRemainingBackqqqqX=0, eventsCutBackqqqqX=0, eventsRemainingBackqqHX=0, eventsCutBackqqHX=0;
	//double cutStep=0.001;
	//cout<<endl<<"cutStep: "<<cutStep<<endl;
	for(double cut=bottomHistLimit; cut<topHistLimit; cut+=cutStep)
	{
		while(cut>elementNHH && indexHH<sizeHH)
		{
			indexHH++;
			if(indexHH<sizeHH) elementNHH = NNOutputHH[indexHH];
		}
		while(cut>elementNBackqqX && indexBackqqX<sizeBackqqX)
		{
			indexBackqqX++;
			if(indexBackqqX<sizeBackqqX) elementNBackqqX = NNOutputBackqqX[indexBackqqX];
		}
		while(cut>elementNBackqqqqX && indexBackqqqqX<sizeBackqqqqX)
		{
			indexBackqqqqX++;
			if(indexBackqqqqX<sizeBackqqqqX) elementNBackqqqqX = NNOutputBackqqqqX[indexBackqqqqX];
		}
		while(cut>elementNBackqqHX && indexBackqqHX<sizeBackqqHX)
		{
			indexBackqqHX++;
			if(indexBackqqHX<sizeBackqqHX) elementNBackqqHX = NNOutputBackqqHX[indexBackqqHX];
		}

		eventsRemainingHH = sizeHH-indexHH;
		eventsRemainingBackqqX = sizeBackqqX-indexBackqqX;
		eventsRemainingBackqqqqX = sizeBackqqqqX-indexBackqqqqX;
		eventsRemainingBackqqHX = sizeBackqqHX-indexBackqqHX;

		eventsRemainingBack = eventsRemainingBackqqX+eventsRemainingBackqqqqX+eventsRemainingBackqqHX;
		eventsRemainingBackW = eventsRemainingBackqqX*weightqqX+eventsRemainingBackqqqqX*weightqqqqX+eventsRemainingBackqqHX*weightqqHX;

		significance=(eventsRemainingHH*weightHH)/(sqrt((eventsRemainingHH*weightHH)+(eventsRemainingBackW)));
		/*cout<<"CUT: "<<cut<<endl<<endl;
		cout<<"IndexHH: "<<indexHH<<"    elem: "<<elementNHH<<endl<<"Signal remaining: "<<eventsRemainingHH<<"    signal cut: "<<sizeHH-eventsRemainingHH<<endl<<endl;
		cout<<"IndexBackqq: "<<indexBackqq<<"    elemqq: "<<elementNBackqq<<endl<<"Backqq remaining: "<<eventsRemainingBackqq<<"    backqq cut: "<<sizeBackqq-eventsRemainingBackqq<<endl;
		cout<<"IndexBackttbar: "<<indexBackttbar<<"    elemttbar: "<<elementNBackttbar<<endl<<"Backttbar remaining: "<<eventsRemainingBackttbar<<"    backttbar cut: "<<sizeBackttbar-eventsRemainingBackttbar<<endl;
		cout<<"IndexBackZZ: "<<indexBackZZ<<"    elemZZ: "<<elementNBackZZ<<endl<<"BackZZ remaining: "<<eventsRemainingBackZZ<<"    backZZ cut: "<<sizeBackZZ-eventsRemainingBackZZ<<endl;
		cout<<"IndexBackWW: "<<indexBackWW<<"    elemWW: "<<elementNBackWW<<endl<<"BackWW remaining: "<<eventsRemainingBackWW<<"    backWW cut: "<<sizeBackWW-eventsRemainingBackWW<<endl;
		cout<<"IndexBackqqX: "<<indexBackqqX<<"    elemqqX: "<<elementNBackqqX<<endl<<"BackqqX remaining: "<<eventsRemainingBackqqX<<"    backqqX cut: "<<sizeBackqqX-eventsRemainingBackqqX<<endl;
		cout<<"IndexBackqqqqX: "<<indexBackqqqqX<<"    elemqqqqX: "<<elementNBackqqqqX<<endl<<"BackqqqqX remaining: "<<eventsRemainingBackqqqqX<<"    backqqqqX cut: "<<sizeBackqqqqX-eventsRemainingBackqqqqX<<endl;
		cout<<"IndexBackqqHX: "<<indexBackqqHX<<"    elemqqHX: "<<elementNBackqqHX<<endl<<"BackqqHX remaining: "<<eventsRemainingBackqqHX<<"    backqqHX cut: "<<sizeBackqqHX-eventsRemainingBackqqHX<<endl;
		cout<<"Significance: "<<significance<<"         maxSignificnace: "<<maxSignificance<<endl<<endl;*/
	   	if(significance>maxSignificance) 
	   	{
	   		maxSignificance=significance;
	   		defCut=cut;
	   		backRemaining=eventsRemainingBack;
	   		HHRemaining=eventsRemainingHH*weightHH;
			indexHHMaxSignificance=indexHH;
			indexBackqqXMaxSignificance=indexBackqqX;
			indexBackqqqqXMaxSignificance=indexBackqqqqX;
			indexBackqqHXMaxSignificance=indexBackqqHX;
	   	}

		//cout<<"cut: "<<cut<<endl;
		//cout<<"Significance: "<<significance<<"         maxSignificnace: "<<maxSignificance<<endl<<endl;
		//cout<<"--------"<<endl;
	}
		NNOutputHH.erase(NNOutputHH.begin(), NNOutputHH.begin()+indexHHMaxSignificance);
		NNOutputBackqqX.erase(NNOutputBackqqX.begin(), NNOutputBackqqX.begin()+indexBackqqXMaxSignificance);
		NNOutputBackqqqqX.erase(NNOutputBackqqqqX.begin(), NNOutputBackqqqqX.begin()+indexBackqqqqXMaxSignificance);
		NNOutputBackqqHX.erase(NNOutputBackqqHX.begin(), NNOutputBackqqHX.begin()+indexBackqqHXMaxSignificance);

	
}


////Function that finds the optimal cuts for the eGamma indivivual NNs considering eGamma backs
void findDefEGammaCuts(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, vector<double>& BDTqqXOutput, vector<double>& BDTqqqqXOutput, vector<double>& BDTqqHXOutput, double& defCutqq, double& defCutttbar, double& defCutZZ, double& defCutWW, double& defCutqqX, double& defCutqqqqX, double& defCutqqHX, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, double weightHH, double weightqq, double weightttbar, double weightZZ, double weightWW, double weightqqX, double weightqqqqX, double weightqqHX)
{
	vector<double> BDTqqXOutputCopy(BDTqqXOutput.begin(), BDTqqXOutput.end());
	vector<double> BDTqqqqXOutputCopy(BDTqqqqXOutput.begin(), BDTqqqqXOutput.end());
	vector<double> BDTqqHXOutputCopy(BDTqqHXOutput.begin(), BDTqqHXOutput.end());

	vector<vector<double>> NNOutputs = {BDTqqXOutputCopy, BDTqqqqXOutputCopy, BDTqqHXOutputCopy};
	double defCutqqXSen=bottomHistLimit, defCutqqqqXSen=bottomHistLimit, defCutqqHXSen=bottomHistLimit;
	vector<double> NNCuts = {defCutqqXSen, defCutqqqqXSen, defCutqqHXSen};
	vector<double> maxSignificances = {-999, -999, -999};
	vector<double> nBacks = {4, 5, 6};
	double numberOfBacks=nBacks.size();
	double cutPrev=0, cutStep=0;
	int contEventsRounds=0;
	for(int i=0; i<numberOfBacks; i++)
	{
		//cout<<endl<<endl<<endl<<endl<<"Starting round "<<i<<": "<<endl;

		for(int j=0; j<NNOutputs.size(); j++)
		{
			//cout<<"In position "<<j<<": "<<endl;
			//cout<<"Back: "<<nBacks[j]<<endl<<"number of events staring the round: ";
			for(int k=0; k<NNOutputs[j].size(); k++) if(NNOutputs[j][k] != bottomHistLimit) contEventsRounds++;
			//cout<<contEventsRounds<<endl;
			contEventsRounds=0;
			//cout<<"Before findSignificanceEGammaBacks: "<<endl<<"maxSignificance: "<<maxSignificances[j]<<endl<<"Cut in NN"<<NNCuts[j]<<endl<<endl;
			cutStep = 1.0 / pow(10, i + 1);
			findSignificanceEGammaBacks(bottomHistLimit, topHistLimit, nbin, NNOutputs[j], NNCuts[j], maxSignificances[j], nBacks[j], sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, cutStep);

			//cout<<"After findSignificanceEGammaBacks: "<<endl<<"maxSignificance: "<<maxSignificances[j]<<endl<<"Cut in NN: "<<NNCuts[j]<<endl<<endl;
			//cout<<"maxSignificance: "<<maxSignificances[j]<<endl<<"Cut in NN: "<<NNCuts[j]<<endl<<endl;
		}

		//cout<<endl<<endl<<"Ends round findSignificanceEGammaBacks and starts sortVectorQuadruple: "<<endl<<endl;
		sortVectorQuadruple(maxSignificances, NNCuts, NNOutputs, nBacks);
		//cout<<"last after sortVectorQuadruple: "<<nBacks[nBacks.size()-1]<<endl<<"first after sortVectorQuadruple: "<<nBacks[0]<<endl;
		//cout<<endl<<endl<<"Ends sortVectorQuadruple and starts event elimination: "<<endl<<endl;
		for(int j=0; j<NNOutputs[NNOutputs.size()-1].size(); j++)
		{
			//if(j==0) cout<<"range: "<<NNOutputs[NNOutputs.size()-1].size()<<endl;
			//if(j<30) cout<<"NNOutputs[NNOutputs.size()-1][j]: "<<NNOutputs[NNOutputs.size()-1][j]<<endl;
			//if(j<30) cout<<"NNCuts[NNCuts.size()-1]: "<<NNCuts[NNCuts.size()-1]<<endl;
			if(NNOutputs[NNOutputs.size()-1][j] != bottomHistLimit && NNOutputs[NNOutputs.size()-1][j] < NNCuts[(NNCuts.size()-1)])
			{
				//if(j<30) cout<<"Event eliminated"<<endl;
				for(int k=1; k<=NNOutputs.size(); k++)
				{
					//if(j<30) cout<<"Before: "<<NNOutputs[NNOutputs.size()-k][j]<<"     ";
					NNOutputs[NNOutputs.size()-k][j] = bottomHistLimit;
					//if(j<30) cout<<"After: "<<NNOutputs[NNOutputs.size()-k][j]<<endl;
				}
			} 
		}
		//cout<<"number of events finishing the round: ";
		for(int j=0; j<NNOutputs[0].size(); j++) if(NNOutputs[0][j] != bottomHistLimit) contEventsRounds++;
		//cout<<contEventsRounds<<endl;
		contEventsRounds=0;

		//cout<<endl<<endl<<"Ends event elimination and starts cut assignment: "<<endl<<endl;
		if(nBacks[nBacks.size()-1] == 4) defCutqqX = NNCuts[NNCuts.size()-1];
		else if(nBacks[nBacks.size()-1] == 5) defCutqqqqX = NNCuts[NNCuts.size()-1];
		else if(nBacks[nBacks.size()-1] == 6) defCutqqHX = NNCuts[NNCuts.size()-1];
		maxSignificances.pop_back();
		NNCuts.pop_back();
		NNOutputs.pop_back();
		nBacks.pop_back();
		//cout<<endl<<endl<<"Ends cut assignment and ends round: "<<endl<<endl;
	}
}

void cutEGamma(vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, vector<double>& BDTqqXOutput, vector<double>& BDTqqqqXOutput, vector<double>& BDTqqHXOutput, double& defCutqqX, double& defCutqqqqX, double& defCutqqHX)
{
	for(int i=0; i<BDTqqOutput.size(); i++)
	{
		if((BDTqqXOutput[i] < defCutqqX) || (BDTqqqqXOutput[i] < defCutqqqqX) || (BDTqqHXOutput[i] < defCutqqHX))
		{
			BDTqqOutput[i] = -999;
			BDTttbarOutput[i] = -999;
			BDTZZOutput[i] = -999;
			BDTWWOutput[i] = -999;
			BDTqqXOutput[i] = -999;
			BDTqqqqXOutput[i] = -999;
			BDTqqHXOutput[i] = -999;
		}
	}
}

///Function that calculates the significance of the data after applying multiple NNs.
void findSignificanceCombined(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& NN1Output, vector<double>& NN2Output, int sizeHH, int sizeBack1, int sizeBack2, int topology1, int topology2, double& defCutNN1, double& defCutNN2, double& maxSignificanceCombined, double weightHH, double weightBack1, double weightBack2, TH2F& histROCCombined, TH2F& histROCRejCombined, TH2F& histSignificanceCombined, double& totalRemainingCombined, double& HHRemainingCombined, double& backRemainingCombined)
{
	defCutNN1=bottomHistLimit;
	defCutNN2=bottomHistLimit;
	double significance;
	
	vector<double> NN1OutputHH(NN1Output.begin(), NN1Output.begin()+sizeHH);
	vector<double> NN2OutputHH(NN2Output.begin(), NN2Output.begin()+sizeHH);
	vector<double> NN1OutputBack1(NN1Output.begin()+sizeHH, NN1Output.begin()+sizeHH+sizeBack1);
	vector<double> NN2OutputBack1(NN2Output.begin()+sizeHH, NN2Output.begin()+sizeHH+sizeBack1);
	vector<double> NN1OutputBack2(NN1Output.begin()+sizeHH+sizeBack1, NN1Output.begin()+sizeHH+sizeBack1+sizeBack2);
	vector<double> NN2OutputBack2(NN2Output.begin()+sizeHH+sizeBack1, NN2Output.begin()+sizeHH+sizeBack1+sizeBack2);
	
	
	/*cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<"    NN2OutputHH size: "<<NN2OutputHH.size()<<endl;
	cout<<"NN1OutputBack1 size: "<<NN1OutputBack1.size()<<"    NN2OutputBack1 size: "<<NN2OutputBack1.size()<<endl;
	cout<<"NN1OutputBack2 size: "<<NN1OutputBack2.size()<<"    NN2OutputBack2 size: "<<NN2OutputBack2.size()<<endl<<endl<<endl;
	
	for(int i=0; i<6; i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<6; i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<4; i++) cout<<NN1OutputBack1[i]<<", ";
	cout<<endl;
	for(int i=0; i<4; i++) cout<<NN2OutputBack1[i]<<", ";
	cout<<endl;
	for(int i=0; i<3; i++) cout<<NN1OutputBack2[i]<<", ";
	cout<<endl;
	for(int i=0; i<3; i++) cout<<NN2OutputBack2[i]<<", ";
	cout<<endl;*/
	
	sortVectorPair(NN1OutputHH, NN2OutputHH);
	sortVectorPair(NN1OutputBack1, NN2OutputBack1);
	sortVectorPair(NN1OutputBack2, NN2OutputBack2);
	
	/*cout<<"Now sorted: "<<endl;
	for(int i=0; i<6; i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<6; i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<4; i++) cout<<NN1OutputBack1[i]<<", ";
	cout<<endl;
	for(int i=0; i<4; i++) cout<<NN2OutputBack1[i]<<", ";
	cout<<endl;
	for(int i=0; i<3; i++) cout<<NN1OutputBack2[i]<<", ";
	cout<<endl;
	for(int i=0; i<3; i++) cout<<NN2OutputBack2[i]<<", ";
	cout<<endl;*/
	
	int index1HH=0, index2HH=0, index1Back1=0, index2Back1=0, index1Back2=0, index2Back2=0;
	double elementN1HH=NN1OutputHH[index1HH], elementN2HH=NN2OutputHH[index2HH], elementN1Back1=NN1OutputBack1[index1Back1], elementN2Back1=NN2OutputBack1[index2Back1], elementN1Back2=NN1OutputBack2[index1Back2], elementN2Back2=NN2OutputBack2[index2Back2];
	double eventsRemainingHH=0, eventsCutHH=0, eventsRemainingBack1=0, eventsCutBack1=0, eventsRemainingBack2=0, eventsCutBack2=0, eventsRemainingBack=0, eventsCutBack=0;
	double sizeBack=sizeBack1+sizeBack2;
	
	for(double cut1=bottomHistLimit; cut1<=topHistLimit; cut1+=0.001)
	{
		while(cut1>elementN1HH && index1HH<sizeHH)
		{
			index1HH++;
			if(index1HH<sizeHH) elementN1HH = NN1OutputHH[index1HH];
		}
		while(cut1>elementN1Back1 && index1Back1<sizeBack1)
		{
			index1Back1++;
			if(index1Back1<sizeBack1) elementN1Back1 = NN1OutputBack1[index1Back1];
		}
		while(cut1>elementN1Back2 && index1Back2<sizeBack2)
		{
			index1Back2++;
			if(index1Back2<sizeBack2) elementN1Back2 = NN1OutputBack2[index1Back2];
		}
		
		vector<double> NNOutputSurvivorsHH(NN2OutputHH.begin()+index1HH, NN2OutputHH.end());
		vector<double> NNOutputSurvivorsBack1(NN2OutputBack1.begin()+index1Back1, NN2OutputBack1.end());
		vector<double> NNOutputSurvivorsBack2(NN2OutputBack2.begin()+index1Back2, NN2OutputBack2.end());
		
		sort(NNOutputSurvivorsHH.begin(), NNOutputSurvivorsHH.end());
		sort(NNOutputSurvivorsBack1.begin(), NNOutputSurvivorsBack1.end());
		sort(NNOutputSurvivorsBack2.begin(), NNOutputSurvivorsBack2.end());
		
		index2HH=0;
		index2Back1=0;
		index2Back2=0;
		if(NNOutputSurvivorsHH.size()>0) elementN2HH=NNOutputSurvivorsHH[index2HH];
		if(NNOutputSurvivorsBack1.size()>0) elementN2Back1=NNOutputSurvivorsBack1[index2Back1];
		if(NNOutputSurvivorsBack2.size()>0) elementN2Back2=NNOutputSurvivorsBack2[index2Back2];
		
		for(double cut2=bottomHistLimit; cut2<=topHistLimit; cut2+=0.001)
		{
			while(cut2>elementN2HH && index2HH<NNOutputSurvivorsHH.size())
			{
				index2HH++;
				if(index2HH<NNOutputSurvivorsHH.size()) elementN2HH = NNOutputSurvivorsHH[index2HH];
			}
			while(cut2>elementN2Back1 && index2Back1<NNOutputSurvivorsBack1.size())
			{
				index2Back1++;
				if(index2Back1<NNOutputSurvivorsBack1.size()) elementN2Back1 = NNOutputSurvivorsBack1[index2Back1];
			}
			while(cut2>elementN2Back2 && index2Back2<NNOutputSurvivorsBack2.size())
			{
				index2Back2++;
				if(index2Back2<NNOutputSurvivorsBack2.size()) elementN2Back2 = NNOutputSurvivorsBack2[index2Back2];
			}
			
			eventsRemainingHH = (sizeHH-index1HH)-(NNOutputSurvivorsHH.size()-(NNOutputSurvivorsHH.size()-index2HH));
			eventsRemainingBack1 = (sizeBack1-index1Back1)-(NNOutputSurvivorsBack1.size()-(NNOutputSurvivorsBack1.size()-index2Back1));
			eventsRemainingBack2 = (sizeBack2-index1Back2)-(NNOutputSurvivorsBack2.size()-(NNOutputSurvivorsBack2.size()-index2Back2));
			
			
			/*weightHH=1;
			weightBack1=1;
			weightBack2=1;*/
			
			eventsRemainingBack = eventsRemainingBack1*weightBack1 + eventsRemainingBack2*weightBack2;
			
			significance=(eventsRemainingHH*weightHH)/(sqrt((eventsRemainingHH*weightHH)+(eventsRemainingBack1*weightBack1+eventsRemainingBack2*weightBack2)));
			
			if(significance>maxSignificanceCombined) 
		   	{
		   		maxSignificanceCombined=significance;
		   		defCutNN1=cut1;
		   		defCutNN2=cut2;
		   		totalRemainingCombined=eventsRemainingHH*weightHH+eventsRemainingBack;
		   		backRemainingCombined=eventsRemainingBack;
		   		HHRemainingCombined=eventsRemainingHH*weightHH;
		   	}
		   	
		   	
		   	/*cout<<"CUT1: "<<cut1<<endl;
			cout<<"Index1HH: "<<index1HH<<"    elementN1HH: "<<elementN1HH<<endl;
			cout<<"Index1Back1: "<<index1Back1<<"    elementN1Back1: "<<elementN1Back1<<endl;
			cout<<"Index1Back2: "<<index1Back2<<"    elementN1Back2: "<<elementN1Back2<<endl;
			cout<<"CUT2: "<<cut2<<endl;
			cout<<"Index2HH: "<<index2HH<<"    elementN2HH: "<<elementN2HH<<endl;
			cout<<"Index2Back1: "<<index2Back1<<"    elementN2Back1: "<<elementN2Back1<<endl;
			cout<<"Index2Back2: "<<index2Back2<<"    elementN2Back2: "<<elementN2Back2<<endl;
			cout<<"Signal remaining: "<<eventsRemainingHH<<"    signal cut: "<<sizeHH-eventsRemainingHH<<endl;
			cout<<"Back1 remaining: "<<eventsRemainingBack1<<"    back1 cut: "<<sizeBack1-eventsRemainingBack1<<endl;
			cout<<"Back2 remaining: "<<eventsRemainingBack2<<"    back2 cut: "<<sizeBack2-eventsRemainingBack2<<endl;
			for(int i=0; i<NNOutputSurvivorsHH.size(); i++) cout<<NNOutputSurvivorsHH[i]<<", ";
			cout<<endl;
			for(int i=0; i<NNOutputSurvivorsBack1.size(); i++) cout<<NNOutputSurvivorsBack1[i]<<", ";
			cout<<endl;
			for(int i=0; i<NNOutputSurvivorsBack2.size(); i++) cout<<NNOutputSurvivorsBack2[i]<<", ";
			cout<<endl;
			cout<<"Significance: "<<significance<<"         maxSignificance: "<<maxSignificanceCombined<<endl<<endl;
			cout<<"--------"<<endl<<endl;*/
		
		
		}
	}
	
	/*cout<<"maxSignificanceCombined: "<<maxSignificanceCombined<<endl;
	cout<<"For cutNN1: "<<defCutNN1<<"    and cutNN2: "<<defCutNN2<<endl;
	cout<<"totalRemainingCombined: "<<totalRemainingCombined<<endl;
	cout<<"HHRemainingCombined: "<<HHRemainingCombined<<endl;
	cout<<"backRemainingCombined: "<<backRemainingCombined<<endl;*/

}

void findSignificanceCutsCombined(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& NN1Output, vector<double>& NN2Output, vector<double>& NN3Output, vector<double>& NN4Output, vector<double>& NN5Output, vector<double>& NN6Output, vector<double>& NN7Output, int sizeHH, int sizeBack1, int sizeBack2, int sizeBack3, int sizeBack4, int sizeBack5, int sizeBack6, int sizeBack7, double defCutNN1, double defCutNN2, double defCutNN3, double defCutNN4, double defCutNN5, double defCutNN6, double defCutNN7, double& maxSignificanceCombined, double weightHH, double weightBack1, double weightBack2, double weightBack3, double weightBack4, double weightBack5, double weightBack6, double weightBack7, TH2F& histROCCombined, TH2F& histROCRejCombined, TH2F& histSignificanceCombined, double& totalRemainingCombined, double& HHRemainingCombined, double& backRemainingCombined)
{
	vector<double> topologyTracker;
	for(int i=0; i<NN1Output.size(); i++)
	{
		if(i<sizeHH) topologyTracker.push_back(0);
		else if(i<sizeHH+sizeBack1) topologyTracker.push_back(1); 
		else if(i<sizeHH+sizeBack1+sizeBack2) topologyTracker.push_back(2);
		else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) topologyTracker.push_back(3);  
		else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) topologyTracker.push_back(4);
		else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) topologyTracker.push_back(5);
		else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) topologyTracker.push_back(6);
		else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) topologyTracker.push_back(7);  
	}
	
	/*cout<<endl<<"topologyTracker: ";
	for(int i=0;i<topologyTracker.size();i++) cout<<topologyTracker[i]<<", ";
	
	cout<<endl<<endl<<endl<<"NN1Output size: "<<NN1Output.size()<<"    NN2Output size: "<<NN2Output.size()<<"    NN3Output size: "<<NN3Output.size()<<"    NN4Output size: "<<NN3Output.size()<<endl<<endl;
	
	cout<<"NN1Output: ";
	for(int i=0;i<NN1Output.size();i++) cout<<NN1Output[i]<<", ";
	cout<<endl<<"NN2Output: ";
	for(int i=0;i<NN2Output.size();i++) cout<<NN2Output[i]<<", ";
	cout<<endl<<"NN3Output: ";
	for(int i=0;i<NN3Output.size();i++) cout<<NN3Output[i]<<", ";
	cout<<endl<<"NN4Output: ";
	for(int i=0;i<NN4Output.size();i++) cout<<NN4Output[i]<<", ";
	cout<<endl<<"topologyTracker: ";
	for(int i=0;i<topologyTracker.size();i++) cout<<topologyTracker[i]<<", ";*/
	
	vector<double> NN1OutputHH(NN1Output.begin(), NN1Output.begin()+sizeHH);
	vector<double> NN1OutputBack(NN1Output.begin()+sizeHH, NN1Output.end());
	vector<double> NN2OutputHH(NN2Output.begin(), NN2Output.begin()+sizeHH);
	vector<double> NN2OutputBack(NN2Output.begin()+sizeHH, NN2Output.end());
	vector<double> NN3OutputHH(NN3Output.begin(), NN3Output.begin()+sizeHH);
	vector<double> NN3OutputBack(NN3Output.begin()+sizeHH, NN3Output.end());
	vector<double> NN4OutputHH(NN4Output.begin(), NN4Output.begin()+sizeHH);
	vector<double> NN4OutputBack(NN4Output.begin()+sizeHH, NN4Output.end());
	vector<double> NN5OutputHH(NN5Output.begin(), NN5Output.begin()+sizeHH);
	vector<double> NN5OutputBack(NN5Output.begin()+sizeHH, NN5Output.end());
	vector<double> NN6OutputHH(NN6Output.begin(), NN6Output.begin()+sizeHH);
	vector<double> NN6OutputBack(NN6Output.begin()+sizeHH, NN6Output.end());
	vector<double> NN7OutputHH(NN7Output.begin(), NN7Output.begin()+sizeHH);
	vector<double> NN7OutputBack(NN7Output.begin()+sizeHH, NN7Output.end());
	vector<double> topologyTrackerHH(topologyTracker.begin(), topologyTracker.begin()+sizeHH);
	vector<double> topologyTrackerBack(topologyTracker.begin()+sizeHH, topologyTracker.end());
	
	/*cout<<endl<<endl<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<"topologyTrackerHH size: "<<topologyTrackerHH.size()<<endl<<"topologyTrackerBack size: "<<topologyTrackerBack.size()<<endl<<endl;
	
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	sortVectorEightuple(NN1OutputHH, NN2OutputHH, NN3OutputHH, NN4OutputHH, NN5OutputHH, NN6OutputHH, NN7OutputHH, topologyTrackerHH);
	sortVectorEightuple(NN1OutputBack, NN2OutputBack, NN3OutputBack, NN4OutputBack, NN5OutputBack, NN6OutputBack, NN7OutputBack, topologyTrackerBack);
	
	/*cout<<endl<<endl<<"Now after qq sorting: ";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	//cout<<"HH starting: "<<topologyTrackerHH.size()<<endl;
	
	double event;
	
	int indexNN1OutputHH = 0;
	event = NN1OutputHH[indexNN1OutputHH];
	while(event<defCutNN1 && indexNN1OutputHH<NN1OutputHH.size())
	{
		indexNN1OutputHH++;
		event = NN1OutputHH[indexNN1OutputHH];
	}
	
	int indexNN1OutputBack = 0;
	event = NN1OutputBack[indexNN1OutputBack];
	while(event<defCutNN1 && indexNN1OutputBack<NN1OutputBack.size())
	{
		indexNN1OutputBack++;
		event = NN1OutputBack[indexNN1OutputBack];
	}
	
	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN1OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN1OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN1OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN1OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN1OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN1OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN1OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN1OutputBack);
	NN5OutputHH.erase(NN5OutputHH.begin(), NN5OutputHH.begin()+indexNN1OutputHH);
	NN5OutputBack.erase(NN5OutputBack.begin(), NN5OutputBack.begin()+indexNN1OutputBack);
	NN6OutputHH.erase(NN6OutputHH.begin(), NN6OutputHH.begin()+indexNN1OutputHH);
	NN6OutputBack.erase(NN6OutputBack.begin(), NN6OutputBack.begin()+indexNN1OutputBack);
	NN7OutputHH.erase(NN7OutputHH.begin(), NN7OutputHH.begin()+indexNN1OutputHH);
	NN7OutputBack.erase(NN7OutputBack.begin(), NN7OutputBack.begin()+indexNN1OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN1OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN1OutputBack);
	

	/*cout<<endl<<endl<<"After qq cutting: "<<endl;
	cout<<"indexNN1OutputHH: "<<indexNN1OutputHH<<"    indexNN1OutputBack: "<<indexNN1OutputBack<<endl;
	cout<<"HH remaining: "<<topologyTrackerHH.size()<<endl;*/
	/*cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<endl;
	
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	
	
	sortVectorEightuple(NN2OutputHH, NN3OutputHH, NN4OutputHH, NN5OutputHH, NN6OutputHH, NN7OutputHH, NN1OutputHH, topologyTrackerHH);
	sortVectorEightuple(NN2OutputBack, NN3OutputBack, NN4OutputBack, NN5OutputBack, NN6OutputBack, NN7OutputBack, NN1OutputBack, topologyTrackerBack);
	
	/*cout<<endl<<endl<<"Now after ttbar sorting: ";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	int indexNN2OutputHH = 0;
	event = NN2OutputHH[indexNN2OutputHH];
	while(event<defCutNN2 && indexNN2OutputHH<NN2OutputHH.size())
	{
		indexNN2OutputHH++;
		event = NN2OutputHH[indexNN2OutputHH];
	}
	
	int indexNN2OutputBack = 0;
	event = NN2OutputBack[indexNN2OutputBack];
	while(event<defCutNN2 && indexNN2OutputBack<NN2OutputBack.size())
	{
		indexNN2OutputBack++;
		event = NN2OutputBack[indexNN2OutputBack];
	}
	
	//cout<<endl<<endl<<"indexNN2OutputHH: "<<indexNN2OutputHH<<"    indexNN2OutputBack: "<<indexNN2OutputBack<<endl;
	
	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN2OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN2OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN2OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN2OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN2OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN2OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN2OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN2OutputBack);
	NN5OutputHH.erase(NN5OutputHH.begin(), NN5OutputHH.begin()+indexNN2OutputHH);
	NN5OutputBack.erase(NN5OutputBack.begin(), NN5OutputBack.begin()+indexNN2OutputBack);
	NN6OutputHH.erase(NN6OutputHH.begin(), NN6OutputHH.begin()+indexNN2OutputHH);
	NN6OutputBack.erase(NN6OutputBack.begin(), NN6OutputBack.begin()+indexNN2OutputBack);
	NN7OutputHH.erase(NN7OutputHH.begin(), NN7OutputHH.begin()+indexNN2OutputHH);
	NN7OutputBack.erase(NN7OutputBack.begin(), NN7OutputBack.begin()+indexNN2OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN2OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN2OutputBack);
	
	/*cout<<endl<<endl<<"After ttbar cutting: "<<endl;
	cout<<"indexNN2OutputHH: "<<indexNN2OutputHH<<"    indexNN2OutputBack: "<<indexNN2OutputBack<<endl;
	cout<<"HH remaining: "<<topologyTrackerHH.size()<<endl;*/
	
	//cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<endl;
	
	
	sortVectorEightuple(NN3OutputHH, NN4OutputHH, NN5OutputHH, NN6OutputHH, NN7OutputHH, NN1OutputHH, NN2OutputHH, topologyTrackerHH);
	sortVectorEightuple(NN3OutputBack, NN4OutputBack, NN5OutputBack, NN6OutputBack, NN7OutputBack, NN1OutputBack, NN2OutputBack, topologyTrackerBack);
	
	/*cout<<endl<<endl<<"Now after ZZ sorting: ";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	
	
	int indexNN3OutputHH = 0;
	event = NN3OutputHH[indexNN3OutputHH];
	while(event<defCutNN3 && indexNN3OutputHH<NN3OutputHH.size())
	{
		indexNN3OutputHH++;
		event = NN3OutputHH[indexNN3OutputHH];
	}
	
	int indexNN3OutputBack = 0;
	event = NN3OutputBack[indexNN3OutputBack];
	while(event<defCutNN3 && indexNN3OutputBack<NN3OutputBack.size())
	{
		indexNN3OutputBack++;
		event = NN3OutputBack[indexNN3OutputBack];
	}
	
	//cout<<endl<<endl<<"indexNN3OutputHH: "<<indexNN3OutputHH<<"    indexNN3OutputBack: "<<indexNN3OutputBack<<endl;
	
	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN3OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN3OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN3OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN3OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN3OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN3OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN3OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN3OutputBack);
	NN5OutputHH.erase(NN5OutputHH.begin(), NN5OutputHH.begin()+indexNN3OutputHH);
	NN5OutputBack.erase(NN5OutputBack.begin(), NN5OutputBack.begin()+indexNN3OutputBack);
	NN6OutputHH.erase(NN6OutputHH.begin(), NN6OutputHH.begin()+indexNN3OutputHH);
	NN6OutputBack.erase(NN6OutputBack.begin(), NN6OutputBack.begin()+indexNN3OutputBack);
	NN7OutputHH.erase(NN7OutputHH.begin(), NN7OutputHH.begin()+indexNN3OutputHH);
	NN7OutputBack.erase(NN7OutputBack.begin(), NN7OutputBack.begin()+indexNN3OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN3OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN3OutputBack);
	
	/*cout<<endl<<endl<<"After ZZ cutting: "<<endl;
	cout<<"indexNN3OutputHH: "<<indexNN3OutputHH<<"    indexNN3OutputBack: "<<indexNN3OutputBack<<endl;
	cout<<"HH remaining: "<<topologyTrackerHH.size()<<endl;*/
	
	/*cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<endl;*/
	
	sortVectorEightuple(NN4OutputHH, NN5OutputHH, NN6OutputHH, NN7OutputHH, NN1OutputHH, NN2OutputHH, NN3OutputHH, topologyTrackerHH);
	sortVectorEightuple(NN4OutputBack, NN5OutputBack, NN6OutputBack, NN7OutputBack, NN1OutputBack, NN2OutputBack, NN3OutputBack, topologyTrackerBack);
	
	/*cout<<endl<<endl<<"Now after WW sorting: ";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	
	
	int indexNN4OutputHH = 0;
	event = NN4OutputHH[indexNN4OutputHH];
	while(event<defCutNN4 && indexNN4OutputHH<NN4OutputHH.size())
	{
		indexNN4OutputHH++;
		event = NN4OutputHH[indexNN4OutputHH];
	}
	
	int indexNN4OutputBack = 0;
	event = NN4OutputBack[indexNN4OutputBack];
	while(event<defCutNN4 && indexNN4OutputBack<NN4OutputBack.size())
	{
		indexNN4OutputBack++;
		event = NN4OutputBack[indexNN4OutputBack];
	}
	
	//cout<<endl<<endl<<"indexNN4OutputHH: "<<indexNN4OutputHH<<"    indexNN4OutputBack: "<<indexNN4OutputBack<<endl;
	
	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN4OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN4OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN4OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN4OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN4OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN4OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN4OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN4OutputBack);
	NN5OutputHH.erase(NN5OutputHH.begin(), NN5OutputHH.begin()+indexNN4OutputHH);
	NN5OutputBack.erase(NN5OutputBack.begin(), NN5OutputBack.begin()+indexNN4OutputBack);
	NN6OutputHH.erase(NN6OutputHH.begin(), NN6OutputHH.begin()+indexNN4OutputHH);
	NN6OutputBack.erase(NN6OutputBack.begin(), NN6OutputBack.begin()+indexNN4OutputBack);
	NN7OutputHH.erase(NN7OutputHH.begin(), NN7OutputHH.begin()+indexNN4OutputHH);
	NN7OutputBack.erase(NN7OutputBack.begin(), NN7OutputBack.begin()+indexNN4OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN4OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN4OutputBack);

	////UPDATES EGAMMA BACKGROUNDS
	sortVectorEightuple(NN5OutputHH, NN6OutputHH, NN7OutputHH, NN1OutputHH, NN2OutputHH, NN3OutputHH, NN4OutputHH, topologyTrackerHH);
	sortVectorEightuple(NN5OutputBack, NN6OutputBack, NN7OutputBack, NN1OutputBack, NN2OutputBack, NN3OutputBack, NN4OutputBack, topologyTrackerBack);

	int indexNN5OutputHH = 0;
	event = NN5OutputHH[indexNN5OutputHH];
	while(event<defCutNN5 && indexNN5OutputHH<NN5OutputHH.size())
	{
		indexNN5OutputHH++;
		event = NN5OutputHH[indexNN5OutputHH];
	}
	
	int indexNN5OutputBack = 0;
	event = NN5OutputBack[indexNN5OutputBack];
	while(event<defCutNN5 && indexNN5OutputBack<NN5OutputBack.size())
	{
		indexNN5OutputBack++;
		event = NN5OutputBack[indexNN5OutputBack];
	}

	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN5OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN5OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN5OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN5OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN5OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN5OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN5OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN5OutputBack);
	NN5OutputHH.erase(NN5OutputHH.begin(), NN5OutputHH.begin()+indexNN5OutputHH);
	NN5OutputBack.erase(NN5OutputBack.begin(), NN5OutputBack.begin()+indexNN5OutputBack);
	NN6OutputHH.erase(NN6OutputHH.begin(), NN6OutputHH.begin()+indexNN5OutputHH);
	NN6OutputBack.erase(NN6OutputBack.begin(), NN6OutputBack.begin()+indexNN5OutputBack);
	NN7OutputHH.erase(NN7OutputHH.begin(), NN7OutputHH.begin()+indexNN5OutputHH);
	NN7OutputBack.erase(NN7OutputBack.begin(), NN7OutputBack.begin()+indexNN5OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN5OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN5OutputBack);


	sortVectorEightuple(NN6OutputHH, NN7OutputHH, NN1OutputHH, NN2OutputHH, NN3OutputHH, NN4OutputHH, NN5OutputHH, topologyTrackerHH);
	sortVectorEightuple(NN6OutputBack, NN7OutputBack, NN1OutputBack, NN2OutputBack, NN3OutputBack, NN4OutputBack, NN5OutputBack, topologyTrackerBack);

	int indexNN6OutputHH = 0;
	event = NN6OutputHH[indexNN6OutputHH];
	while(event<defCutNN6 && indexNN6OutputHH<NN6OutputHH.size())
	{
		indexNN6OutputHH++;
		event = NN6OutputHH[indexNN6OutputHH];
	}
	
	int indexNN6OutputBack = 0;
	event = NN6OutputBack[indexNN6OutputBack];
	while(event<defCutNN6 && indexNN6OutputBack<NN6OutputBack.size())
	{
		indexNN6OutputBack++;
		event = NN6OutputBack[indexNN6OutputBack];
	}

	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN6OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN6OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN6OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN6OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN6OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN6OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN6OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN6OutputBack);
	NN5OutputHH.erase(NN5OutputHH.begin(), NN5OutputHH.begin()+indexNN6OutputHH);
	NN5OutputBack.erase(NN5OutputBack.begin(), NN5OutputBack.begin()+indexNN6OutputBack);
	NN6OutputHH.erase(NN6OutputHH.begin(), NN6OutputHH.begin()+indexNN6OutputHH);
	NN6OutputBack.erase(NN6OutputBack.begin(), NN6OutputBack.begin()+indexNN6OutputBack);
	NN7OutputHH.erase(NN7OutputHH.begin(), NN7OutputHH.begin()+indexNN6OutputHH);
	NN7OutputBack.erase(NN7OutputBack.begin(), NN7OutputBack.begin()+indexNN6OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN6OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN6OutputBack);


	sortVectorEightuple(NN7OutputHH, NN1OutputHH, NN2OutputHH, NN3OutputHH, NN4OutputHH, NN5OutputHH, NN6OutputHH, topologyTrackerHH);
	sortVectorEightuple(NN7OutputBack, NN1OutputBack, NN2OutputBack, NN3OutputBack, NN4OutputBack, NN5OutputBack, NN6OutputBack, topologyTrackerBack);

	int indexNN7OutputHH = 0;
	event = NN7OutputHH[indexNN7OutputHH];
	while(event<defCutNN7 && indexNN7OutputHH<NN7OutputHH.size())
	{
		indexNN7OutputHH++;
		event = NN7OutputHH[indexNN7OutputHH];
	}
	
	int indexNN7OutputBack = 0;
	event = NN7OutputBack[indexNN7OutputBack];
	while(event<defCutNN7 && indexNN7OutputBack<NN7OutputBack.size())
	{
		indexNN7OutputBack++;
		event = NN7OutputBack[indexNN7OutputBack];
	}

	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN7OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN7OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN7OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN7OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN7OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN7OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN7OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN7OutputBack);
	NN5OutputHH.erase(NN5OutputHH.begin(), NN5OutputHH.begin()+indexNN7OutputHH);
	NN5OutputBack.erase(NN5OutputBack.begin(), NN5OutputBack.begin()+indexNN7OutputBack);
	NN6OutputHH.erase(NN6OutputHH.begin(), NN6OutputHH.begin()+indexNN7OutputHH);
	NN6OutputBack.erase(NN6OutputBack.begin(), NN6OutputBack.begin()+indexNN7OutputBack);
	NN7OutputHH.erase(NN7OutputHH.begin(), NN7OutputHH.begin()+indexNN7OutputHH);
	NN7OutputBack.erase(NN7OutputBack.begin(), NN7OutputBack.begin()+indexNN7OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN7OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN7OutputBack);
	
	/*cout<<endl<<endl<<"After WW cutting: "<<endl;
	cout<<"indexNN4OutputHH: "<<indexNN4OutputHH<<"    indexNN4OutputBack: "<<indexNN4OutputBack<<endl;
	cout<<"HH remaining: "<<topologyTrackerHH.size()<<endl;*/
	
	/*cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<endl;*/

	/*cout<<endl<<endl<<"Final vectors";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	int contBack1=0, contBack2=0, contBack3=0, contBack4=0, contBack5=0, contBack6=0, contBack7=0;
	for(int i=0; i<topologyTrackerBack.size(); i++)
	{
		if(topologyTrackerBack[i]==1) contBack1++;
		else if(topologyTrackerBack[i]==2) contBack2++;
		else if(topologyTrackerBack[i]==3) contBack3++;
		else if(topologyTrackerBack[i]==4) contBack4++;
		else if(topologyTrackerBack[i]==5) contBack5++;
		else if(topologyTrackerBack[i]==6) contBack6++;
		else if(topologyTrackerBack[i]==7) contBack7++;
	}
	
	HHRemainingCombined = NN1OutputHH.size()*weightHH;
	backRemainingCombined = contBack1*weightBack1+contBack2*weightBack2+contBack3*weightBack3+contBack4*weightBack4+contBack5*weightBack5+contBack6*weightBack6+contBack7*weightBack7;
	int HHRemainingCombinedUnweighted = NN1OutputHH.size();
	int backRemainingCombinedUnweighted = contBack1+contBack2+contBack3+contBack4+contBack5+contBack6+contBack7;
	
	maxSignificanceCombined = HHRemainingCombined/sqrt(HHRemainingCombined+backRemainingCombined);
	
	cout<<"remainingHH unweighted: "<<HHRemainingCombinedUnweighted<<endl<<"remainingBackUnweighted: "<<backRemainingCombinedUnweighted<<endl;
	cout<<"remainingHH: "<<HHRemainingCombined<<endl<<"remainingBack: "<<backRemainingCombined<<endl;
	
	
	return;
}

void findSignificanceCutsCombinedCheck(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& NN1Output, vector<double>& NN2Output, vector<double>& NN3Output, vector<double>& NN4Output, vector<double>& NN5Output, vector<double>& NN6Output, vector<double>& NN7Output, int sizeHH, int sizeBack1, int sizeBack2, int sizeBack3, int sizeBack4, int sizeBack5, int sizeBack6, int sizeBack7, double defCutNN1, double defCutNN2, double defCutNN3, double defCutNN4, double defCutNN5, double defCutNN6, double defCutNN7, double& maxSignificanceCombined, double weightHH, double weightBack1, double weightBack2, double weightBack3, double weightBack4, double weightBack5, double weightBack6, double weightBack7, TH2F& histROCCombined, TH2F& histROCRejCombined, TH2F& histSignificanceCombined, double& totalRemainingCombined, double& HHRemainingCombined, double& backRemainingCombined)
{
	vector<double> NN1OutputCopy(NN1Output.begin(), NN1Output.end());
	vector<double> NN2OutputCopy(NN2Output.begin(), NN2Output.end());
	vector<double> NN3OutputCopy(NN3Output.begin(), NN3Output.end());
	vector<double> NN4OutputCopy(NN4Output.begin(), NN4Output.end());
	vector<double> NN5OutputCopy(NN5Output.begin(), NN5Output.end());
	vector<double> NN6OutputCopy(NN6Output.begin(), NN6Output.end());
	vector<double> NN7OutputCopy(NN7Output.begin(), NN7Output.end());
	
	double contHHSurviving=0, contBack1Surviving=0, contBack2Surviving=0, contBack3Surviving=0, contBack4Surviving=0, contBack5Surviving=0, contBack6Surviving=0, contBack7Surviving=0, contHHInitial=0, contBack1Initial=0, contBack2Initial=0, contBack3Initial=0, contBack4Initial=0, contBack5Initial=0, contBack6Initial=0, contBack7Initial=0;
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999)
		{
			if(i<sizeHH) contHHInitial++;
			else if(i<sizeHH+sizeBack1) contBack1Initial++;
			else if(i<sizeHH+sizeBack1+sizeBack2) contBack2Initial++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3Initial++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4Initial++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5Initial++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6Initial++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7Initial++;
		}
	}
	
	cout<<"Initial HH events (before any NN): "<<contHHInitial*weightHH<<endl<<"Initial qq events (before any NN): "<<contBack1Initial*weightBack1<<endl<<"Initial ttbar events (before any NN): "<<contBack2Initial*weightBack2<<endl<<"Initial ZZ events (before any NN): "<<contBack3Initial*weightBack3<<endl<<"Initial WW events (before any NN): "<<contBack4Initial*weightBack4<<endl<<"Initial qqX events (before any NN): "<<contBack5Initial*weightBack5<<endl<<"Initial qqqqX events (before any NN): "<<contBack6Initial*weightBack6<<endl<<"Initial qqHX events (before any NN): "<<contBack7Initial*weightBack7<<endl<<endl;
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]<defCutNN1 && NN1OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
			NN5OutputCopy[i]=-999;
			NN6OutputCopy[i]=-999;
			NN7OutputCopy[i]=-999;
		}
	}
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999)
		{
			if(i<sizeHH) contHHSurviving++;
			else if(i<sizeHH+sizeBack1) contBack1Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2) contBack2Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7Surviving++;
		}
	}
	
	cout<<"Initial HH events (before NNqq): "<<contHHInitial*weightHH<<endl<<"Surviving HH events (after NNqq): "<<contHHSurviving*weightHH<<endl<<"Initial qq events (before NNqq): "<<contBack1Initial*weightBack1<<endl<<"Surviving qq events (after NNqq): "<<contBack1Surviving*weightBack1<<endl<<"Initial ttbar events (before NNqq): "<<contBack2Initial*weightBack2<<endl<<"Surviving ttbar events (after NNqq): "<<contBack2Surviving*weightBack2<<endl<<"Initial ZZ events (before NNqq): "<<contBack3Initial*weightBack3<<endl<<"Surviving ZZ events (after NNqq): "<<contBack3Surviving*weightBack3<<endl<<"Initial WW events (before NNqq): "<<contBack4Initial*weightBack4<<endl<<"Surviving WW events (after NNqq): "<<contBack4Surviving*weightBack4<<endl<<"Initial qqX events (before NNqq): "<<contBack5Initial*weightBack5<<endl<<"Surviving qqX events (after NNqq): "<<contBack5Surviving*weightBack5<<endl<<"Initial qqqqX events (before NNqq): "<<contBack6Initial*weightBack6<<endl<<"Surviving qqqqX events (after NNqq): "<<contBack6Surviving*weightBack6<<endl<<"Initial qqHX events (before NNqq): "<<contBack7Initial*weightBack7<<endl<<"Surviving qqHX events (after NNqq): "<<contBack7Surviving*weightBack7<<endl<<endl;
	
	contHHInitial=contHHSurviving;
	contBack1Initial=contBack1Surviving;
	contBack2Initial=contBack2Surviving;
	contBack3Initial=contBack3Surviving;
	contBack4Initial=contBack4Surviving;
	contBack5Initial=contBack5Surviving;
	contBack6Initial=contBack6Surviving;
	contBack7Initial=contBack7Surviving;
	
	contHHSurviving=0;
	contBack1Surviving=0;
	contBack2Surviving=0;
	contBack3Surviving=0;
	contBack4Surviving=0;
	contBack5Surviving=0;
	contBack6Surviving=0;
	contBack7Surviving=0;
	
	for(int i=0; i<NN2OutputCopy.size(); i++)
	{
		if(NN2OutputCopy[i]<defCutNN2 && NN2OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
			NN5OutputCopy[i]=-999;
			NN6OutputCopy[i]=-999;
			NN7OutputCopy[i]=-999;
		}
	}
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999)
		{
			if(i<sizeHH) contHHSurviving++;
			else if(i<sizeHH+sizeBack1) contBack1Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2) contBack2Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7Surviving++;
		}
	}
	
	cout<<"Initial HH events (before NNttbar): "<<contHHInitial*weightHH<<endl<<"Surviving HH events (after NNttbar): "<<contHHSurviving*weightHH<<endl<<"Initial qq events (before NNttbar): "<<contBack1Initial*weightBack1<<endl<<"Surviving qq events (after NNttbar): "<<contBack1Surviving*weightBack1<<endl<<"Initial ttbar events (before NNttbar): "<<contBack2Initial*weightBack2<<endl<<"Surviving ttbar events (after NNttbar): "<<contBack2Surviving*weightBack2<<endl<<"Initial ZZ events (before NNttbar): "<<contBack3Initial*weightBack3<<endl<<"Surviving ZZ events (after NNttbar): "<<contBack3Surviving*weightBack3<<endl<<"Initial WW events (before NNttbar): "<<contBack4Initial*weightBack4<<endl<<"Surviving WW events (after NNttbar): "<<contBack4Surviving*weightBack4<<endl<<"Initial qqX events (before NNttbar): "<<contBack5Initial*weightBack5<<endl<<"Surviving qqX events (after NNttbar): "<<contBack5Surviving*weightBack5<<endl<<"Initial qqqqX events (before NNttbar): "<<contBack6Initial*weightBack6<<endl<<"Surviving qqqqX events (after NNttbar): "<<contBack6Surviving*weightBack6<<endl<<"Initial qqHX events (before NNttbar): "<<contBack7Initial*weightBack7<<endl<<"Surviving qqHX events (after NNttbar): "<<contBack7Surviving*weightBack7<<endl<<endl;
	
	contHHInitial=contHHSurviving;
	contBack1Initial=contBack1Surviving;
	contBack2Initial=contBack2Surviving;
	contBack3Initial=contBack3Surviving;
	contBack4Initial=contBack4Surviving;
	contBack5Initial=contBack5Surviving;
	contBack6Initial=contBack6Surviving;
	contBack7Initial=contBack7Surviving;
	
	contHHSurviving=0;
	contBack1Surviving=0;
	contBack2Surviving=0;
	contBack3Surviving=0;
	contBack4Surviving=0;
	contBack5Surviving=0;
	contBack6Surviving=0;
	contBack7Surviving=0;
	
	for(int i=0; i<NN3OutputCopy.size(); i++)
	{
		if(NN3OutputCopy[i]<defCutNN3 && NN3OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
			NN5OutputCopy[i]=-999;
			NN6OutputCopy[i]=-999;
			NN7OutputCopy[i]=-999;
		}
	}
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999)
		{
			if(i<sizeHH) contHHSurviving++;
			else if(i<sizeHH+sizeBack1) contBack1Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2) contBack2Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7Surviving++;
		}
	}
	
	cout<<"Initial HH events (before NNZZ): "<<contHHInitial*weightHH<<endl<<"Surviving HH events (after NNZZ): "<<contHHSurviving*weightHH<<endl<<"Initial qq events (before NNZZ): "<<contBack1Initial*weightBack1<<endl<<"Surviving qq events (after NNZZ): "<<contBack1Surviving*weightBack1<<endl<<"Initial ttbar events (before NNZZ): "<<contBack2Initial*weightBack2<<endl<<"Surviving ttbar events (after NNZZ): "<<contBack2Surviving*weightBack2<<endl<<"Initial ZZ events (before NNZZ): "<<contBack3Initial*weightBack3<<endl<<"Surviving ZZ events (after NNZZ): "<<contBack3Surviving*weightBack3<<endl<<"Initial WW events (before NNZZ): "<<contBack4Initial*weightBack4<<endl<<"Surviving WW events (after NNZZ): "<<contBack4Surviving*weightBack4<<endl<<"Initial qqX events (before NNZZ): "<<contBack5Initial*weightBack5<<endl<<"Surviving qqX events (after NNZZ): "<<contBack5Surviving*weightBack5<<endl<<"Initial qqqqX events (before NNZZ): "<<contBack6Initial*weightBack6<<endl<<"Surviving qqqqX events (after NNZZ): "<<contBack6Surviving*weightBack6<<endl<<"Initial qqHX events (before NNZZ): "<<contBack7Initial*weightBack7<<endl<<"Surviving qqHX events (after NNZZ): "<<contBack7Surviving*weightBack7<<endl<<endl;
	
	contHHInitial=contHHSurviving;
	contBack1Initial=contBack1Surviving;
	contBack2Initial=contBack2Surviving;
	contBack3Initial=contBack3Surviving;
	contBack4Initial=contBack4Surviving;
	contBack5Initial=contBack5Surviving;
	contBack6Initial=contBack6Surviving;
	contBack7Initial=contBack7Surviving;
	
	contHHSurviving=0;
	contBack1Surviving=0;
	contBack2Surviving=0;
	contBack3Surviving=0;
	contBack4Surviving=0;
	contBack5Surviving=0;
	contBack6Surviving=0;
	contBack7Surviving=0;
	
	for(int i=0; i<NN4OutputCopy.size(); i++)
	{
		if(NN4OutputCopy[i]<defCutNN4 && NN4OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
			NN5OutputCopy[i]=-999;
			NN6OutputCopy[i]=-999;
			NN7OutputCopy[i]=-999;
		}
	}
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999)
		{
			if(i<sizeHH) contHHSurviving++;
			else if(i<sizeHH+sizeBack1) contBack1Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2) contBack2Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7Surviving++;
		}
	}
	
	cout<<"Initial HH events (before NNWW): "<<contHHInitial*weightHH<<endl<<"Surviving HH events (after NNWW): "<<contHHSurviving*weightHH<<endl<<"Initial qq events (before NNWW): "<<contBack1Initial*weightBack1<<endl<<"Surviving qq events (after NNWW): "<<contBack1Surviving*weightBack1<<endl<<"Initial ttbar events (before NNWW): "<<contBack2Initial*weightBack2<<endl<<"Surviving ttbar events (after NNWW): "<<contBack2Surviving*weightBack2<<endl<<"Initial ZZ events (before NNWW): "<<contBack3Initial*weightBack3<<endl<<"Surviving ZZ events (after NNWW): "<<contBack3Surviving*weightBack3<<endl<<"Initial WW events (before NNWW): "<<contBack4Initial*weightBack4<<endl<<"Surviving WW events (after NNWW): "<<contBack4Surviving*weightBack4<<endl<<"Initial qqX events (before NNWW): "<<contBack5Initial*weightBack5<<endl<<"Surviving qqX events (after NNWW): "<<contBack5Surviving*weightBack5<<endl<<"Initial qqqqX events (before NNWW): "<<contBack6Initial*weightBack6<<endl<<"Surviving qqqqX events (after NNWW): "<<contBack6Surviving*weightBack6<<endl<<"Initial qqHX events (before NNWW): "<<contBack7Initial*weightBack7<<endl<<"Surviving qqHX events (after NNWW): "<<contBack7Surviving*weightBack7<<endl<<endl;
	
	contHHInitial=contHHSurviving;
	contBack1Initial=contBack1Surviving;
	contBack2Initial=contBack2Surviving;
	contBack3Initial=contBack3Surviving;
	contBack4Initial=contBack4Surviving;
	contBack5Initial=contBack5Surviving;
	contBack6Initial=contBack6Surviving;
	contBack7Initial=contBack7Surviving;
	
	contHHSurviving=0;
	contBack1Surviving=0;
	contBack2Surviving=0;
	contBack3Surviving=0;
	contBack4Surviving=0;
	contBack5Surviving=0;
	contBack6Surviving=0;
	contBack7Surviving=0;

	for(int i=0; i<NN5OutputCopy.size(); i++)
	{
		if(NN5OutputCopy[i]<defCutNN5 && NN5OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
			NN5OutputCopy[i]=-999;
			NN6OutputCopy[i]=-999;
			NN7OutputCopy[i]=-999;
		}
	}
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999)
		{
			if(i<sizeHH) contHHSurviving++;
			else if(i<sizeHH+sizeBack1) contBack1Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2) contBack2Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7Surviving++;
		}
	}

	cout<<"Initial HH events (before NNqqX): "<<contHHInitial*weightHH<<endl<<"Surviving HH events (after NNqqX): "<<contHHSurviving*weightHH<<endl<<"Initial qq events (before NNqqX): "<<contBack1Initial*weightBack1<<endl<<"Surviving qq events (after NNqqX): "<<contBack1Surviving*weightBack1<<endl<<"Initial ttbar events (before NNqqX): "<<contBack2Initial*weightBack2<<endl<<"Surviving ttbar events (after NNqqX): "<<contBack2Surviving*weightBack2<<endl<<"Initial ZZ events (before NNqqX): "<<contBack3Initial*weightBack3<<endl<<"Surviving ZZ events (after NNqqX): "<<contBack3Surviving*weightBack3<<endl<<"Initial WW events (before NNqqX): "<<contBack4Initial*weightBack4<<endl<<"Surviving WW events (after NNqqX): "<<contBack4Surviving*weightBack4<<endl<<"Initial qqX events (before NNqqX): "<<contBack5Initial*weightBack5<<endl<<"Surviving qqX events (after NNqqX): "<<contBack5Surviving*weightBack5<<endl<<"Initial qqqqX events (before NNqqX): "<<contBack6Initial*weightBack6<<endl<<"Surviving qqqqX events (after NNqqX): "<<contBack6Surviving*weightBack6<<endl<<"Initial qqHX events (before NNqqX): "<<contBack7Initial*weightBack7<<endl<<"Surviving qqHX events (after NNqqX): "<<contBack7Surviving*weightBack7<<endl<<endl;

	contHHInitial=contHHSurviving;
	contBack1Initial=contBack1Surviving;
	contBack2Initial=contBack2Surviving;
	contBack3Initial=contBack3Surviving;
	contBack4Initial=contBack4Surviving;
	contBack5Initial=contBack5Surviving;
	contBack6Initial=contBack6Surviving;
	contBack7Initial=contBack7Surviving;
	
	contHHSurviving=0;
	contBack1Surviving=0;
	contBack2Surviving=0;
	contBack3Surviving=0;
	contBack4Surviving=0;
	contBack5Surviving=0;
	contBack6Surviving=0;
	contBack7Surviving=0;

	for(int i=0; i<NN6OutputCopy.size(); i++)
	{
		if(NN6OutputCopy[i]<defCutNN6 && NN6OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
			NN5OutputCopy[i]=-999;
			NN6OutputCopy[i]=-999;
			NN7OutputCopy[i]=-999;
		}
	}
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999)
		{
			if(i<sizeHH) contHHSurviving++;
			else if(i<sizeHH+sizeBack1) contBack1Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2) contBack2Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7Surviving++;
		}
	}

	cout<<"Initial HH events (before NNqqqqX): "<<contHHInitial*weightHH<<endl<<"Surviving HH events (after NNqqqqX): "<<contHHSurviving*weightHH<<endl<<"Initial qq events (before NNqqqqX): "<<contBack1Initial*weightBack1<<endl<<"Surviving qq events (after NNqqqqX): "<<contBack1Surviving*weightBack1<<endl<<"Initial ttbar events (before NNqqqqX): "<<contBack2Initial*weightBack2<<endl<<"Surviving ttbar events (after NNqqqqX): "<<contBack2Surviving*weightBack2<<endl<<"Initial ZZ events (before NNqqqqX): "<<contBack3Initial*weightBack3<<endl<<"Surviving ZZ events (after NNqqqqX): "<<contBack3Surviving*weightBack3<<endl<<"Initial WW events (before NNqqqqX): "<<contBack4Initial*weightBack4<<endl<<"Surviving WW events (after NNqqqqX): "<<contBack4Surviving*weightBack4<<endl<<"Initial qqX events (before NNqqqqX): "<<contBack5Initial*weightBack5<<endl<<"Surviving qqX events (after NNqqqqX): "<<contBack5Surviving*weightBack5<<endl<<"Initial qqqqX events (before NNqqqqX): "<<contBack6Initial*weightBack6<<endl<<"Surviving qqqqX events (after NNqqqqX): "<<contBack6Surviving*weightBack6<<endl<<"Initial qqHX events (before NNqqqqX): "<<contBack7Initial*weightBack7<<endl<<"Surviving qqHX events (after NNqqqqX): "<<contBack7Surviving*weightBack7<<endl<<endl;

	contHHInitial=contHHSurviving;
	contBack1Initial=contBack1Surviving;
	contBack2Initial=contBack2Surviving;
	contBack3Initial=contBack3Surviving;
	contBack4Initial=contBack4Surviving;
	contBack5Initial=contBack5Surviving;
	contBack6Initial=contBack6Surviving;
	contBack7Initial=contBack7Surviving;
	
	contHHSurviving=0;
	contBack1Surviving=0;
	contBack2Surviving=0;
	contBack3Surviving=0;
	contBack4Surviving=0;
	contBack5Surviving=0;
	contBack6Surviving=0;
	contBack7Surviving=0;

	for(int i=0; i<NN7OutputCopy.size(); i++)
	{
		if(NN7OutputCopy[i]<defCutNN7 && NN7OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
			NN5OutputCopy[i]=-999;
			NN6OutputCopy[i]=-999;
			NN7OutputCopy[i]=-999;
		}
	}
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999)
		{
			if(i<sizeHH) contHHSurviving++;
			else if(i<sizeHH+sizeBack1) contBack1Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2) contBack2Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6Surviving++;
			else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7Surviving++;
		}
	}

	cout<<"Initial HH events (before NNqqHX): "<<contHHInitial*weightHH<<endl<<"Surviving HH events (after NNqqHX): "<<contHHSurviving*weightHH<<endl<<"Initial qq events (before NNqqHX): "<<contBack1Initial*weightBack1<<endl<<"Surviving qq events (after NNqqHX): "<<contBack1Surviving*weightBack1<<endl<<"Initial ttbar events (before NNqqHX): "<<contBack2Initial*weightBack2<<endl<<"Surviving ttbar events (after NNqqHX): "<<contBack2Surviving*weightBack2<<endl<<"Initial ZZ events (before NNqqHX): "<<contBack3Initial*weightBack3<<endl<<"Surviving ZZ events (after NNqqHX): "<<contBack3Surviving*weightBack3<<endl<<"Initial WW events (before NNqqHX): "<<contBack4Initial*weightBack4<<endl<<"Surviving WW events (after NNqqHX): "<<contBack4Surviving*weightBack4<<endl<<"Initial qqX events (before NNqqHX): "<<contBack5Initial*weightBack5<<endl<<"Surviving qqX events (after NNqqHX): "<<contBack5Surviving*weightBack5<<endl<<"Initial qqqqX events (before NNqqHX): "<<contBack6Initial*weightBack6<<endl<<"Surviving qqqqX events (after NNqqHX): "<<contBack6Surviving*weightBack6<<endl<<"Initial qqHX events (before NNqqHX): "<<contBack7Initial*weightBack7<<endl<<"Surviving qqHX events (after NNqqHX): "<<contBack7Surviving*weightBack7<<endl<<endl;

	contHHInitial=contHHSurviving;
	contBack1Initial=contBack1Surviving;
	contBack2Initial=contBack2Surviving;
	contBack3Initial=contBack3Surviving;
	contBack4Initial=contBack4Surviving;
	contBack5Initial=contBack5Surviving;
	contBack6Initial=contBack6Surviving;
	contBack7Initial=contBack7Surviving;
	
	contHHSurviving=0;
	contBack1Surviving=0;
	contBack2Surviving=0;
	contBack3Surviving=0;
	contBack4Surviving=0;
	contBack5Surviving=0;
	contBack6Surviving=0;
	contBack7Surviving=0;

	int contHH=0, contBack1=0, contBack2=0, contBack3=0, contBack4=0, contBack5=0, contBack6=0, contBack7=0;
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999 && i<sizeHH) contHH++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1) contBack1++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2) contBack2++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5) contBack5++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6) contBack6++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4+sizeBack5+sizeBack6+sizeBack7) contBack7++;
	}
	
	HHRemainingCombined = contHH*weightHH;
	backRemainingCombined = contBack1*weightBack1+contBack2*weightBack2+contBack3*weightBack3+contBack4*weightBack4+contBack5*weightBack5+contBack6*weightBack6+contBack7*weightBack7;
	int HHRemainingCombinedUnweighted = contHH;
	int backRemainingCombinedUnweighted = contBack1+contBack2+contBack3+contBack4+contBack5+contBack6+contBack7;
	
	maxSignificanceCombined = HHRemainingCombined/sqrt(HHRemainingCombined+backRemainingCombined);
	
	cout<<endl<<endl<<"From check: "<<endl<<"remainingHH unweighted: "<<HHRemainingCombinedUnweighted<<endl<<"remainingBackUnweighted: "<<backRemainingCombinedUnweighted<<endl<<endl;
	cout<<"remainingHH: "<<HHRemainingCombined<<endl<<"remainingBack: "<<backRemainingCombined<<endl<<endl<<endl;
}

 
 void splitHistogram(TH1F& hist, double bottomHistLimit, double topHistLimit, double cut, string title)
 {
 	TH1F *histRemaining = new TH1F("histRemaining", title.c_str(), hist.GetNbinsX(), bottomHistLimit, topHistLimit);
 	histRemaining->SetFillColor(hist.GetFillColor());
 	
 	
 	int binLow = hist.FindBin(cut);
 	int binHigh = hist.FindBin(topHistLimit);
 	for(int i=binLow; i<=binHigh; i++)
 	{
 	 	double binContent = hist.GetBinContent(i);
        	histRemaining->SetBinContent(i, binContent);
 	} 
 	
 	TH1F *histRemainingZoommed = (TH1F*)histRemaining->Clone("histRemainingZoommed");
	histRemainingZoommed->GetYaxis()->SetRangeUser(0, histRemaining->GetMaximum()*1.05);
 	
 	TLine* linex = new TLine(cut,0.0,cut,(hist.GetMaximum())*1.05);
	linex->SetLineColor(kOrange);
	linex->SetLineWidth(2);
	linex->SetLineStyle(7);
 	
 	TCanvas *csplitHistogram1 = new TCanvas();
 	TCanvas *csplitHistogram2 = new TCanvas();
 	TCanvas *csplitHistogram3 = new TCanvas();
 	
 	csplitHistogram1->cd();
 	hist.Draw("HIST");
	linex->Draw("same");
	
	histRemaining->GetYaxis()->SetRangeUser(0, hist.GetMaximum()*1.05);
	csplitHistogram2->cd();
 	histRemaining->Draw("HIST");
	linex->Draw("same");
	
	csplitHistogram3->cd();
 	histRemainingZoommed->Draw("HIST");
 	TLine* linex2 = new TLine(cut,0.0,cut,(histRemainingZoommed->GetMaximum())*1.05);
	linex2->SetLineColor(kOrange);
	linex2->SetLineWidth(2);
	linex2->SetLineStyle(7);
	linex2->Draw("same");
	binLow = hist.FindBin(bottomHistLimit);
	int entries = histRemainingZoommed->Integral(binLow, binHigh);
	string entriesText = "Entries: " + to_string(entries);
	TPaveText *pave = new TPaveText(0.3, 100, 0.5, 150, "");	
        pave->AddText(entriesText.c_str());
        pave->SetBorderSize(0);  // optionally set border size
        pave->SetFillColor(0);   // set fill color (0 is transparent)
        pave->SetTextFont(1);   // optional: set the font
        pave->SetTextSize(0.03); // optional: set the text size
        pave->Draw("same");
	
	
 }
 
void pairMassCut(int topology, double targetFractionHH, double& effi)
{
    	TFile* input = new TFile("outputTreesHHbbbbESpreadDurham.root");
    
 	std::cout << "--- pairMassCut   : Using input file: " << input->GetName() << std::endl;
 
    	TTree* theTree;
    
    	if(topology == 1) theTree = (TTree*)input->Get("TreeS");
    	else if(topology == 2) theTree = (TTree*)input->Get("TreeBqq");
    	else if(topology == 3) theTree = (TTree*)input->Get("TreeBtt");
    
    	Float_t aplanarity, cosThetaB1, cosThetaB2, cosThetaB3, cosThetaB4, invMassB1, invMassB2, jetB1M, jetB2M, jetB3M, jetB4M, jetB1Pt, jetB2Pt, jetB3Pt, jetB4Pt, minJetM, sphericity, sumPt;
    	theTree->SetBranchAddress( "aplanarity", &aplanarity );
    	theTree->SetBranchAddress( "cosThetaB1", &cosThetaB1 );
    	theTree->SetBranchAddress( "cosThetaB2", &cosThetaB2 );
    	theTree->SetBranchAddress( "cosThetaB3", &cosThetaB3 );
    	theTree->SetBranchAddress( "cosThetaB4", &cosThetaB4 );
    	theTree->SetBranchAddress( "invMassB1", &invMassB1 );
    	theTree->SetBranchAddress( "invMassB2", &invMassB2 );
    	theTree->SetBranchAddress( "jetB1M", &jetB1M );
    	theTree->SetBranchAddress( "jetB2M", &jetB2M );
    	theTree->SetBranchAddress( "jetB3M", &jetB3M );
    	theTree->SetBranchAddress( "jetB4M", &jetB4M );
    	theTree->SetBranchAddress( "jetB1Pt", &jetB1Pt );
    	theTree->SetBranchAddress( "jetB2Pt", &jetB2Pt );
    	theTree->SetBranchAddress( "jetB3Pt", &jetB3Pt );
    	theTree->SetBranchAddress( "jetB4Pt", &jetB4Pt );
    	theTree->SetBranchAddress( "minJetM", &minJetM );
    	theTree->SetBranchAddress( "sphericity", &sphericity );
    	theTree->SetBranchAddress( "sumPt", &sumPt );
    	
    	double nEntries = theTree->GetEntries();
    	double invMassB1Remaining=0, invMassB2Remaining=0, eventsRemaining=0;
    	for (Long64_t entry = 0; entry < nEntries; entry++) 
    	{
		theTree->GetEntry(entry);
		//std::cout << "Entry " << entry << ": invMassB = " << invMassB << std::endl;
		
		if(50<invMassB1 && invMassB1<180 && 50<invMassB2 && invMassB2<180)
		{ 
			invMassB1Remaining++;
			invMassB2Remaining++;
			eventsRemaining++;
		}
		
		/*if(73<invMassB2 && invMassB2<140)
		{ 
			invMassB2Remaining++;
		}*/
    	}
    	
    	cout<<"invMassB1Remaining: "<<eventsRemaining<<"      fraction: "<<eventsRemaining/nEntries<<endl<<"invMassB2Remaining: "<<invMassB2Remaining<<"      fraction: "<<invMassB2Remaining/nEntries<<endl;
    	
    	effi = eventsRemaining/nEntries;
    	
    	
    	return;
}

void addPointTo2DHistogram(TH2F& hist, double x, double y, int color, double radiusX, double radiusY)
{
	TEllipse* point = new TEllipse(x,y,radiusX,radiusY); 
        point->SetLineStyle(1);
        point->SetFillColor(color);
        point->SetLineColor(color);
   	point->SetLineWidth(1);
   	
   	TCanvas *caddPointTo2DHistogram1 = new TCanvas();
   	
   	//gStyle->SetOptStat(0);
   	caddPointTo2DHistogram1->cd();
	caddPointTo2DHistogram1->SetLogy();
   	hist.Draw("BOX");
	point->Draw("same");
}

void addPointToCanvas(double x, double y, int color, double radiusX, double radiusY)
{
	TEllipse* point = new TEllipse(x,y,radiusX,radiusY); 
        point->SetLineStyle(1);
        point->SetFillColor(color);
        point->SetLineColor(color);
   	point->SetLineWidth(1);
	point->Draw("same");
}

void add2DHistToFile(TH2F& hist, string fileName, string histName)
{
	int color = hist.GetLineColor();
	TH2F *histClone = (TH2F*)hist.Clone(histName.c_str());
	histClone->SetName(histName.c_str());
	histClone->SetLineColor(color);
	TFile *file = TFile::Open(fileName.c_str(), "UPDATE");
	histClone->Write();
	file->Close();
	return;
}

void draw2DHistsFromFile(string fileName, string canvasScale, bool point, double x, double y, int pointColor, double radiusX, double radiusY)
{
	TFile *file = TFile::Open(fileName.c_str(), "READ");
	//gStyle->SetOptStat(0);
	TCanvas *cdraw2DHistsFromFile1 = new TCanvas();
	cdraw2DHistsFromFile1->cd();
	if(canvasScale == "log") cdraw2DHistsFromFile1->SetLogy();
	TIter next(file->GetListOfKeys());
    	TKey *key;
    	bool isFirstHistogram = true;
    	while ((key = (TKey*)next()))
    	{
    		if (strstr(key->GetClassName(), "TH2F"))
    		{
    			TH2F *hist = (TH2F*)key->ReadObj();
    			string str(hist->GetName());
    			if(str == "ROCRejBDT") hist->SetFillColor(kRed);
    			else if(str == "ROCRejBDTInvMass") hist->SetFillColor(kRed-9);
    			else if(str == "ROCRejBDTG") hist->SetFillColor(kBlue);
    			else if(str == "ROCRejBDTGInvMass") hist->SetFillColor(kAzure+5);
    			else if(str == "ROCRejMLP") hist->SetFillColor(kViolet);
    			else if(str == "ROCRejMLPInvMass") hist->SetFillColor(kMagenta-4);
    			/*hist->SetFillColor(kRed);
    			hist->SetLineColor(kRed);*/
    			hist->Draw("BOX same");
    			/*hist->Draw(isFirstHistogram ? "" : "BOX");
    			isFirstHistogram = false;
    			int color = hist->GetFillColor();
    			cout<<endl<<endl<<endl<<"Hist name: "<<hist->GetName()<<"       Color: "<<color<<endl<<endl<<endl;*/
    		}
    	}
    	if(point == true) addPointToCanvas(x, y, pointColor, radiusX, radiusY);
}


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

	errorBottom = errorLeft;
	errorTop = errorRight;
	////////
	
}

//////Function that generates the files to run genetic algorithm
void generateFilesGA(vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, vector<double>& BDTqqXOutput, vector<double>& BDTqqqqXOutput, vector<double>& BDTqqHXOutput, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, string rtdCut, string preselection, string sampleName)
{
	string outputTreeSGAText = "analysis/outputTreeSGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
	TFile *outputTreeSGA = new TFile(outputTreeSGAText.c_str(), "recreate");	
	TTree TreeSGA("TreeSGA","a simple Tree with simple variables (GA)");
	string outputTreeBqqGAText = "analysis/outputTreeBqqGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
	TFile *outputTreeBqqGA = new TFile(outputTreeBqqGAText.c_str(), "recreate");
	TTree TreeBqqGA("TreeBqqGA","a bqqimple Tree with bqqimple variables (GA)");
	string outputTreeBttbarGAText = "analysis/outputTreeBttbarGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
	TFile *outputTreeBttbarGA = new TFile(outputTreeBttbarGAText.c_str(), "recreate");
	TTree TreeBttbarGA("TreeBttbarGA","a bttbarimple Tree with bttbarimple variables (GA)");
	string outputTreeBZZGAText = "analysis/outputTreeBZZGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
	TFile *outputTreeBZZGA = new TFile(outputTreeBZZGAText.c_str(), "recreate");
	TTree TreeBZZGA("TreeBZZGA","a bZZimple Tree with bZZimple variables (GA)");
	string outputTreeBWWGAText = "analysis/outputTreeBWWGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
	TFile *outputTreeBWWGA = new TFile(outputTreeBWWGAText.c_str(), "recreate");
	TTree TreeBWWGA("TreeBWWGA","a bWWimple Tree with bWWimple variables (GA)");
	string outputTreeBqqXGAText = "analysis/outputTreeBqqXGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
	TFile *outputTreeBqqXGA = new TFile(outputTreeBqqXGAText.c_str(), "recreate");
	TTree TreeBqqXGA("TreeBqqXGA","a bqqXimple Tree with bqqXimple variables (GA)");
	string outputTreeBqqqqXGAText = "analysis/outputTreeBqqqqXGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
	TFile *outputTreeBqqqqXGA = new TFile(outputTreeBqqqqXGAText.c_str(), "recreate");
	TTree TreeBqqqqXGA("TreeBqqqqXGA","a bqqqqXimple Tree with bqqqqXimple variables (GA)");
	string outputTreeBqqHXGAText = "analysis/outputTreeBqqHXGAESpreadDurham"+rtdCut+preselection+sampleName+".root";
	TFile *outputTreeBqqHXGA = new TFile(outputTreeBqqHXGAText.c_str(), "recreate");
	TTree TreeBqqHXGA("TreeBqqHXGA","a bqqHXimple Tree with bqqHXimple variables (GA)");

	float entryIndex, NN1Output, NN2Output, NN3Output, NN4Output, NN5Output, NN6Output, NN7Output;
  	TreeSGA.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeSGA.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeSGA.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeSGA.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeSGA.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	TreeSGA.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeSGA.Branch("NN6Output",&NN6Output,"NN6Output/F");
 	TreeSGA.Branch("NN7Output",&NN7Output,"NN7Output/F");	
  	TreeBqqGA.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBqqGA.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBqqGA.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBqqGA.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBqqGA.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqGA.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqGA.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqGA.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBttbarGA.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBttbarGA.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBttbarGA.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBttbarGA.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBttbarGA.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBttbarGA.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBttbarGA.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBttbarGA.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBZZGA.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBZZGA.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBZZGA.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBZZGA.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBZZGA.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBZZGA.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBZZGA.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBZZGA.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBWWGA.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBWWGA.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBWWGA.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBWWGA.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBWWGA.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBWWGA.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBWWGA.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBWWGA.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqXGA.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqXGA.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqXGA.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqXGA.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqXGA.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqXGA.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqXGA.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqXGA.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqqqXGA.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqqqXGA.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqqqXGA.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqqqXGA.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqqqXGA.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqqqXGA.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqqqXGA.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqqqXGA.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqHXGA.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqHXGA.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqHXGA.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqHXGA.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqHXGA.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqHXGA.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqHXGA.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqHXGA.Branch("NN7Output",&NN7Output,"NN7Output/F");

	for(int i=0; i<BDTqqOutput.size(); i++)
  	{
  		entryIndex = i;
  		NN1Output = BDTqqOutput[i];
  		NN2Output = BDTttbarOutput[i];
  		NN3Output = BDTZZOutput[i];
  		NN4Output = BDTWWOutput[i];
		NN5Output = BDTqqXOutput[i];
		NN6Output = BDTqqqqXOutput[i];
		NN7Output = BDTqqHXOutput[i];

		if(i<sizeHH) TreeSGA.Fill();
		else if(i<sizeHH+sizeqq) TreeBqqGA.Fill();
		else if(i<sizeHH+sizeqq+sizettbar) TreeBttbarGA.Fill();
		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ) TreeBZZGA.Fill();
		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW) TreeBWWGA.Fill();
		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX) TreeBqqXGA.Fill();
		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX) TreeBqqqqXGA.Fill();
		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX+sizeqqHX) TreeBqqHXGA.Fill();
		else throw std::runtime_error("ERROR: event out of bounds in generate!");	
	}	

	outputTreeSGA->cd();
   	TreeSGA.Write();
   	outputTreeBqqGA->cd();
   	TreeBqqGA.Write();
   	outputTreeBttbarGA->cd();
   	TreeBttbarGA.Write();
   	outputTreeBZZGA->cd();
   	TreeBZZGA.Write();
   	outputTreeBWWGA->cd();
   	TreeBWWGA.Write();
	outputTreeBqqXGA->cd();
	TreeBqqXGA.Write();
	outputTreeBqqqqXGA->cd();
	TreeBqqqqXGA.Write();
	outputTreeBqqHXGA->cd();
	TreeBqqHXGA.Write();
}

/////Function that generates the files to train NN to find optimal cut for NNs outputs
void generateFilesOutputNN(vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, vector<double>& BDTqqXOutput, vector<double>& BDTqqqqXOutput, vector<double>& BDTqqHXOutput, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, string rtdCut, string preselection, string sampleName)
{
	string outputTreeSNNTrainText = "analysis/outputTreeSNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	TFile *outputTreeSNNTrain = new TFile(outputTreeSNNTrainText.c_str(), "recreate");
 	TTree TreeSNNTrain("TreeSNNTrain","a simple Tree with simple variables (Train)");
 	string outputTreeSNNTestText = "analysis/outputTreeSNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeSNNTest = new TFile(outputTreeSNNTestText.c_str(), "recreate");
  	TTree TreeSNNTest("TreeSNNTest","a simple Tree with simple variables (Test)");
  	string outputTreeBqqNNTrainText = "analysis/outputTreeBqqNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	TFile *outputTreeBqqNNTrain = new TFile(outputTreeBqqNNTrainText.c_str(), "recreate");
 	TTree TreeBqqNNTrain("TreeBqqNNTrain","a bqqimple Tree with simple variables (Train)");
 	string outputTreeBqqNNTestText = "analysis/outputTreeBqqNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeBqqNNTest = new TFile(outputTreeBqqNNTestText.c_str(), "recreate");
  	TTree TreeBqqNNTest("TreeBqqNNTest","a bqqimple Tree with simple variables (Test)");
  	string outputTreeBttbarNNTrainText = "analysis/outputTreeBttbarNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	TFile *outputTreeBttbarNNTrain = new TFile(outputTreeBttbarNNTrainText.c_str(), "recreate");
 	TTree TreeBttbarNNTrain("TreeBttbarNNTrain","a bttbarimple Tree with simple variables (Train)");
 	string outputTreeBttbarNNTestText = "analysis/outputTreeBttbarNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeBttbarNNTest = new TFile(outputTreeBttbarNNTestText.c_str(), "recreate");
  	TTree TreeBttbarNNTest("TreeBttbarNNTest","a bttbarimple Tree with simple variables (Test)");
  	string outputTreeBZZNNTrainText = "analysis/outputTreeBZZNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	TFile *outputTreeBZZNNTrain = new TFile(outputTreeBZZNNTrainText.c_str(), "recreate");
 	TTree TreeBZZNNTrain("TreeBZZNNTrain","a bZZimple Tree with simple variables (Train)");
 	string outputTreeBZZNNTestText = "analysis/outputTreeBZZNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeBZZNNTest = new TFile(outputTreeBZZNNTestText.c_str(), "recreate");
  	TTree TreeBZZNNTest("TreeBZZNNTest","a bZZimple Tree with simple variables (Test)");
  	string outputTreeBWWNNTrainText = "analysis/outputTreeBWWNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	TFile *outputTreeBWWNNTrain = new TFile(outputTreeBWWNNTrainText.c_str(), "recreate");
 	TTree TreeBWWNNTrain("TreeBWWNNTrain","a bWWimple Tree with simple variables (Train)");
 	string outputTreeBWWNNTestText = "analysis/outputTreeBWWNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeBWWNNTest = new TFile(outputTreeBWWNNTestText.c_str(), "recreate");
  	TTree TreeBWWNNTest("TreeBWWNNTest","a bWWimple Tree with simple variables (Test)");
	string outputTreeBqqXNNTrainText = "analysis/outputTreeBqqXNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
   	TFile *outputTreeBqqXNNTrain = new TFile(outputTreeBqqXNNTrainText.c_str(), "recreate");
   	TTree TreeBqqXNNTrain("TreeBqqXNNTrain","a bqqXimple Tree with simple variables (Train)");
   	string outputTreeBqqXNNTestText = "analysis/outputTreeBqqXNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
   	TFile *outputTreeBqqXNNTest = new TFile(outputTreeBqqXNNTestText.c_str(), "recreate");
   	TTree TreeBqqXNNTest("TreeBqqXNNTest","a bqqXimple Tree with simple variables (Test)");
	string outputTreeBqqqqXNNTrainText = "analysis/outputTreeBqqqqXNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	TFile *outputTreeBqqqqXNNTrain = new TFile(outputTreeBqqqqXNNTrainText.c_str(), "recreate");
	TTree TreeBqqqqXNNTrain("TreeBqqqqXNNTrain","a bqqqqXimple Tree with simple variables (Train)");
	string outputTreeBqqqqXNNTestText = "analysis/outputTreeBqqqqXNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile *outputTreeBqqqqXNNTest = new TFile(outputTreeBqqqqXNNTestText.c_str(), "recreate");
	TTree TreeBqqqqXNNTest("TreeBqqqqXNNTest","a bqqqqXimple Tree with simple variables (Test)");
	string outputTreeBqqHXNNTrainText = "analysis/outputTreeBqqHXNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	TFile *outputTreeBqqHXNNTrain = new TFile(outputTreeBqqHXNNTrainText.c_str(), "recreate");
	TTree TreeBqqHXNNTrain("TreeBqqHXNNTrain","a bqqHXimple Tree with simple variables (Train)");
	string outputTreeBqqHXNNTestText = "analysis/outputTreeBqqHXNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile *outputTreeBqqHXNNTest = new TFile(outputTreeBqqHXNNTestText.c_str(), "recreate");
	TTree TreeBqqHXNNTest("TreeBqqHXNNTest","a bqqHXimple Tree with simple variables (Test)");


  	float entryIndex, NN1Output, NN2Output, NN3Output, NN4Output, NN5Output, NN6Output, NN7Output;
  	TreeSNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeSNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeSNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeSNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeSNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	TreeSNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeSNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
 	TreeSNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");	
  	TreeBqqNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBqqNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBqqNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBqqNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBqqNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBttbarNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBttbarNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBttbarNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBttbarNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBttbarNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBttbarNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBttbarNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBttbarNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBZZNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBZZNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBZZNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBZZNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBZZNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBZZNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBZZNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBZZNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBWWNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBWWNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBWWNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBWWNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBWWNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBWWNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBWWNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBWWNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqXNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqXNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqXNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqXNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqXNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqXNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqXNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqXNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqqqXNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqqqXNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqqqXNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqqqXNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqqqXNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqqqXNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqqqXNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqqqXNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqHXNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqHXNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqHXNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqHXNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqHXNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqHXNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqHXNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqHXNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");	
  	
  	TreeSNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
   TreeSNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
   TreeSNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
   TreeSNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
   TreeSNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
   TreeSNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
   TreeSNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
   TreeSNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F"); 
   TreeBqqNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
   TreeBqqNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
   TreeBqqNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
   TreeBqqNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
   TreeBqqNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
   TreeBqqNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
   TreeBqqNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
   TreeBqqNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
   TreeBttbarNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
   TreeBttbarNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
   TreeBttbarNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
   TreeBttbarNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
   TreeBttbarNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
   TreeBttbarNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
   TreeBttbarNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
   TreeBttbarNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
   TreeBZZNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
   TreeBZZNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
   TreeBZZNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
   TreeBZZNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
   TreeBZZNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
   TreeBZZNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
   TreeBZZNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
   TreeBZZNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
   TreeBWWNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
   TreeBWWNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
   TreeBWWNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
   TreeBWWNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
   TreeBWWNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
   TreeBWWNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
   TreeBWWNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
   TreeBWWNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
   TreeBqqXNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
   TreeBqqXNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
   TreeBqqXNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
   TreeBqqXNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
   TreeBqqXNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
   TreeBqqXNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
   TreeBqqXNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
   TreeBqqXNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
   TreeBqqqqXNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
   TreeBqqqqXNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
   TreeBqqqqXNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
   TreeBqqqqXNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
   TreeBqqqqXNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
   TreeBqqqqXNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
   TreeBqqqqXNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
   TreeBqqqqXNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
   TreeBqqHXNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
   TreeBqqHXNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
   TreeBqqHXNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
   TreeBqqHXNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
   TreeBqqHXNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
   TreeBqqHXNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
   TreeBqqHXNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
   TreeBqqHXNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");

   int contEventsSurvived=0;
   int contEventsEliminated=0;
   //cout<<"size of NNOutputs: "<<BDTqqOutput.size()<<endl;
  	
  	for(int i=0; i<BDTqqOutput.size(); i++)
  	{
  		entryIndex = i;
  		NN1Output = BDTqqOutput[i];
  		NN2Output = BDTttbarOutput[i];
  		NN3Output = BDTZZOutput[i];
  		NN4Output = BDTWWOutput[i];
		NN5Output = BDTqqXOutput[i];
		NN6Output = BDTqqqqXOutput[i];
		NN7Output = BDTqqHXOutput[i];

		if(NN1Output != -999) contEventsSurvived++;
		else contEventsEliminated++;

		if(NN1Output == -999) continue;

  		random_device rd;  // Will be used to obtain a seed for the random number engine
		mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
		uniform_int_distribution<> distrib(1, 4); // Uniform distribution between 1 and 4
		if(distrib(gen) != 1)
		{ 
			if(i<sizeHH) TreeSNNTrain.Fill();
			else if(i<sizeHH+sizeqq) TreeBqqNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar) TreeBttbarNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ) TreeBZZNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW) TreeBWWNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX) TreeBqqXNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX) TreeBqqqqXNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX+sizeqqHX) TreeBqqHXNNTrain.Fill();
			else throw std::runtime_error("ERROR: event out of bounds in generate!");	
		}
		else
		{ 
			if(i<sizeHH) TreeSNNTest.Fill();
			else if(i<sizeHH+sizeqq) TreeBqqNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar) TreeBttbarNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ) TreeBZZNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW) TreeBWWNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX) TreeBqqXNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX) TreeBqqqqXNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX+sizeqqHX) TreeBqqHXNNTest.Fill();
			else throw std::runtime_error("ERROR: event out of bounds in generate!");	
		}
						
  	}
	//cout<<"contEventsSurvived: "<<contEventsSurvived<<endl;
	//cout<<"contEventsEliminated: "<<contEventsEliminated<<endl;
  	
  	outputTreeSNNTrain->cd();
   	TreeSNNTrain.Write();
   	outputTreeBqqNNTrain->cd();
   	TreeBqqNNTrain.Write();
   	outputTreeBttbarNNTrain->cd();
   	TreeBttbarNNTrain.Write();
   	outputTreeBZZNNTrain->cd();
   	TreeBZZNNTrain.Write();
   	outputTreeBWWNNTrain->cd();
   	TreeBWWNNTrain.Write();
	outputTreeBqqXNNTrain->cd();
	TreeBqqXNNTrain.Write();
	outputTreeBqqqqXNNTrain->cd();
	TreeBqqqqXNNTrain.Write();
	outputTreeBqqHXNNTrain->cd();
	TreeBqqHXNNTrain.Write();	
   	
   	outputTreeSNNTest->cd();
   	TreeSNNTest.Write();
   	outputTreeBqqNNTest->cd();
   	TreeBqqNNTest.Write();
   	outputTreeBttbarNNTest->cd();
   	TreeBttbarNNTest.Write();
   	outputTreeBZZNNTest->cd();
   	TreeBZZNNTest.Write();
   	outputTreeBWWNNTest->cd();
   	TreeBWWNNTest.Write();
	outputTreeBqqXNNTest->cd();
	TreeBqqXNNTest.Write();
	outputTreeBqqqqXNNTest->cd();
	TreeBqqqqXNNTest.Write();
	outputTreeBqqHXNNTest->cd();
	TreeBqqHXNNTest.Write();   	
}

/////Function that REgenerates the files to train NN to find optimal cut for NNs outputs to have the same sampling as some previously given files.
void reGenerateFilesOutputNN(vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, vector<double>& BDTqqXOutput, vector<double>& BDTqqqqXOutput, vector<double>& BDTqqHXOutput,  int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, string rtdCut, string preselection, string sampleName)
{
	TFile *fileTrainHH, *fileTestHH, *fileTrainqq, *fileTestqq, *fileTrainttbar, *fileTestttbar, *fileTrainZZ, *fileTestZZ, *fileTrainWW, *fileTestWW, *fileTrainqqX, *fileTestqqX, *fileTrainqqqqX, *fileTestqqqqX, *fileTrainqqHX, *fileTestqqHX;
	TTree *TreeTrainHH, *TreeTestHH, *TreeTrainqq, *TreeTestqq, *TreeTrainttbar, *TreeTestttbar, *TreeTrainZZ, *TreeTestZZ, *TreeTrainWW, *TreeTestWW, *TreeTrainqqX, *TreeTestqqX, *TreeTrainqqqqX, *TreeTestqqqqX, *TreeTrainqqHX, *TreeTestqqHX;
	fileTrainHH = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeSNNTrain.root");
	fileTestHH = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeSNNTest.root");
	fileTrainHH->GetObject("TreeSNNTrain", TreeTrainHH);
	fileTestHH->GetObject("TreeSNNTest", TreeTestHH);
	fileTrainqq = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBqqNNTrain.root");
	fileTestqq = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBqqNNTest.root");
	fileTrainqq->GetObject("TreeBqqNNTrain", TreeTrainqq);
	fileTestqq->GetObject("TreeBqqNNTest", TreeTestqq);
	fileTrainttbar = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBttbarNNTrain.root");
	fileTestttbar = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBttbarNNTest.root");
	fileTrainttbar->GetObject("TreeBttbarNNTrain", TreeTrainttbar);
	fileTestttbar->GetObject("TreeBttbarNNTest", TreeTestttbar);
	fileTrainZZ = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBZZNNTrain.root");
	fileTestZZ = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBZZNNTest.root");
	fileTrainZZ->GetObject("TreeBZZNNTrain", TreeTrainZZ);
	fileTestZZ->GetObject("TreeBZZNNTest", TreeTestZZ);
	fileTrainWW = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBWWNNTrain.root");
	fileTestWW = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBWWNNTest.root");
	fileTrainWW->GetObject("TreeBWWNNTrain", TreeTrainWW);
	fileTestWW->GetObject("TreeBWWNNTest", TreeTestWW);
	fileTrainqqX = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBqqXNNTrain.root");
	fileTestqqX = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBqqXNNTest.root");
	fileTrainqqX->GetObject("TreeBqqXNNTrain", TreeTrainqqX);
	fileTestqqX->GetObject("TreeBqqXNNTest", TreeTestqqX);
	fileTrainqqqqX = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBqqqqXNNTrain.root");
	fileTestqqqqX = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBqqqqXNNTest.root");
	fileTrainqqqqX->GetObject("TreeBqqqqXNNTrain", TreeTrainqqqqX);
	fileTestqqqqX->GetObject("TreeBqqqqXNNTest", TreeTestqqqqX);
	fileTrainqqHX = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBqqHXNNTrain.root");
	fileTestqqHX = TFile::Open("analysis/SampleOGNNOfNNs/outputTreeBqqHXNNTest.root");
	fileTrainqqHX->GetObject("TreeBqqHXNNTrain", TreeTrainqqHX);
	fileTestqqHX->GetObject("TreeBqqHXNNTest", TreeTestqqHX);
	
	float NN1OutputTrainHH, NN2OutputTrainHH, NN3OutputTrainHH, NN4OutputTrainHH, NN5OutputTrainHH, NN6OutputTrainHH, NN7OutputTrainHH, NN1OutputTestHH, NN2OutputTestHH, NN3OutputTestHH, NN4OutputTestHH, NN5OutputTestHH, NN6OutputTestHH, NN7OutputTestHH, NN1OutputTrainqq, NN2OutputTrainqq, NN3OutputTrainqq, NN4OutputTrainqq, NN5OutputTrainqq, NN6OutputTrainqq, NN7OutputTrainqq, NN1OutputTestqq, NN2OutputTestqq, NN3OutputTestqq, NN4OutputTestqq, NN5OutputTestqq, NN6OutputTestqq, NN7OutputTestqq, NN1OutputTrainttbar, NN2OutputTrainttbar, NN3OutputTrainttbar, NN4OutputTrainttbar, NN5OutputTrainttbar, NN6OutputTrainttbar, NN7OutputTrainttbar, NN1OutputTestttbar, NN2OutputTestttbar, NN3OutputTestttbar, NN4OutputTestttbar, NN5OutputTestttbar, NN6OutputTestttbar, NN7OutputTestttbar, NN1OutputTrainZZ, NN2OutputTrainZZ, NN3OutputTrainZZ, NN4OutputTrainZZ, NN5OutputTrainZZ, NN6OutputTrainZZ, NN7OutputTrainZZ, NN1OutputTestZZ, NN2OutputTestZZ, NN3OutputTestZZ, NN4OutputTestZZ, NN5OutputTestZZ, NN6OutputTestZZ, NN7OutputTestZZ, NN1OutputTrainWW, NN2OutputTrainWW, NN3OutputTrainWW, NN4OutputTrainWW, NN5OutputTrainWW, NN6OutputTrainWW, NN7OutputTrainWW, NN1OutputTestWW, NN2OutputTestWW, NN3OutputTestWW, NN4OutputTestWW, NN5OutputTestWW, NN6OutputTestWW, NN7OutputTestWW, NN1OutputTrainqqX, NN2OutputTrainqqX, NN3OutputTrainqqX, NN4OutputTrainqqX, NN5OutputTrainqqX, NN6OutputTrainqqX, NN7OutputTrainqqX, NN1OutputTestqqX, NN2OutputTestqqX, NN3OutputTestqqX, NN4OutputTestqqX, NN5OutputTestqqX, NN6OutputTestqqX, NN7OutputTestqqX, NN1OutputTrainqqqqX, NN2OutputTrainqqqqX, NN3OutputTrainqqqqX, NN4OutputTrainqqqqX, NN5OutputTrainqqqqX, NN6OutputTrainqqqqX, NN7OutputTrainqqqqX, NN1OutputTestqqqqX, NN2OutputTestqqqqX, NN3OutputTestqqqqX, NN4OutputTestqqqqX, NN5OutputTestqqqqX, NN6OutputTestqqqqX, NN7OutputTestqqqqX, NN1OutputTrainqqHX, NN2OutputTrainqqHX, NN3OutputTrainqqHX, NN4OutputTrainqqHX, NN5OutputTrainqqHX, NN6OutputTrainqqHX, NN7OutputTrainqqHX, NN1OutputTestqqHX, NN2OutputTestqqHX, NN3OutputTestqqHX, NN4OutputTestqqHX, NN5OutputTestqqHX, NN6OutputTestqqHX, NN7OutputTestqqHX, entryIndexTrain, entryIndexTest;

	TreeTrainHH->SetBranchAddress("entryIndex", &entryIndexTrain);
	TreeTestHH->SetBranchAddress("entryIndex", &entryIndexTest);
	TreeTrainHH->SetBranchAddress("NN1Output", &NN1OutputTrainHH);
	TreeTestHH->SetBranchAddress("NN1Output", &NN1OutputTestHH);
	TreeTrainHH->SetBranchAddress("NN2Output", &NN2OutputTrainHH);
	TreeTestHH->SetBranchAddress("NN2Output", &NN2OutputTestHH);
	TreeTrainHH->SetBranchAddress("NN3Output", &NN3OutputTrainHH);
	TreeTestHH->SetBranchAddress("NN3Output", &NN3OutputTestHH);
	TreeTrainHH->SetBranchAddress("NN4Output", &NN4OutputTrainHH);
	TreeTestHH->SetBranchAddress("NN4Output", &NN4OutputTestHH);
	TreeTrainHH->SetBranchAddress("NN5Output", &NN5OutputTrainHH);
	TreeTestHH->SetBranchAddress("NN5Output", &NN5OutputTestHH);
	TreeTrainHH->SetBranchAddress("NN6Output", &NN6OutputTrainHH);
	TreeTestHH->SetBranchAddress("NN6Output", &NN6OutputTestHH);
	TreeTrainHH->SetBranchAddress("NN7Output", &NN7OutputTrainHH);
	TreeTestHH->SetBranchAddress("NN7Output", &NN7OutputTestHH);
	TreeTrainqq->SetBranchAddress("entryIndex", &entryIndexTrain);
	TreeTestqq->SetBranchAddress("entryIndex", &entryIndexTest);
	TreeTrainqq->SetBranchAddress("NN1Output", &NN1OutputTrainqq);
	TreeTestqq->SetBranchAddress("NN1Output", &NN1OutputTestqq);
	TreeTrainqq->SetBranchAddress("NN2Output", &NN2OutputTrainqq);
	TreeTestqq->SetBranchAddress("NN2Output", &NN2OutputTestqq);
	TreeTrainqq->SetBranchAddress("NN3Output", &NN3OutputTrainqq);
	TreeTestqq->SetBranchAddress("NN3Output", &NN3OutputTestqq);
	TreeTrainqq->SetBranchAddress("NN4Output", &NN4OutputTrainqq);
	TreeTestqq->SetBranchAddress("NN4Output", &NN4OutputTestqq);
	TreeTrainqq->SetBranchAddress("NN5Output", &NN5OutputTrainqq);
	TreeTestqq->SetBranchAddress("NN5Output", &NN5OutputTestqq);
	TreeTrainqq->SetBranchAddress("NN6Output", &NN6OutputTrainqq);
	TreeTestqq->SetBranchAddress("NN6Output", &NN6OutputTestqq);
	TreeTrainqq->SetBranchAddress("NN7Output", &NN7OutputTrainqq);
	TreeTestqq->SetBranchAddress("NN7Output", &NN7OutputTestqq);
	TreeTrainttbar->SetBranchAddress("entryIndex", &entryIndexTrain);
	TreeTestttbar->SetBranchAddress("entryIndex", &entryIndexTest);
	TreeTrainttbar->SetBranchAddress("NN1Output", &NN1OutputTrainttbar);
	TreeTestttbar->SetBranchAddress("NN1Output", &NN1OutputTestttbar);
	TreeTrainttbar->SetBranchAddress("NN2Output", &NN2OutputTrainttbar);
	TreeTestttbar->SetBranchAddress("NN2Output", &NN2OutputTestttbar);
	TreeTrainttbar->SetBranchAddress("NN3Output", &NN3OutputTrainttbar);
	TreeTestttbar->SetBranchAddress("NN3Output", &NN3OutputTestttbar);
	TreeTrainttbar->SetBranchAddress("NN4Output", &NN4OutputTrainttbar);
	TreeTestttbar->SetBranchAddress("NN4Output", &NN4OutputTestttbar);
	TreeTrainttbar->SetBranchAddress("NN5Output", &NN5OutputTrainttbar);
	TreeTestttbar->SetBranchAddress("NN5Output", &NN5OutputTestttbar);
	TreeTrainttbar->SetBranchAddress("NN6Output", &NN6OutputTrainttbar);
	TreeTestttbar->SetBranchAddress("NN6Output", &NN6OutputTestttbar);
	TreeTrainttbar->SetBranchAddress("NN7Output", &NN7OutputTrainttbar);
	TreeTestttbar->SetBranchAddress("NN7Output", &NN7OutputTestttbar);
	TreeTrainZZ->SetBranchAddress("entryIndex", &entryIndexTrain);
	TreeTestZZ->SetBranchAddress("entryIndex", &entryIndexTest);
	TreeTrainZZ->SetBranchAddress("NN1Output", &NN1OutputTrainZZ);
	TreeTestZZ->SetBranchAddress("NN1Output", &NN1OutputTestZZ);
	TreeTrainZZ->SetBranchAddress("NN2Output", &NN2OutputTrainZZ);
	TreeTestZZ->SetBranchAddress("NN2Output", &NN2OutputTestZZ);
	TreeTrainZZ->SetBranchAddress("NN3Output", &NN3OutputTrainZZ);
	TreeTestZZ->SetBranchAddress("NN3Output", &NN3OutputTestZZ);
	TreeTrainZZ->SetBranchAddress("NN4Output", &NN4OutputTrainZZ);
	TreeTestZZ->SetBranchAddress("NN4Output", &NN4OutputTestZZ);
	TreeTrainZZ->SetBranchAddress("NN5Output", &NN5OutputTrainZZ);
	TreeTestZZ->SetBranchAddress("NN5Output", &NN5OutputTestZZ);
	TreeTrainZZ->SetBranchAddress("NN6Output", &NN6OutputTrainZZ);
	TreeTestZZ->SetBranchAddress("NN6Output", &NN6OutputTestZZ);
	TreeTrainZZ->SetBranchAddress("NN7Output", &NN7OutputTrainZZ);
	TreeTestZZ->SetBranchAddress("NN7Output", &NN7OutputTestZZ);
	TreeTrainWW->SetBranchAddress("entryIndex", &entryIndexTrain);
	TreeTestWW->SetBranchAddress("entryIndex", &entryIndexTest);
	TreeTrainWW->SetBranchAddress("NN1Output", &NN1OutputTrainWW);
	TreeTestWW->SetBranchAddress("NN1Output", &NN1OutputTestWW);
	TreeTrainWW->SetBranchAddress("NN2Output", &NN2OutputTrainWW);
	TreeTestWW->SetBranchAddress("NN2Output", &NN2OutputTestWW);
	TreeTrainWW->SetBranchAddress("NN3Output", &NN3OutputTrainWW);
	TreeTestWW->SetBranchAddress("NN3Output", &NN3OutputTestWW);
	TreeTrainWW->SetBranchAddress("NN4Output", &NN4OutputTrainWW);
	TreeTestWW->SetBranchAddress("NN4Output", &NN4OutputTestWW);
	TreeTrainWW->SetBranchAddress("NN5Output", &NN5OutputTrainWW);
	TreeTestWW->SetBranchAddress("NN5Output", &NN5OutputTestWW);
	TreeTrainWW->SetBranchAddress("NN6Output", &NN6OutputTrainWW);
	TreeTestWW->SetBranchAddress("NN6Output", &NN6OutputTestWW);
	TreeTrainWW->SetBranchAddress("NN7Output", &NN7OutputTrainWW);
	TreeTestWW->SetBranchAddress("NN7Output", &NN7OutputTestWW);
	TreeTrainqqX->SetBranchAddress("entryIndex", &entryIndexTrain);
	TreeTestqqX->SetBranchAddress("entryIndex", &entryIndexTest);
	TreeTrainqqX->SetBranchAddress("NN1Output", &NN1OutputTrainqqX);
	TreeTestqqX->SetBranchAddress("NN1Output", &NN1OutputTestqqX);
	TreeTrainqqX->SetBranchAddress("NN2Output", &NN2OutputTrainqqX);
	TreeTestqqX->SetBranchAddress("NN2Output", &NN2OutputTestqqX);
	TreeTrainqqX->SetBranchAddress("NN3Output", &NN3OutputTrainqqX);
	TreeTestqqX->SetBranchAddress("NN3Output", &NN3OutputTestqqX);
	TreeTrainqqX->SetBranchAddress("NN4Output", &NN4OutputTrainqqX);
	TreeTestqqX->SetBranchAddress("NN4Output", &NN4OutputTestqqX);
	TreeTrainqqX->SetBranchAddress("NN5Output", &NN5OutputTrainqqX);
	TreeTestqqX->SetBranchAddress("NN5Output", &NN5OutputTestqqX);
	TreeTrainqqX->SetBranchAddress("NN6Output", &NN6OutputTrainqqX);
	TreeTestqqX->SetBranchAddress("NN6Output", &NN6OutputTestqqX);
	TreeTrainqqX->SetBranchAddress("NN7Output", &NN7OutputTrainqqX);
	TreeTestqqX->SetBranchAddress("NN7Output", &NN7OutputTestqqX);
	TreeTrainqqqqX->SetBranchAddress("entryIndex", &entryIndexTrain);
	TreeTestqqqqX->SetBranchAddress("entryIndex", &entryIndexTest);
	TreeTrainqqqqX->SetBranchAddress("NN1Output", &NN1OutputTrainqqqqX);
	TreeTestqqqqX->SetBranchAddress("NN1Output", &NN1OutputTestqqqqX);
	TreeTrainqqqqX->SetBranchAddress("NN2Output", &NN2OutputTrainqqqqX);
	TreeTestqqqqX->SetBranchAddress("NN2Output", &NN2OutputTestqqqqX);
	TreeTrainqqqqX->SetBranchAddress("NN3Output", &NN3OutputTrainqqqqX);
	TreeTestqqqqX->SetBranchAddress("NN3Output", &NN3OutputTestqqqqX);
	TreeTrainqqqqX->SetBranchAddress("NN4Output", &NN4OutputTrainqqqqX);
	TreeTestqqqqX->SetBranchAddress("NN4Output", &NN4OutputTestqqqqX);
	TreeTrainqqqqX->SetBranchAddress("NN5Output", &NN5OutputTrainqqqqX);
	TreeTestqqqqX->SetBranchAddress("NN5Output", &NN5OutputTestqqqqX);
	TreeTrainqqqqX->SetBranchAddress("NN6Output", &NN6OutputTrainqqqqX);
	TreeTestqqqqX->SetBranchAddress("NN6Output", &NN6OutputTestqqqqX);
	TreeTrainqqqqX->SetBranchAddress("NN7Output", &NN7OutputTrainqqqqX);
	TreeTestqqqqX->SetBranchAddress("NN7Output", &NN7OutputTestqqqqX);
	TreeTrainqqHX->SetBranchAddress("entryIndex", &entryIndexTrain);
	TreeTestqqHX->SetBranchAddress("entryIndex", &entryIndexTest);
	TreeTrainqqHX->SetBranchAddress("NN1Output", &NN1OutputTrainqqHX);
	TreeTestqqHX->SetBranchAddress("NN1Output", &NN1OutputTestqqHX);
	TreeTrainqqHX->SetBranchAddress("NN2Output", &NN2OutputTrainqqHX);
	TreeTestqqHX->SetBranchAddress("NN2Output", &NN2OutputTestqqHX);
	TreeTrainqqHX->SetBranchAddress("NN3Output", &NN3OutputTrainqqHX);
	TreeTestqqHX->SetBranchAddress("NN3Output", &NN3OutputTestqqHX);
	TreeTrainqqHX->SetBranchAddress("NN4Output", &NN4OutputTrainqqHX);
	TreeTestqqHX->SetBranchAddress("NN4Output", &NN4OutputTestqqHX);
	TreeTrainqqHX->SetBranchAddress("NN5Output", &NN5OutputTrainqqHX);
	TreeTestqqHX->SetBranchAddress("NN5Output", &NN5OutputTestqqHX);
	TreeTrainqqHX->SetBranchAddress("NN6Output", &NN6OutputTrainqqHX);
	TreeTestqqHX->SetBranchAddress("NN6Output", &NN6OutputTestqqHX);
	TreeTrainqqHX->SetBranchAddress("NN7Output", &NN7OutputTrainqqHX);
	TreeTestqqHX->SetBranchAddress("NN7Output", &NN7OutputTestqqHX);


	string outputTreeSNNTrainText = "analysis/outputTreeSNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	TFile *outputTreeSNNTrain = new TFile(outputTreeSNNTrainText.c_str(), "recreate");
 	TTree TreeSNNTrain("TreeSNNTrain","a simple Tree with simple variables (Train)");
 	string outputTreeSNNTestText = "analysis/outputTreeSNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeSNNTest = new TFile(outputTreeSNNTestText.c_str(), "recreate");
  	TTree TreeSNNTest("TreeSNNTest","a simple Tree with simple variables (Test)");
  	string outputTreeBqqNNTrainText = "analysis/outputTreeBqqNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	TFile *outputTreeBqqNNTrain = new TFile(outputTreeBqqNNTrainText.c_str(), "recreate");
 	TTree TreeBqqNNTrain("TreeBqqNNTrain","a bqqimple Tree with simple variables (Train)");
 	string outputTreeBqqNNTestText = "analysis/outputTreeBqqNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeBqqNNTest = new TFile(outputTreeBqqNNTestText.c_str(), "recreate");
  	TTree TreeBqqNNTest("TreeBqqNNTest","a bqqimple Tree with simple variables (Test)");
  	string outputTreeBttbarNNTrainText = "analysis/outputTreeBttbarNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	TFile *outputTreeBttbarNNTrain = new TFile(outputTreeBttbarNNTrainText.c_str(), "recreate");
 	TTree TreeBttbarNNTrain("TreeBttbarNNTrain","a bttbarimple Tree with simple variables (Train)");
 	string outputTreeBttbarNNTestText = "analysis/outputTreeBttbarNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeBttbarNNTest = new TFile(outputTreeBttbarNNTestText.c_str(), "recreate");
  	TTree TreeBttbarNNTest("TreeBttbarNNTest","a bttbarimple Tree with simple variables (Test)");
  	string outputTreeBZZNNTrainText = "analysis/outputTreeBZZNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	TFile *outputTreeBZZNNTrain = new TFile(outputTreeBZZNNTrainText.c_str(), "recreate");
 	TTree TreeBZZNNTrain("TreeBZZNNTrain","a bZZimple Tree with simple variables (Train)");
 	string outputTreeBZZNNTestText = "analysis/outputTreeBZZNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeBZZNNTest = new TFile(outputTreeBZZNNTestText.c_str(), "recreate");
  	TTree TreeBZZNNTest("TreeBZZNNTest","a bZZimple Tree with simple variables (Test)");
  	string outputTreeBWWNNTrainText = "analysis/outputTreeBWWNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	TFile *outputTreeBWWNNTrain = new TFile(outputTreeBWWNNTrainText.c_str(), "recreate");
 	TTree TreeBWWNNTrain("TreeBWWNNTrain","a bWWimple Tree with simple variables (Train)");
 	string outputTreeBWWNNTestText = "analysis/outputTreeBWWNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	TFile *outputTreeBWWNNTest = new TFile(outputTreeBWWNNTestText.c_str(), "recreate");
  	TTree TreeBWWNNTest("TreeBWWNNTest","a bWWimple Tree with simple variables (Test)");
	string outputTreeBqqXNNTrainText = "analysis/outputTreeBqqXNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	TFile *outputTreeBqqXNNTrain = new TFile(outputTreeBqqXNNTrainText.c_str(), "recreate");
	TTree TreeBqqXNNTrain("TreeBqqXNNTrain","a bqqXimple Tree with simple variables (Train)");
	string outputTreeBqqXNNTestText = "analysis/outputTreeBqqXNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile *outputTreeBqqXNNTest = new TFile(outputTreeBqqXNNTestText.c_str(), "recreate");
	TTree TreeBqqXNNTest("TreeBqqXNNTest","a bqqXimple Tree with simple variables (Test)");
	string outputTreeBqqqqXNNTrainText = "analysis/outputTreeBqqqqXNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	TFile *outputTreeBqqqqXNNTrain = new TFile(outputTreeBqqqqXNNTrainText.c_str(), "recreate");
	TTree TreeBqqqqXNNTrain("TreeBqqqqXNNTrain","a bqqqqXimple Tree with simple variables (Train)");
	string outputTreeBqqqqXNNTestText = "analysis/outputTreeBqqqqXNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile *outputTreeBqqqqXNNTest = new TFile(outputTreeBqqqqXNNTestText.c_str(), "recreate");
	TTree TreeBqqqqXNNTest("TreeBqqqqXNNTest","a bqqqqXimple Tree with simple variables (Test)");
	string outputTreeBqqHXNNTrainText = "analysis/outputTreeBqqHXNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	TFile *outputTreeBqqHXNNTrain = new TFile(outputTreeBqqHXNNTrainText.c_str(), "recreate");
	TTree TreeBqqHXNNTrain("TreeBqqHXNNTrain","a bqqHXimple Tree with simple variables (Train)");
	string outputTreeBqqHXNNTestText = "analysis/outputTreeBqqHXNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	TFile *outputTreeBqqHXNNTest = new TFile(outputTreeBqqHXNNTestText.c_str(), "recreate");
	TTree TreeBqqHXNNTest("TreeBqqHXNNTest","a bqqHXimple Tree with simple variables (Test)");


  	float entryIndex, NN1Output, NN2Output, NN3Output, NN4Output, NN5Output, NN6Output, NN7Output;
  	TreeSNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeSNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeSNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeSNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeSNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeSNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
  	TreeSNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
  	TreeSNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");	
  	TreeBqqNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBqqNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBqqNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBqqNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBqqNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
  	TreeBqqNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBttbarNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBttbarNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBttbarNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBttbarNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBttbarNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBttbarNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBttbarNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBttbarNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBZZNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBZZNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBZZNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBZZNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBZZNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBZZNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBZZNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBZZNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	TreeBWWNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
  	TreeBWWNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBWWNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBWWNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBWWNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBWWNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBWWNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBWWNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqXNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqXNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqXNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqXNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqXNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqXNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqXNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqXNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqqqXNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqqqXNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqqqXNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqqqXNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqqqXNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqqqXNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqqqXNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqqqXNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqHXNNTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqHXNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqHXNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqHXNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqHXNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqHXNNTrain.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqHXNNTrain.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqHXNNTrain.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	
  	TreeSNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeSNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeSNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeSNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeSNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeSNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeSNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeSNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F"); 
	TreeBqqNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBttbarNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBttbarNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBttbarNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBttbarNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBttbarNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBttbarNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBttbarNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBttbarNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBZZNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBZZNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBZZNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBZZNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBZZNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBZZNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBZZNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBZZNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBWWNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBWWNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBWWNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBWWNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBWWNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBWWNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBWWNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBWWNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqXNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqXNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqXNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqXNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqXNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqXNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqXNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqXNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqqqXNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqqqXNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqqqXNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqqqXNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqqqXNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqqqXNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqqqXNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqqqXNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
	TreeBqqHXNNTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeBqqHXNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
	TreeBqqHXNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
	TreeBqqHXNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
	TreeBqqHXNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
	TreeBqqHXNNTest.Branch("NN5Output",&NN5Output,"NN5Output/F");
	TreeBqqHXNNTest.Branch("NN6Output",&NN6Output,"NN6Output/F");
	TreeBqqHXNNTest.Branch("NN7Output",&NN7Output,"NN7Output/F");
  	
  	int contTrainHH=0, contTestHH=0, contTrainqq=0, contTestqq=0, contTrainttbar=0, contTestttbar=0, contTrainZZ=0, contTestZZ=0, contTrainWW=0, contTestWW=0, contTrainqqX=0, contTestqqX=0, contTrainqqqqX=0, contTestqqqqX=0, contTrainqqHX=0, contTestqqHX=0;
  	for(int i=0; i<BDTqqOutput.size(); i++)
  	{
  		entryIndex = i;
  		NN1Output = BDTqqOutput[i];
  		NN2Output = BDTttbarOutput[i];
  		NN3Output = BDTZZOutput[i];
  		NN4Output = BDTWWOutput[i];
		NN5Output = BDTqqXOutput[i];
  		NN6Output = BDTqqqqXOutput[i];
		NN7Output = BDTqqHXOutput[i];
  		
  		if(i<sizeHH)
  		{
  			TreeTrainHH->GetEntry(contTrainHH);
  			TreeTestHH->GetEntry(contTestHH);
  			if(entryIndex == entryIndexTrain)
  			{
  				TreeSNNTrain.Fill();
  				contTrainHH++;
  			}
  			else if(entryIndex == entryIndexTest)
  			{
  				TreeSNNTest.Fill();
  				contTestHH++;
  			}
			else throw std::runtime_error("ERROR: no match for HH entry in regenerate!");		
  		}
  		else if(i<sizeHH+sizeqq)
  		{
  			TreeTrainqq->GetEntry(contTrainqq);
  			TreeTestqq->GetEntry(contTestqq);
  			if(entryIndex == entryIndexTrain)
  			{
  				TreeBqqNNTrain.Fill();
  				contTrainqq++;
  			}
  			else if(entryIndex == entryIndexTest)
  			{
  				TreeBqqNNTest.Fill();
  				contTestqq++;
  			}
			else throw std::runtime_error("ERROR: no match for qq entry in regenerate!");			
  		}
  		else if(i<sizeHH+sizeqq+sizettbar)
  		{
  			TreeTrainttbar->GetEntry(contTrainttbar);
  			TreeTestttbar->GetEntry(contTestttbar);
  			if(entryIndex == entryIndexTrain)
  			{
  				TreeBttbarNNTrain.Fill();
  				contTrainttbar++;
  			}
  			else if(entryIndex == entryIndexTest)
  			{
  				TreeBttbarNNTest.Fill();
  				contTestttbar++;
  			}
			else throw std::runtime_error("ERROR: no match for ttbar entry in regenerate!");			
  		}
  		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ)
  		{
  			TreeTrainZZ->GetEntry(contTrainZZ);
  			TreeTestZZ->GetEntry(contTestZZ);
  			if(entryIndex == entryIndexTrain)
  			{
  				TreeBZZNNTrain.Fill();
  				contTrainZZ++;
  			}
  			else if(entryIndex == entryIndexTest)
  			{
  				TreeBZZNNTest.Fill();
  				contTestZZ++;
  			}
			else throw std::runtime_error("ERROR: no match for ZZ entry in regenerate!");			
  		}
  		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW)
  		{
  			TreeTrainWW->GetEntry(contTrainWW);
  			TreeTestWW->GetEntry(contTestWW);
  			if(entryIndex == entryIndexTrain)
  			{
  				TreeBWWNNTrain.Fill();
  				contTrainWW++;
  			}
  			else if(entryIndex == entryIndexTest)
  			{
  				TreeBWWNNTest.Fill();
  				contTestWW++;
  			}
			else throw std::runtime_error("ERROR: no match for WW entry in regenerate!");				
  		}
		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX)
  		{
  			TreeTrainqqX->GetEntry(contTrainqqX);
  			TreeTestqqX->GetEntry(contTestqqX);
  			if(entryIndex == entryIndexTrain)
  			{
  				TreeBqqXNNTrain.Fill();
  				contTrainqqX++;
  			}
  			else if(entryIndex == entryIndexTest)
  			{
  				TreeBqqXNNTest.Fill();
  				contTestqqX++;
  			}
			else throw std::runtime_error("ERROR: no match for qqX entry in regenerate!");				
  		}			
		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX)
  		{
  			TreeTrainqqqqX->GetEntry(contTrainqqqqX);
  			TreeTestqqqqX->GetEntry(contTestqqqqX);
  			if(entryIndex == entryIndexTrain)
  			{
  				TreeBqqqqXNNTrain.Fill();
  				contTrainqqqqX++;
  			}
  			else if(entryIndex == entryIndexTest)
  			{
  				TreeBqqqqXNNTest.Fill();
  				contTestqqqqX++;
  			}
			else throw std::runtime_error("ERROR: no match for qqqqX entry in regenerate!");	
		}
		else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX+sizeqqHX)
  		{
  			TreeTrainqqHX->GetEntry(contTrainqqHX);
  			TreeTestqqHX->GetEntry(contTestqqHX);
  			if(entryIndex == entryIndexTrain)
  			{
  				TreeBqqHXNNTrain.Fill();
  				contTrainqqHX++;
  			}
  			else if(entryIndex == entryIndexTest)
  			{
  				TreeBqqHXNNTest.Fill();
  				contTestqqHX++;
  			}
			else throw std::runtime_error("ERROR: no match for qqHX entry in regenerate!");
		}
  		else throw std::runtime_error("ERROR: event out of bounds in regenerate!");	
  	}
  	
  	outputTreeSNNTrain->cd();
   	TreeSNNTrain.Write();
   	outputTreeBqqNNTrain->cd();
   	TreeBqqNNTrain.Write();
   	outputTreeBttbarNNTrain->cd();
   	TreeBttbarNNTrain.Write();
   	outputTreeBZZNNTrain->cd();
   	TreeBZZNNTrain.Write();
   	outputTreeBWWNNTrain->cd();
   	TreeBWWNNTrain.Write();
  	outputTreeBqqXNNTrain->cd();
	TreeBqqXNNTrain.Write();
	outputTreeBqqqqXNNTrain->cd();
	TreeBqqqqXNNTrain.Write();
	outputTreeBqqHXNNTrain->cd();
	TreeBqqHXNNTrain.Write();
   	
   	outputTreeSNNTest->cd();
   	TreeSNNTest.Write();
   	outputTreeBqqNNTest->cd();
   	TreeBqqNNTest.Write();
   	outputTreeBttbarNNTest->cd();
   	TreeBttbarNNTest.Write();
   	outputTreeBZZNNTest->cd();
   	TreeBZZNNTest.Write();
   	outputTreeBWWNNTest->cd();
   	TreeBWWNNTest.Write();
	outputTreeBqqXNNTest->cd();
	TreeBqqXNNTest.Write();
	outputTreeBqqqqXNNTest->cd();
	TreeBqqqqXNNTest.Write();
	outputTreeBqqHXNNTest->cd();
	TreeBqqHXNNTest.Write();
}


void checkSigma(vector<double>& BDTqqOutput, int sizeHH, int sizeqq, double weightHH, double weightqq)
{
	double signalPass=0, signalCut=0, backPass=0, backCut=0;
	for(int i=0; i<BDTqqOutput.size(); i++)
  	{
  		double NN1Output = BDTqqOutput[i];
  		
  		
		if(i<sizeHH)
		{
			if(NN1Output > .6977) signalPass++;
			else signalCut++;
		}
		else if(i<sizeHH+sizeqq)
		{
			if(NN1Output > .6977) backPass++;
			else backCut++;
		}

						
  	}
  	
  	cout<<endl<<endl<<"From checkSigma: "<<endl<<"remainingHH unweighted: "<<signalPass<<endl<<"remainingBackUnweighted: "<<backPass<<endl<<endl;
	cout<<"remainingHH: "<<signalPass*weightHH<<endl<<"remainingBack: "<<backPass*weightqq<<endl<<endl<<endl;
	cout<<"significance: "<<(signalPass*weightHH)/sqrt((signalPass*weightHH)+(backPass*weightqq))<<endl<<endl;
  	
}

//////function to plot the outputs of the NNs
void plotNNOutputs(double bottomHistLimit, double topHistLimit, int nbinNN, vector<double>& BDTOutput, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int sizeqqX, int sizeqqqqX, int sizeqqHX, string title)
{
	string completeTitle = title + " on HH";
	TH1F *histNNHH = new TH1F("NNHH", completeTitle.c_str(), nbinNN, bottomHistLimit-.05, topHistLimit+.05);
	histNNHH->SetLineColor(kBlue);
	histNNHH->SetFillColor(kBlue);
	for(int i=0; i<sizeHH; i++)
	{
		histNNHH->Fill(BDTOutput[i]);
	}
	completeTitle = title + " on qq";
	TH1F *histNNqq = new TH1F("NNqq", completeTitle.c_str(), nbinNN, bottomHistLimit-.05, topHistLimit+.05);
	histNNqq->SetLineColor(kRed+2);
	histNNqq->SetFillColor(kRed+2);
	for(int i=sizeHH; i<sizeHH+sizeqq; i++)
	{
		histNNqq->Fill(BDTOutput[i]);
	}
	completeTitle = title + " on ttbar";
	TH1F *histNNttbar = new TH1F("NNttbar", completeTitle.c_str(), nbinNN, bottomHistLimit-.05, topHistLimit+.05);
	histNNttbar->SetLineColor(kPink+6);
	histNNttbar->SetFillColor(kPink+6);
	for(int i=sizeHH+sizeqq; i<sizeHH+sizeqq+sizettbar; i++)
	{
		histNNttbar->Fill(BDTOutput[i]);
	}
	completeTitle = title + " on ZZ";
	TH1F *histNNZZ = new TH1F("NNZZ", completeTitle.c_str(), nbinNN, bottomHistLimit-.05, topHistLimit+.05);
	histNNZZ->SetLineColor(kGreen+1);
	histNNZZ->SetFillColor(kGreen+1);
	for(int i=sizeHH+sizeqq+sizettbar; i<sizeHH+sizeqq+sizettbar+sizeZZ; i++)
	{
		histNNZZ->Fill(BDTOutput[i]);
	}
	completeTitle = title + " on WW";
	TH1F *histNNWW = new TH1F("NNWW", completeTitle.c_str(), nbinNN, bottomHistLimit-.05, topHistLimit+.05);
	histNNWW->SetLineColor(kYellow-2);
	histNNWW->SetFillColor(kYellow-2);
	for(int i=sizeHH+sizeqq+sizettbar+sizeZZ; i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW; i++)
	{
		histNNWW->Fill(BDTOutput[i]);
	}
	completeTitle = title + " on qqX";
	TH1F *histNNqqX = new TH1F("NNqqX", completeTitle.c_str(), nbinNN, bottomHistLimit-.05, topHistLimit+.05);
	histNNqqX->SetLineColor(kTeal-6);
	histNNqqX->SetFillColor(kTeal-6);
	for(int i=sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW; i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX; i++)
	{
		histNNqqX->Fill(BDTOutput[i]);
	}
	completeTitle = title + " on qqqqX";
	TH1F *histNNqqqqX = new TH1F("NNqqqqX", completeTitle.c_str(), nbinNN, bottomHistLimit-.05, topHistLimit+.05);
	histNNqqqqX->SetLineColor(kMagenta+2);
	histNNqqqqX->SetFillColor(kMagenta+2);
	for(int i=sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX; i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX; i++)
	{
		histNNqqqqX->Fill(BDTOutput[i]);
	}
	completeTitle = title + " on qqHX";
	TH1F *histNNqqHX = new TH1F("NNqqHX", completeTitle.c_str(), nbinNN, bottomHistLimit-.05, topHistLimit+.05);
	histNNqqHX->SetLineColor(kOrange+7);
	histNNqqHX->SetFillColor(kOrange+7);
	for(int i=sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX; i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW+sizeqqX+sizeqqqqX+sizeqqHX; i++)
	{
		histNNqqHX->Fill(BDTOutput[i]);
	}

	string saveTitle = "analysis/generatePlots/analysisPlots/NNOutputs/" + title + "HH" + ".png";
	TCanvas *cPlotNNOutputs1 = new TCanvas();
   	cPlotNNOutputs1->cd();
	histNNHH->Draw();
	cPlotNNOutputs1->Update();
	cPlotNNOutputs1->SaveAs(saveTitle.c_str());

	saveTitle = "analysis/generatePlots/analysisPlots/NNOutputs/" + title + "qq" + ".png";
	TCanvas *cPlotNNOutputs2 = new TCanvas();
   	cPlotNNOutputs2->cd();
	histNNqq->Draw();
	cPlotNNOutputs2->Update();
	cPlotNNOutputs2->SaveAs(saveTitle.c_str());

	saveTitle = "analysis/generatePlots/analysisPlots/NNOutputs/" + title + "ttbar" + ".png";
	TCanvas *cPlotNNOutputs3 = new TCanvas();
   	cPlotNNOutputs3->cd();
	histNNttbar->Draw();
	cPlotNNOutputs3->Update();
	cPlotNNOutputs3->SaveAs(saveTitle.c_str());

	saveTitle = "analysis/generatePlots/analysisPlots/NNOutputs/" + title + "ZZ" + ".png";
	TCanvas *cPlotNNOutputs4 = new TCanvas();
   	cPlotNNOutputs4->cd();
	histNNZZ->Draw();
	cPlotNNOutputs4->Update();
	cPlotNNOutputs4->SaveAs(saveTitle.c_str());

	saveTitle = "analysis/generatePlots/analysisPlots/NNOutputs/" + title + "WW" + ".png";
	TCanvas *cPlotNNOutputs5 = new TCanvas();
   	cPlotNNOutputs5->cd();
	histNNWW->Draw();
	cPlotNNOutputs5->Update();
	cPlotNNOutputs5->SaveAs(saveTitle.c_str());

	saveTitle = "analysis/generatePlots/analysisPlots/NNOutputs/" + title + "qqX" + ".png";
	TCanvas *cPlotNNOutputs6 = new TCanvas();
   	cPlotNNOutputs6->cd();
		histNNqqX->Draw();
	cPlotNNOutputs6->Update();
	cPlotNNOutputs6->SaveAs(saveTitle.c_str());

	saveTitle = "analysis/generatePlots/analysisPlots/NNOutputs/" + title + "qqqqX" + ".png";
	TCanvas *cPlotNNOutputs7 = new TCanvas();
   	cPlotNNOutputs7->cd();
	histNNqqqqX->Draw();
	cPlotNNOutputs7->Update();
	cPlotNNOutputs7->SaveAs(saveTitle.c_str());

	saveTitle = "analysis/generatePlots/analysisPlots/NNOutputs/" + title + "qqHX" + ".png";
	TCanvas *cPlotNNOutputs8 = new TCanvas();
   	cPlotNNOutputs8->cd();
	histNNqqHX->Draw();
	cPlotNNOutputs8->Update();
	cPlotNNOutputs8->SaveAs(saveTitle.c_str());
}

 ///Function that finds the limits of the NN output
 void findHistLimits(string method, double& bottomHistLimit, double& topHistLimit, int& methodColor)
 {
 	if(method == "Likelihood")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "LikelihoodD")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 0.9999;
 	}
 	if(method == "LikelihoodPCA")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "LikelihoodKDE")
 	{
 		bottomHistLimit = -0.00001;
 		topHistLimit = 0.99999;
 	}
 	if(method == "LikelihoodMIX")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "PDERS")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "PDERSD")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "PDERSPCA")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "KNN")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "HMatrix")
 	{
 		bottomHistLimit = -0.95;
 		topHistLimit = 1.55;
 	}
 	if(method == "Fisher")
 	{
 		bottomHistLimit = -4.0;
 		topHistLimit = 4.0;
 	}
 	else if(method == "FisherG")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "BoostedFisher")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 2.0;
 	}
 	else if(method == "LD")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 2.0;
 	}
 	else if(method == "MLP")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 		methodColor = kViolet;
 	}
 	else if(method == "MLPBFGS")
 	{
 		bottomHistLimit = -1.25;
 		topHistLimit = 1.5;
 	}
 	else if(method == "MLPBNN")
 	{
 		bottomHistLimit = -1.25;
 		topHistLimit = 1.5;
 	}
 	else if(method == "CFMlpANN")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "TMlpANN")
 	{
 		bottomHistLimit = -1.3;
 		topHistLimit = 1.3;
 	}
 	else if(method == "DNN_GPU")
 	{
 		bottomHistLimit = -0.1;
 		topHistLimit = 1.1;
 	}
 	else if(method == "DNN_CPU")
 	{
 		bottomHistLimit = -0.1;
 		topHistLimit = 1.1;
 	}
 	else if(method == "BDT")
 	{
 		bottomHistLimit = -0.8;
 		topHistLimit = 0.8;
 		methodColor = kRed;
 	}
 	else if(method == "BDTG")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 		methodColor = kBlue;
 	}
 	else if(method == "BDTB")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "BDTD")
 	{
 		bottomHistLimit = -0.8;
 		topHistLimit = 0.8;
 	}
 	else if(method == "BDTF")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "RuleFit")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 2.0;
 	}
 	else if(method == "SVM_Gauss")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "SVM_Poly")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "SVM_Lin")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "FDA_MT")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 3.0;
 	}
 	else if(method == "FDA_GA")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 3.0;
 	}
 	else if(method == "Category")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 2.0;
 	}
 	else if(method == "Plugin")
 	{
 		bottomHistLimit = -0.8;
 		topHistLimit = 0.8;
 	}
 	else cout<<"ERROOOOOOOOOOOOOOOOOOR UNKNOWN METHOD"<<endl;
 	
 	return;
 }
 
 void FSRTMVAClassificationApplicationHHbbbbGeneratesNNs(string fileFunction, string rtdCut, string preselection, string varVersion, string sampleName)
 {
 	cout<<"Aquì empieza main()"<<endl;
 	gStyle->SetOptStat(0);
 	
 	string method = "BDTG";
 	string variables = varVersion + preselection;
 	//string rtdCut = "10";
 	
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

    cout<<"WeightFactorHH: "<<weightFactorHH<<endl;
    cout<<"WeightFactorqq: "<<weightFactorqq<<endl;
    cout<<"WeightFactorttbar: "<<weightFactorttbar<<endl;
    cout<<"WeightFactorZZ: "<<weightFactorZZ<<endl;
    cout<<"WeightFactorWW: "<<weightFactorWW<<endl;
	cout<<"WeightFactorqqX: "<<weightFactorqqX<<endl;
	cout<<"WeightFactorqqqqX: "<<weightFactorqqqqX<<endl;
	cout<<"WeightFactorqqHX: "<<weightFactorqqHX<<endl;	
    cout<<endl<<endl<<endl;
    	
 	
 	double bottomHistLimit=0.0, topHistLimit=0.0;
 	int methodColor;
 	findHistLimits(method, bottomHistLimit, topHistLimit, methodColor);
 	int sizeHH=0, sizeqq=0, sizettbar=0, sizeZZ=0, sizeWW=0, sizeqqX=0, sizeqqqqX=0, sizeqqHX=0;
 	int nbin=10000, nbinNN=100;
 	double weightHH=0.001225*weightFactorHH, weightqq=0.0349*weightFactorqq, weightttbar=0.503*weightFactorttbar, weightZZ=0.8167*weightFactorZZ, weightWW=0.5149*weightFactorWW, weightqqX=0.04347826*weightFactorqqX, weightqqqqX=0.04*weightFactorqqqqX, weightqqHX=0.001*weightFactorqqHX;
 	vector<double> BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput;
 	
 	
 	TH1F *histBdtqqHH     = new TH1F( "MVA_BDTqqHH",           "MVA_BDTqq on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtqqqq     = new TH1F( "MVA_BDTqqqq",           "MVA_BDTqq on qq",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtttbarHH     = new TH1F( "MVA_BDTttbarHH",           "MVA_BDTttbar on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtttbarttbar     = new TH1F( "MVA_BDTttbarttbar",           "MVA_BDTttbar on ttbar",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtZZHH     = new TH1F( "MVA_BDTZZHH",           "MVA_BDTZZ on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtZZZZ     = new TH1F( "MVA_BDTZZZZ",           "MVA_BDTZZ on ZZ",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtWWHH     = new TH1F( "MVA_BDTWWHH",           "MVA_BDTWW on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtWWWW     = new TH1F( "MVA_BDTWWWW",           "MVA_BDTWW on WW",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
	TH1F *histBdtqqXHH     = new TH1F( "MVA_BDTqqXHH",           "MVA_BDTqqX on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
   	TH1F *histBdtqqXqqX     = new TH1F( "MVA_BDTqqXqqX",           "MVA_BDTqqX on qqX",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
	TH1F *histBdtqqqqXHH     = new TH1F( "MVA_BDTqqqqXHH",           "MVA_BDTqqqqX on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
   	TH1F *histBdtqqqqXqqqqX     = new TH1F( "MVA_BDTqqqqXqqqqX",           "MVA_BDTqqqqX on qqqqX",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
	TH1F *histBdtqqHXHH     = new TH1F( "MVA_BDTqqHXHH",           "MVA_BDTqqHX on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
  	TH1F *histBdtqqHXqqHX     = new TH1F( "MVA_BDTqqHXqqHX",           "MVA_BDTqqHX on qqHX",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtqqHHW    = new TH1F( "MVA_BDTqqHHW",           "MVA_BDTqq on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtqqqqW     = new TH1F( "MVA_BDTqqqqW",           "MVA_BDTqq on qq",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtttbarHHW     = new TH1F( "MVA_BDTttbarHHW",           "MVA_BDTttbar on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtttbarttbarW     = new TH1F( "MVA_BDTttbarttbarW",           "MVA_BDTttbar on ttbar",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtZZHHW     = new TH1F( "MVA_BDTZZHHW",           "MVA_BDTZZ on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtZZZZW     = new TH1F( "MVA_BDTZZZZW",           "MVA_BDTZZ on ZZ",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtWWHHW     = new TH1F( "MVA_BDTWWHHW",           "MVA_BDTWW on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtWWWWW     = new TH1F( "MVA_BDTWWWWW",           "MVA_BDTWW on WW",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
	TH1F *histBdtqqXHHW     = new TH1F( "MVA_BDTqqXHHW",           "MVA_BDTqqX on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
   	TH1F *histBdtqqXqqXW     = new TH1F( "MVA_BDTqqXqqXW",           "MVA_BDTqqX on qqX",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
	TH1F *histBdtqqqqXHHW     = new TH1F( "MVA_BDTqqqqXHHW",           "MVA_BDTqqqqX on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
   	TH1F *histBdtqqqqXqqqqXW     = new TH1F( "MVA_BDTqqqqXqqqqXW",           "MVA_BDTqqqqX on qqqqX",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
	TH1F *histBdtqqHXHHW     = new TH1F( "MVA_BDTqqHXHHW",           "MVA_BDTqqHX on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
   	TH1F *histBdtqqHXqqHXW     = new TH1F( "MVA_BDTqqHXqqHXW",           "MVA_BDTqqHX on qqHX",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histDummy     = new TH1F( "histDummy",           "histDummy",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	
 	
 	histBdtqqHH->SetLineColor(kBlue);
 	histBdtqqqq->SetLineColor(kRed+2);
 	histBdtttbarHH->SetLineColor(kBlue);
 	histBdtttbarttbar->SetLineColor(kPink+6);
 	histBdtZZHH->SetLineColor(kBlue);
 	histBdtZZZZ->SetLineColor(kGreen+1);
 	histBdtWWHH->SetLineColor(kBlue);
 	histBdtWWWW->SetLineColor(kYellow-2);
	histBdtqqXHH->SetLineColor(kBlue);
	histBdtqqXqqX->SetLineColor(kTeal-6);
	histBdtqqqqXHH->SetLineColor(kBlue);
	histBdtqqqqXqqqqX->SetLineColor(kMagenta+2);
	histBdtqqHXHH->SetLineColor(kBlue);
	histBdtqqHXqqHX->SetLineColor(kOrange+7);	
 	histBdtqqHHW->SetLineColor(kBlue);
 	histBdtqqqqW->SetLineColor(kRed+2);
 	histBdtttbarHHW->SetLineColor(kBlue);
 	histBdtttbarttbarW->SetLineColor(kPink+6);
 	histBdtZZHHW->SetLineColor(kBlue);
 	histBdtZZZZW->SetLineColor(kGreen+1);
 	histBdtWWHHW->SetLineColor(kBlue);
 	histBdtWWWWW->SetLineColor(kYellow-2);
	histBdtqqXHHW->SetLineColor(kBlue);
	histBdtqqXqqXW->SetLineColor(kTeal-6);
	histBdtqqqqXHHW->SetLineColor(kBlue);
	histBdtqqqqXqqqqXW->SetLineColor(kMagenta+2);
	histBdtqqHXHHW->SetLineColor(kBlue);
	histBdtqqHXqqHXW->SetLineColor(kOrange+7);
 
 	histBdtqqHH->SetFillColor(kBlue);
 	histBdtqqqq->SetFillColor(kRed+2);
 	histBdtttbarHH->SetFillColor(kBlue);
 	histBdtttbarttbar->SetFillColor(kPink+6);
 	histBdtZZHH->SetFillColor(kBlue);
 	histBdtZZZZ->SetFillColor(kGreen+1);
 	histBdtWWHH->SetFillColor(kBlue);
 	histBdtWWWW->SetFillColor(kYellow-2);
	histBdtqqXHH->SetFillColor(kBlue);
	histBdtqqXqqX->SetFillColor(kTeal-6);
	histBdtqqqqXHH->SetFillColor(kBlue);
	histBdtqqqqXqqqqX->SetFillColor(kMagenta+2);
	histBdtqqHXHH->SetFillColor(kBlue);
	histBdtqqHXqqHX->SetFillColor(kOrange+7);
 	histBdtqqHHW->SetFillColor(kBlue);
 	histBdtqqqqW->SetFillColor(kRed+2);
 	histBdtttbarHHW->SetFillColor(kBlue);
 	histBdtttbarttbarW->SetFillColor(kPink+6);
 	histBdtZZHHW->SetFillColor(kBlue);
 	histBdtZZZZW->SetFillColor(kGreen+1);
 	histBdtWWHHW->SetFillColor(kBlue);
 	histBdtWWWWW->SetFillColor(kYellow-2);
	histBdtqqXHHW->SetFillColor(kBlue);
	histBdtqqXqqXW->SetFillColor(kTeal-6);
	histBdtqqqqXHHW->SetFillColor(kBlue);
	histBdtqqqqXqqqqXW->SetFillColor(kMagenta+2);
	histBdtqqHXHHW->SetFillColor(kBlue);
	histBdtqqHXqqHXW->SetFillColor(kOrange+7);
 	
 	double targetFractionHH=0.8;
 	double effiHH=0, effiqq=0, effittbar=0, effiZZ=0, effiWW=0, effiqqX=0, effiqqqqX=0, effiqqHX=0;
 	//pairMassCut(1, targetFractionHH, effiHH);
 	//pairMassCut(2, targetFractionHH, effiqq);
 	//pairMassCut(3, targetFractionHH, effittbar);
 	
 	/////For arguments: empty string, topology (1: HH, 2: qq, 3:ttbar, 4: ZZ, 5: WW, 6: qqX, 7: qqqqX, 8: qqHX), nBack (0: qq, 1: ttbar, 2: ZZ, 3: WW, 4: qqX, 5: qqqqX, 6: qqHX), vector for BDTqq output, vector for BDTttbar output, var for size of tree, nbin for hists, histogram for BDT output, weighted histogram for BDT output, weight, method, NNVars (e.g. All if using all variables, PairMasses if only using the pair masses, etc), self descriptive the rest up until preselection (from analysis)
 	//////////
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, nbinNN, *histBdtqqHH, *histBdtqqHHW, weightHH, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, nbinNN, *histBdtttbarHH, *histBdtttbarHHW, weightHH, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, nbinNN, *histBdtZZHH, *histBdtZZHHW, weightHH, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, nbinNN, *histBdtWWHH, *histBdtWWHHW, weightHH, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 4, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, nbinNN, *histBdtqqXHH, *histBdtqqXHHW, weightHH, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 5, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, nbinNN, *histBdtqqqqXHH, *histBdtqqqqXHHW, weightHH, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 6, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, nbinNN, *histBdtqqHXHH, *histBdtqqHXHHW, weightHH, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqq, nbinNN, *histBdtqqqq, *histBdtqqqqW, weightqq, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 4, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 5, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 6, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizettbar, nbinNN, *histBdtttbarttbar, *histBdtttbarttbarW, weightttbar, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 4, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 5, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 6, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightZZ, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightZZ, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeZZ, nbinNN, *histBdtZZZZ, *histBdtZZZZW, weightZZ, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 4, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 5, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 6, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeWW, nbinNN, *histBdtWWWW, *histBdtWWWWW, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 4, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 5, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 6, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 6, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqX, nbinNN, *histDummy, *histDummy, weightqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 6, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqX, nbinNN, *histDummy, *histDummy, weightqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 6, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqX, nbinNN, *histDummy, *histDummy, weightqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 6, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqX, nbinNN, *histDummy, *histDummy, weightqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 6, 4, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqX, nbinNN, *histBdtqqXqqX, *histBdtqqXqqXW, weightqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 6, 5, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqX, nbinNN, *histDummy, *histDummy, weightqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 6, 6, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqX, nbinNN, *histDummy, *histDummy, weightqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 7, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqqqX, nbinNN, *histDummy, *histDummy, weightqqqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 7, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqqqX, nbinNN, *histDummy, *histDummy, weightqqqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 7, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqqqX, nbinNN, *histDummy, *histDummy, weightqqqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 7, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqqqX, nbinNN, *histDummy, *histDummy, weightqqqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 7, 4, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqqqX, nbinNN, *histDummy, *histDummy, weightqqqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 7, 5, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqqqX, nbinNN, *histBdtqqqqXqqqqX, *histBdtqqqqXqqqqXW, weightqqqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 7, 6, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqqqX, nbinNN, *histDummy, *histDummy, weightqqqqX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 8, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqHX, nbinNN, *histDummy, *histDummy, weightqqHX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 8, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqHX, nbinNN, *histDummy, *histDummy, weightqqHX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 8, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqHX, nbinNN, *histDummy, *histDummy, weightqqHX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 8, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqHX, nbinNN, *histDummy, *histDummy, weightqqHX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 8, 4, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqHX, nbinNN, *histDummy, *histDummy, weightqqHX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 8, 5, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqHX, nbinNN, *histDummy, *histDummy, weightqqHX, method, variables, rtdCut, sampleName, varVersion, preselection);
	FSRTMVAClassificationApplicationHHbbbbHelper("", 8, 6, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeqqHX, nbinNN, *histBdtqqHXqqHX, *histBdtqqHXqqHXW, weightqqHX, method, variables, rtdCut, sampleName, varVersion, preselection);	
 	//////////////
 	
 	cout<<"BDTqqOutput size: "<<BDTqqOutput.size()<<endl;
 	cout<<"BDTttbarOutput size: "<<BDTttbarOutput.size()<<endl;
 	cout<<"BDTZZOutput size: "<<BDTZZOutput.size()<<endl;
 	cout<<"BDTWWOutput size: "<<BDTWWOutput.size()<<endl;
 	cout<<"BDTqqXOutput size: "<<BDTqqXOutput.size()<<endl;
 	cout<<"BDTqqqqXOutput size: "<<BDTqqqqXOutput.size()<<endl;
 	cout<<"BDTqqHXOutput size: "<<BDTqqHXOutput.size()<<endl;
 	
 	TH2F *histBDTsOutputHH = new TH2F("btsoutputHH", "NNqq vs NNttbar", 100, bottomHistLimit-0.1, topHistLimit+0.1, 100, bottomHistLimit-0.1, topHistLimit+0.1);
 	histBDTsOutputHH->GetXaxis()->SetTitle("NNqq output");
   	histBDTsOutputHH->GetYaxis()->SetTitle("NNttbar output");
 	//BDTOutputPlane(bottomHistLimit, topHistLimit, BDTqqOutput, BDTttbarOutput, sizeHH, sizeqq, sizettbar, *histBDTsOutputHH);
 	
 	///////////////////////
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqHH, kBlue, *histBdtqqqq, kRed+2, method + "qq on HH and qq (normalized)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqHHW, kBlue, *histBdtqqqqW, kRed+2, method + "qq on HH and qq (weighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtttbarHH, kBlue, *histBdtttbarttbar, kPink+6, method + "ttbar on HH and ttbar (normalized)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtttbarHHW, kBlue, *histBdtttbarttbarW, kPink+6, method + "ttbar on HH and ttbar (weighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtZZHH, kBlue, *histBdtZZZZ, kGreen+1, method + "ZZ on HH and ZZ (normalized)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtZZHHW, kBlue, *histBdtZZZZW, kGreen+1, method + "ZZ on HH and ZZ (weighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtWWHH, kBlue, *histBdtWWWW, kYellow-2, method + "WW on HH and WW (normalized)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtWWHHW, kBlue, *histBdtWWWWW, kYellow-2, method + "WW on HH and WW (weighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqXHH, kBlue, *histBdtqqXqqX, kTeal-6, method + "qqX on HH and qqX (normalized)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqXHHW, kBlue, *histBdtqqXqqXW, kTeal-6, method + "qqX on HH and qqX (weighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqqqXHH, kBlue, *histBdtqqqqXqqqqX, kMagenta+2, method + "qqqqX on HH and qqqqX (normalized)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqqqXHHW, kBlue, *histBdtqqqqXqqqqXW, kMagenta+2, method + "qqqqX on HH and qqqqX (weighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqHXHH, kBlue, *histBdtqqHXqqHX, kOrange+7, method + "qqHX on HH and qqHX (normalized)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqHXHHW, kBlue, *histBdtqqHXqqHXW, kOrange+7, method + "qqHX on HH and qqHX (weighted)");
 	/////////////////////////

	/*////////plotting outputs of NNs
	plotNNOutputs(bottomHistLimit, topHistLimit, nbinNN, BDTqqOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, "NNqq");
	plotNNOutputs(bottomHistLimit, topHistLimit, nbinNN, BDTttbarOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, "NNttbar");
	plotNNOutputs(bottomHistLimit, topHistLimit, nbinNN, BDTZZOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, "NNZZ");
	plotNNOutputs(bottomHistLimit, topHistLimit, nbinNN, BDTWWOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, "NNWW");
	plotNNOutputs(bottomHistLimit, topHistLimit, nbinNN, BDTqqXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, "NNqqX");
	plotNNOutputs(bottomHistLimit, topHistLimit, nbinNN, BDTqqqqXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, "NNqqqqX");
	plotNNOutputs(bottomHistLimit, topHistLimit, nbinNN, BDTqqHXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, "NNqqHX");
	/////plotting outputs of NNs*/

	
 	
 	//double n, b, optCut;
 	//double maxSignificance=0;
 	//maxSignificance = findMaxSignificanceAndROC(bottomHistLimit, topHistLimit, nbin, *histBdtqqHHW, *histBdtqqqqW, n, b, optCut);
 	
 	
 	double maxSignificanceqq=0, defCutqq, totalRemainingqq=0, HHRemainingqq=0, backRemainingqq=0;
 	double maxSignificancettbar=0, defCutttbar, totalRemainingttbar=0, HHRemainingttbar=0, backRemainingttbar=0;
 	double maxSignificanceZZ=0, defCutZZ, totalRemainingZZ=0, HHRemainingZZ=0, backRemainingZZ=0;
 	double maxSignificanceWW=0, defCutWW, totalRemainingWW=0, HHRemainingWW=0, backRemainingWW=0;
	double maxSignificanceqqX=0, defCutqqX, totalRemainingqqX=0, HHRemainingqqX=0, backRemainingqqX=0;
	double maxSignificanceqqqqX=0, defCutqqqqX, totalRemainingqqqqX=0, HHRemainingqqqqX=0, backRemainingqqqqX=0;
	double maxSignificanceqqHX=0, defCutqqHX, totalRemainingqqHX=0, HHRemainingqqHX=0, backRemainingqqHX=0;
 	double maxSignificanceCombined=0, defCutqqCombined, defCutttbarCombined=0, totalRemainingCombined=0, HHRemainingCombined=0, backRemainingCombined=0;
 	
 	TH2F *histROCqq = new TH2F("hROCqq", "Signal and Background Acceptance (NNqq)", nbin, 0.0, 1.0, nbin, 0, 1.0);
 	TH2F *histROCttbar = new TH2F("hROCttbar", "Signal and Background Acceptance (NNttbar)", nbin, 0.0, 1.0, nbin, 0, 1.0);
 	TH2F *histROCZZ = new TH2F("hROCZZ", "Signal and Background Acceptance (NNZZ)", nbin, 0.0, 1.0, nbin, 0, 1.0);
 	TH2F *histROCWW = new TH2F("hROCWW", "Signal and Background Acceptance (NNWW)", nbin, 0.0, 1.0, nbin, 0, 1.0);
	TH2F *histROCqqX = new TH2F("hROCqqX", "Signal and Background Acceptance (NNqqX)", nbin, 0.0, 1.0, nbin, 0, 1.0);
	TH2F *histROCqqqqX = new TH2F("hROCqqqqX", "Signal and Background Acceptance (NNqqqqX)", nbin, 0.0, 1.0, nbin, 0, 1.0);
	TH2F *histROCqqHX = new TH2F("hROCqqHX", "Signal and Background Acceptance (NNqqHX)", nbin, 0.0, 1.0, nbin, 0, 1.0);
 	TH2F *histROCqqttbar = new TH2F("hROCqqttbar", "Signal and Background Acceptance (NNqq and NNttbar)", nbin, 0.0, 1.0, nbin, 0, 1.0);
  	TH2F *histROCRejqq = new TH2F("hROCRejqq", "Signal Acceptance and Background Rejection (NNqq)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	TH2F *histROCRejttbar = new TH2F("hROCRejttbar", "Signal Acceptance and Background Rejection (NNttbar)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	TH2F *histROCRejZZ = new TH2F("hROCRejZZ", "Signal Acceptance and Background Rejection (NNZZ)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	TH2F *histROCRejWW = new TH2F("hROCRejWW", "Signal Acceptance and Background Rejection (NNWW)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
	TH2F *histROCRejqqX = new TH2F("hROCRejqqX", "Signal Acceptance and Background Rejection (NNqqX)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
	TH2F *histROCRejqqqqX = new TH2F("hROCRejqqqqX", "Signal Acceptance and Background Rejection (NNqqqqX)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
	TH2F *histROCRejqqHX = new TH2F("hROCRejqqHX", "Signal Acceptance and Background Rejection (NNqqHX)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	TH2F *histROCRejqqttbar = new TH2F("hROCRejqqttbar", "Signal Acceptance and Background Rejection (NNqq and NNttbar)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	
  	TH2F *histSignificanceqq     = new TH2F( "hist_significanceqq", "Significance s/sqrt(s+b) (NNqq)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceqq->GetXaxis()->SetTitle("Cut value for NNqq");
   	histSignificanceqq->GetYaxis()->SetTitle("Significance for NNqq");
   	TH2F *histSignificancettbar     = new TH2F( "hist_significancettbar", "Significance s/sqrt(s+b) (NNttbar)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificancettbar->GetXaxis()->SetTitle("Cut value for NNttbar)");
   	histSignificancettbar->GetYaxis()->SetTitle("Significance for NNttbar");
   	TH2F *histSignificanceZZ     = new TH2F( "hist_significanceZZ", "Significance s/sqrt(s+b) (NNZZ)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceZZ->GetXaxis()->SetTitle("Cut value for NNZZ");
   	histSignificanceZZ->GetYaxis()->SetTitle("Significance for NNZZ");
   	TH2F *histSignificanceWW     = new TH2F( "hist_significanceWW", "Significance s/sqrt(s+b) (NNWW)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceWW->GetXaxis()->SetTitle("Cut value for NNWW");
   	histSignificanceWW->GetYaxis()->SetTitle("Significance for NNWW");
	TH2F *histSignificanceqqX     = new TH2F( "hist_significanceqqX", "Significance s/sqrt(s+b) (NNqqX)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceqqX->GetXaxis()->SetTitle("Cut value for NNqqX");
   	histSignificanceqqX->GetYaxis()->SetTitle("Significance for NNqqX");
   	TH2F *histSignificanceqqqqX     = new TH2F( "hist_significanceqqqqX", "Significance s/sqrt(s+b) (NNqqqqX)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceqqqqX->GetXaxis()->SetTitle("Cut value for NNqqqqX");
   	histSignificanceqqqqX->GetYaxis()->SetTitle("Significance for NNqqqqX");
	TH2F *histSignificanceqqHX     = new TH2F( "hist_significanceqqHX", "Significance s/sqrt(s+b) (NNqqHX)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
	histSignificanceqqHX->GetXaxis()->SetTitle("Cut value for NNqqHX");
  	histSignificanceqqHX->GetYaxis()->SetTitle("Significance for NNqqHX");
   	TH2F *histSignificanceqqttbar     = new TH2F( "hist_significanceqqttbar", "Significance s/sqrt(s+b) (NNqq and NNttbar)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceqqttbar->GetXaxis()->SetTitle("Cut value for NNqq");
   	histSignificanceqqttbar->GetYaxis()->SetTitle("Significance for NNqq"); //////////CHECAR COMO DIMENSIONALIZAR
   	
   	///////////////////
   	cout<<endl<<"Significance for "<<method<<"qq applied only to HH and qq: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, 0, *histBdtqqHH, *histBdtqqqq, defCutqq, maxSignificanceqq, weightHH, weightqq, targetFractionHH, *histROCqq, *histROCRejqq, *histSignificanceqq, totalRemainingqq, HHRemainingqq, backRemainingqq);
 	cout<<"maxSignificance (qq): "<<maxSignificanceqq<<endl<<"Cut in NNqq for max significance (qq): "<<defCutqq<<endl;
 	
 	cout<<endl<<"Significance for "<<method<<"ttbar applied only to HH and ttbar: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTttbarOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, 1, *histBdtttbarHH, *histBdtttbarttbar, defCutttbar, maxSignificancettbar, weightHH, weightttbar, targetFractionHH, *histROCttbar, *histROCRejttbar, *histSignificancettbar, totalRemainingttbar, HHRemainingttbar, backRemainingttbar);
 	cout<<"maxSignificance (ttbar): "<<maxSignificancettbar<<endl<<"Cut in NNttbar for max significance (ttbar): "<<defCutttbar<<endl;
 	
 	cout<<endl<<"Significance for "<<method<<"ZZ applied only to HH and ZZ: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTZZOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, 2, *histBdtZZHH, *histBdtZZZZ, defCutZZ, maxSignificanceZZ, weightHH, weightZZ, targetFractionHH, *histROCZZ, *histROCRejZZ, *histSignificanceZZ, totalRemainingZZ, HHRemainingZZ, backRemainingZZ);
 	cout<<"maxSignificance (ZZ): "<<maxSignificanceZZ<<endl<<"Cut in NNZZ for max significance (ZZ): "<<defCutZZ<<endl;
 	
 	cout<<endl<<"Significance for "<<method<<"WW applied only to HH and WW: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTWWOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, 3, *histBdtWWHH, *histBdtWWWW, defCutWW, maxSignificanceWW, weightHH, weightWW, targetFractionHH, *histROCWW, *histROCRejWW, *histSignificanceWW, totalRemainingWW, HHRemainingWW, backRemainingWW);
 	cout<<"maxSignificance (WW): "<<maxSignificanceWW<<endl<<"Cut in NNWW for max significance (WW): "<<defCutWW<<endl;

	cout<<endl<<"Significance for "<<method<<"qqX applied only to HH and qqX: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTqqXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, 4, *histBdtqqXHH, *histBdtqqXqqX, defCutqqX, maxSignificanceqqX, weightHH, weightqqX, targetFractionHH, *histROCqqX, *histROCRejqqX, *histSignificanceqqX, totalRemainingqqX, HHRemainingqqX, backRemainingqqX);
 	cout<<"maxSignificance (qqX): "<<maxSignificanceqqX<<endl<<"Cut in NNqqX for max significance (qqX): "<<defCutqqX<<endl;

	cout<<endl<<"Significance for "<<method<<"qqqqX applied only to HH and qqqqX: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTqqqqXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, 5, *histBdtqqqqXHH, *histBdtqqqqXqqqqX, defCutqqqqX, maxSignificanceqqqqX, weightHH, weightqqqqX, targetFractionHH, *histROCqqqqX, *histROCRejqqqqX, *histSignificanceqqqqX, totalRemainingqqqqX, HHRemainingqqqqX, backRemainingqqqqX);
 	cout<<"maxSignificance (qqqqX): "<<maxSignificanceqqqqX<<endl<<"Cut in NNqqqqX for max significance (qqqqX): "<<defCutqqqqX<<endl;

	cout<<endl<<"Significance for "<<method<<"qqHX applied only to HH and qqHX: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTqqHXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, 6, *histBdtqqHXHH, *histBdtqqHXqqHX, defCutqqHX, maxSignificanceqqHX, weightHH, weightqqHX, targetFractionHH, *histROCqqHX, *histROCRejqqHX, *histSignificanceqqHX, totalRemainingqqHX, HHRemainingqqHX, backRemainingqqHX);
 	cout<<"maxSignificance (qqHX): "<<maxSignificanceqqHX<<endl<<"Cut in NNqqHX for max significance (qqHX): "<<defCutqqHX<<endl;
 	///////////////////////
 	
 	//cout<<endl<<"(using crazy for loops) Significance for "<<method<<"qq and "<<method<<"ttbar combined applied to HH, qq and, ttbar: "<<endl;
 	//findSignificanceCombined(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, BDTttbarOutput, sizeHH, sizeqq, sizettbar, 0, 1, defCutqqCombined, defCutttbarCombined, maxSignificanceCombined, weightHH, weightqq, weightttbar, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined);
 	//cout<<"maxSignificance (combined): "<<maxSignificanceCombined<<endl<<"Cut in NNqq for max significance (combined): "<<defCutqqCombined<<endl<<"Cut in NNttbar for max significance (combined): "<<defCutttbarCombined<<endl;
 	
 	/*//////////////
 	//defCutqq=0.636;
 	//defCutttbar=0.802;
 	cout<<endl<<"Significance for "<<method<<"qq, "<<method<<"ttbar, "<<method<<"ZZ, and "<<method<<"WW, "<<method<<"qqX, "<<method<<"qqqqX, and "<<method<<"qqHX combined applied to HH, qq, ttbar, ZZ, WW, qqX, qqqqX, and qqHX (trained and applied independently): "<<endl;
 	findSignificanceCutsCombined(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, defCutqq, defCutttbar, defCutZZ, defCutWW, defCutqqX, defCutqqqqX, defCutqqHX, maxSignificanceCombined, weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined);
 	cout<<"maxSignificance (combined): "<<maxSignificanceCombined<<endl<<endl<<endl;
 	//////////////*/
 	
 	///////Find significance by max. each NNs significance against its own back -- check / brute force
	cout<<endl<<"CHECK Significance for "<<method<<"qq, "<<method<<"ttbar, "<<method<<"ZZ, and "<<method<<"WW, "<<method<<"qqX, "<<method<<"qqqqX, and "<<method<<"qqHX combined applied to HH, qq, ttbar, ZZ, WW, qqX, qqqqX, and qqHX (trained and applied independently): "<<endl;
	findSignificanceCutsCombinedCheck (bottomHistLimit, topHistLimit, nbin, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, defCutqq, defCutttbar, defCutZZ, defCutWW, defCutqqX, defCutqqqqX, defCutqqHX, maxSignificanceCombined, weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined);
	cout<<"CHECK maxSignificance (combined): "<<maxSignificanceCombined<<endl<<endl<<endl;	
	double preliminarySignificance = maxSignificanceCombined;
	//////////////
 	
 	/*//////////DEBUGGING PURPOSES
 	vector<double> testNN1Output = {0.01, 0.61, 0.11, 0.81, 0.91, 0.98, 0.03, 0.933, 0.13, 0.23, 0.98, 0.85, 0.95, 0.92, 0.92};
 	vector<double> testNN2Output = {0.02, 0.03, 0.12, 0.82, 0.92, 0.98, 0.94, 0.95, 0.84, 0.94, 0.06, 0.066, 0.96, 0.93, 0.91};
 	vector<double> testNN3Output = {0.08, 0.16, 0.65, 0.89, 0.96, 0.98, 0.95, 0.96, 0.93, 0.82, 0.95, 0.84, 0.97, 0.91, 0.904};
 	//findSignificanceCombined(bottomHistLimit, topHistLimit, nbin, testNN1Output, testNN2Output, 6, 4, 3, 0, 1, defCutqqCombined, defCutttbarCombined, maxSignificanceCombined, weightHH, weightqq, weightttbar, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined); 
 	findSignificanceCutsCombined(bottomHistLimit, topHistLimit, nbin, testNN1Output, testNN2Output, testNN3Output,  6, 4, 3, 2, 0.6, 0.7, 0.8, maxSignificanceCombined, weightHH, weightqq, weightttbar, weightZZ, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined);
 	/////////DEBUGGING PURPOSES*/



	//////////////find significance by finding cuts for each NNs significance against all backs
	findDefCuts(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, defCutqq, defCutttbar, defCutZZ, defCutWW, defCutqqX, defCutqqqqX, defCutqqHX, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX);
	cout<<endl<<"(optimizing NNs for all backs.) CHECK Significance for "<<method<<"qq, "<<method<<"ttbar, "<<method<<"ZZ, and "<<method<<"WW, "<<method<<"qqX, "<<method<<"qqqqX, and "<<method<<"qqHX combined applied to HH, qq, ttbar, ZZ, WW, qqX, qqqqX, and qqHX (trained and applied independently): "<<endl;
	findSignificanceCutsCombinedCheck(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, defCutqq, defCutttbar, defCutZZ, defCutWW, defCutqqX, defCutqqqqX, defCutqqHX, maxSignificanceCombined, weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined);
	cout<<"(optimizing NNs for all backs.) CHECK maxSignificance (combined): "<<maxSignificanceCombined<<endl<<endl<<endl;	
	preliminarySignificance = maxSignificanceCombined;
 	/////////////find significance by finding cuts for each NNs significance against all backs

	////////////find def cuts for each eGamma NNs against eGamma backs and eliminate events with eGamma cuts
	findDefEGammaCuts(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, defCutqq, defCutttbar, defCutZZ, defCutWW, defCutqqX, defCutqqqqX, defCutqqHX, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, weightHH, weightqq, weightttbar, weightZZ, weightWW, weightqqX, weightqqqqX, weightqqHX);
	cout<<endl<<"defCutqqX: "<<defCutqqX<<endl<<"defCutqqqqX: "<<defCutqqqqX<<endl<<"defCutqqHX: "<<defCutqqHX<<endl;
	cutEGamma(BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, defCutqqX, defCutqqqqX, defCutqqHX);
	///////////////find def cuts for each eGamma NNs against eGamma backs and eliminate events with eGamma cuts 	
 	
 	//add2DHistToFile(*histROCRej, "ROCHists.root", "ROCRej"+method);
 	
 	//splitHistogram(*histBdtqqqq, bottomHistLimit, topHistLimit, defCut, "BDTqq on qq remaining after cut");
 	
 	//addPointTo2DHistogram(*histROCRej, effiHH, 1/effiqq, kRed, .005, 1);
 	
 	///For arguments: (name of file with hists, scale (log or linear), bool for whether a point will be graffed, x coord, y coord, point color, radius x, radius y 
 	//draw2DHistsFromFile("ROCHists.root", "log", true, effiHH, 1/effiqq, kGreen, .005, 1);
 	
 	
 	
 	
 	
 	
 	double luminosity = 4900, crossSectionqq, errorTopqq, errorBottomqq;
 	double crossSectionttbar, errorTopttbar, errorBottomttbar;
 	double crossSectionZZ, errorTopZZ, errorBottomZZ;
 	double crossSectionWW, errorTopWW, errorBottomWW;
	double crossSectionqqX, errorTopqqX, errorBottomqqX;
	double crossSectionqqqqX, errorTopqqqqX, errorBottomqqqqX;
	double crossSectionqqHX, errorTopqqHX, errorBottomqqHX;
 	double crossSectionCombined, errorTopCombined, errorBottomCombined;
 	
 	cout<<endl<<endl;
 	
 	/*cout<<"Cross-section and error for the output of "<<method<<"qq applied only to HH and qq: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingqq+backRemainingqq, HHRemainingqq, backRemainingqq, luminosity, crossSectionqq, errorTopqq, errorBottomqq, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"ttbar applied only to HH and ttbar: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingttbar+backRemainingttbar, HHRemainingttbar, backRemainingttbar, luminosity, crossSectionttbar, errorTopttbar, errorBottomttbar, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"ZZ applied only to HH and ZZ: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingZZ+backRemainingZZ, HHRemainingZZ, backRemainingZZ, luminosity, crossSectionttbar, errorTopZZ, errorBottomZZ, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"WW applied only to HH and WW: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingWW+backRemainingWW, HHRemainingWW, backRemainingWW, luminosity, crossSectionttbar, errorTopWW, errorBottomWW, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"qqX applied only to HH and qqX: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingqqX+backRemainingqqX, HHRemainingqqX, backRemainingqqX, luminosity, crossSectionqqX, errorTopqqX, errorBottomqqX, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"qqqqX applied only to HH and qqqqX: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingqqqqX+backRemainingqqqqX, HHRemainingqqqqX, backRemainingqqqqX, luminosity, crossSectionqqqqX, errorTopqqqqX, errorBottomqqqqX, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"qqHX applied only to HH and qqHX: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingqqHX+backRemainingqqHX, HHRemainingqqHX, backRemainingqqHX, luminosity, crossSectionqqHX, errorTopqqHX, errorBottomqqHX, nbin);*/
 	cout<<"Cross-section and error for the output of "<<method<<" combined applied to HH, qq, ttbar, ZZ, WW, qqX, qqqqX, qqHX: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingCombined+backRemainingCombined, HHRemainingCombined, backRemainingCombined, luminosity, crossSectionCombined, errorTopCombined, errorBottomCombined, nbin);
 	double preliminaryErrorLeft = errorBottomCombined;
	double preliminaryErrorRight = errorTopCombined;
 	

	//////saving preliminarySignificance in a txt
	ofstream outFile("analysis/preliminarySignificanceAndErrorsFile.txt");
	if (outFile.is_open()) {
	    outFile << preliminarySignificance << endl;
		outFile << preliminaryErrorLeft << endl;
		outFile << preliminaryErrorRight << endl;
	    outFile.close();
	} else {
	    cerr << "Unable to open file for writing preliminary significance" << endl;
	}
	/////////saving preliminarySignificance in a txt	
 	/////// 
 	
 	
 	
 	
 	
 	/*TCanvas *cMain2 = new TCanvas();
 	TCanvas *cMain3 = new TCanvas();
 	TCanvas *cMain4 = new TCanvas();
 	TCanvas *cMain5 = new TCanvas();
 	
 	cMain2->cd();
 	histBdtqqHH->Draw("HIST");
 	
 	cMain3->cd();
 	histBdtqqqq->Draw("HIST");
 	
 	cMain4->cd();
 	histBdtttbarHH->Draw("HIST");
 	
 	cMain5->cd();
 	histBdtttbarttbar->Draw("HIST");*/
 	
 	
 	cout<<endl<<endl<<endl<<"General info of input files: "<<endl;
 	cout<<"w HH size: "<<sizeHH*weightHH<<endl<<"uw HH size: "<<sizeHH<<endl<<endl;
 	cout<<"w qq size: "<<sizeqq*weightqq<<endl<<"uw qq size: "<<sizeqq<<endl<<endl;
 	cout<<"w ttbar size: "<<sizettbar*weightttbar<<endl<<"uw ttbar size: "<<sizettbar<<endl<<endl;
 	cout<<"w ZZ size: "<<sizeZZ*weightZZ<<endl<<"uw ZZ size: "<<sizeZZ<<endl<<endl;
 	cout<<"w WW size: "<<sizeWW*weightWW<<endl<<"uw WW size: "<<sizeWW<<endl<<endl;
 	cout<<"w qqX size: "<<sizeqqX*weightqqX<<endl<<"uw qqX size: "<<sizeqqX<<endl<<endl;
 	cout<<"w qqqqX size: "<<sizeqqqqX*weightqqqqX<<endl<<"uw qqqqX size: "<<sizeqqqqX<<endl<<endl;
 	cout<<"w qqHX size: "<<sizeqqHX*weightqqHX<<endl<<"uw qqHX size: "<<sizeqqHX<<endl<<endl;
 	
 	cout<<"BDTqqOutput size: "<<BDTqqOutput.size()<<endl;
 	cout<<"BDTttbarOutput size: "<<BDTttbarOutput.size()<<endl;
 	cout<<"BDTZZOutput size: "<<BDTZZOutput.size()<<endl;
 	cout<<"BDTWWOutput size: "<<BDTWWOutput.size()<<endl;
 	cout<<"BDTqqXOutput size: "<<BDTqqXOutput.size()<<endl;
 	cout<<"BDTqqqqXOutput size: "<<BDTqqqqXOutput.size()<<endl;
 	cout<<"BDTqqHXOutput size: "<<BDTqqHXOutput.size()<<endl;	
 	
 	///////Filling out tree for output NN
	if(fileFunction == "generateGA") generateFilesGA(BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, rtdCut, preselection, sampleName);
 	if(fileFunction == "generate") generateFilesOutputNN(BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, rtdCut, preselection, sampleName);
  	else if(fileFunction == "merge") reGenerateFilesOutputNN(BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, BDTqqXOutput, BDTqqqqXOutput, BDTqqHXOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, sizeqqX, sizeqqqqX, sizeqqHX, rtdCut, preselection, sampleName);
  	///////Filling out tree for output NN
  	
  	//checkSigma(BDTqqOutput, sizeHH, sizeqq, weightHH, weightqq);
 	
 	
 	
 	return;
 }
 
 /*int main( int argc, char** argv )
 {
    cout<<"Aquì empieza main()"<<endl;
    TString methodList;
    for (int i=1; i<argc; i++) {
       TString regMethod(argv[i]);
       if(regMethod=="-b" || regMethod=="--batch") continue;
       if (!methodList.IsNull()) methodList += TString(",");
       methodList += regMethod;
    }
    FSRTMVAClassificationApplicationHHbbbb(methodList);
    return 0;
 }*/


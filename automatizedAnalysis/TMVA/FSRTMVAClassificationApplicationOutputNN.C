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
 
 void FSRTMVAClassificationApplicationOutputNNHelper( TString myMethodList, vector<double>& NNOutput, int& sizeTree, int nbin, TH1F& histOutput, string inputMethod, TString NNVars, int topology, string rtdCut, string preselection, string varVersion, string sampleName)
 {
 	
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
 
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplicationOutputNN" << std::endl;
 
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
    Float_t NN1Output, NN2Output, NN3Output, NN4Output;
    reader->AddVariable( "NN1Output", &NN1Output );
    reader->AddVariable( "NN2Output", &NN2Output );
    reader->AddVariable( "NN3Output", &NN3Output );
    reader->AddVariable( "NN4Output", &NN4Output );
 
 
    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //
    string inputSText = "analysis/outputTreeSNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    string inputBqqText = "analysis/outputTreeBqqNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    string inputBttbarText = "analysis/outputTreeBttbarNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    string inputBZZText = "analysis/outputTreeBZZNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    string inputBWWText = "analysis/outputTreeBWWNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    TFile* inputS = new TFile(inputSText.c_str());
    TFile* inputBqq = new TFile(inputBqqText.c_str());
    TFile* inputBttbar = new TFile(inputBttbarText.c_str());
    TFile* inputBZZ = new TFile(inputBZZText.c_str());
    TFile* inputBWW = new TFile(inputBWWText.c_str());
    
    
    if(topology==0) std::cout << "--- TMVAClassificationApplicationOutputNN   : Using input file for signal: " << inputS->GetName() << std::endl;
    else if(topology==1) std::cout << "--- TMVAClassificationApplicationOutputNN   : Using input file for qq background: " << inputBqq->GetName() << std::endl;
    else if(topology==2) std::cout << "--- TMVAClassificationApplicationOutputNN   : Using input file for ttbar background: " << inputBttbar->GetName() << std::endl;
    else if(topology==3) std::cout << "--- TMVAClassificationApplicationOutputNN   : Using input file for ZZ background: " << inputBZZ->GetName() << std::endl;
    else if(topology==4) std::cout << "--- TMVAClassificationApplicationOutputNN   : Using input file for WW background: " << inputBWW->GetName() << std::endl;
    
    // Event loop
 
    // Prepare the event tree
    // - Here the variable names have to corresponds to your tree
    // - You can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //
    
    TTree* theTree;
    
    if(topology == 0) theTree = (TTree*)inputS->Get("TreeSNNTest");
    else if(topology == 1) theTree = (TTree*)inputBqq->Get("TreeBqqNNTest");
    else if(topology == 2) theTree = (TTree*)inputBttbar->Get("TreeBttbarNNTest");
    else if(topology == 3) theTree = (TTree*)inputBZZ->Get("TreeBZZNNTest");
    else if(topology == 4) theTree = (TTree*)inputBWW->Get("TreeBWWNNTest");
    
    
    theTree->SetBranchAddress( "NN1Output", &NN1Output );
    theTree->SetBranchAddress( "NN2Output", &NN2Output );
    theTree->SetBranchAddress( "NN3Output", &NN3Output );
    theTree->SetBranchAddress( "NN4Output", &NN4Output );

    
    
    
    // Efficiency calculator for cut method
    Int_t    nSelCutsGA = 0;
    Double_t effS       = 0.7;
 
    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
 
    if(topology == 0) std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
    else if(topology == 1) std::cout << "--- Processing: " << theTree->GetEntries() << " qq background events" << std::endl;
    else if(topology == 2) std::cout << "--- Processing: " << theTree->GetEntries() << " ttbar background events" << std::endl;
    else if(topology == 3) std::cout << "--- Processing: " << theTree->GetEntries() << " ZZ background events" << std::endl;
    else if(topology == 4) std::cout << "--- Processing: " << theTree->GetEntries() << " WW background events" << std::endl;
    
    TStopwatch sw;
    sw.Start();
   
	    // Book the MVA methods
	    TString dir;
	    dir    = "analysis/datasetOutputNNqqttbarZZWW"+rtdCut+varVersion+preselection+sampleName+"/weights/";
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
	       
	       TString methodName = inputMethod + TString(" method");
	       
	        sum+=reader->EvaluateMVA( methodName );
	       
	       ///Printing for qualitative feeling
	       if(topology == 0 && ievt == 0) cout<<"POS 0: "<<reader->EvaluateMVA( methodName )<<endl;
	       if(topology == 0 && ievt == 70) cout<<"POS 70: "<<reader->EvaluateMVA( methodName )<<endl;
	       
	       ////Filling out arrays of NN outputs
	      NNOutput.push_back(reader->EvaluateMVA(methodName)); 
	       
	       /////Filling out histograms 
	       histOutput.Fill(reader->EvaluateMVA(methodName));
	       
	       sizeTree = ievt+1;
	       
	    }
    
    
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();
    
    cout<<endl<<endl<<"Average of total: "<<sum/(theTree->GetEntries())<<endl<<endl;
 
    // Get efficiency for cuts classifier
    if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                 << " (for a required signal efficiency of " << effS << ")" << std::endl;
 
    if (Use["CutsGA"]) {
 
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
    }
 
    // Write histograms
 
    TFile *target  = new TFile( "analysis/TMVAppOutputNN.root","RECREATE" );
    
    target->Close();
 
    //std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
 
    delete reader;
 
    std::cout << "==> TMVAClassificationApplicationOutputNN is done!" << std::endl << std::endl;
    
 	
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
   		histBDTHHClone->Draw("HIST");
   		histBDTBackClone->Draw("HIST same");
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
 void findSignificance(double bottomHistLimit, double topHistLimit, double nbin, vector<double> NNOutput, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, double& defCut, double& maxSignificance, double weightHH, double weightBackqq, double weightBackttbar, double weightBackZZ, double& weightBackWW, TH2F& histROC, TH2F& histROCRej, TH2F& histSignificance, double& totalRemaining, double& HHRemaining, double& backRemaining, double& backqqRemaining, double& backttbarRemaining, double& backZZRemaining, double& backWWRemaining)
 {
 	defCut=bottomHistLimit;
	double significance; 	
  	 
	double size = NNOutput.size();
	
	//cout<<endl<<endl<<endl<<endl<<endl<<"NNOutputSize: "<<size<<endl<<endl<<endl<<endl;
	
	
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
	topLimitBack = size;
	vector<double> NNOutputBackWW(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	
	
	cout<<"HH size in significance function: "<<NNOutputHH.size()<<endl<<"qq size in significance function: "<<NNOutputBackqq.size()<<endl<<"ttbar size in significance function: "<<NNOutputBackttbar.size()<<endl<<"ZZ size in significance function: "<<NNOutputBackZZ.size()<<endl<<"WW size in significance function: "<<NNOutputBackWW.size()<<endl<<endl;
	//cout<<"bottomLimitBack: "<<bottomLimitBack<<"      topLimitBack: "<<topLimitBack<<endl<<endl<<endl<<endl;
	//cout<<"No crashea en la creacion del vector de back ouputs"<<endl<<endl<<endl;
	sort(NNOutputHH.begin(), NNOutputHH.end());
	sort(NNOutputBackqq.begin(), NNOutputBackqq.end());
	sort(NNOutputBackttbar.begin(), NNOutputBackttbar.end());
	sort(NNOutputBackZZ.begin(), NNOutputBackZZ.end());
	sort(NNOutputBackWW.begin(), NNOutputBackWW.end());
	
	int sizeBackqq = NNOutputBackqq.size();
	int sizeBackttbar = NNOutputBackttbar.size();
	int sizeBackZZ = NNOutputBackZZ.size();
	int sizeBackWW = NNOutputBackWW.size();
	int sizeBack = sizeqq + sizettbar + sizeZZ + sizeWW;
	
	cout<<"sizeBack: "<<sizeBack<<endl<<"sizeBackqq: "<<sizeBackqq<<endl<<endl;
	
	//cout<<"Size: "<<size<<endl<<"SizeHH: "<<NNOutputHH.size()<<endl<<"SizeBack: "<<NNOutputBack.size()<<endl<<endl;
	/*for(int i=0; i<100; i++) cout<<NNOutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<100; i++) cout<<NNOutputBack[i]<<", ";
	cout<<endl;*/
	
	int indexHH=0, indexBackqq=0, indexBackttbar=0, indexBackZZ=0, indexBackWW=0;
	double elementNHH=NNOutputHH[indexHH],elementN1HH, elementNBackqq=NNOutputBackqq[indexBackqq], elementN1Backqq;
	double elementNBackttbar=NNOutputBackttbar[indexBackttbar], elementN1Backttbar, elementNBackZZ=NNOutputBackZZ[indexBackZZ], elementN1BackZZ, elementNBackWW=NNOutputBackWW[indexBackWW], elementN1BackWW;
	double eventsRemainingHH=0, eventsCutHH=0, eventsRemainingBack=0, eventsCutBack=0, eventsRemainingBackqq=0, eventsCutBackqq=0, eventsRemainingBackttbar=0, eventsCutBackttbar=0, eventsRemainingBackZZ=0, eventsCutBackZZ=0, eventsRemainingBackWW=0, eventsCutBackWW=0;
	bool flagFraction=false;
	//for(double cut=bottomHistLimit; cut<=topHistLimit; cut+=0.000001)
	for(double cut=bottomHistLimit; cut<=topHistLimit; cut+=0.0001)
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
		
		eventsRemainingHH = sizeHH-indexHH;
		eventsRemainingBackqq = sizeBackqq-indexBackqq;
		eventsRemainingBackttbar = sizeBackttbar-indexBackttbar;
		eventsRemainingBackZZ = sizeBackZZ-indexBackZZ;
		eventsRemainingBackWW = sizeBackWW-indexBackWW;
		eventsRemainingBack = eventsRemainingBackqq + eventsRemainingBackttbar + eventsRemainingBackZZ + eventsRemainingBackWW;
		
		double fractionHH = eventsRemainingHH/sizeHH;
		double fractionBack = eventsRemainingBack/sizeBack;
		histROCRej.Fill(fractionHH, 1/(fractionBack));
		histROC.Fill(eventsRemainingHH/sizeHH, eventsRemainingBack/sizeBack);
		
		significance=(eventsRemainingHH*weightHH)/(sqrt((eventsRemainingHH*weightHH)+(eventsRemainingBackqq*weightBackqq+eventsRemainingBackttbar*weightBackttbar+eventsRemainingBackZZ*weightBackZZ+eventsRemainingBackWW*weightBackWW)));
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
	   		backRemaining=eventsRemainingBackqq*weightBackqq+eventsRemainingBackttbar*weightBackttbar+eventsRemainingBackZZ*weightBackZZ+eventsRemainingBackWW*weightBackWW;
	   		backqqRemaining=eventsRemainingBackqq*weightBackqq;
	   		backttbarRemaining=eventsRemainingBackttbar*weightBackttbar;
	   		backZZRemaining=eventsRemainingBackZZ*weightBackZZ;
	   		backWWRemaining=eventsRemainingBackWW*weightBackWW;
	   		HHRemaining=eventsRemainingHH*weightHH;
	   		totalRemaining=eventsRemainingHH+eventsRemainingBack;
	   	}
	}
	cout<<"HHRemainingUnweighted: "<<HHRemaining/weightHH<<endl<<"backqqRemainingUnweighted: "<<backqqRemaining/weightBackqq<<endl<<"backttbarRemainingUnweighted: "<<backttbarRemaining/weightBackttbar<<endl<<"backZZRemainingUnweighted: "<<backZZRemaining/weightBackZZ<<endl<<"backWWRemainingUnweighted: "<<backWWRemaining/weightBackWW<<endl<<endl; 
	cout<<"HHRemaining: "<<HHRemaining<<endl<<"backRemaining: "<<backRemaining<<endl<<"backqqRemaining: "<<backqqRemaining<<endl<<"backttbarRemaining: "<<backttbarRemaining<<endl<<"backZZRemaining: "<<backZZRemaining<<endl<<"backWWRemaining: "<<backWWRemaining<<endl<<endl; 
	cout<<"Significance from manual formula: "<<HHRemaining/sqrt(HHRemaining+backRemaining)<<endl;
	
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

void checkSigmaOGStyle(double bottomHistLimit, double topHistLimit, int nBack, double weightHH, double weightBack)
 {
	cout<<endl<<"From checkSIgmaOGStyle: "<<endl;
	double signalPass=0, signalCut=0, backPass=0, backCut=0;
	vector<double> NNOutput;
	double totalRemaining=0, HHRemaining=0, backRemaining=0, defCut=0, maxSignificance=0;
	
	TFile *fileHH = TFile::Open("outputTreeSNNTest.root");
	TFile *fileqq = TFile::Open("outputTreeBqqNNTest.root");
	
	TTree *treeHH, *treeqq;
	fileHH->GetObject("TreeSNNTest", treeHH);
	fileqq->GetObject("TreeBqqNNTest", treeqq);
	
	float NN1Output;
	treeHH->SetBranchAddress("NN1Output", &NN1Output);
	
	 Long64_t sizeHH = treeHH->GetEntries();
	 for (Long64_t i = 0; i < sizeHH; ++i) {
	 	treeHH->GetEntry(i);
	 	NNOutput.push_back(NN1Output);
	 }
	 treeqq->SetBranchAddress("NN1Output", &NN1Output);
	 Long64_t sizeqq = treeqq->GetEntries();
	 for (Long64_t i = 0; i < sizeqq; ++i) {
	 	treeqq->GetEntry(i);
	 	NNOutput.push_back(NN1Output);
	 }
	
	double significance; 	
	
	int sizettbar, sizeZZ, sizeWW;
  	 
	double size = NNOutput.size();
	
	cout<<endl<<endl<<endl<<endl<<endl<<"NNOutputSize: "<<size<<endl<<endl<<endl<<endl;
	
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
		topLimitBack = size;
	}
	//cout<<"bottomLimitBack: "<<bottomLimitBack<<"      topLimitBack: "<<topLimitBack<<endl<<endl<<endl<<endl;
	vector<double> NNOutputBack(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	//cout<<"No crashea en la creacion del vector de back ouputs"<<endl<<endl<<endl;
	sort(NNOutputHH.begin(), NNOutputHH.end());
	sort(NNOutputBack.begin(), NNOutputBack.end());
	
	int sizeBack = NNOutputBack.size();
	
	cout<<"sizeBack: "<<sizeBack<<endl<<endl<<endl;
	
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
	for(double cut=bottomHistLimit; cut<=topHistLimit; cut+=0.0001)
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
		
		significance=(eventsRemainingHH*weightHH)/(sqrt((eventsRemainingHH*weightHH)+(eventsRemainingBack*weightBack)));
		/*cout<<"CUT: "<<cut<<endl<<endl;
		cout<<"IndexHH: "<<indexHH<<"    elem: "<<elementNHH<<endl<<"Signal remaining: "<<eventsRemainingHH<<"    signal cut: "<<sizeHH-eventsRemainingHH<<endl<<endl;
		cout<<"IndexBack: "<<indexBack<<"    elem: "<<elementNBack<<endl<<"Back remaining: "<<eventsRemainingBack<<"    back cut: "<<sizeBack-eventsRemainingBack<<endl;*/
		//cout<<"Significance: "<<significance<<"         maxSignificnace: "<<maxSignificance<<endl<<endl;
		//cout<<"--------"<<endl;
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
	cout<<endl<<"Max significance: "<<maxSignificance<<endl<<endl;
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
 
 void FSRTMVAClassificationApplicationOutputNN(string rtdCut, string preselection, string varVersion, string sampleName)
 {
 	cout<<"Aquì empieza main()"<<endl;
 	gStyle->SetOptStat(0);
 	
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
    	
    	string inputTrainHHNNText = "analysis/outputTreeSNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
    	string inputTestHHNNText = "analysis/outputTreeSNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    	TFile* inputTrainHHNN = new TFile(inputTrainHHNNText.c_str());
    	TTree* theTreeTrainHHNN = (TTree*)inputTrainHHNN->Get("TreeSNNTrain");
 	TFile* inputTestHHNN = new TFile(inputTestHHNNText.c_str());
    	TTree* theTreeTestHHNN = (TTree*)inputTestHHNN->Get("TreeSNNTest");
    	double weightFactorHHNN = (theTreeTestHHNN->GetEntries()+theTreeTrainHHNN->GetEntries());
    	weightFactorHHNN = weightFactorHHNN/(theTreeTestHHNN->GetEntries());
    	string inputTrainqqNNText = "analysis/outputTreeBqqNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
    	string inputTestqqNNText = "analysis/outputTreeBqqNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    	TFile* inputTrainqqNN = new TFile(inputTrainqqNNText.c_str());
    	TTree* theTreeTrainqqNN = (TTree*)inputTrainqqNN->Get("TreeBqqNNTrain");
 	TFile* inputTestqqNN = new TFile(inputTestqqNNText.c_str());
    	TTree* theTreeTestqqNN = (TTree*)inputTestqqNN->Get("TreeBqqNNTest");
    	double weightFactorqqNN = (theTreeTestqqNN->GetEntries()+theTreeTrainqqNN->GetEntries());
    	weightFactorqqNN = weightFactorqqNN/(theTreeTestqqNN->GetEntries());
    	string inputTrainttbarNNText = "analysis/outputTreeBttbarNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
    	string inputTestttbarNNText = "analysis/outputTreeBttbarNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    	TFile* inputTrainttbarNN = new TFile(inputTrainttbarNNText.c_str());
    	TTree* theTreeTrainttbarNN = (TTree*)inputTrainttbarNN->Get("TreeBttbarNNTrain");
 	TFile* inputTestttbarNN = new TFile(inputTestttbarNNText.c_str());
    	TTree* theTreeTestttbarNN = (TTree*)inputTestttbarNN->Get("TreeBttbarNNTest");
    	double weightFactorttbarNN = (theTreeTestttbarNN->GetEntries()+theTreeTrainttbarNN->GetEntries());
    	weightFactorttbarNN = weightFactorttbarNN/(theTreeTestttbarNN->GetEntries());
    	string inputTrainZZNNText = "analysis/outputTreeBZZNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
    	string inputTestZZNNText = "analysis/outputTreeBZZNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    	TFile* inputTrainZZNN = new TFile(inputTrainZZNNText.c_str());
    	TTree* theTreeTrainZZNN = (TTree*)inputTrainZZNN->Get("TreeBZZNNTrain");
 	TFile* inputTestZZNN = new TFile(inputTestZZNNText.c_str());
    	TTree* theTreeTestZZNN = (TTree*)inputTestZZNN->Get("TreeBZZNNTest");
    	double weightFactorZZNN = (theTreeTestZZNN->GetEntries()+theTreeTrainZZNN->GetEntries());
    	weightFactorZZNN = weightFactorZZNN/(theTreeTestZZNN->GetEntries());
    	string inputTrainWWNNText = "analysis/outputTreeBWWNNESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
    	string inputTestWWNNText = "analysis/outputTreeBWWNNESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
    	TFile* inputTrainWWNN = new TFile(inputTrainWWNNText.c_str());
    	TTree* theTreeTrainWWNN = (TTree*)inputTrainWWNN->Get("TreeBWWNNTrain");
 	TFile* inputTestWWNN = new TFile(inputTestWWNNText.c_str());
    	TTree* theTreeTestWWNN = (TTree*)inputTestWWNN->Get("TreeBWWNNTest");
    	double weightFactorWWNN = (theTreeTestWWNN->GetEntries()+theTreeTrainWWNN->GetEntries());
    	weightFactorWWNN = weightFactorWWNN/(theTreeTestWWNN->GetEntries());
    	
    	cout<<"WeightFactorHH: "<<weightFactorHH<<endl<<"WeightFactorHHNN: "<<weightFactorHHNN<<endl;
    	cout<<"WeightFactorqq: "<<weightFactorqq<<endl<<"WeightFactorqqNN: "<<weightFactorqqNN<<endl;
    	cout<<"WeightFactorttbar: "<<weightFactorttbar<<endl<<"WeightFactorttbarNN: "<<weightFactorttbarNN<<endl;
    	cout<<"WeightFactorZZ: "<<weightFactorZZ<<endl<<"WeightFactorZZNN: "<<weightFactorZZNN<<endl;
    	cout<<"WeightFactorWW: "<<weightFactorWW<<endl<<"WeightFactorWWNN: "<<weightFactorWWNN<<endl;
    	
    	//double weightFactorHH=1, weightFactorqq=1, weightFactorttbar=1, weightFactorZZ=1, weightFactorWW=1; 
    	//cout<<"Issue before"<<endl;
 	
 	string method = "BDTG";
 	string variables = "qqttbarZZWW";
 	double bottomHistLimit=0.0, topHistLimit=0.0;
 	int methodColor;
 	findHistLimits(method, bottomHistLimit, topHistLimit, methodColor);
 
 	int sizeHH=0, sizeqq=0, sizettbar=0, sizeZZ=0, sizeWW=0;
 	int nbin=10000, nbinNN=100;
 	double weightHH=0.001225*weightFactorHH*weightFactorHHNN, weightqq=0.0349*weightFactorqq*weightFactorqqNN, weightttbar=0.503*weightFactorttbar*weightFactorttbarNN, weightZZ=0.8167*weightFactorZZ*weightFactorZZNN, weightWW=0.5149*weightFactorWW*weightFactorWWNN;
 	vector<double> NNOutput;
 	
 	
 	TH1F *histNNHH     = new TH1F( "NNHH",           "NNoutput on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histNNB     = new TH1F( "NNBZZ",            "NNoutput on all back (qq+ttbar+ZZ)",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histNNBqq     = new TH1F( "NNBqq",            "NNoutput on qq",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histNNBttbar     = new TH1F( "NNBttbar",      "NNoutput on ttbar",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histNNBZZ     = new TH1F( "NNBZZ",            "NNoutput on ZZ",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histNNBWW     = new TH1F( "NNBWW",            "NNoutput on WW",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	
 	histNNHH->SetLineColor(kBlue);
 	histNNB->SetLineColor(kRed+2);
 	histNNBqq->SetLineColor(kRed+2);
 	histNNBttbar->SetLineColor(kRed+2);
 	histNNBZZ->SetLineColor(kRed+2);
 	histNNBWW->SetLineColor(kRed+2);
 	
 	histNNHH->SetFillColor(kBlue);
 	histNNB->SetFillColor(kRed+2);
 	histNNBqq->SetFillColor(kRed+2);
 	histNNBttbar->SetFillColor(kRed+2);
 	histNNBZZ->SetFillColor(kRed+2);
 	histNNBWW->SetFillColor(kRed+2);
 	
 	/////For arguments: empty string, vector for output of NN on all events, size of sample, number of bins, hist for output, string with type of NN, string with vars in the dataset, topology (i.e. 0 for HH, 1 for qq, 2 for ttbar, 3 for ZZ).
 	//////////
 	FSRTMVAClassificationApplicationOutputNNHelper("", NNOutput, sizeHH, nbinNN, *histNNHH, method, variables, 0, rtdCut, preselection, varVersion, sampleName);
 	FSRTMVAClassificationApplicationOutputNNHelper("", NNOutput, sizeqq, nbinNN, *histNNBqq, method, variables, 1, rtdCut, preselection, varVersion, sampleName);
 	FSRTMVAClassificationApplicationOutputNNHelper("", NNOutput, sizettbar, nbinNN, *histNNBttbar, method, variables, 2, rtdCut, preselection, varVersion, sampleName);
 	FSRTMVAClassificationApplicationOutputNNHelper("", NNOutput, sizeZZ, nbinNN, *histNNBZZ, method, variables, 3, rtdCut, preselection, varVersion, sampleName);
 	FSRTMVAClassificationApplicationOutputNNHelper("", NNOutput, sizeWW, nbinNN, *histNNBWW, method, variables, 4, rtdCut, preselection, varVersion, sampleName);
 	//////////////
 	
 	cout<<"NNOutput size: "<<NNOutput.size()<<endl;
 	cout<<"sizeHH: "<<sizeHH<<endl<<"sizeqq: "<<sizeqq<<endl<<"sizettbar: "<<sizettbar<<endl<<"sizeZZ: "<<sizeZZ<<endl<<"sizeWW: "<<sizeWW<<endl;
 	
 	
 	///////////////////////
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histNNHH, kBlue, *histNNB, kRed+2, method + "NN of NN outputs on HH and all backs (qq+ttbar+ZZ) (unweighted)");
 	/////////////////////////
 	
 	///////////////////
 	double maxSignificance=0, defCut, totalRemaining=0, HHRemaining=0, backqqRemaining=0, backttbarRemaining=0, backZZRemaining=0, backWWRemaining=0, backRemaining=0;
 	
 	TH2F *histROC = new TH2F("hROC", "Signal and Background Acceptance (NN of NN outputs)", nbin, 0.0, 1.0, nbin, 0, 1.0);
  	TH2F *histROCRej = new TH2F("hROCRej", "Signal Acceptance and Background Rejection (NN of NN outputs)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	
  	TH2F *histSignificance     = new TH2F( "hist_significance", "Significance s/sqrt(s+b) (NN of NN outputs)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificance->GetXaxis()->SetTitle("Cut value for NN of NN outputs");
   	histSignificance->GetYaxis()->SetTitle("Significance for NN of NN outputs");
   	
   	///////////////TERMINAR SIGNIFICANCE
   	
   	cout<<endl<<"Significance for "<<method<<" of NN outputs applied to HH and qq"<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, NNOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, defCut, maxSignificance, weightHH, weightqq, weightttbar, weightZZ, weightWW, *histROC, *histROCRej, *histSignificance, totalRemaining, HHRemaining, backRemaining, backqqRemaining, backttbarRemaining, backZZRemaining, backWWRemaining);
 	cout<<"maxSignificance: "<<maxSignificance<<endl<<"Cut in NN of NN outputs for max significance: "<<defCut<<endl;
 	double finalSignificance = maxSignificance;
 	///////////////////////

 	

 	
 	
 	
 	
 	
 	
 	double luminosity = 4900, crossSectionqq, errorTopqq, errorBottomqq;
 	double crossSectionttbar, errorTopttbar, errorBottomttbar;
 	double crossSectionZZ, errorTopZZ, errorBottomZZ;
 	double crossSectionCombined, errorTopCombined, errorBottomCombined;
 	
 	cout<<endl<<endl;
 	
 	/*cout<<"Cross-section and error for the output of "<<method<<"qq applied only to HH and qq: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingqq+backRemainingqq, HHRemainingqq, backRemainingqq, luminosity, crossSectionqq, errorTopqq, errorBottomqq, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"ttbar applied only to HH and ttbar: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingttbar+backRemainingttbar, HHRemainingttbar, backRemainingttbar, luminosity, crossSectionttbar, errorTopttbar, errorBottomttbar, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"ZZ applied only to HH and ZZ: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingZZ+backRemainingZZ, HHRemainingZZ, backRemainingZZ, luminosity, crossSectionttbar, errorTopZZ, errorBottomZZ, nbin);*/
 	cout<<"Cross-section and error for the output of "<<method<<" combined applied to HH, qq, ttbar, ZZ, and WW: "<<endl;
 	cout<<"HHRemaining: "<<HHRemaining<<endl<<"backRemaining: "<<backRemaining<<endl;
 	findCrossSectionHHbbbb(HHRemaining+backRemaining, HHRemaining, backRemaining, luminosity, crossSectionCombined, errorTopCombined, errorBottomCombined, nbin);
 	double finalErrorLeft = errorBottomCombined;
 	double finalErrorRight = errorTopCombined;
 	
 	/////saving significance and errors in txt
 	ofstream outFile("analysis/significanceAndErrorsFile.txt");
	if (outFile.is_open()) {
	    outFile << finalSignificance << " " << finalErrorLeft << " " << finalErrorRight << endl;
	    outFile.close();
	} else {
	    cerr << "Unable to open file for writing significance and errors" << endl;
	}
	///////saving significance and errors in txt
 	
 	
 	
 	cout<<endl<<endl;
 	
 	//checkSigmaOGStyle(bottomHistLimit, topHistLimit, 0, weightHH, weightqq);
 	
 	
 	
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
    FSRTMVAClassificationApplicationOutputNN(methodList);
    return 0;
 }*/

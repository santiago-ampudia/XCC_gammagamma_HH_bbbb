#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <algorithm>

using namespace std;

void createEntryIndexFiles(TString nthSample)
{
	 TFile *inputFileSTrain = TFile::Open("analysis/outputTreeSHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileSTest = TFile::Open("analysis/outputTreeSHHbbbbESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileSNNTrain = TFile::Open("analysis/outputTreeSNNESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileSNNTest = TFile::Open("analysis/outputTreeSNNESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileBqqTrain = TFile::Open("analysis/outputTreeBqqHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileBqqTest = TFile::Open("analysis/outputTreeBqqHHbbbbESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileBqqNNTrain = TFile::Open("analysis/outputTreeBqqNNESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileBqqNNTest = TFile::Open("analysis/outputTreeBqqNNESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileBttTrain = TFile::Open("analysis/outputTreeBttHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileBttTest = TFile::Open("analysis/outputTreeBttHHbbbbESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileBttNNTrain = TFile::Open("analysis/outputTreeBttbarNNESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileBttNNTest = TFile::Open("analysis/outputTreeBttbarNNESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileBZZTrain = TFile::Open("analysis/outputTreeBZZHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileBZZTest = TFile::Open("analysis/outputTreeBZZHHbbbbESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileBZZNNTrain = TFile::Open("analysis/outputTreeBZZNNESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileBZZNNTest = TFile::Open("analysis/outputTreeBZZNNESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileBWWTrain = TFile::Open("analysis/outputTreeBWWHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileBWWTest = TFile::Open("analysis/outputTreeBWWHHbbbbESpreadDurham1034BSplitTestSampleN.root");
	 TFile *inputFileBWWNNTrain = TFile::Open("analysis/outputTreeBWWNNESpreadDurham1034BSplitTrainSampleN.root");
	 TFile *inputFileBWWNNTest = TFile::Open("analysis/outputTreeBWWNNESpreadDurham1034BSplitTestSampleN.root");
	 
	 TTree *originalTreeSTrain = nullptr;
    	 inputFileSTrain->GetObject("TreeSTrain", originalTreeSTrain);
    	 TTree *originalTreeSTest = nullptr;
    	 inputFileSTest->GetObject("TreeSTest", originalTreeSTest);
    	 TTree *originalTreeSNNTrain = nullptr;
    	 inputFileSNNTrain->GetObject("TreeSNNTrain", originalTreeSNNTrain);
    	 TTree *originalTreeSNNTest = nullptr;
    	 inputFileSNNTest->GetObject("TreeSNNTest", originalTreeSNNTest);
    	 TTree *originalTreeBqqTrain = nullptr;
    	 inputFileBqqTrain->GetObject("TreeBqqTrain", originalTreeBqqTrain);
    	 TTree *originalTreeBqqTest = nullptr;
    	 inputFileBqqTest->GetObject("TreeBqqTest", originalTreeBqqTest);
    	 TTree *originalTreeBqqNNTrain = nullptr;
    	 inputFileBqqNNTrain->GetObject("TreeBqqNNTrain", originalTreeBqqNNTrain);
    	 TTree *originalTreeBqqNNTest = nullptr;
    	 inputFileBqqNNTest->GetObject("TreeBqqNNTest", originalTreeBqqNNTest);
    	 TTree *originalTreeBttTrain = nullptr;
    	 inputFileBttTrain->GetObject("TreeBttTrain", originalTreeBttTrain);
    	 TTree *originalTreeBttTest = nullptr;
    	 inputFileBttTest->GetObject("TreeBttTest", originalTreeBttTest);
    	 TTree *originalTreeBttNNTrain = nullptr;
    	 inputFileBttNNTrain->GetObject("TreeBttbarNNTrain", originalTreeBttNNTrain);
    	 TTree *originalTreeBttNNTest = nullptr;
    	 inputFileBttNNTest->GetObject("TreeBttbarNNTest", originalTreeBttNNTest);
    	 TTree *originalTreeBZZTrain = nullptr;
    	 inputFileBZZTrain->GetObject("TreeBZZTrain", originalTreeBZZTrain);
    	 TTree *originalTreeBZZTest = nullptr;
    	 inputFileBZZTest->GetObject("TreeBZZTest", originalTreeBZZTest);
    	 TTree *originalTreeBZZNNTrain = nullptr;
    	 inputFileBZZNNTrain->GetObject("TreeBZZNNTrain", originalTreeBZZNNTrain);
    	 TTree *originalTreeBZZNNTest = nullptr;
    	 inputFileBZZNNTest->GetObject("TreeBZZNNTest", originalTreeBZZNNTest);
    	 TTree *originalTreeBWWTrain = nullptr;
    	 inputFileBWWTrain->GetObject("TreeBWWTrain", originalTreeBWWTrain);
    	 TTree *originalTreeBWWTest = nullptr;
    	 inputFileBWWTest->GetObject("TreeBWWTest", originalTreeBWWTest);
    	 TTree *originalTreeBWWNNTrain = nullptr;
    	 inputFileBWWNNTrain->GetObject("TreeBWWNNTrain", originalTreeBWWNNTrain);
    	 TTree *originalTreeBWWNNTest = nullptr;
    	 inputFileBWWNNTest->GetObject("TreeBWWNNTest", originalTreeBWWNNTest);
    	 
    	 TString sampleNameEntryIndexFile = "Sample" + TString::Format("%d", nthSample);
    	 TString outputFileText = "analysis/entryIndexFile"+sampleNameEntryIndexFile+".root";
    	 TFile *outputFile = new TFile(outputFileText.Data(), "RECREATE");
    	 
    	 TTree *TreeEntryIndex = new TTree("TreeEntryIndex", "Tree with only entryIndex branch"); // Give a new name to the new tree
    	 float entryIndexSTrain, entryIndexSTest, entryIndexSNNTrain, entryIndexSNNTest, entryIndexqqTrain, entryIndexqqTest, entryIndexqqNNTrain, entryIndexqqNNTest, entryIndexttTrain, entryIndexttTest, entryIndexttNNTrain, entryIndexttNNTest, entryIndexZZTrain, entryIndexZZTest, entryIndexZZNNTrain, entryIndexZZNNTest, entryIndexWWTrain, entryIndexWWTest, entryIndexWWNNTrain, entryIndexWWNNTest;
    	 
    	 originalTreeSTrain->SetBranchAddress("entryIndex", &entryIndexSTrain);
    	 originalTreeSTest->SetBranchAddress("entryIndex", &entryIndexSTest);
    	 originalTreeSNNTrain->SetBranchAddress("entryIndex", &entryIndexSNNTrain);
    	 originalTreeSNNTest->SetBranchAddress("entryIndex", &entryIndexSNNTest);
    	 originalTreeBqqTrain->SetBranchAddress("entryIndex", &entryIndexqqTrain);
    	 originalTreeBqqTest->SetBranchAddress("entryIndex", &entryIndexqqTest);
    	 originalTreeBqqNNTrain->SetBranchAddress("entryIndex", &entryIndexqqNNTrain);
    	 originalTreeBqqNNTest->SetBranchAddress("entryIndex", &entryIndexqqNNTest);
    	 originalTreeBttTrain->SetBranchAddress("entryIndex", &entryIndexttTrain);
    	 originalTreeBttTest->SetBranchAddress("entryIndex", &entryIndexttTest);
    	 originalTreeBttNNTrain->SetBranchAddress("entryIndex", &entryIndexttNNTrain);
    	 originalTreeBttNNTest->SetBranchAddress("entryIndex", &entryIndexttNNTest);
    	 originalTreeBZZTrain->SetBranchAddress("entryIndex", &entryIndexZZTrain);
    	 originalTreeBZZTest->SetBranchAddress("entryIndex", &entryIndexZZTest);
    	 originalTreeBZZNNTrain->SetBranchAddress("entryIndex", &entryIndexZZNNTrain);
    	 originalTreeBZZNNTest->SetBranchAddress("entryIndex", &entryIndexZZNNTest);
    	 originalTreeBWWTrain->SetBranchAddress("entryIndex", &entryIndexWWTrain);
    	 originalTreeBWWTest->SetBranchAddress("entryIndex", &entryIndexWWTest);
    	 originalTreeBWWNNTrain->SetBranchAddress("entryIndex", &entryIndexWWNNTrain);
    	 originalTreeBWWNNTest->SetBranchAddress("entryIndex", &entryIndexWWNNTest);
    	 
    	 TreeEntryIndex->Branch("entryIndexSTrain", &entryIndexSTrain);
    	 TreeEntryIndex->Branch("entryIndexSTest", &entryIndexSTest);
    	 TreeEntryIndex->Branch("entryIndexSNNTrain", &entryIndexSNNTrain);
    	 TreeEntryIndex->Branch("entryIndexSNNTest", &entryIndexSNNTest);
    	 TreeEntryIndex->Branch("entryIndexqqTrain", &entryIndexqqTrain);
	 TreeEntryIndex->Branch("entryIndexqqTest", &entryIndexqqTest);
	 TreeEntryIndex->Branch("entryIndexqqNNTrain", &entryIndexqqNNTrain);
	 TreeEntryIndex->Branch("entryIndexqqNNTest", &entryIndexqqNNTest);
	 TreeEntryIndex->Branch("entryIndexttTrain", &entryIndexttTrain);
	 TreeEntryIndex->Branch("entryIndexttTest", &entryIndexttTest);
	 TreeEntryIndex->Branch("entryIndexttNNTrain", &entryIndexttNNTrain);
	 TreeEntryIndex->Branch("entryIndexttNNTest", &entryIndexttNNTest);
	 TreeEntryIndex->Branch("entryIndexZZTrain", &entryIndexZZTrain);
	 TreeEntryIndex->Branch("entryIndexZZTest", &entryIndexZZTest);
	 TreeEntryIndex->Branch("entryIndexZZNNTrain", &entryIndexZZNNTrain);
	 TreeEntryIndex->Branch("entryIndexZZNNTest", &entryIndexZZNNTest);
	 TreeEntryIndex->Branch("entryIndexWWTrain", &entryIndexWWTrain);
	 TreeEntryIndex->Branch("entryIndexWWTest", &entryIndexWWTest);
	 TreeEntryIndex->Branch("entryIndexWWNNTrain", &entryIndexWWNNTrain);
	 TreeEntryIndex->Branch("entryIndexWWNNTest", &entryIndexWWNNTest);
	 
	 Long64_t nEntriesSTrain = originalTreeSTrain->GetEntries();
    	Long64_t nEntriesSTest = originalTreeSTest->GetEntries();
	 Long64_t nEntriesSNNTrain = originalTreeSNNTrain->GetEntries();
	    Long64_t nEntriesSNNTest = originalTreeSNNTest->GetEntries();
	    Long64_t nEntriesBqqTrain = originalTreeBqqTrain->GetEntries();
	    Long64_t nEntriesBqqTest = originalTreeBqqTest->GetEntries();
	    Long64_t nEntriesBqqNNTrain = originalTreeBqqNNTrain->GetEntries();
	    Long64_t nEntriesBqqNNTest = originalTreeBqqNNTest->GetEntries();
	    Long64_t nEntriesBttTrain = originalTreeBttTrain->GetEntries();
	    Long64_t nEntriesBttTest = originalTreeBttTest->GetEntries();
	    Long64_t nEntriesBttNNTrain = originalTreeBttNNTrain->GetEntries();
	    Long64_t nEntriesBttNNTest = originalTreeBttNNTest->GetEntries();
	    Long64_t nEntriesBZZTrain = originalTreeBZZTrain->GetEntries();
	    Long64_t nEntriesBZZTest = originalTreeBZZTest->GetEntries();
	    Long64_t nEntriesBZZNNTrain = originalTreeBZZNNTrain->GetEntries();
	    Long64_t nEntriesBZZNNTest = originalTreeBZZNNTest->GetEntries();
	    Long64_t nEntriesBWWTrain = originalTreeBWWTrain->GetEntries();
	    Long64_t nEntriesBWWTest = originalTreeBWWTest->GetEntries();
	    Long64_t nEntriesBWWNNTrain = originalTreeBWWNNTrain->GetEntries();
	    Long64_t nEntriesBWWNNTest = originalTreeBWWNNTest->GetEntries();

	    for (Long64_t i = 0; i < max({nEntriesSTrain, nEntriesSTest, nEntriesSNNTrain, nEntriesSNNTest,
		                          nEntriesBqqTrain, nEntriesBqqTest, nEntriesBqqNNTrain, nEntriesBqqNNTest,
		                          nEntriesBttTrain, nEntriesBttTest, nEntriesBttNNTrain, nEntriesBttNNTest,
		                          nEntriesBZZTrain, nEntriesBZZTest, nEntriesBZZNNTrain, nEntriesBZZNNTest,
		                          nEntriesBWWTrain, nEntriesBWWTest, nEntriesBWWNNTrain, nEntriesBWWNNTest}); ++i) {
		if (i < nEntriesSTrain) originalTreeSTrain->GetEntry(i);
		if (i < nEntriesSTest) originalTreeSTest->GetEntry(i);
		if (i < nEntriesSNNTrain) originalTreeSNNTrain->GetEntry(i);
		if (i < nEntriesSNNTest) originalTreeSNNTest->GetEntry(i);
		if (i < nEntriesBqqTrain) originalTreeBqqTrain->GetEntry(i);
		if (i < nEntriesBqqTest) originalTreeBqqTest->GetEntry(i);
		if (i < nEntriesBqqNNTrain) originalTreeBqqNNTrain->GetEntry(i);
		if (i < nEntriesBqqNNTest) originalTreeBqqNNTest->GetEntry(i);
		if (i < nEntriesBttTrain) originalTreeBttTrain->GetEntry(i);
		if (i < nEntriesBttTest) originalTreeBttTest->GetEntry(i);
		if (i < nEntriesBttNNTrain) originalTreeBttNNTrain->GetEntry(i);
		if (i < nEntriesBttNNTest) originalTreeBttNNTest->GetEntry(i);
		if (i < nEntriesBZZTrain) originalTreeBZZTrain->GetEntry(i);
		if (i < nEntriesBZZTest) originalTreeBZZTest->GetEntry(i);
		if (i < nEntriesBZZNNTrain) originalTreeBZZNNTrain->GetEntry(i);
		if (i < nEntriesBZZNNTest) originalTreeBZZNNTest->GetEntry(i);
		if (i < nEntriesBWWTrain) originalTreeBWWTrain->GetEntry(i);
		if (i < nEntriesBWWTest) originalTreeBWWTest->GetEntry(i);
		if (i < nEntriesBWWNNTrain) originalTreeBWWNNTrain->GetEntry(i);
		if (i < nEntriesBWWNNTest) originalTreeBWWNNTest->GetEntry(i);
		TreeEntryIndex->Fill();
    }

    outputFile->cd();
    TreeEntryIndex->Write();
    outputFile->Close();
    	 
}

double getMean(const vector<double>& dataVector)
{
	double sum=0;
	for(int i=0; i<dataVector.size(); i++) sum += dataVector[i];
	return sum/(dataVector.size());
}

void getStats(int samples, const vector<double>& preliminarySignificances, const vector<double>& significances, const vector<double>& errorsLeft, const vector<double>& errorsRight, double& meanPreliminarySignificances, double& meanSignificances, double& meanErrorsLeft, double& meanErrorsRight, int& bestSample, int& worstSample, double& bestPreliminarySignificances, double& bestSignificances, double& bestErrorsLeft, double& bestErrorsRight, double& worstPreliminarySignificances, double& worstSignificances, double& worstErrorsLeft, double& worstErrorsRight)
{
	//if(samples != significances.size()) throw std::runtime_error("ERROR: size of signifcances vector != n of samples");
	
	meanPreliminarySignificances = getMean(preliminarySignificances);
	meanSignificances = getMean(significances);
	meanErrorsLeft = getMean(errorsLeft);
	meanErrorsRight = getMean(errorsRight);
	auto max_iter = max_element(significances.begin(), significances.end()); 
	bestSample = distance(significances.begin(), max_iter); 
	bestPreliminarySignificances = preliminarySignificances[bestSample];	
	bestSignificances = significances[bestSample];
	bestErrorsRight = errorsRight[bestSample];
	bestErrorsLeft = errorsLeft[bestSample];
	auto min_iter = min_element(significances.begin(), significances.end()); 
	worstSample = distance(significances.begin(), min_iter);
	worstPreliminarySignificances = preliminarySignificances[worstSample];	
	worstSignificances = significances[worstSample];
	worstErrorsRight = errorsRight[worstSample];
	worstErrorsLeft = errorsLeft[worstSample];

	
	cout<<"For "<<samples<<" samples: "<<endl;
	cout<<"preliminarySignificances: ";
	for(int i=0; i<significances.size(); i++) cout<<preliminarySignificances[i]<<", ";
	cout<<endl<<"significances: ";
	for(int i=0; i<significances.size(); i++) cout<<significances[i]<<", ";
	cout<<endl<<"errorsLeft: ";
	for(int i=0; i<significances.size(); i++) cout<<errorsLeft[i]<<", ";
	cout<<endl<<"errorsRight: ";
	for(int i=0; i<significances.size(); i++) cout<<errorsRight[i]<<", ";
	cout<<endl<<"meanPreliminarySignificances: "<<meanPreliminarySignificances<<endl;
	cout<<"meanSignificances: "<<meanSignificances<<endl;
	cout<<"meanErrorsLeft: "<<meanErrorsLeft<<endl;
	cout<<"meanErrorsRight: "<<meanErrorsRight<<endl;
	cout<<"sample that generates the best final significance: "<<bestSample<<endl;
	cout<<"bestPreliminarySignificances: "<<bestPreliminarySignificances<<endl;
	cout<<"bestSignificances: "<<bestSignificances<<endl;
	cout<<"bestErrorsLeft: "<<bestErrorsLeft<<endl;
	cout<<"bestErrorsRight: "<<bestErrorsRight<<endl;
	cout<<"sample that generates the worst final significance: "<<worstSample<<endl;
	cout<<"worstPreliminarySignificances: "<<worstPreliminarySignificances<<endl;
	cout<<"worstSignificances: "<<worstSignificances<<endl;
	cout<<"worstErrorsLeft: "<<worstErrorsLeft<<endl;
	cout<<"worstErrorsRight: "<<worstErrorsRight<<endl;
}

void trainAllBacks(TString rtdCut, TString preselection, TString vars, TString sampleName)
{
	TString back = "qq";
	gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationHHbbbb.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", back.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));

	back = "ttbar";
	gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationHHbbbb.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", back.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
	back = "ZZ";
	gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationHHbbbb.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", back.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
	back = "WW";
	gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationHHbbbb.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", back.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
}


void runAll(int samples, TString fileFunction, TString preselection, TString vars, TString rtdCut, vector<double>& preliminarySignificances, vector<double>& significances, vector<double>& errorsLeft, vector<double>& errorsRight, TString outFileText)
{
	ofstream outFile(outFileText.Data(), ios::app);
	double prevPreliminarySignificance=-999, prevSignificance=-999, prevErrorRight=-999, prevErrorLeft=-999;
	for(int nthSample=0; nthSample<samples; nthSample++)
    	{
    		TString sampleName = "Sample" + TString::Format("%d", nthSample);
    		sampleName = "SampleN";
    		//gSystem->Exec("root -l -b -q 'analysis/FSRGammaGammaHHbbbbAnalysis.C'");
    		//trainAllBacks(rtdCut, preselection, vars, sampleName);
    		gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationApplicationHHbbbbGeneratesNNs.C+(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", fileFunction.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    		gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationOutputNN.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\")'", rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    		gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationApplicationOutputNN.C+(\"%s\", \"%s\", \"%s\", \"%s\")'", rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    		
    		double preliminarySignificance=-999, significance=-999, errorLeft=-999, errorRight=-999;
    		ifstream file1("analysis/preliminarySignificanceFile.txt");
		file1 >> preliminarySignificance;
		file1.close();
		ifstream file2("analysis/significanceAndErrorsFile.txt");
		file2 >> significance >> errorLeft >> errorRight;
		file2.close();
    		
    		if(preliminarySignificance == prevPreliminarySignificance && significance == prevSignificance && errorLeft == prevErrorLeft && errorRight == prevErrorRight)  outFile << "sample: " << nthSample << endl << "CRASH SOMEWHERE" << endl;
    		else 
    		{
    			preliminarySignificances.push_back(preliminarySignificance);
			significances.push_back(significance);
		    	errorsLeft.push_back(errorLeft);
		    	errorsRight.push_back(errorRight);
    			
    			outFile << "sample: " << nthSample << endl;
			outFile << "preliminary significance: " << preliminarySignificance << endl;
			outFile << "significance: " << significance << endl;
			outFile << "errorsLeft: " << errorLeft << endl;
			outFile << "errorsRight: " << errorRight << endl;
			outFile << endl;
    		} 
    		prevPreliminarySignificance = preliminarySignificance;
    		prevSignificance = significance;
    		prevErrorLeft = errorLeft;
    		prevErrorRight = errorRight; 
    		
    		cout<<"Sample: "<<nthSample;
    		cout<<endl<<"preliminarySignificance: "<<preliminarySignificance<<endl;
    		cout<<endl<<"significance: "<<significance<<endl<<endl;
    		
    		createEntryIndexFiles(nthSample);
    	}
    	outFile.close();
}

/////For fileFunction, "generate" for new sampling or "merge" to use new data for same sampling using entryIndex
void FSRWholeClassification(TString fileFunction, TString preselection, TString vars, int samples)
{
	TString rtdCut="10";
	vector<double> preliminarySignificances, significances, errorsLeft, errorsRight;
	TString outFileText = "analysis/sampleStats"+rtdCut+vars+preselection+samples+"Samples"+".txt"; 
	ofstream outFileCreate(outFileText.Data());
	runAll(samples, fileFunction, preselection, vars, rtdCut, preliminarySignificances, significances, errorsLeft, errorsRight, outFileText);
	
	double meanPreliminarySignificances=-999, meanSignificances=-999, meanErrorsLeft=-999, meanErrorsRight=-999, bestPreliminarySignificances=-999, bestSignificances=-999, bestErrorsLeft=-999, bestErrorsRight=-999, worstPreliminarySignificances=-999, worstSignificances=-999, worstErrorsLeft=-999, worstErrorsRight=-999;
	int bestSample = -999, worstSample=-999;
	getStats(samples, preliminarySignificances, significances, errorsLeft, errorsRight, meanPreliminarySignificances, meanSignificances, meanErrorsLeft, meanErrorsRight, bestSample, worstSample, bestPreliminarySignificances, bestSignificances, bestErrorsLeft, bestErrorsRight, worstPreliminarySignificances, worstSignificances, worstErrorsLeft, worstErrorsRight);
	ofstream outFile(outFileText.Data(), ios::app); 
    	outFile << endl << endl << "Summary:" << endl;
    	outFile << "meanPreliminarySignificances: " << meanPreliminarySignificances << endl;
    	outFile << "meanSignificances: " << meanSignificances << endl;
    	outFile << "meanErrorsLeft: " << meanErrorsLeft << endl;
    	outFile << "meanErrorsRight: " << meanErrorsRight << endl;
    	outFile << "Best sample: " << bestSample << endl;
    	outFile << "bestPreliminarySignificances: " << bestPreliminarySignificances << endl;
    	outFile << "bestSignificances: " << bestSignificances << endl;
    	outFile << "bestErrorsLeft: " << bestErrorsLeft << endl;
    	outFile << "bestErrorsRight: " << bestErrorsRight << endl;
    	outFile << "Worst sample: " << worstSample << endl;
    	outFile << "worstPreliminarySignificances: " << worstPreliminarySignificances << endl;
    	outFile << "worstSignificances: " << worstSignificances << endl;
    	outFile << "worstErrorsLeft: " << worstErrorsLeft << endl;
    	outFile << "worstErrorsRight: " << worstErrorsRight << endl;
    	outFile << endl << endl << "fileFunction: " << fileFunction << endl;
    	outFile << "preselection: " << preselection << endl;
    	outFile << "vars: " << vars << endl;
    	outFile << "samples: " << samples << endl;
    	outFile << "rtdCut: " << rtdCut << endl;
    	outFile.close();
    	
    	return 0;
}

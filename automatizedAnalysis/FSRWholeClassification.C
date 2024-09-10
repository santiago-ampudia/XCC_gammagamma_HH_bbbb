#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <algorithm>

using namespace std;

void createEntryIndexFiles(int nthSample)
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
	TFile *inputFileBqqXTrain = TFile::Open("analysis/outputTreeBqqXHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	TFile *inputFileBqqXTest = TFile::Open("analysis/outputTreeBqqXHHbbbbESpreadDurham1034BSplitTestSampleN.root");
	TFile *inputFileBqqXNNTrain = TFile::Open("analysis/outputTreeBqqXNNESpreadDurham1034BSplitTrainSampleN.root");
	TFile *inputFileBqqXNNTest = TFile::Open("analysis/outputTreeBqqXNNESpreadDurham1034BSplitTestSampleN.root");
	TFile *inputFileBqqqqXTrain = TFile::Open("analysis/outputTreeBqqqqXHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	TFile *inputFileBqqqqXTest = TFile::Open("analysis/outputTreeBqqqqXHHbbbbESpreadDurham1034BSplitTestSampleN.root");
	TFile *inputFileBqqqqXNNTrain = TFile::Open("analysis/outputTreeBqqqqXNNESpreadDurham1034BSplitTrainSampleN.root");
	TFile *inputFileBqqqqXNNTest = TFile::Open("analysis/outputTreeBqqqqXNNESpreadDurham1034BSplitTestSampleN.root");
	TFile *inputFileBqqHXTrain = TFile::Open("analysis/outputTreeBqqHXHHbbbbESpreadDurham1034BSplitTrainSampleN.root");
	TFile *inputFileBqqHXTest = TFile::Open("analysis/outputTreeBqqHXHHbbbbESpreadDurham1034BSplitTestSampleN.root");
	TFile *inputFileBqqHXNNTrain = TFile::Open("analysis/outputTreeBqqHXNNESpreadDurham1034BSplitTrainSampleN.root");
	TFile *inputFileBqqHXNNTest = TFile::Open("analysis/outputTreeBqqHXNNESpreadDurham1034BSplitTestSampleN.root");
	 
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
	TTree *originalTreeBqqXTrain = nullptr;
    inputFileBqqXTrain->GetObject("TreeBqqXTrain", originalTreeBqqXTrain);
	TTree *originalTreeBqqXTest = nullptr;
    inputFileBqqXTest->GetObject("TreeBqqXTest", originalTreeBqqXTest);
    TTree *originalTreeBqqXNNTrain = nullptr;
	inputFileBqqXNNTrain->GetObject("TreeBqqXNNTrain", originalTreeBqqXNNTrain);
	TTree *originalTreeBqqXNNTest = nullptr;
	inputFileBqqXNNTest->GetObject("TreeBqqXNNTest", originalTreeBqqXNNTest);
	TTree *originalTreeBqqqqXTrain = nullptr;
    inputFileBqqqqXTrain->GetObject("TreeBqqqqXTrain", originalTreeBqqqqXTrain);
    TTree *originalTreeBqqqqXTest = nullptr;
    inputFileBqqqqXTest->GetObject("TreeBqqqqXTest", originalTreeBqqqqXTest);
    TTree *originalTreeBqqqqXNNTrain = nullptr;
    inputFileBqqqqXNNTrain->GetObject("TreeBqqqqXNNTrain", originalTreeBqqqqXNNTrain);
    TTree *originalTreeBqqqqXNNTest = nullptr;
    inputFileBqqqqXNNTest->GetObject("TreeBqqqqXNNTest", originalTreeBqqqqXNNTest);
	TTree *originalTreeBqqHXTrain = nullptr;
    inputFileBqqHXTrain->GetObject("TreeBqqHXTrain", originalTreeBqqHXTrain);
    TTree *originalTreeBqqHXTest = nullptr;
    inputFileBqqHXTest->GetObject("TreeBqqHXTest", originalTreeBqqHXTest);
    TTree *originalTreeBqqHXNNTrain = nullptr;
    inputFileBqqHXNNTrain->GetObject("TreeBqqHXNNTrain", originalTreeBqqHXNNTrain);
    TTree *originalTreeBqqHXNNTest = nullptr;
    inputFileBqqHXNNTest->GetObject("TreeBqqHXNNTest", originalTreeBqqHXNNTest);


    TString sampleNameEntryIndexFile = "Sample" + TString::Format("%d", nthSample);
    TString outputFileText = "analysis/entryIndexFile"+sampleNameEntryIndexFile+".root";
    TFile *outputFile = new TFile(outputFileText.Data(), "RECREATE");
    	 
    TTree *TreeSTrainEntryIndex = new TTree("TreeSTrainEntryIndex", "TreeSTrain with only entryIndex branch"); 
    TTree *TreeSTestEntryIndex = new TTree("TreeSTestEntryIndex", "TreeSTest with only entryIndex branch"); 
    TTree *TreeSNNTrainEntryIndex = new TTree("TreeSNNTrainEntryIndex", "TreeSNNTrain with only entryIndex branch");
	TTree *TreeSNNTestEntryIndex = new TTree("TreeSNNTestEntryIndex", "TreeSNNTest with only entryIndex branch");
	TTree *TreeBqqTrainEntryIndex = new TTree("TreeBqqTrainEntryIndex", "TreeBqqTrain with only entryIndex branch");
	TTree *TreeBqqTestEntryIndex = new TTree("TreeBqqTestEntryIndex", "TreeBqqTest with only entryIndex branch");
	TTree *TreeBqqNNTrainEntryIndex = new TTree("TreeBqqNNTrainEntryIndex", "TreeBqqNNTrain with only entryIndex branch");
	TTree *TreeBqqNNTestEntryIndex = new TTree("TreeBqqNNTestEntryIndex", "TreeBqqNNTest with only entryIndex branch");
	TTree *TreeBttTrainEntryIndex = new TTree("TreeBttTrainEntryIndex", "TreeBttTrain with only entryIndex branch");
	TTree *TreeBttTestEntryIndex = new TTree("TreeBttTestEntryIndex", "TreeBttTest with only entryIndex branch");
	TTree *TreeBttNNTrainEntryIndex = new TTree("TreeBttNNTrainEntryIndex", "TreeBttNNTrain with only entryIndex branch");
	TTree *TreeBttNNTestEntryIndex = new TTree("TreeBttNNTestEntryIndex", "TreeBttNNTest with only entryIndex branch");
	TTree *TreeBZZTrainEntryIndex = new TTree("TreeBZZTrainEntryIndex", "TreeBZZTrain with only entryIndex branch");
	TTree *TreeBZZTestEntryIndex = new TTree("TreeBZZTestEntryIndex", "TreeBZZTest with only entryIndex branch");
	TTree *TreeBZZNNTrainEntryIndex = new TTree("TreeBZZNNTrainEntryIndex", "TreeBZZNNTrain with only entryIndex branch");
	TTree *TreeBZZNNTestEntryIndex = new TTree("TreeBZZNNTestEntryIndex", "TreeBZZNNTest with only entryIndex branch");
	TTree *TreeBWWTrainEntryIndex = new TTree("TreeBWWTrainEntryIndex", "TreeBWWTrain with only entryIndex branch");
	TTree *TreeBWWTestEntryIndex = new TTree("TreeBWWTestEntryIndex", "TreeBWWTest with only entryIndex branch");
	TTree *TreeBWWNNTrainEntryIndex = new TTree("TreeBWWNNTrainEntryIndex", "TreeBWWNNTrain with only entryIndex branch");
	TTree *TreeBWWNNTestEntryIndex = new TTree("TreeBWWNNTestEntryIndex", "TreeBWWNNTest with only entryIndex branch");
	TTree *TreeBqqXTrainEntryIndex = new TTree("TreeBqqXTrainEntryIndex", "TreeBqqXTrain with only entryIndex branch");
	TTree *TreeBqqXTestEntryIndex = new TTree("TreeBqqXTestEntryIndex", "TreeBqqXTest with only entryIndex branch");
	TTree *TreeBqqXNNTrainEntryIndex = new TTree("TreeBqqXNNTrainEntryIndex", "TreeBqqXNNTrain with only entryIndex branch");
	TTree *TreeBqqXNNTestEntryIndex = new TTree("TreeBqqXNNTestEntryIndex", "TreeBqqXNNTest with only entryIndex branch");
	TTree *TreeBqqqqXTrainEntryIndex = new TTree("TreeBqqqqXTrainEntryIndex", "TreeBqqqqXTrain with only entryIndex branch");
	TTree *TreeBqqqqXTestEntryIndex = new TTree("TreeBqqqqXTestEntryIndex", "TreeBqqqqXTest with only entryIndex branch");
	TTree *TreeBqqqqXNNTrainEntryIndex = new TTree("TreeBqqqqXNNTrainEntryIndex", "TreeBqqqqXNNTrain with only entryIndex branch");
	TTree *TreeBqqqqXNNTestEntryIndex = new TTree("TreeBqqqqXNNTestEntryIndex", "TreeBqqqqXNNTest with only entryIndex branch");
	TTree *TreeBqqHXTrainEntryIndex = new TTree("TreeBqqHXTrainEntryIndex", "TreeBqqHXTrain with only entryIndex branch");
	TTree *TreeBqqHXTestEntryIndex = new TTree("TreeBqqHXTestEntryIndex", "TreeBqqHXTest with only entryIndex branch");
	TTree *TreeBqqHXNNTrainEntryIndex = new TTree("TreeBqqHXNNTrainEntryIndex", "TreeBqqHXNNTrain with only entryIndex branch");
	TTree *TreeBqqHXNNTestEntryIndex = new TTree("TreeBqqHXNNTestEntryIndex", "TreeBqqHXNNTest with only entryIndex branch");

    float entryIndex, entryIndexSTrain, entryIndexSTest, entryIndexSNNTrain, entryIndexSNNTest, entryIndexqqTrain, entryIndexqqTest, entryIndexqqNNTrain, entryIndexqqNNTest, entryIndexttTrain, entryIndexttTest, entryIndexttNNTrain, entryIndexttNNTest, entryIndexZZTrain, entryIndexZZTest, entryIndexZZNNTrain, entryIndexZZNNTest, entryIndexWWTrain, entryIndexWWTest, entryIndexWWNNTrain, entryIndexWWNNTest, entryIndexqqXTrain, entryIndexqqXTest, entryIndexqqXNNTrain, entryIndexqqXNNTest, entryIndexqqqqXTrain, entryIndexqqqqXTest, entryIndexqqqqXNNTrain, entryIndexqqqqXNNTest, entryIndexqqHXTrain, entryIndexqqHXTest, entryIndexqqHXNNTrain, entryIndexqqHXNNTest;
    	 
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
	originalTreeBqqXTrain->SetBranchAddress("entryIndex", &entryIndexqqXTrain);
	originalTreeBqqXTest->SetBranchAddress("entryIndex", &entryIndexqqXTest);
	originalTreeBqqXNNTrain->SetBranchAddress("entryIndex", &entryIndexqqXNNTrain);
	originalTreeBqqXNNTest->SetBranchAddress("entryIndex", &entryIndexqqXNNTest);
	originalTreeBqqqqXTrain->SetBranchAddress("entryIndex", &entryIndexqqqqXTrain);
	originalTreeBqqqqXTest->SetBranchAddress("entryIndex", &entryIndexqqqqXTest);
	originalTreeBqqqqXNNTrain->SetBranchAddress("entryIndex", &entryIndexqqqqXNNTrain);
	originalTreeBqqqqXNNTest->SetBranchAddress("entryIndex", &entryIndexqqqqXNNTest);
	originalTreeBqqHXTrain->SetBranchAddress("entryIndex", &entryIndexqqHXTrain);
	originalTreeBqqHXTest->SetBranchAddress("entryIndex", &entryIndexqqHXTest);
	originalTreeBqqHXNNTrain->SetBranchAddress("entryIndex", &entryIndexqqHXNNTrain);
	originalTreeBqqHXNNTest->SetBranchAddress("entryIndex", &entryIndexqqHXNNTest);
    	 
    TreeSTrainEntryIndex->Branch("entryIndex", &entryIndexSTrain);
    TreeSTestEntryIndex->Branch("entryIndex", &entryIndexSTest);
    TreeSTrainEntryIndex->Branch("entryIndex", &entryIndexSTrain);
	TreeSTestEntryIndex->Branch("entryIndex", &entryIndexSTest);
	TreeSNNTrainEntryIndex->Branch("entryIndex", &entryIndexSNNTrain);
	TreeSNNTestEntryIndex->Branch("entryIndex", &entryIndexSNNTest);
	TreeBqqTrainEntryIndex->Branch("entryIndex", &entryIndexqqTrain);
	TreeBqqTestEntryIndex->Branch("entryIndex", &entryIndexqqTest);
	TreeBqqNNTrainEntryIndex->Branch("entryIndex", &entryIndexqqNNTrain);
	TreeBqqNNTestEntryIndex->Branch("entryIndex", &entryIndexqqNNTest);
	TreeBttTrainEntryIndex->Branch("entryIndex", &entryIndexttTrain);
	TreeBttTestEntryIndex->Branch("entryIndex", &entryIndexttTest);
	TreeBttNNTrainEntryIndex->Branch("entryIndex", &entryIndexttNNTrain);
	TreeBttNNTestEntryIndex->Branch("entryIndex", &entryIndexttNNTest);
	TreeBZZTrainEntryIndex->Branch("entryIndex", &entryIndexZZTrain);
	TreeBZZTestEntryIndex->Branch("entryIndex", &entryIndexZZTest);
	TreeBZZNNTrainEntryIndex->Branch("entryIndex", &entryIndexZZNNTrain);
	TreeBZZNNTestEntryIndex->Branch("entryIndex", &entryIndexZZNNTest);
	TreeBWWTrainEntryIndex->Branch("entryIndex", &entryIndexWWTrain);
	TreeBWWTestEntryIndex->Branch("entryIndex", &entryIndexWWTest);
	TreeBWWNNTrainEntryIndex->Branch("entryIndex", &entryIndexWWNNTrain);
	TreeBWWNNTestEntryIndex->Branch("entryIndex", &entryIndexWWNNTest);
	TreeBqqXTrainEntryIndex->Branch("entryIndex", &entryIndexqqXTrain);
	TreeBqqXTestEntryIndex->Branch("entryIndex", &entryIndexqqXTest);
	TreeBqqXNNTrainEntryIndex->Branch("entryIndex", &entryIndexqqXNNTrain);
	TreeBqqXNNTestEntryIndex->Branch("entryIndex", &entryIndexqqXNNTest);
	TreeBqqqqXTrainEntryIndex->Branch("entryIndex", &entryIndexqqqqXTrain);
	TreeBqqqqXTestEntryIndex->Branch("entryIndex", &entryIndexqqqqXTest);
	TreeBqqqqXNNTrainEntryIndex->Branch("entryIndex", &entryIndexqqqqXNNTrain);
	TreeBqqqqXNNTestEntryIndex->Branch("entryIndex", &entryIndexqqqqXNNTest);
	TreeBqqHXTrainEntryIndex->Branch("entryIndex", &entryIndexqqHXTrain);
	TreeBqqHXTestEntryIndex->Branch("entryIndex", &entryIndexqqHXTest);
	TreeBqqHXNNTrainEntryIndex->Branch("entryIndex", &entryIndexqqHXNNTrain);
	TreeBqqHXNNTestEntryIndex->Branch("entryIndex", &entryIndexqqHXNNTest);
	 
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
	Long64_t nEntriesBqqXTrain = originalTreeBqqXTrain->GetEntries();
	Long64_t nEntriesBqqXTest = originalTreeBqqXTest->GetEntries();
	Long64_t nEntriesBqqXNNTrain = originalTreeBqqXNNTrain->GetEntries();
	Long64_t nEntriesBqqXNNTest = originalTreeBqqXNNTest->GetEntries();
	Long64_t nEntriesBqqqqXTrain = originalTreeBqqqqXTrain->GetEntries();
	Long64_t nEntriesBqqqqXTest = originalTreeBqqqqXTest->GetEntries();
	Long64_t nEntriesBqqqqXNNTrain = originalTreeBqqqqXNNTrain->GetEntries();
	Long64_t nEntriesBqqqqXNNTest = originalTreeBqqqqXNNTest->GetEntries();
	Long64_t nEntriesBqqHXTrain = originalTreeBqqHXTrain->GetEntries();
	Long64_t nEntriesBqqHXTest = originalTreeBqqHXTest->GetEntries();
	Long64_t nEntriesBqqHXNNTrain = originalTreeBqqHXNNTrain->GetEntries();
	Long64_t nEntriesBqqHXNNTest = originalTreeBqqHXNNTest->GetEntries();

	for (Long64_t i = 0; i < max({nEntriesSTrain, nEntriesSTest, nEntriesSNNTrain, nEntriesSNNTest,
		                          nEntriesBqqTrain, nEntriesBqqTest, nEntriesBqqNNTrain, nEntriesBqqNNTest,
		                          nEntriesBttTrain, nEntriesBttTest, nEntriesBttNNTrain, nEntriesBttNNTest,
		                          nEntriesBZZTrain, nEntriesBZZTest, nEntriesBZZNNTrain, nEntriesBZZNNTest,
		                          nEntriesBWWTrain, nEntriesBWWTest, nEntriesBWWNNTrain, nEntriesBWWNNTest,
		                          nEntriesBqqXTrain, nEntriesBqqXTest, nEntriesBqqXNNTrain, nEntriesBqqXNNTest,
		                          nEntriesBqqqqXTrain, nEntriesBqqqqXTest, nEntriesBqqqqXNNTrain, nEntriesBqqqqXNNTest,
		                          nEntriesBqqHXTrain, nEntriesBqqHXTest, nEntriesBqqHXNNTrain, nEntriesBqqHXNNTest}); ++i) {
		if (i < nEntriesSTrain) {
			originalTreeSTrain->GetEntry(i);
			TreeSTrainEntryIndex->Fill();
		}
		if (i < nEntriesSTest) {
			originalTreeSTest->GetEntry(i);
			TreeSTestEntryIndex->Fill();
	    }
	    if (i < nEntriesSNNTrain) {
		    originalTreeSNNTrain->GetEntry(i);
		    TreeSNNTrainEntryIndex->Fill();
		}
		if (i < nEntriesSNNTest) {
		    originalTreeSNNTest->GetEntry(i);
		    TreeSNNTestEntryIndex->Fill();
		}
		if (i < nEntriesBqqTrain) {
		    originalTreeBqqTrain->GetEntry(i);
		    TreeBqqTrainEntryIndex->Fill();
		}
		if (i < nEntriesBqqTest) {
		    originalTreeBqqTest->GetEntry(i);
		    TreeBqqTestEntryIndex->Fill();
		}
		if (i < nEntriesBqqNNTrain) {
		    originalTreeBqqNNTrain->GetEntry(i);
		    TreeBqqNNTrainEntryIndex->Fill();
		}
		if (i < nEntriesBqqNNTest) {
		    originalTreeBqqNNTest->GetEntry(i);
		    TreeBqqNNTestEntryIndex->Fill();
		}
		if (i < nEntriesBttTrain) {
		    originalTreeBttTrain->GetEntry(i);
		    TreeBttTrainEntryIndex->Fill();
		}
		if (i < nEntriesBttTest) {
		    originalTreeBttTest->GetEntry(i);
		    TreeBttTestEntryIndex->Fill();
		}
		if (i < nEntriesBttNNTrain) {
		    originalTreeBttNNTrain->GetEntry(i);
		    TreeBttNNTrainEntryIndex->Fill();
		}
		if (i < nEntriesBttNNTest) {
		    originalTreeBttNNTest->GetEntry(i);
		    TreeBttNNTestEntryIndex->Fill();
		}
		if (i < nEntriesBZZTrain) {
		    originalTreeBZZTrain->GetEntry(i);
		    TreeBZZTrainEntryIndex->Fill();
		}
		if (i < nEntriesBZZTest) {
		    originalTreeBZZTest->GetEntry(i);
		    TreeBZZTestEntryIndex->Fill();
		}
		if (i < nEntriesBZZNNTrain) {
		    originalTreeBZZNNTrain->GetEntry(i);
		    TreeBZZNNTrainEntryIndex->Fill();
		}
		if (i < nEntriesBZZNNTest) {
		    originalTreeBZZNNTest->GetEntry(i);
		    TreeBZZNNTestEntryIndex->Fill();
		}
		if (i < nEntriesBWWTrain) {
		    originalTreeBWWTrain->GetEntry(i);
		    TreeBWWTrainEntryIndex->Fill();
		}
		if (i < nEntriesBWWTest) {
		    originalTreeBWWTest->GetEntry(i);
		    TreeBWWTestEntryIndex->Fill();
		}
		if (i < nEntriesBWWNNTrain) {
		    originalTreeBWWNNTrain->GetEntry(i);
		    TreeBWWNNTrainEntryIndex->Fill();
		}
		if (i < nEntriesBWWNNTest) {
		    originalTreeBWWNNTest->GetEntry(i);
		    TreeBWWNNTestEntryIndex->Fill();
		}
		if (i < nEntriesBqqXTrain) {
		    originalTreeBqqXTrain->GetEntry(i);
		    TreeBqqXTrainEntryIndex->Fill();
		}
		if (i < nEntriesBqqXTest) {
		    originalTreeBqqXTest->GetEntry(i);
		    TreeBqqXTestEntryIndex->Fill();
		}
		if (i < nEntriesBqqXNNTrain) {
		    originalTreeBqqXNNTrain->GetEntry(i);
		    TreeBqqXNNTrainEntryIndex->Fill();
		}
		if (i < nEntriesBqqXNNTest) {
		    originalTreeBqqXNNTest->GetEntry(i);
		    TreeBqqXNNTestEntryIndex->Fill();
		}
		if (i < nEntriesBqqqqXTrain) {
		    originalTreeBqqqqXTrain->GetEntry(i);
		    TreeBqqqqXTrainEntryIndex->Fill();
		}
		if (i < nEntriesBqqqqXTest) {
		    originalTreeBqqqqXTest->GetEntry(i);
		    TreeBqqqqXTestEntryIndex->Fill();
		}
		if (i < nEntriesBqqqqXNNTrain) {
		    originalTreeBqqqqXNNTrain->GetEntry(i);
		    TreeBqqqqXNNTrainEntryIndex->Fill();
		}
		if (i < nEntriesBqqqqXNNTest) {
		    originalTreeBqqqqXNNTest->GetEntry(i);
		    TreeBqqqqXNNTestEntryIndex->Fill();
		}
		if (i < nEntriesBqqHXTrain) {
		    originalTreeBqqHXTrain->GetEntry(i);
		    TreeBqqHXTrainEntryIndex->Fill();
		}
		if (i < nEntriesBqqHXTest) {
		    originalTreeBqqHXTest->GetEntry(i);
		    TreeBqqHXTestEntryIndex->Fill();
		}
		if (i < nEntriesBqqHXNNTrain) {
		    originalTreeBqqHXNNTrain->GetEntry(i);
		    TreeBqqHXNNTrainEntryIndex->Fill();
		}
		if (i < nEntriesBqqHXNNTest) {
		    originalTreeBqqHXNNTest->GetEntry(i);
		    TreeBqqHXNNTestEntryIndex->Fill();
		}
    }

    outputFile->cd();
    TreeSTrainEntryIndex->Write();
    TreeSTestEntryIndex->Write();
    TreeSNNTrainEntryIndex->Write();
    TreeSNNTestEntryIndex->Write();
    TreeBqqTrainEntryIndex->Write();
    TreeBqqTestEntryIndex->Write();
    TreeBqqNNTrainEntryIndex->Write();
    TreeBqqNNTestEntryIndex->Write();
    TreeBttTrainEntryIndex->Write();
    TreeBttTestEntryIndex->Write();
    TreeBttNNTrainEntryIndex->Write();
    TreeBttNNTestEntryIndex->Write();
    TreeBZZTrainEntryIndex->Write();
    TreeBZZTestEntryIndex->Write();
    TreeBZZNNTrainEntryIndex->Write();
    TreeBZZNNTestEntryIndex->Write();
    TreeBWWTrainEntryIndex->Write();
    TreeBWWTestEntryIndex->Write();
    TreeBWWNNTrainEntryIndex->Write();
    TreeBWWNNTestEntryIndex->Write();
	TreeBqqXTrainEntryIndex->Write();
	TreeBqqXTestEntryIndex->Write();
	TreeBqqXNNTrainEntryIndex->Write();
	TreeBqqXNNTestEntryIndex->Write();
	TreeBqqqqXTrainEntryIndex->Write();
	TreeBqqqqXTestEntryIndex->Write();
	TreeBqqqqXNNTrainEntryIndex->Write();
	TreeBqqqqXNNTestEntryIndex->Write();
	TreeBqqHXTrainEntryIndex->Write();
	TreeBqqHXTestEntryIndex->Write();
	TreeBqqHXNNTrainEntryIndex->Write();
	TreeBqqHXNNTestEntryIndex->Write();
    outputFile->Close();	 
}

double getMean(const vector<double>& dataVector)
{
	double sum=0;
	for(int i=0; i<dataVector.size(); i++) sum += dataVector[i];
	return sum/(dataVector.size());
}

void getStats(int samples, const vector<double>& preliminarySignificances, const vector<double>& preliminaryErrorsLeft, const vector<double>& preliminaryErrorsRight, const vector<double>& significances, const vector<double>& errorsLeft, const vector<double>& errorsRight, double& meanPreliminarySignificances, double& meanPreliminaryErrorsLeft, double& meanPreliminaryErrorsRight, double& meanSignificances, double& meanErrorsLeft, double& meanErrorsRight, int& bestSample, int& worstSample, double& bestPreliminarySignificances, double& bestPreliminaryErrorsLeft, double& bestPreliminaryErrorsRight, double& bestSignificances, double& bestErrorsLeft, double& bestErrorsRight, double& worstPreliminarySignificances, double& worstPreliminaryErrorsLeft, double& worstPreliminaryErrorsRight, double& worstSignificances, double& worstErrorsLeft, double& worstErrorsRight)
{
	//if(samples != significances.size()) throw std::runtime_error("ERROR: size of signifcances vector != n of samples");
	
	meanPreliminarySignificances = getMean(preliminarySignificances);
	meanPreliminaryErrorsLeft = getMean(preliminaryErrorsLeft);
	meanPreliminaryErrorsRight = getMean(preliminaryErrorsRight);
	meanErrorsLeft = getMean(errorsLeft);
	meanErrorsRight = getMean(errorsRight);
	meanSignificances = getMean(significances);
	meanErrorsLeft = getMean(errorsLeft);
	meanErrorsRight = getMean(errorsRight);
	auto max_iter = max_element(significances.begin(), significances.end()); 
	bestSample = distance(significances.begin(), max_iter); 
	bestPreliminarySignificances = preliminarySignificances[bestSample];
	bestPreliminaryErrorsLeft = preliminaryErrorsLeft[bestSample];
	bestPreliminaryErrorsRight = preliminaryErrorsRight[bestSample];	
	bestSignificances = significances[bestSample];
	bestErrorsRight = errorsRight[bestSample];
	bestErrorsLeft = errorsLeft[bestSample];
	auto min_iter = min_element(significances.begin(), significances.end()); 
	worstSample = distance(significances.begin(), min_iter);
	worstPreliminarySignificances = preliminarySignificances[worstSample];	
	worstPreliminaryErrorsLeft = preliminaryErrorsLeft[worstSample];
	worstPreliminaryErrorsRight = preliminaryErrorsRight[worstSample];
	worstSignificances = significances[worstSample];
	worstErrorsRight = errorsRight[worstSample];
	worstErrorsLeft = errorsLeft[worstSample];

	
	cout<<"For "<<samples<<" samples: "<<endl;
	cout<<"preliminarySignificances: ";
	for(int i=0; i<significances.size(); i++) cout<<preliminarySignificances[i]<<", ";
	cout<<endl<<"preliminaryErrorsLeft: ";
	for(int i=0; i<significances.size(); i++) cout<<preliminaryErrorsLeft[i]<<", ";
	cout<<endl<<"preliminaryErrorsRight: ";
	for(int i=0; i<significances.size(); i++) cout<<preliminaryErrorsRight[i]<<", ";
	cout<<endl<<"significances: ";
	for(int i=0; i<significances.size(); i++) cout<<significances[i]<<", ";
	cout<<endl<<"errorsLeft: ";
	for(int i=0; i<significances.size(); i++) cout<<errorsLeft[i]<<", ";
	cout<<endl<<"errorsRight: ";
	for(int i=0; i<significances.size(); i++) cout<<errorsRight[i]<<", ";
	cout<<endl<<"meanPreliminarySignificances: "<<meanPreliminarySignificances<<endl;
	cout<<"meanPreliminaryErrorsLeft: "<<meanPreliminaryErrorsLeft<<endl;
	cout<<"meanPreliminaryErrorsRight: "<<meanPreliminaryErrorsRight<<endl;
	cout<<"meanSignificances: "<<meanSignificances<<endl;
	cout<<"meanErrorsLeft: "<<meanErrorsLeft<<endl;
	cout<<"meanErrorsRight: "<<meanErrorsRight<<endl;
	cout<<"sample that generates the best final significance: "<<bestSample<<endl;
	cout<<"bestPreliminarySignificances: "<<bestPreliminarySignificances<<endl;
	cout<<"bestPreliminaryErrorsLeft: "<<bestPreliminaryErrorsLeft<<endl;
	cout<<"bestPreliminaryErrorsRight: "<<bestPreliminaryErrorsRight<<endl;
	cout<<"bestSignificances: "<<bestSignificances<<endl;
	cout<<"bestErrorsLeft: "<<bestErrorsLeft<<endl;
	cout<<"bestErrorsRight: "<<bestErrorsRight<<endl;
	cout<<"sample that generates the worst final significance: "<<worstSample<<endl;
	cout<<"worstPreliminarySignificances: "<<worstPreliminarySignificances<<endl;
	cout<<"worstPreliminaryErrorsLeft: "<<worstPreliminaryErrorsLeft<<endl;
	cout<<"worstPreliminaryErrorsRight: "<<worstPreliminaryErrorsRight<<endl;
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
	back = "qqX";
	gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationHHbbbb.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", back.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
	back = "qqqqX";
	gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationHHbbbb.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", back.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
	back = "qqHX";
	gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationHHbbbb.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", back.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
}


void runAll(int samples, TString fileFunction, TString preselection, TString vars, TString rtdCut, vector<double>& preliminarySignificances, vector<double>& preliminaryErrorsLeft, vector<double>& preliminaryErrorsRight, vector<double>& significances, vector<double>& errorsLeft, vector<double>& errorsRight, TString outFileText)
{
	ofstream outFile(outFileText.Data(), ios::app);
	double prevPreliminarySignificance=-999, prevPreliminaryErrorRight=-999, prevPreliminaryErrorLeft=-999, prevSignificance=-999, prevErrorRight=-999, prevErrorLeft=-999;
	for(int nthSample=0; nthSample<samples; nthSample++)
    {
    	TString sampleName = "Sample" + TString::Format("%d", nthSample);
    	sampleName = "SampleN";
    	//gSystem->Exec("root -l -b -q 'analysis/FSRGammaGammaHHbbbbAnalysis.C'");
    	//trainAllBacks(rtdCut, preselection, vars, sampleName);
    	//gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationApplicationHHbbbbGeneratesNNs.C+(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", fileFunction.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    	//gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationOutputNN.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\")'", rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    	gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationApplicationOutputNN.C+(\"%s\", \"%s\", \"%s\", \"%s\")'", rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    
    	double preliminarySignificance=-999, preliminaryErrorLeft=-999, preliminaryErrorRight=-999, significance=-999, errorLeft=-999, errorRight=-999;
    	ifstream file1("analysis/preliminarySignificanceAndErrorsFile.txt");
		file1 >> preliminarySignificance >> preliminaryErrorLeft >> preliminaryErrorRight;
		file1.close();
		ifstream file2("analysis/significanceAndErrorsFile.txt");
		file2 >> significance >> errorLeft >> errorRight;
		file2.close();
    		
    	if(preliminarySignificance == prevPreliminarySignificance && preliminaryErrorLeft == prevPreliminaryErrorLeft && preliminaryErrorRight == prevPreliminaryErrorRight && significance == prevSignificance && errorLeft == prevErrorLeft && errorRight == prevErrorRight)  outFile << "sample: " << nthSample << endl << "CRASH SOMEWHERE" << endl;
    	else 
    	{
    		preliminarySignificances.push_back(preliminarySignificance);
			preliminaryErrorsLeft.push_back(preliminaryErrorLeft);
			preliminaryErrorsRight.push_back(preliminaryErrorRight);
			significances.push_back(significance);
			errorsLeft.push_back(errorLeft);
	    	errorsRight.push_back(errorRight);
    			
    		outFile << "sample: " << nthSample << endl;
			outFile << "preliminary significance: " << preliminarySignificance << endl;
			outFile << "preliminary errorsLeft: " << preliminaryErrorLeft << endl;
			outFile << "preliminary errorsRight: " << preliminaryErrorRight << endl;
			outFile << "significance: " << significance << endl;
			outFile << "errorsLeft: " << errorLeft << endl;
			outFile << "errorsRight: " << errorRight << endl;
			outFile << endl;
    	} 
    		prevPreliminarySignificance = preliminarySignificance;
			prevPreliminaryErrorLeft = preliminaryErrorLeft;
			prevPreliminaryErrorRight = preliminaryErrorRight;
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
	vector<double> preliminarySignificances, preliminaryErrorsLeft, preliminaryErrorsRight, significances, errorsLeft, errorsRight;
	TString outFileText = "analysis/sampleStats"+rtdCut+vars+preselection+samples+"Samples"+".txt"; 
	ofstream outFileCreate(outFileText.Data());
	runAll(samples, fileFunction, preselection, vars, rtdCut, preliminarySignificances, preliminaryErrorsLeft, preliminaryErrorsRight, significances, errorsLeft, errorsRight, outFileText);
	
	double meanPreliminarySignificances=-999, meanPreliminaryErrorsLeft=-999, meanPreliminaryErrorsRight=-999, meanSignificances=-999, meanErrorsLeft=-999, meanErrorsRight=-999, bestPreliminarySignificances=-999, bestPreliminaryErrorsLeft=-999, bestPreliminaryErrorsRight=-999, bestSignificances=-999, bestErrorsLeft=-999, bestErrorsRight=-999, worstPreliminarySignificances=-999, worstPreliminaryErrorsLeft=-999, worstPreliminaryErrorsRight=-999, worstSignificances=-999, worstErrorsLeft=-999, worstErrorsRight=-999;
	int bestSample = -999, worstSample=-999;
	getStats(samples, preliminarySignificances, preliminaryErrorsLeft, preliminaryErrorsRight, significances, errorsLeft, errorsRight, meanPreliminarySignificances, meanPreliminaryErrorsLeft, meanPreliminaryErrorsRight, meanSignificances, meanErrorsLeft, meanErrorsRight, bestSample, worstSample, bestPreliminarySignificances, bestPreliminaryErrorsLeft, bestPreliminaryErrorsRight, bestSignificances, bestErrorsLeft, bestErrorsRight, worstPreliminarySignificances, worstPreliminaryErrorsLeft, worstPreliminaryErrorsRight, worstSignificances, worstErrorsLeft, worstErrorsRight);
	ofstream outFile(outFileText.Data(), ios::app); 
    outFile << endl << endl << "Summary:" << endl;
    outFile << "meanPreliminarySignificances: " << meanPreliminarySignificances << endl;
	outFile << "meanPreliminaryErrorsLeft: " << meanPreliminaryErrorsLeft << endl;
	outFile << "meanPreliminaryErrorsRight: " << meanPreliminaryErrorsRight << endl;
    outFile << "meanSignificances: " << meanSignificances << endl;
    outFile << "meanErrorsLeft: " << meanErrorsLeft << endl;
    outFile << "meanErrorsRight: " << meanErrorsRight << endl;
	outFile << "stdPreliminarySignificances: " << TMath::StdDev(preliminarySignificances.size(), &preliminarySignificances[0]) << endl;
    outFile << "variancePreliminarySignificances: " << (TMath::StdDev(preliminarySignificances.size(), &preliminarySignificances[0]))*(TMath::StdDev(preliminarySignificances.size(), &preliminarySignificances[0])) << endl;
    outFile << "stdSignificances: " << TMath::StdDev(significances.size(), &significances[0]) << endl;
    outFile << "varianceSignificances: " << (TMath::StdDev(significances.size(), &significances[0]))*(TMath::StdDev(significances.size(), &significances[0])) << endl;
    outFile << "Best sample: " << bestSample << endl;
    outFile << "bestPreliminarySignificances: " << bestPreliminarySignificances << endl;
	outFile << "bestPreliminaryErrorsLeft: " << bestPreliminaryErrorsLeft << endl;
	outFile << "bestPreliminaryErrorsRight: " << bestPreliminaryErrorsRight << endl;
    outFile << "bestSignificances: " << bestSignificances << endl;
    outFile << "bestErrorsLeft: " << bestErrorsLeft << endl;
    outFile << "bestErrorsRight: " << bestErrorsRight << endl;
    outFile << "Worst sample: " << worstSample << endl;
    outFile << "worstPreliminarySignificances: " << worstPreliminarySignificances << endl;
	outFile << "worstPreliminaryErrorsLeft: " << worstPreliminaryErrorsLeft << endl;
	outFile << "worstPreliminaryErrorsRight: " << worstPreliminaryErrorsRight << endl;
    outFile << "worstSignificances: " << worstSignificances << endl;
    outFile << "worstErrorsLeft: " << worstErrorsLeft << endl;
    outFile << "worstErrorsRight: " << worstErrorsRight << endl;
    outFile << endl << endl << "fileFunction: " << fileFunction << endl;
    outFile << "preselection: " << preselection << endl;
    outFile << "vars: " << vars << endl;
    outFile << "samples: " << samples << endl;
    outFile << "rtdCut: " << rtdCut << endl;
    outFile.close();
    	
}

#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <algorithm>

using namespace std;

double getMean(const vector<double>& dataVector)
{
	double sum=0;
	for(int i=0; i<dataVector.size(); i++) sum += dataVector[i];
	cout<<"size: "<<dataVector.size()<<endl;
	return sum/(dataVector.size());
}

void getStats(int samples, const vector<double>& preliminarySignificances, const vector<double>& significances, const vector<double>& errorsLeft, const vector<double>& errorsRight, double& meanPreliminarySignificances, double& meanSignificances, double& meanErrorsLeft, double& meanErrorsRight, int& bestSample)
{
	//if(samples != significances.size()) throw std::runtime_error("ERROR: size of signifcances vector != n of samples");
	
	meanPreliminarySignificances = getMean(preliminarySignificances);
	meanSignificances = getMean(significances);
	meanErrorsLeft = getMean(errorsLeft);
	meanErrorsRight = getMean(errorsRight);
	auto max_iter = max_element(significances.begin(), significances.end()); 
	bestSample = distance(significances.begin(), max_iter); 
	
	
	cout<<"For "<<samples<<" samples: "<<endl;
	cout<<"preliminarySignificances: ";
	for(int i=0; i<significances.size(); i++) cout<<preliminarySignificances[i]<<", ";
	cout<<endl<<"significances: ";
	for(int i=0; i<significances.size(); i++) cout<<significances[i]<<", ";
	cout<<endl<<"errorsLeft: ";
	for(int i=0; i<significances.size(); i++) cout<<errorsLeft[i]<<", ";
	cout<<endl<<"errorsRight: ";
	for(int i=0; i<significances.size(); i++) cout<<errorsRight[i]<<", ";
	cout<<"meanPreliminarySignificances: "<<meanPreliminarySignificances<<endl;
	cout<<"meanSignificances: "<<meanSignificances<<endl;
	cout<<"meanErrorsLeft: "<<meanErrorsLeft<<endl;
	cout<<"meanErrorsRight: "<<meanErrorsRight<<endl<<endl;
	cout<<"sample that generates the best final significance: "<<bestSample<<endl<<endl;
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


void runAll(int samples, TString fileFunction, TString preselection, TString vars, TString rtdCut, vector<double>& preliminarySignificances, vector<double>& significances, vector<double>& errorsLeft, vector<double>& errorsRight)
{
	ofstream outFile("analysis/samplesStats.txt", ios::app);
	double prevPreliminarySignificance=-999, prevSignificance=-999, prevErrorRight=-999, prevErrorLeft=-999;
	for(int nthSample=0; nthSample<samples; nthSample++)
    	{
    		TString sampleName = "Sample" + TString::Format("%d", nthSample);
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
    	}
    	outFile.close();
}

/////For fileFunction, "generate" for new sampling or "merge" to use new data for same sampling using entryIndex
void FSRWholeClassification(TString fileFunction, TString preselection, TString vars, int samples)
{
	TString rtdCut="10";
	vector<double> preliminarySignificances, significances, errorsLeft, errorsRight;
	ofstream outFileCreate("analysis/samplesStats.txt");
	runAll(samples, fileFunction, preselection, vars, rtdCut, preliminarySignificances, significances, errorsLeft, errorsRight);
	
	double meanPreliminarySignificances=-999, meanSignificances=-999, meanErrorsLeft=-999, meanErrorsRight=-999;
	int bestSample = -999;
	getStats(samples, preliminarySignificances, significances, errorsLeft, errorsRight, meanPreliminarySignificances, meanSignificances, meanErrorsLeft, meanErrorsRight, bestSample);
	ofstream outFile("analysis/samplesStats.txt", ios::app); 
    	outFile << endl << endl << "Summary:" << endl;
    	outFile << "meanPreliminarySignificances: " << meanPreliminarySignificances << endl;
    	outFile << "meanSignificances: " << meanSignificances << endl;
    	outFile << "meanErrorsLeft: " << meanErrorsLeft << endl;
    	outFile << "meanErrorsRight: " << meanErrorsRight << endl;
    	outFile << "Best sample: " << bestSample << endl;
    	outFile << endl << endl << "fileFunction: " << fileFunction << endl;
    	outFile << "preselection: " << preselection << endl;
    	outFile << "vars: " << vars << endl;
    	outFile << "samples: " << samples << endl;
    	outFile << "rtdCut: " << rtdCut << endl;
    	outFile.close();
    	
    	return 0;
}

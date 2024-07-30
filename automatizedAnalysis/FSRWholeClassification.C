#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>

using namespace std;

double getMean(const vector<double>& dataVector)
{
	double sum=0;
	for(int i=0; i<dataVector.size(); i++) sum += dataVector[i];
	return sum/(dataVector.size());
}

void getStats(int samples, const vector<double>& preliminarySignificances, const vector<double>& significances, const vector<double>& errorsLeft, const vector<double>& errorsRight, double& meanPreliminarySignificances, double& meanSignificances, double& meanErrorsLeft, double& meanErrorsRight)
{
	if(samples != significances.size()) throw std::runtime_error("ERROR: size of signifcances vector != n of samples");
	
	meanPreliminarySignificances = getMean(preliminarySignificances);
	meanSignificances = getMean(significances);
	meanErrorsLeft = getMean(errorsLeft);
	meanErrorsRight = getMean(errorsRight);
	
	cout<<"For "<<samples<<" samples: "<<endl;
	cout<<"meanPreliminarySignificances: "<<meanPreliminarySignificances<<endl;
	cout<<"meanSignificances: "<<meanSignificances<<endl;
	cout<<"meanErrorsLeft: "<<meanErrorsLeft<<endl;
	cout<<"meanErrorsRight: "<<meanErrorsRight<<endl<<endl;
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
	for(int nthSample=0; nthSample<samples; nthSample++)
    	{
    		TString sampleName = "Sample" + TString::Format("%d", nthSample);
    		//gSystem->Exec("root -l -b -q 'analysis/FSRGammaGammaHHbbbbAnalysis.C'");
    		//trainAllBacks(rtdCut, preselection, vars, sampleName);
    		gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationApplicationHHbbbbGeneratesNNs.C+(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\")'", fileFunction.Data(), rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    		//gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationOutputNN.C+(\"\", \"%s\", \"%s\", \"%s\", \"%s\")'", rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    		gSystem->Exec(Form("root -l -b -q 'analysis/FSRTMVAClassificationApplicationOutputNN.C+(\"%s\", \"%s\", \"%s\", \"%s\")'", rtdCut.Data(), preselection.Data(), vars.Data(), sampleName.Data()));
    		
    		ifstream file1("analysis/preliminarySignificanceFile.txt");
    		double preliminarySignificance;
		if (file1 >> preliminarySignificance) preliminarySignificances.push_back(preliminarySignificance);
		else preliminarySignificances.push_back(-999); 
		file1.close();
		ifstream file2("analysis/significanceAndErrorsFile.txt");
		double significance, errorLeft, errorRight;
		if (file2 >> significance >> errorLeft >> errorRight) {
		    significances.push_back(significance);
		    errorsLeft.push_back(errorLeft);
		    errorsRight.push_back(errorRight);
		} else {
		    significances.push_back(-999);
		    errorsLeft.push_back(-999);
		    errorsRight.push_back(-999); 
		}
		file2.close();
    		
    		cout<<endl<<"preliminarySignificance: "<<preliminarySignificance<<endl<<endl;
    		cout<<endl<<"significance: "<<significance<<endl<<endl;
    		
    	}
}

/////For fileFunction, "generate" for new sampling or "merge" to use new data for same sampling using entryIndex
void FSRWholeClassification(TString fileFunction, TString preselection, TString vars, int samples)
{
	////34BSplit
	TString rtdCut="10";
	vector<double> preliminarySignificances, significances, errorsLeft, errorsRight;
	runAll(samples, fileFunction, preselection, vars, rtdCut, preliminarySignificances, significances, errorsLeft, errorsRight);
	double meanPreliminarySignificances=-999, meanSignificances=-999, meanErrorsLeft=-999, meanErrorsRight=-999;
	getStats(samples, preliminarySignificances, significances, errorsLeft, errorsRight, meanPreliminarySignificances, meanSignificances, meanErrorsLeft, meanErrorsRight);
    	return 0;
}

If using analysis code for backgrounds post WW, skip step 2 and download all the root files for said backgrounds in same folders as the first 4.




Instructions to run the whole automatized analyisis:
   1) Setting up the environment:

1.1) All of these files are in github.com/santiago-ampudia/XCC_gammagamma_HH_bbbb, in the branch "fittedAutomatizedAnalysisAfterqqttbarZZWWqqXqqqqXqqHXZH".

  1.2) In the main delphes folder, create a folder called "analysis".
  
  1.3) In the analysis folder, create a folder called "filesPostDelphes". In that folder, download from https://drive.google.com/drive/folders/1J3-QTyudBDuOUP8xE-UI6KjnA1BFSu4H?usp=drive_link the following files (conserve the names):
  
    1.3.1) GammaGammaHHESpreadAll.root (https://drive.google.com/file/d/1RYew-EAp7osoZN1_TDRtElfalj0ktf4W/view?usp=drive_link)
    
    1.3.2) GammaGammabbbbqqESpreadAll.root (https://drive.google.com/file/d/1Itd33G35XDI5zQxDRCihQEDOuea-a7wb/view?usp=drive_link)
    
    1.3.3) GammaGammattAll.root (https://drive.google.com/file/d/1sgGuIcBX_9xLRiCRS6HF2gT5Ae1IWh-T/view?usp=drive_link)
    
    1.3.4) GammaGammaZZESpreadAll.root (https://drive.google.com/file/d/13L4Ahgb8dTT4okx2W5KHMd1_T2NXC8T7/view?usp=drive_link)
    
    1.3.5) GammaGammaWWESpreadAll.root (https://drive.google.com/file/d/1x6FT7avVS4VxrHtToImq2Jf7CJbYl5FZ/view?usp=drive_link)

    1.3.6) eGammaqqXESpreadAll.root (https://drive.google.com/file/d/133gpNNG6VWHIT79wuP2paYfFXOAKV5ue/view?usp=share_link)

    1.3.7) eGammaqqqqXESpreadAll.root (https://drive.google.com/file/d/1dtZ882B9XKreyCj8earzoii1QL6K_N1f/view?usp=share_link)

    1.3.8) eGammaqqHXESpreadAll.root (https://drive.google.com/file/d/1_-nl4B7McsB0a4b_c_thcVJJ-MDVUXLe/view?usp=share_link)


SKIPPPPPPPPP
    
  1.4) In the analysis folder, create a folder called "SampleOG" and download there the files post-initial preselection used for benchmark significance; found in automatizedAnalysis/TMVA/postPreselectionTrainingAndTestingFiles (link to drive: https://drive.google.com/drive/folders/1mVbKmLJLOGIIjmmfpWWev66vmpt5V6Zu?usp=drive_link); conserve the names of the files

SKIPPPPPPP
  1.5) In the analysis folder, create a folder called "SampleOGNNOfNNs" and download there the files for NNofNNs training and testing; files found in automatizedAnalysis/TMVA/filesForTrainignAndTestingNNOfNNs (link to drive: https://drive.google.com/drive/folders/1TplUy9B_DsVpDBlkfLFu29ZGILggFN_B?usp=sharing); conserve the names.
  
  1.6) In the main delphes folder, download FSRWholeClassification.C, found in automatizedAnalysis, as well as all the files in automatizedAnalysis/kinematicFit.
  
  1.7) In the analysis folder, download FSRGammaGammaHHbbbbAnalysis.C and samplingTrainTest.C,found in automatizedAnalysis/analysis.
  
  1.8) In the analysis folder, download the following files found in automatizedAnalysis/TMVA:
  
    1.8.1) FSRTMVAClassificationHHbbbb.C (https://github.com/santiago-ampudia/XCC_gammagamma_HH_bbbb/blob/fittedAutomatizedAnalysisAfterqqttbarZZWWqqXqqqqXqqHXZH/automatizedAnalysis/TMVA/FSRTMVAClassificationHHbbbb.C)
    
    1.8.2) FSRTMVAClassificationApplicationHHbbbbGeneratesNNs.C (https://github.com/santiago-ampudia/XCC_gammagamma_HH_bbbb/blob/fittedAutomatizedAnalysisAfterqqttbarZZWWqqXqqqqXqqHXZH/automatizedAnalysis/TMVA/FSRTMVAClassificationApplicationHHbbbbGeneratesNNs.C)

    1.8.3) TMVACutsGA.C (https://github.com/santiago-ampudia/XCC_gammagamma_HH_bbbb/blob/fittedAutomatizedAnalysisAfterqqttbarZZWWqqXqqqqXqqHXZH/automatizedAnalysis/TMVA/TMVACutsGA.C)
    
    1.8.4) FSRTMVAClassificationOutputNN.C (https://github.com/santiago-ampudia/XCC_gammagamma_HH_bbbb/blob/fittedAutomatizedAnalysisAfterqqttbarZZWWqqXqqqqXqqHXZH/automatizedAnalysis/TMVA/FSRTMVAClassificationOutputNN.C)
    
    1.8.5) FSRTMVAClassificationApplicationOutputNN.C (https://github.com/santiago-ampudia/XCC_gammagamma_HH_bbbb/blob/fittedAutomatizedAnalysisAfterqqttbarZZWWqqXqqqqXqqHXZH/automatizedAnalysis/TMVA/FSRTMVAClassificationApplicationOutputNN.C)
    
/////skip if post WW back
2) Replicating result and confirming proper set up:
   
  2.1) In FSRGammaGammaHHbbbbAnalysis.C, change in line 3935, "string fileFunction = "generate";" for "string fileFunction = "merge";"
  
  2.2) cd into the main delphes folder and run "root -l -b -q FSRWholeClassification.C'( "merge", "34BSplit", "All4V", 1)' "
  
    2.2.1) First argument indicates that we are using the same sampling as the one used for the benchmark result (alternatively, we could use "generate" instead of "merge" for new sampling; use "merge" to replicate result).
    
    2.2.2) Second argument indicates the preselection ("34BSplit" means that we are accepting events with 4b jets or 3b 1 non-b jets, and splitting the data for training and testing 75%-25%).
    
    2.2.3) Third argument indicates the batch of variables used (check file in variablesNN/variablesVersionDescription.odt for the different batches of variables (link: https://github.com/santiago-ampudia/XCC_gammagamma_HH_bbbb/blob/automatizedAnalysisAfterqqttbarZZWW/variablesNN/variablesVersionDescription.odt))
    
    2.2.4) Fourth argument is how many trials (samples) you want for the analysis stats. To replicate a result, just run it once.
    
  2.3) Compare results (for benchmark result, these are the stats: 
  
      preliminarySignificance: 5.84281
      
      significance: 7.2256
      
      For 1 samples: 
      
      meanPreliminarySignificances: 5.84281
      
      meanSignificances: 7.2256
      
      meanErrorsLeft: 0.135458
      
      meanErrorsRight: 0.141374
      
      sample that generates the best final significance: 0)
   //////skip up until here if post WW back

3) Getting new sampling-independent statistics:

  //////Skip if skipped step 2  3.1) In the analysis folder, delete the files generated by the previous analysisis (i.e. outputTreeSHHbbbbESpreadDurham1034BSplitTrainSample0.root, ..., outputTreeBWWHHbbbbESpreadDurham1034BSplitTestSample0.root, as well as outputTreeSNNESpreadDurham1034BSplitTrainSample0.root, ..., outputTreeBWWNNESpreadDurham1034BSplitTestSample0.root)///////
  
  /////skip if skipped step 2 3.2) In FSRGammaGammaHHbbbbAnalysis.C, change in line 4193, "string fileFunction = "merge";" for "string fileFunction = "generate";"
  
  3.3) cd into the main delphes folder and run rroot -l FSRGammaGammaHHbbbbAnalysis.C. make sure all analysis files post preselction are generated for all topologies
  
  3.3) run "root -l -b -q FSRWholeClassification.C'( "generateGA", "34BSplit", "All9V", 1)'  ".

  3.4) Now in FSRWholeClassification comment out line 533 ("gSystem->Exec("root -l -b -q 'analysis/FSRGammaGammaHHbbbbAnalysis.C'");")
  
  3.5) Run "root -l -b -q FSRWholeClassification.C'( "generateGA", "34BSplit", "All9V", n)'  ". for high n
  
  3.6) Get the results and save them, as well as the files associated to the sample with the highest significance.

  Benchmark results: 
      Summary:
      preliminaryErrorsLeft: 0.0952726, 0.0944913, 0.0928169, 0.094268, 
      preliminaryErrorsRight: 0.098733, 0.0981191, 0.0961098, 0.0977842, 
      significances: 10.6277, 10.8693, 11.0237, 10.8426, 
      errorsLeft: 0.0922029, 0.090082, 0.0889658, 0.0904727, 
      errorsRight: 0.096054, 0.0939331, 0.0924262, 0.0939889, 
      meanPreliminarySignificances: 10.4243
      meanPreliminaryErrorsLeft: 0.0942122
      meanPreliminaryErrorsRight: 0.0976865
      meanSignificances: 10.8408
      meanErrorsLeft: 0.0904309
      meanErrorsRight: 0.0941006
      stdPreliminarySignificances: 0.116846
      variancePreliminarySignificances: 0.0136529
      stdSignificances: 0.162972
      varianceSignificances: 0.0265598
      sample that generates the best final significance: 2
      bestPreliminarySignificances: 10.5869
      bestPreliminaryErrorsLeft: 0.0928169
      bestPreliminaryErrorsRight: 0.0961098
      bestSignificances: 11.0237
      bestErrorsLeft: 0.0889658
      bestErrorsRight: 0.0924262
      sample that generates the worst final significance: 0
      worstPreliminarySignificances: 10.3115
      worstPreliminaryErrorsLeft: 0.0952726
      worstPreliminaryErrorsRight: 0.098733
      worstSignificances: 10.6277
      worstErrorsLeft: 0.0922029
      worstErrorsRight: 0.096054
      ampudia@DNa8127d3 delphes % 
      ampudia@DN51udmt noFit % 
      
      
      
      
      **4) Get groundbreaking results and save physics and the future of academia.**
      

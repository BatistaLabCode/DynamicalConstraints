# Code Package for Dynamical Constraints (Nature Neuroscience, 2025)

## Overview
This repository contains the data and code used for the analysis and figure generation in our recently published paper, [Dynamical constraints on neural population activity]. The data provided will include processed trial data from example session (E20190719), processed trial data of the two target task from a single experiment session(E20181004), flow field data from all valid sessions, and summary source data for all statistical analysis in the paper. For the code, we provide all the standard functions we used process, analyze, and plot the data as well simplified scripts and functions that produce the values of the summary data structures and all the plots found in the main figures.

*Notes: *
+*The dataHigh function "orthogonalize" is called in several functions in this code pack. Orthogonalize is performed using the MATLAB built-in function svd which can report different, but valid singlar vectors depending on MATLAB version and machine use. This may cause visual descripancies in when Figure 3 panels or plotting of the online separationMaximizingProjection. See MATLAB [svd documentation page](https://www.mathworks.com/help/releases/R2024b/matlab/ref/double.svd.html?overload=svd+false#bu2_0hq-3) for more information.*
+*There are slight discrepancies in the trajectory arrows shown in the paper and the arrows plotted using* createAllFigures *script. This is because the trajectory arrows shown in the main figures of the paper were added using Adobe Illustrator, which calculates direction of the arrow head differently than the matlab code.*

## Basic File Organization 
### 1) Data Location: DynamicalConstraints\DynamicalConstraints_NatNeuro_2024_data\
All data is found in this folder. See below for a brief summary of each:
  + *exampleDatasetCatalog*: Structure that containes file details about the two example sessions provided.
  + *dPrime_workspace+neuralSpace_fig3*: Summary structures of the discriminability metric (d') for sessions. This metric is calculated the MoveInt trials projected into both the MoveInt workspace and the SepMax workspace. (See Figure 3 for plots related to the data)
  + *dPrime_EarlyLate+ShuffleData*: Structure storing each sessions d' values calculated over the early two target trials, the late two target trials, and with random shuffling. (See figure 4)
  + *initialAng_tt_vs_uncon+con_Simplified*: Stores AD_compare, the raw initial angle values used in figures 6 and 7. This includes the individual trial comparisons of IT task and the instructed path task to the direct IT target path and the initial angles for the two target task performed with SepMax visual feedback. Note that these metrics are calculated both for the online SepMax control trials and for MoveInt trials projected into the SepMax space. (See figure 6 and 7).
  + *centerOut_AngleAnalatsisData_Signed*: Contains all sessions' raw initial angle values that are used to calculate the "full change control". Stores the structure as AD_CO (See figure 6 and 7).
  + Folder summary:
	+ *Example Sessions\...: Folder that contains the processed data, latent space parameters, and BCI decoders for the two example session (monkeyE20181004 and monkeyE20190719). The neural and behavior data is saved as TrajectoryData object (lab-custom design object included in the code pact) and everything else is saved as structures. Each unique session task or parameter change is save as a unique folder (ie monkeyE20190719_04_condGridTask_01_SI_trajectoryData or monkeyE20190719_13_twoTargetABBA_intTarg_startTarg_T1_tube_04_SI_trajectoryData). **Example session: monkeyE20190719 is the input example session used in the example analysis scripts described below.**
   	+ *flowAnalysis\mat\...*: Includes FlowResult structures for each valid experiment with each file using the following naming structure: [animalID][session date]_[MoveInt two target folder number] _[SepMax two target folder number]_FlowResults. See the code flow.main_analysis for greater description of the structure.
	+ *ConstrainedPath\mat\int_targ_data\...*: Stores the Intermediate Target (intTarg) structure of the example session monkeyE20190719. This data is used in figure 6 and 7 and stores the IT task analysis and the instructed path task analysis.

### 2) Analysis examples location: DynamicalConstraints\exampleAnalysisScripts\
There are six example analysis script that we provided. These a simplified versions many of the "batch" analysis scripts provided. They are also base code in case user would like to implement any of our methods on their data. The scripts are as follow:
+ identify_separationMaximizingProjection
+ dPrime_calc
+ dPrime_calc_EL_shuffle
+ flowField_calc
+ calc_initialAngle_data
+ calc_centerOutAngle_control

### 3) Recreate manuscript figures: DynamicalConstraints\
In the main folder of DynamicalConstraints, there are six functions that will produce figure and one wrap script that will call all figure functions and produce their figures and statistical test results. This wrap script is called *createAllPaperFigures* and will automatically plot all figures as it calls the following figure functions:
 
1) fig_2_two_target
2) fig_3
3) fig_4_early_vs_late_two_target_trajectories_comparisons
4) fig_4_plot_View_vs_Projections
5) fig_5_plot_flow_field_analysis
6) fig_6_fig_7

Users have the option to save the figures and even subselect the desired plots (ie only plot figure 3 and 5). Furthermore there are additional variables that the user can adjust for each figure function. Review the header information in each figure code of to determine the options that can be modified.

## Getting started (Create the paper results)

To run the main analysis, you will primarily interact with the following functions:

- **DynamicalConstraints\serverPath.m**: location of data folder on the computer. This is a function that will attempt to identify the location of the data information, the example session information, and where the user wants to save the material. If the code is unable to find either the automated path or the user input path, than the code will prompt the user to identify the correct folders using the file GUI.

- **DynamicalConstraints\createAllFigures**: produces all the main paper figures and their related statistical tests. 

### Steps
1) Open *serverPath.m* and confirm that is using the correct locations for a) dataLoc, b) example session (ie exSessDataLoc), and c) save figure location (ie saveFigLoc)
	+ Default settings assume that the data folder <<DynamicalConstraints_NatNeuro_2024_data>> in still maintained in the github folder and that the folder with the processed example session data is stored in the <<DynamicalConstraints_NatNeuro_2024_data\Example Sessions>> folder.
2) Update *serverPath.m* variables dataLoc, exSessDataLoc, and saveFigLoc are different from where the user has saved the data.
3) Confirm *serverPath* runs without error.
	+ Code will prompt the user to give a new location if the defined location does not exist.
4) Run either *createAllFigures.m* to create the main paper figures or any script in the *exampleAnalysisScripts* folder to test the analysis

It is possible to also run the individual figure codes separately, e.g. fig_2_two_target(dataLoc), or adjust optional inputs/variables in the functions.

## Additional Folder Organization
The repository is organized into several folders, each with a specific purpose and even anaylsis focus.

- @ folders: These are custom MATLAB classes built specifically for projects in the Batista lab. The advantage of these classes is that they enable users to create oject specific anaylsis functions (An example of this is the TubeOject.plot and TrajectoryData.plot are difference functions that will produce difference results). While they could be replaced with structures in the future, they are necessary for several minor codes and will remain until we finalize our shared data structure. Here is a brief summary of each object
	-  @TrajectoryData
  	-  @KinematicData
  	-  @VRCursorData
  	-  @GPFAData
  	-  @EL_ExperimentInfo
  	-  @IntTargExp
  	-  @TubeObject
- +db: Database code for organizing and subselecting data.
- +flow: Code related to flow field analysis.
- +opt: Optimization code that identifies SepMax projections.
- +plt: General plotting functions that are generalized and used across a number of projects.
- +SDT: Functions for computing d' and ROC curves.
- +tube: Contains code for intermediate target (IT) task, instructed path, and initial angle calculations.
- +util: Utility functions built specifically for this project.
- figure_saving: Contains code for saving figures as PDFs.
- gfap_v0203: This contains the modified version of DataHigh([website](https://users.ece.cmu.edu/~byronyu/software/DataHigh/datahigh.html) and [git](https://github.com/BenjoCowley/DataHigh)) data we used for the analysis.

## Requirements
- [Signal Processing Toolbox]
- [Statistics and Machine Learning Toolbox]
- [invToeplitzFastZohar.c]
- [makePautoSum.c]
MATLAB 2024a 

## Installation
Clone this repository:
bash
git clone [https://github.com/BatistaLabCode/DynamicalConstraints](https://github.com/BatistaLabCode/DynamicalConstraints)

## Authorship and Contributions
This code package was developed by [Erinn Grigsby](mailto:erinn.grigsby@gmail.com) and [Alan Degenhart](mailto:alan.degenhart@gmail.com)

## Acknowledgments
We would like to thank the [Batista Lab](https://smile.pitt.edu) for their support in this project.

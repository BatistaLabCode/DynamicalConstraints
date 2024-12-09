
# **Dynamical Constraints Code Package**  
*(Published in Nature Neuroscience, 2025)*

## **Overview**
This repository contains the data and code used for the analysis and figure generation in our recently published paper, [Dynamical constraints on neural population activity]. 

1. The **data** provided will include:
	-Processed trial data from example session (E20190719)
	-Processed trial data of the two-target task from a single experiment session(E20181004)
	-Flow field data from all valid sessions
	-Summary source data for all statistical analysis in the paper. 
		
2. The **code** provided will include:
	- All the standard functions we used to process, analyze, and plot the data (*See folder organization description of function types*)
	- Simplified scripts that will reproduce the summary data values and information for the example session, E20190719.
	- Scripts and functions that produce all the plots found in the main figures.

---

## **Requirements**
- MATLAB 2024a
- Toolboxes:
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - DataHigh (*see note below*)

> **Note**: 
> - We provide a slightly modified version of DataHigh within the code pack, including all its external dependencies. Please keep in mind if you already have DataHigh installed on your computer.
> - The `dataHigh` function `orthogonalize` relies on MATLAB's `svd` function, which may produce slight visual discrepancies in Figure 3 panels or SepMax projections, depending on MATLAB version and machine. See MATLAB’s [SVD documentation](https://www.mathworks.com/help/releases/R2024b/matlab/ref/double.svd.html?overload=svd+false#bu2_0hq-3) for details. 
> - This `svd` function variability will not cause any issue with identifying the Separation Maximizing projection from scratch or new data. It simply means that the orthogonal vector may have sign flips that change matrix **M** weights, but not the relationship of the trajectories plotted together.

---

## **Installation**
Clone the repository:  
```bash
git clone https://github.com/BatistaLabCode/DynamicalConstraints
```

---

## **Getting Started**

To run the main analysis, you will primarily interact with the following functions:

- **`DynamicalConstraints\serverPath.m`**: location of data folder on the computer. This is a function that will attempt to identify the location of the data information, the example session information, and where the user wants to save the material. If the code is unable to find either the automated path or the user input path, then the code will prompt the user to identify the correct folders using the file GUI.

- **`DynamicalConstraints\createAllFigures.m`**: produces all the main paper figures and their related statistical tests. Note that the serverPath and createAllFigures will both the GitHub path and all subfolders to MATLAB, but not all provided codes do this.

> **Note**:   
> - There are slight discrepancies in the trajectory arrows shown in the paper and the arrows plotted using `createAllFigures.m` script. This is because the trajectory arrows shown in the main figures of the paper were added using Adobe Illustrator, which calculates the direction of the arrowhead differently than the Matlab code.

### **Step 1: Configure Paths**  
Edit `serverPath.m` to update or confirm the following:
- **`dataLoc`**: Data folder location.
- **`exSessDataLoc`**: Example session folder.
- **`saveFigLoc`**: Save location for figures.

*The default settings of this code assume that the data folder is still stored in the same location as the main GitHub folder.*

### **Step 2: Verify `serverPath.m`**  
Run `serverPath.m` to ensure paths are correctly configured. The code will prompt for updates if paths are invalid.

### **Step 3: Run Scripts**  
- To generate all figures:  
  ```matlab
  createAllFigures
  ```
- To test single example analysis scripts:  
  ```matlab
  exampleAnalysisScripts\identify_separationMaximizingProjection
  exampleAnalysisScripts\dPrime_calc
  exampleAnalysisScripts\dPrime_calc_EL_shuffle
  exampleAnalysisScripts\flowField_calc
  exampleAnalysisScripts\calc_initialAngle_data
  exampleAnalysisScripts\calc_centerOutAngle_control
  ```
  
> **Note**: 
> - It is possible to also run the individual figure codes separately, e.g. `fig_2_two_target(dataLoc)`, or adjust optional inputs/variables in the functions.  
 ---

## **File Organization**

### **1. Data**:  
Located in: `DynamicalConstraints\DynamicalConstraints_NatNeuro_2024_data\`  
Key files and folders:
- **`exampleDatasetCatalog`**: Structure with details of the provided example sessions.
- **`dPrime_workspace+neuralSpace_fig3`**: Structure with the discriminability metric (d') for all sessions. This metric is calculated with the MoveInt trials projected into both the MoveInt workspace and the SepMax workspace. (used in Figure 3).
- **`dPrime_EarlyLate+ShuffleData`**: d' values of all sessions for early/late trials and shuffled data (used in Figure 4).
- **`initialAng_tt_vs_uncon+con_Simplified`**: Initial angle data for all sessions where the animal performed IT task, instructed path task, and the two-target task (used in Figures 6 and 7). The structure *AD_compare* includes both the individual trials' values and averages for the initial angle. These values were calculated in the SepMax space (used in Figures 6 and 7).
- **`centerOut_AngleAnalatsisData_Signed`**: Initial angle data for full change control. Calculated for the same session as the other initial angle structure. The structure is *AD_CO* (used in Figures 6 and 7).
>*Note: these were calculated using only the last 32 trials of the center out block since those were the only trials where the animals had full control of the cursor*

Subfolders:
- **`Example Sessions\`**: Folder contains the processed data, latent space parameters, and BCI decoders for two example sessions (`monkeyE20181004` and `monkeyE20190719`). The neural and behavior data is saved as a TrajectoryData object (lab-custom design object included in the code pact) and everything else is saved as structures. Each unique session task or parameter change is save as a unique folder (ie monkeyE20190719_04_condGridTask_01_SI_trajectoryData or monkeyE20190719_13_twoTargetABBA_intTarg_startTarg_T1_tube_04_SI_trajectoryData). Example session `monkeyE20190719` is used for all example analysis scripts (*see below*).  
- **`flowAnalysis\mat\`**: FlowResult structures for each valid experiment with each file using the following naming structure: [animalID][session date]_[MoveInt two target folder number] _[SepMax two target folder number]_FlowResults. See the code `flow.main_analysis` for a greater description of the structure.
- **`ConstrainedPath\mat\int_targ_data\`**: Intermediate target (intTarg) data for `monkeyE20190719`. This data structure is used in Figures 6 and 7 and stores all the data relevant to the IT task analysis and the instructed path task analysis.

---

### **2. Example Analysis Scripts**:  
There are six example analysis scripts that we provided. These a simplified versions of many of the "batch" analysis scripts provided. They are also base code in case the user would like to implement any of our methods on their data.
Located in: `DynamicalConstraints\exampleAnalysisScripts\`  
Scripts provided:
- `identify_separationMaximizingProjection`
- `dPrime_calc`
- `dPrime_calc_EL_shuffle`
- `flowField_calc`
- `calc_initialAngle_data`
- `calc_centerOutAngle_control`

---

### **3. Manuscript Figures**: 
Six functions will produce the main figures of the paper.
Located in: `DynamicalConstraints\`  
Functions to recreate figures:
1. `fig_2_two_target`
2. `fig_3`
3. `fig_4_early_vs_late_two_target_trajectories_comparisons`
4. `fig_4_plot_View_vs_Projections`
5. `fig_5_plot_flow_field_analysis`
6. `fig_6_fig_7`

To generate all figures and statistical tests, run the wrapper script:  
```matlab
createAllPaperFigures
```
In the wrapper script, there is the option to save the figures and even subselect the desired plots (ie only plot figures 3 and 5). For other functions, there are options to adjust and save figures are included in each function's header.

---

## **Folder Structure**
The repository is organized into several folders, each with a specific purpose or anaylsis focus. Every function has a header with a detailed description of its use, required inputs, optional inputs, and outputs.

| Folder          | Description                                                                                       |
|----------------------|-------------------------------------------------------------------------------------------------|
| `@`             | Custom MATLAB classes built specifically for projects in the Batista lab for object-specific analysis. The advantage of these classes is that they enable users to create object specific analysis functions (An example of this is the TubeOject.plot and TrajectoryData.plot are different functions that will produce different results). While they could be replaced with structures in the future, they are necessary for several minor codes.|
| `@TrajectoryData`		| Object with data for an individual trial. These are often combined to create structures of whole experiment conditions.|
| `@KinematicData`		| Object for the kinematic source data that control the computer cursor. The default for this experiment is BCI control unless stated explicitly. |
| `@VRCursorData`		| Object with transformation information to convert kinematic data into the visual workspace. |
| `@GPFAData`		| Object of the GPFA and neural space data. This is the data that is used to calculate the *MoveInt* and the *SepMax* workspace.|
| `@EL_ExperimentInfo`	| Object stores the experiment session information and specifically the information related to the decoder for the *SepMax* information.
| `@IntTargExp`		| Object the information and data relevant to the intermediate target and initial angle analysis.
| `@TubeObject`	| Object related to the instructed path and its display on the screen.
| `+db`           | Database management and data sub-selection scripts.                                               |
| `+flow`         | Flow field analysis scripts and functions.                                                                      |
| `+opt`          | Code for identifying and optimizing the SepMax projection.                                                          |
| `+plt`          | Generalized plotting functions.                                                                   |
| `+SDT`          | Functions for computing d' and ROC curves.                                                        |
| `+tube`         | Contains analysis and plotting code for Intermediate target (IT) task and instructed path task as well as initial angle calculations.                        |
| `+util`         | Utility functions built specifically for this project.                                                                  |
| `figure_saving` | Functions for saving figures as PDFs.                                                             |
| `gpfa_v0203`    | Modified version of [DataHigh](https://github.com/BenjoCowley/DataHigh) that was used for all experiments and analysis.                          |

---

## **Authorship**
This code package was developed by:  
- **Erinn Grigsby** ([erinn.grigsby@gmail.com](mailto:erinn.grigsby@gmail.com))  
- **Alan Degenhart** ([alan.degenhart@gmail.com](mailto:alan.degenhart@gmail.com))  

---

## **Acknowledgments**
We thank the [Batista Lab](https://smile.pitt.edu) for their invaluable support and collaboration

---

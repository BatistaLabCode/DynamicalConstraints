# Code Package for Dynamical Constraints (Nature Neuroscience, 2025)

## Overview
This repository contains the code used for the analysis and figure generation in our recently published paper, [Paper Title]. The primary functions users will interact with are:

`serverPath.m` function: This is a function that will attempt to identify the location of the data information, the example session information, and where the user wants to save the material. If the code is unable to find either the automated path or the user input path, than the code will prompt the user to identify the correct folders using the file GUI.

`createAllFigures` function: This script produces all the main paper figures and their related statistical tests. The user can adjust the code to only produce a subselection of the main figure plots. This will also add paths for all the codes into MATLAB.

Additionally, the repository is organized into several folders, each with a specific purpose. Please refer to the "File Organization" section below for a quick overview.

## File Organization
- @ folders: These are custom MATLAB classes built specifically for projects in the Batista lab. The advantage of these classes is that they enable users to create While they could be replaced with structures in the future, they are necessary for several minor codes and will remain until we finalize our shared data structure.
- gfap_v0203: This contains the modified version of DataHigh data we used for the analysis.
- figure_saving: Contains code for saving figures as PDFs.
- +plt: General plotting functions that generalized used across a number of projects.
- +util: Utility functions built specifically for this project.
- +tube: Contains code for constrained path, instructed path, and initial angle calculations.
- +SDT: Functions for computing d' and ROC curves.
- +opt: Optimization code that identifies SepMax projections.
- +flow: Code related to flow field analysis.
- +db: Database code for organizing and subselecting data.

## Requirements
[List software or libraries required to run the code]
MATLAB 2024a

## Installation
Clone this repository:
bash
git clone https://github.com/BatistaLabCode/EnergyLandscape

## Usage
To run the main analysis, you will primarily interact with the following functions:

serverPath: location of data folder on the computer.
	Code call: [dataLoc, exSessDataLoc, saveFigLoc] = serverPath
createAllFigures: Simple script just hit run.

It is possible to also run the individual figure codes separately.
E.g. fig_2_two_target(dataLoc)

## Authorship and Contributions
This code package was developed by [Erinn Grigsby](mailto:erinn.grigsby@gmail.com) and [Alan Degenhart](mailto:alan.degenhart@gmail.com)

## Acknowledgments
We would like to thank the [Batista Lab](https://smile.pitt.edu) for their support in this project.

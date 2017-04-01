# MicroarrayGeneClustering
Shell scripts and C programs for DNA microarray calibration and analysis with gene clustering.

Author: Jay Battle
Course: BIEN 4290
Refrences: My alma mater, Chu et. al.
Languages: C, Bash
Platform: UNIX
Purpose: Sample C Project
Status: In revisions

Description:
This class inspired project creates a automated command line interface to calibrate, analyze, 
and cluster genes from DNA microarray data files. This project consists of C programs, Shell
scripts, and a Makefile. A sample data set from the paper Chu et. Al is included for testing.
The output of the project is a series of text files containing the expressed, stationary, and 
supressed genes for each data set.

Instructions:
Step 1: Make sure you have a copy of the entire MicroarrayGeneClustering directory
Step 2: Run bash setup.sh to set all file permissions
Step 3: Run bash calibrationAndAnalysis.sh to compile all programs and generate cluster data
	, follow any prompts
Step 4: Run bash dataRetrival.sh to search a gene's status for all times, follow any prompts

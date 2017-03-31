#!/bin/bash
# Jay Battle
# BIEN 4290 -
# Febuary 27, 2017
# Lab 3
# Program Description:
# This program does calibration and analysis on a large dataset. 
# compile with -lm
#
rm MicroarrayData/supressed_genes.txt
rm MicroarrayData/stationary_genes.txt
rm MicroarrayData/expressed_genes.txt
rm MicroarrayData/Summary.txt
make
FileCount=7
SupressedClusterMean=-.5
StationaryClusterMean=0
ExpressedClusterMean=.5
Folder="MicroarrayData/"
for ((FileNum=0;FileNum<FileCount;FileNum++))
{
	#rm log_ratio_input.dat
	./preProcessing ${Folder}red_${FileNum}.dat ${Folder}red_background_${FileNum}.dat ${Folder}green_${FileNum}.dat ${Folder}green_background_${FileNum}.dat ${Folder}log_ratio_${FileNum}.dat 6118
	cp ${Folder}log_ratio_${FileNum}.dat ${Folder}log_ratio_input.dat
	echo "For the ${FileNum} time point dataset" >> ${Folder}Summary.txt
	./clustering ${SupressedClusterMean} ${StationaryClusterMean} ${ExpressedClusterMean} >> ${Folder}Summary.txt
	cp ${Folder}supressed_genes.txt ${Folder}supressed_genes_${FileNum}.txt
	cp ${Folder}stationary_genes.txt ${Folder}stationary_genes_${FileNum}.txt
	cp ${Folder}expressed_genes.txt ${Folder}expressed_genes_${FileNum}.txt
}
cat ${Folder}Summary.txt
echo "All data is located in the directory ~/${Folder}"

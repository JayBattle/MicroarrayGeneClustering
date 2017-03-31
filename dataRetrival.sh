#!/bin/bash
# Jay Battle
# BIEN 4290 -
# Febuary 27, 2017
# Lab 3
# dataRetrival
# Program Description:
# This program show the status of the specified gene at every time point. 
# compile with -lm
#
rm MicroarrayData/supressed_genes.txt
rm MicroarrayData/stationary_genes.txt
rm MicroarrayData/expressed_genes.txt
FileCount=7
echo "Enter the desired Gene: (ex: YPR199C): "
read DesiredGene;
for ((FileNum=0;FileNum<FileCount;FileNum++))
{
	if grep -q ${DesiredGene} MicroarrayData/supressed_genes_${FileNum}.txt; then
		echo "${DesiredGene} is supressed for time point ${FileNum}"
	fi
	if grep -q ${DesiredGene} MicroarrayData/stationary_genes_${FileNum}.txt; then
		echo "${DesiredGene} is stationary for time point ${FileNum}"
	fi
	if grep -q ${DesiredGene} MicroarrayData/expressed_genes_${FileNum}.txt; then
		echo "${DesiredGene} is expressed for time point ${FileNum}"
	fi
}

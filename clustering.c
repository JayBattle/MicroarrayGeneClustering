/*
Jay Battle
BIEN 4290 -
Febuary 14, 2017
Lab 3
clustering.c
Program Description:
This program does clustering for microarray data
remember to compile with -lm 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector_ops.h"
#include "stats.h"

float* LoadDataFile(char*, int**);
void WriteDataFile(char*, float*, int);
void WriteTextFile(char *, char*[], int);

int debug = 0; 

int main(int argc, char *argv[])  {
	
	float *logVector, *geneLocations, *clusterMeans, clusterCriteria;
	float *distancesToSuppressedCluster, *distancesToStationaryCluster, *distancesToExpressedCluster;
	float suppressedGeneClusterMean,stationaryGeneClusterMean,expressedGeneClusterMean;
	float newSuppressedGeneClusterMean,newStationaryGeneClusterMean,newExpressedGeneClusterMean;
	float *suppressedGeneCluster, *stationaryGeneCluster, *expressedGeneCluster;
	int *logVectorCount, dataCount, *geneCount;
	int suppressedCount, stationaryCount, expressedCount;
	FILE *outputFile,*inputFile;
	char *geneNames[6118],*geneNamesPointer;
	char *suppressedGeneClusterByName[6118], *stationaryGeneClusterByName[6118], *expressedGeneClusterByName[6118];

	
	if( argc > 4 ) {
      printf("Too many arguments supplied.\n");
	  return 0;
	} else if (argc < 4 ) {
      printf("Too few arguments supplied.\n");
	  return 0;
	} else if( argc == 4 ) {
      if (debug==1) printf("The mean of the suppressed gene cluster is: %s\n", argv[1]); //for debug
	  if (debug==1) printf("The mean of the stationary gene cluster is: %s\n", argv[2]); //for debug
	  if (debug==1) printf("The mean of the expressed gene cluster is: %s\n", argv[3]); //for debug
	}
	
	suppressedGeneClusterMean = (float) atof(argv[1]);
	stationaryGeneClusterMean = (float) atof(argv[2]);
	expressedGeneClusterMean = (float) atof(argv[3]);
		
	logVector = LoadDataFile("MicroarrayData/log_ratio_input.dat", &logVectorCount);
	dataCount = logVectorCount;
	
	char currLine[12];
	int lineCount, numOfLines;	
	inputFile = fopen("MicroarrayData/gene_list.txt", "r");
	if (!(inputFile)){ 
        printf("MicroarrayData/gene_list.txt does not exist.\n");
        exit(0); 
	}	
	int j = 0;
	while(fgets(currLine, 12, inputFile)) {
		if (debug==2) printf("%s", currLine); //for debug
		geneNames[j] = strdup(currLine);
		if (debug==2) printf("%s", geneNames[j]); //for debug
		j++;
	}
	fclose(inputFile); 

	distancesToSuppressedCluster = (float *) malloc(dataCount*sizeof(float));
	distancesToStationaryCluster = (float *) malloc(dataCount*sizeof(float));	
	distancesToExpressedCluster = (float *) malloc(dataCount*sizeof(float));
	geneLocations = (float *) malloc(dataCount*sizeof(float));
	suppressedGeneCluster = (float *) malloc(dataCount*sizeof(float));
	stationaryGeneCluster = (float *) malloc(dataCount*sizeof(float));
	expressedGeneCluster = (float *) malloc(dataCount*sizeof(float));
	
	clusterCriteria = 1.0;
	while (clusterCriteria > 0.0001) {
		if (debug==1) printf("Cluster Recalculation Required! \n"); //for debug
		for (int i = 0; i < dataCount; i++) { 
			distancesToSuppressedCluster[i] = logVector[i];
			distancesToStationaryCluster[i] = logVector[i];
			distancesToExpressedCluster[i] = logVector[i];
			geneLocations[i] = 0;
			suppressedGeneCluster[i] = 0;
			stationaryGeneCluster[i] = 0;
			expressedGeneCluster[i] = 0;
		}	
		
		SubtractVector(distancesToSuppressedCluster, dataCount, suppressedGeneClusterMean);
		SubtractVector(distancesToStationaryCluster, dataCount, stationaryGeneClusterMean);
		SubtractVector(distancesToExpressedCluster, dataCount, expressedGeneClusterMean);
		
		for (int i = 0; i < dataCount; i++) { 
			float closestCluster = 1;
			float curCluster = fabs(distancesToSuppressedCluster[i]);
			float nextCluster = fabs(distancesToStationaryCluster[i]);
			if (curCluster > nextCluster){
				closestCluster=2;
				curCluster = nextCluster;
			}
			nextCluster =fabs(distancesToExpressedCluster[i]);
			if (curCluster > nextCluster){
				closestCluster=3;
				curCluster = nextCluster;
			}
			geneLocations[i] = closestCluster;
		}
		
		suppressedCount = 0;
		stationaryCount = 0;
		expressedCount = 0;
		for (int i = 0; i < dataCount; i++) {	
			if (geneLocations[i] == 1){
				if (debug==2) printf("Gene 1: %s", geneNames[i]); //for debug
				suppressedGeneCluster[suppressedCount] = logVector[i];
				suppressedGeneClusterByName[suppressedCount] = geneNames[i];
				suppressedCount++;
			} else if (geneLocations[i] == 2){
				if (debug==2) printf("Gene 2: %s", geneNames[i]); //for debug
				stationaryGeneCluster[stationaryCount] = logVector[i];
				stationaryGeneClusterByName[stationaryCount] = geneNames[i];
				stationaryCount++;
			} else {
				if (debug==2) printf("Gene 3: %s", geneNames[i]); //for debug
				expressedGeneCluster[expressedCount] = logVector[i];
				expressedGeneClusterByName[expressedCount] = geneNames[i];
				expressedCount++;
			}
		}	
		
		newSuppressedGeneClusterMean = CalcMean(suppressedGeneCluster, suppressedCount);
		newStationaryGeneClusterMean = CalcMean(stationaryGeneCluster, stationaryCount);
		newExpressedGeneClusterMean = CalcMean(expressedGeneCluster, expressedCount);

		clusterCriteria = fabs(suppressedGeneClusterMean - newSuppressedGeneClusterMean);
		clusterCriteria += fabs(stationaryGeneClusterMean - newStationaryGeneClusterMean);
		clusterCriteria += fabs(expressedGeneClusterMean - newExpressedGeneClusterMean);

		suppressedGeneClusterMean = newSuppressedGeneClusterMean;
		stationaryGeneClusterMean = newStationaryGeneClusterMean;
		expressedGeneClusterMean = newExpressedGeneClusterMean;
		
		if (debug==1) printf("newsuppressedGeneClusterMean: %f \n", suppressedGeneClusterMean); //for debug		
		if (debug==1) printf("newStationaryGeneClusterMean: %f \n", stationaryGeneClusterMean); //for debug		
		if (debug==1) printf("newExpressedGeneClusterMean: %f \n", expressedGeneClusterMean); //for debug				
		if (debug==1) printf("Cluster Criteria: %f \n", clusterCriteria); //for debug		
	}
	
	clusterMeans = (float *) malloc(3*sizeof(float));
	clusterMeans[0] = suppressedGeneClusterMean;
	clusterMeans[1] = stationaryGeneClusterMean;
	clusterMeans[2] = expressedGeneClusterMean;
	
	float minCluster = GetMinimum(clusterMeans, 3);	
	if (minCluster == suppressedGeneClusterMean){
		printf("The suppressed cluster of %d genes has a mean of %f \n",suppressedCount, minCluster);
		WriteTextFile("MicroarrayData/supressed_genes.txt", suppressedGeneClusterByName, suppressedCount);
	} else if (minCluster == stationaryGeneClusterMean){
		printf("The suppressed cluster of %d genes has a mean of %f \n",stationaryCount, minCluster);
		WriteTextFile("MicroarrayData/supressed_genes.txt",stationaryGeneClusterByName, stationaryCount);
	} else {
		printf("The suppressed cluster of %d genes has a mean of %f \n",expressedCount, minCluster);
		WriteTextFile("MicroarrayData/supressed_genes.txt",expressedGeneClusterByName, expressedCount);
	}
	
	float maxCluster = GetMaximum(clusterMeans, 3);	
	if (maxCluster == suppressedGeneClusterMean){
		printf("The expressed cluster of %d genes has a mean of %f \n",suppressedCount, maxCluster);
		WriteTextFile("MicroarrayData/expressed_genes.txt",suppressedGeneClusterByName, suppressedCount);
	} else if (maxCluster == stationaryGeneClusterMean){
		printf("The expressed cluster of %d genes has a mean of %f \n",stationaryCount, maxCluster);
		WriteTextFile("MicroarrayData/expressed_genes.txt",stationaryGeneClusterByName, stationaryCount);
	} else {
		printf("The expressed cluster of %d genes has a mean of %f \n",expressedCount, maxCluster);
		WriteTextFile("MicroarrayData/expressed_genes.txt",expressedGeneClusterByName, expressedCount);
	}
	
	float midCluster = 0;
	for (int i = 0; i < 3; i++) {
		if (clusterMeans[i] != maxCluster && clusterMeans[i] != minCluster){
			midCluster = clusterMeans[i];
		}
	}
	if (midCluster == suppressedGeneClusterMean){
		printf("The stationary cluster of %d genes has a mean of %f \n",suppressedCount, midCluster);
		WriteTextFile("MicroarrayData/stationary_genes.txt",suppressedGeneClusterByName, suppressedCount);
	} else if (midCluster == stationaryGeneClusterMean){
		printf("The stationary cluster of %d genes has a mean of %f \n",stationaryCount, midCluster);
		WriteTextFile("MicroarrayData/stationary_genes.txt",stationaryGeneClusterByName, stationaryCount);
	} else {
		printf("The stationary cluster of %d genes has a mean of %f \n",expressedCount, midCluster);
		WriteTextFile("MicroarrayData/stationary_genes.txt",expressedGeneClusterByName, expressedCount);
	}
	
	free(distancesToSuppressedCluster);
	free(distancesToStationaryCluster);
	free(distancesToExpressedCluster);
	free(geneLocations);
	free(suppressedGeneCluster);
	free(stationaryGeneCluster);
	free(expressedGeneCluster);
	free(clusterMeans);
	return 0;
}

////Method: LoadDataFile
////Description: This method loads the specified data file
////Args: char *inputFileName, int dataCountPointer
////Returns: dataVector, float*, changes the value of dataCountPointer
float* LoadDataFile(char *inputFileName, int **dataCountPointer) {
	
	FILE *inputFile;
	float *dataVector, dataPoint;
	int dataCount, numOfDataPoints;
	
	inputFile = fopen(inputFileName, "r");
	if (!(inputFile)){ 
        printf("%s does not exist.\n", inputFileName);
        exit(0); 
	}
	
	// //User prompt for Number of data points with checks
	// //commented out for command line and shell script use
	// numOfDataPoints = 1;
	// while (numOfDataPoints != dataCount) {
		// printf("Enter the number of data points in the file (Ex: 500) : \n");
		// scanf("%d", &dataCount);
		// numOfDataPoints = 0;
		// while ((fscanf(inputFile, "%f", &dataPoint)) != EOF) {
			// numOfDataPoints++; 
		// }		
		// rewind(inputFile);
			// printf("The number of data points in the file is actually : %d \n", numOfDataPoints); 
	// }
	
	dataCount = 0;
	while ((fscanf(inputFile, "%f", &dataPoint)) != EOF) {
		dataCount++; 
	}		
	rewind(inputFile);
	if (debug==1) printf("The number of data points in the file is: %d \n", dataCount); //for debug

	dataVector = (float *) malloc(dataCount*sizeof(float));	
	*dataCountPointer = dataCount;
	
	int j = 0;
	while ((fscanf(inputFile, "%f", &dataPoint)) != EOF) {
		dataVector[j] = dataPoint;
		j++;
	}
	fclose(inputFile); 
	return dataVector;
}

////Method: WriteDataFile
////Description: This method writes the specified data to the specified data file
////Args: char *outputFileName, float *dataVector, int dataCount)
////Returns: void
void WriteDataFile(char *outputFileName, float *dataVector, int dataCount) {
	
	FILE *outputFile;	
	outputFile = fopen(outputFileName, "w"); 
    for (int i = 0; i < dataCount; i++) { 
        fprintf(outputFile, "%f \n", dataVector[i]);
    }
    fclose(outputFile); //close file
    if (debug==1) printf("Data Saved to %s \n", outputFileName); //for debug
	
}

////Method: WriteTextFile
////Description: This method writes the specified text to the specified text file
////Args: char *outputFileName, float *textVector, int lineCount)
////Returns: void
void WriteTextFile(char *outputFileName, char *textVector[], int lineCount) {
	
	FILE *outputFile;	
	outputFile = fopen(outputFileName, "w"); 
    for (int i = 0; i < lineCount; i++) { 
		//if (debug==1) printf("%s", textVector[i]); //for debug
        fprintf(outputFile, "%s", textVector[i]);
    }
    fclose(outputFile); //close file
    if (debug==1) printf("Data Saved to %s \n", outputFileName); //for debug
	
}
/*
Jalen Battle
BIEN 4290 -
Febuary 14, 2017
Lab 3
preProcessing.c
Program Description:
This program does preprocessing for microarray data
remember to compile with -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector_ops.h"
#include "stats.h"

float* LoadDataFile(char*, int**);

int debug = 0; 

int main(int argc, char *argv[])  {
	
	float *redVector, *redBackground, *greenVector, *greenBackground, geneNum;
	int *redVectorCount, *redBackgroundCount, *greenVectorCount, *greenBackgroundCount;
	char *outputFileName;
	float redMean, greenMean;
	FILE *outputFile;
	
	
	if( argc > 7 ) {
      printf("Too many arguments supplied.\n");
	  return 0;
	} else if (argc < 7 ) {
      printf("Too few arguments supplied.\n");
	  return 0;
	} else if( argc == 7 ) {
      if (debug==1) printf("The name of the file containing (red) sporulating cells is: %s\n", argv[1]); //for debug
	  if (debug==1) printf("The name of the file containing (red) background data is: %s\n", argv[2]); //for debug
	  if (debug==1) printf("The name of the file containing (green) non-sporulating cells is: %s\n", argv[3]); //for debug
	  if (debug==1) printf("The name of the file containing (green) background data is: %s\n", argv[4]); //for debug
	  if (debug==1) printf("The name of the file output file for calibrated data is: %s\n", argv[5]); //for debug
	  if (debug==1) printf("The number of genes to be analyzed is: %s\n", argv[6]); //for debug
	}
   
	redVector = LoadDataFile(argv[1], &redVectorCount);
	redBackground = LoadDataFile(argv[2], &redBackgroundCount);
	greenVector = LoadDataFile(argv[3], &greenVectorCount);
	greenBackground = LoadDataFile(argv[4], &greenBackgroundCount);
	outputFileName = argv[5];
	geneNum = atoi(argv[6]);
		
	SubtractVectors(redVector, redBackground, redVectorCount, redBackgroundCount);
	SubtractVectors(greenVector, greenBackground, greenVectorCount, greenBackgroundCount);
	redMean = CalcMean(redVector, redVectorCount);
	if (debug==1) printf("The red mean is: %f\n", redMean);
	greenMean = CalcMean(greenVector, greenVectorCount);
	if (debug==1) printf("The green mean is: %f\n", greenMean);
	DivideVector(redVector, redVectorCount, redMean);
	DivideVector(greenVector, greenVectorCount, greenMean);
	DivideVectors(redVector, greenVector, redVectorCount, greenVectorCount);
	logVector(redVector, redVectorCount, 10);
		
    outputFile = fopen(outputFileName, "w");
    for (int i = 0; i < redVectorCount; i++) { 
        fprintf(outputFile, "%f \n", redVector[i]);
    }
    fclose(outputFile); //close file
    if (debug==1) printf("Data Saved to: %s\n", outputFileName);
    
	return 0;
}

////Method: LoadDataFile
////Description: This method loads the specified data file
////Args: char *inputFileName, int dataCountPointer
////Returns: dataVector, float*, changes the value of dataCountPointer
float* LoadDataFile(char *inputFileName, int **dataCountPointer) {
	
	FILE *inputFile;
	float *dataVector, dataPoint;
	int dataCount;
	
	inputFile = fopen(inputFileName, "r");
	if (!(inputFile)){ 
        printf("%s does not exist.\n", inputFileName);
        exit(0); 
	}
	
	dataCount = 0;
	while ((fscanf(inputFile, "%f", &dataPoint)) != EOF) {
		dataCount++; 
	}		
	rewind(inputFile);
	if (debug==1) printf("The number of data points in the file 1 is: %d \n", dataCount); //for debug

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
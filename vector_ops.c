/*
Jay Battle
BIEN 4290 -
Febuary 14, 2017
Lab 3
vector_ops.c
Program Description:
This program contains some vector math functions
remember to compile with -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector_ops.h"

float SumVector(float*, int);
void SubtractVector(float*, int, float);
void SubtractVectors(float*, float*, int, int);
void DivideVector(float*, int, float);
void DivideVectors(float*, float*, int, int);
void logVector(float*, int, float);

 
////Method: SumVector
////Description: This method adds all of the datapoints in a vector and returns the sum
////Args: float *dataVector, int dataCount
////Returns: sum, float
float SumVector(float *dataVector, int dataCount) {
	float sum = 0;
	for (int i = 0; i < dataCount; i++) {
		sum += dataVector[i]; 
	}
	return sum;
}

////Method: SubtractVector
////Description: This method divides all of the datapoints in a vector by a constant
////Args: float *dataVector, int dataCount, float divisor
////Returns: void, dataVector is altered
void SubtractVector(float *dataVector, int dataCount, float subtractor) {
	for (int i = 0; i < dataCount; i++) {
		dataVector[i] = (dataVector[i] - subtractor); 
	}
}

////Method: SubtractVectors
////Description: This method subrtacts vector A from B element by element
////Args: float *vectorA, float *vectorB, int dataCountA, int dataCountB
////Returns: Void, vectorA is altered
void SubtractVectors(float *vectorA, float *vectorB, int dataCountA, int dataCountB) {
	if (dataCountA<=dataCountB){
		for (int i = 0; i < dataCountA; i++) {
			vectorA[i] = (vectorA[i] - vectorB[i]);
		}
	}
}

////Method: divideVector
////Description: This method divides all of the datapoints in a vector by a constant
////Args: float *dataVector, int dataCount, float divisor
////Returns: void, dataVector is altered
void DivideVector(float *dataVector, int dataCount, float divisor) {
	for (int i = 0; i < dataCount; i++) {
		dataVector[i] = (dataVector[i]/divisor); 
	}
}

////Method: DivideVectors
////Description: This method Divides vector A from B element by element
////Args: float *vectorA, float *vectorB, int dataCountA, int dataCountB
////Returns: Void, vectorA is altered
void DivideVectors(float *vectorA, float *vectorB, int dataCountA, int dataCountB) {
	if (dataCountA<=dataCountB){
		for (int i = 0; i < dataCountA; i++) {
			vectorA[i] = (vectorA[i] / vectorB[i]);
		}
	}
}

////Method: logVector
////Description: This method takes the log all of the datapoints in a vector by specified base
////Args: float *dataVector, int dataCount, float base
////Returns: void, dataVector is altered
void logVector(float *dataVector, int dataCount, float base) {
	for (int i = 0; i < dataCount; i++) {
		dataVector[i] = log(dataVector[i])/log(base);
	}
}



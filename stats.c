/*
Jay Battle
BIEN 4290 -
Febuary 14, 2017
Lab 3
stats.c
Program Description:
This program does vector based statistical analaysis
remember to compile with -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector_ops.h"

float CalcMean(float*, int);
float GetMinimum(float*, int);
float GetMaximum(float*, int);
 
////Method: CalcMean
////Description: This method calculates the mean
////Args: float dataSeries[], int dataCount
////Returns: mean, float
float CalcMean(float *dataVector, int dataCount) {	
	float mean = SumVector(dataVector, dataCount);
	mean = (mean / dataCount);
	return mean;
}

////Method: GetMinimum
////Description: This method finds the minimum
////Args: float* dataSeries, int dataCount
////Returns: min, float
float GetMinimum(float *dataSeries, int dataCount) {
    float min = dataSeries[0];
    for (int i = 0; i < dataCount; i++) { 
        if (min > dataSeries[i]) {
            min = dataSeries[i];
        }
    }
    return min;
}

////Method: GetMaximum
////Description: This method finds the maximum
////Args: float* dataSeries, int dataCount
////Returns: max as a float
float GetMaximum(float *dataSeries, int dataCount) {
    float max = 0;
    for (int i = 0; i < dataCount; i++) {
        if (max < dataSeries[i]) {
            max = dataSeries[i];
        }
    }
    return max;
}

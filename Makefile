#Jay Battle
#BIEN 4290 -
#Febuary 14, 2017
#Lab 3
#Makefile
#this is a Makefile for creating preProcessing and clustering executables

all: preProcessing clustering
	
preProcessing: preProcessing.c stats.c vector_ops.c stats.h vector_ops.h
	gcc -o preProcessing preProcessing.c stats.c vector_ops.c -lm
	
clustering: clustering.c stats.c vector_ops.c stats.h vector_ops.h
	gcc -o clustering clustering.c stats.c vector_ops.c -lm

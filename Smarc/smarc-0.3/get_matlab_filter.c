/*
 * File:   get_matlab_filter.c
 * Author: GabrielSuc
 * Description: affect filter's coefficents created in MATLAB to a new filter usable in c code.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include "get_matlab_filter.h"
 
#define MAXCHAR 10000
int get_matlab_filter(int i) {
   	FILE *fp;
   	char str[MAXCHAR];
   	char buffer[32];
	double bar = 0.0;
	int lines, ch  = 0;

   // char*filename = 
	snprintf(buffer, sizeof(char*) * 32, "/home/gabriel/Documents/End-of-study-Project/SRC_files/PM_Multistage_16_bits/PM_filter%i.txt",i); 
    
	fp = fopen(buffer, "r");
    	if (fp == NULL){
        	printf("Could not open file %s",buffer);
//        	return (double *) 1;
		return 1;
    	}

	while(!feof(fp))
	{
  	ch = fgetc(fp);
  	if(ch == '\n')
  	{
    		lines++;
  	}
	}
    	

	while (fgets(str, MAXCHAR, fp) != NULL)
       // 	printf("%s\n", str);
    	fclose(fp);

	bar = atof(str);

   // 	filter = &bar;

	printf("%d\n",lines);

    	return lines;   // filter;
}



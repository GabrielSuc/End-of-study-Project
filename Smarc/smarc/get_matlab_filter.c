/*
 * File:   get_matlab_filter.c
 * Author: GabrielSuc
 * Description: affect filter's coefficents created in MATLAB to a new filter usable in c code.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "get_matlab_filter.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 1024 

struct getFilter* get_matlab_filter(int i) {
   	FILE *fp;
        
	char buffer[SIZE + 1];
	unsigned int lines = 0;
	char buf[SIZE + 1], lastchar = '\n';
        size_t bytes;
	
	snprintf(buffer, sizeof(char*) * 32, "/home/gabriel/Documents/End-of-study-Project/SRC_files/PM_Multistage_24_bits/PM_filter%i.txt",i+1); 

	fp = fopen(buffer, "r");
    	if (fp == NULL){
        	printf("Could not open file %s",buffer);
        	return NULL;
    	}       
	

	//Counting lines
	while ((bytes = fread(buf, 1, sizeof(buf) - 1, fp))) {
                lastchar = buf[bytes - 1];
                for (char *c = buf; (c = memchr(c, '\n', bytes - (c - buf))); c++) {
                    lines++;
                }
        }

        if (lastchar != '\n') {
                lines++;  /* Count the last line even if it lacks a newline */
        }

	
//	printf("%d\n", lines);  	
//	fflush(stdout);
	fseek(fp, 0, SEEK_SET); /* getting back to first line of fp */

	//Creating filter
	double *filter = (double *) malloc(lines*sizeof(double));


	for (size_t i = 0; i< lines; i++){
		if(fscanf(fp, "%lf\n", filter + i) == 1){
		}//	printf("%.10f\n", *filter+i);	
		else{
			printf("failed to read file.\n");
		}		
	}

                  	
//	printf("%.10f\n", filter[0]);	
	//fflush(stdout);
	fclose(fp);
	
	struct getFilter* getfilter = malloc(sizeof(struct getFilter*));
	
	getfilter->flen = lines;
	getfilter->filter = filter;
	
	return getfilter;
}




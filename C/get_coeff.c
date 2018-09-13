/*
 * File:   get_matlab_filter.c
 * Author: GabrielSuc
 * Description: Read the coefficients made in matlab and stored in a txt file 
 *              Call the function as many times as necessary i.e. in case of multistage
 */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "get_coeff.h"

#define SIZE 1024

//buffer variables
FILE *fp;
char buffer[SIZE + 1]; //= (char *) malloc(SIZE);
unsigned int lines = 0;
char buf[SIZE + 1], lastchar = '\n';
size_t bytes;

struct getFilter* get_matlab_filter(int stage_no){

/* get the coefficients from an external file */
/*if (snprintf(buffer, SIZE, "/home/gabriel/Documents/End-of-study-Project/SRC_files/PM_Multistage_24_bits/PM_filter%d.txt", stage_no+1) >= SIZE) {
  SIZE = SIZE*2;
  printf("Not enough space. Trying %d bytes\n", SIZE);
  free(buffer);
  buffer = malloc(SIZE);

  if (snprintf(buffer, SIZE, "/home/gabriel/Documents/End-of-study-Project/SRC_files/PM_Multistage_24_bits/PM_filter%d.txt", stage_no+1) >= SIZE) {
    printf("Still not enough space. Aborting\n");
    exit(1);
  }
}

else{*/

  snprintf(buffer, sizeof(char*) * SIZE, "/home/gabriel/Documents/End-of-study-Project/SRC_files/PM_Multistage_24_bits/PM_filter%d.txt", stage_no+1);

//}
fp = fopen(buffer, "r");
if (fp == NULL){
  printf("Could not open file %s",buffer);
  exit(1);
}


/* Counting lines */
while ((bytes = fread(buf, 1, sizeof(buf) - 1, fp))) {
  lastchar = buf[bytes - 1];
  for (char *c = buf; (c = memchr(c, '\n', bytes - (c - buf))); c++) {
    lines++;
  }
}

if (lastchar != '\n') {
  lines++;  /* Count the last line even if it lacks a newline */
}


fseek(fp, 0, SEEK_SET); /* getting back to first line of fp */


/* Creating filter */
double *filter = (double *) malloc(lines*sizeof(double));


/* Reading the txt file and storing the coefficients in the temporary taps */
int i = 0;
while(fscanf(fp, "%le\n", &filter[i]) !=EOF){
  //printf("%.20e %d\n", tmp_taps[i], i);
  i++;
}


/* Releasing memory */
fclose(fp);
//free(buffer);

struct getFilter* getfilter = malloc(sizeof(struct getFilter*));

getfilter->flen = lines;
getfilter->filter = filter;

return getfilter;
}



void getFilter_dump(struct getFilter* getfilter){
  
  free(getfilter->filter);
  free(getfilter);

}


/*
 * File:   main.c
 * Author: GabrielSuc
 * Description: Performing Synchronous Sample-Rate Conversion between two known frequencies. 
 *              Choice between Parks-McClellan or Elliptic filters. Multistage and polyphaser 
 *		decomposition for both 
  */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "get_coeff.h"


int main(int argc, char* argv[0]){
  
  struct getFilter* getfilter = malloc(sizeof(struct getFilter*));
  getfilter = get_matlab_filter(1);

  printf("first coefficient %.30e\n", getfilter->filter[0]);

  getFilter_dump(getfilter);	

  return 0;

}


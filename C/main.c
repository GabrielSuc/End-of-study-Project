/*
 * File:        main.c
 * Author:      Gabriel Suc
 * Description: Performaing Synchronous Sample-Rate Conversion between two known frequencies.
 *		Choice between Parks-McClellan or Elliptic filters. Multistage and polyphase
 *              decomposition for both.
 */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "get_coeff.h"
#include "/home/gabriel/Ne10/inc/NE10.h"


int main(int argc, char* argv[0]){

  signed short ss;
  double value;
  FILE* inp = NULL;
  FILE*	wavein;
  FILE* oup = NULL;

  inp = fopen("input.raw","rb");
  oup = fopen("checks.txt","wb");

  ne10_int32_t fftSize = 1024;

/*--------------------- Get Filters Coefficients ---------------------*/

  //Look for begining of the file
  fseek(inp, 0, SEEK_SET);

  struct getFilter* getfilter =  malloc(sizeof(struct getFilter*));
  getfilter = get_matlab_filter(0);

  printf("first coefficient %0.30lf\n", getfilter->filter[0]);

  getFilter_dump(getfilter);


/*-------------------- Process the Input Raw File --------------------*/

  //Look for begining of the file
  fseek(inp, 0, SEEK_SET);

  //Create arrays for left and right channel
  int k, i, j = 0;
  double *raw_data, *raw_data_L, *raw_data_R;
  raw_data = (double *)malloc(i * sizeof(double*));
//  raw_data_L = (double *)malloc(i * sizeof(double*));
//  raw_data_R = (double *)malloc(i * sizeof(double*));
  //memset(raw_data, 0, i * sizeof(char*));

  if (raw_data == NULL){
    printf("Raw_data allocation error\n");
    return 1;
  }

  while (fread(&ss, sizeof(signed short), 1, inp) == 1){
  //Now we have to convert from signed short to double:

    raw_data[i] = 2.0*((double)ss/0xffff);

    i++;

//    value = 2.0*((double)ss/0xffff);

    //count = count + 1;

    //Print the results:
  //  printf("%d\n", count);
  //  printf("%.30lf\n", value);

   // fprintf(oup,"%.30lf", value);
  }

  printf("First raw data value: %.30lf", raw_data[i]);

  free(raw_data);
/*---------------------------- Perform FFT ---------------------------*/

  ne10_fft_r2c_cfg_float32_t cfg = ne10_fft_alloc_r2c_float32_neon(fftSize);

       if (cfg == NULL)
        {
                printf("======ERROR, FFT alloc fails\n");
                return 1;
        }
        printf("cfg allocated\n");

        ne10_fft_cpx_float32_t *fftdstdata1 = (ne10_fft_cpx_float32_t *)malloc((fftSize) * sizeof(ne10_fft_cpx_float32_t));

        //ne10_fft_r2c_1d_float32_c(fftdstdata1, (ne10_float32_t *), cfg);


        printf("fft done\n");


}

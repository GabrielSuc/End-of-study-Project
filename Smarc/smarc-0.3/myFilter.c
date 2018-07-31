//#include <span class="code-string">"stdafx.h"</span>
#include "mat.h" // <span class="code-string">"mat.h"</span>
//#include "matlab.h"//<span class="code-string">"matlab.h"</span>
//#include "matrix.h"
#include "tmwtypes.h"

#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libmatlb.lib")
#pragma comment(lib, "libmmfile.lib")

void main(int argc, char **argv)
{
    MATFile *pmat;
    const char* name=NULL;
    mxArray *pa;
    
    /* open mat file and read it's content */
    pmat = matOpen("~/Document/End-of-study-Project/SRC_files/PM_Multistage_16 bits/PM_filter_stage_1.mat ", "r");
    if (pmat == NULL) 
    {
        printf("Error Opening File: \"%s\"\n", argv[1]);
        return;
    }
    
    /* Read in each array. */
    pa = matGetNextVariable(pmat, &name);
    while (pa!=NULL)
    {
        /*
        * Diagnose array pa
        */
        printf("\nArray %s has %d dimensions.", name, 
               mxGetNumberOfDimensions(pa));
        
        //print matrix elements
        mlfPrintMatrix(pa);
        
        //get next variable
        pa = matGetNextVariable(pmat,&name);
                
        //destroy allocated matrix
        mxDestroyArray(pa);
    }
    
    matClose(pmat);
}

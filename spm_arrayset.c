#include "spm_mex.h"
#include "mex.h"
#include<string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int buflen;
    char *fnc_str;
    const mwSize *d;
    mwSize nd, j, n;
    float *x, *c, *y;

    if (nlhs!=1 || nrhs!=3) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsChar(prhs[0])) mexErrMsgTxt("First argument must be a string.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1])     
       || mxIsSparse(prhs[1]) || !mxIsSingle(prhs[1]))
        mexErrMsgTxt("Second argument must be numeric, real, full and single");

    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2])     
       || mxIsSparse(prhs[2]) || !mxIsSingle(prhs[2]) || !mxIsScalar(prhs[2]))
        mexErrMsgTxt("Third argument must be numeric, real, full, single and scalar");
    
    buflen  = mxGetNumberOfElements(prhs[0]);
    fnc_str = (char *)mxCalloc(buflen+1,sizeof(mxChar));
    mxGetString(prhs[0],fnc_str,buflen+1);

    nd = mxGetNumberOfDimensions(prhs[1]);
    d  = mxGetDimensions(prhs[1]);
    n  = 1;
    for(j=0; j<nd; j++)
        n *= d[j];

    x = (float *)mxGetData(prhs[1]);       
    c = (float *)mxGetData(prhs[2]);
    
    plhs[0] = mxCreateNumericArray(nd,d,mxSINGLE_CLASS,mxREAL);    
    y       = (float *)mxGetData(plhs[0]);
    
    if (!strcmp(fnc_str,"nan"))
    {
        mxFree(fnc_str);
        if (nlhs>0)
        {
            #pragma omp parallel for private(j)
            for(j=0; j<n; j++) 
            {                    
                if (x[j] == x[j]) // this is a bit faster than: if (!mxIsNaN(x[j])), to check for NaN
                    y[j] = x[j];
                else
                    y[j] = c[0];
            }
        }
            
    }    
    else
    {
        mxFree(fnc_str);
        mexErrMsgTxt("Option not recognised.");
    }
}

/* fastranc is a c code MEX-file for the adaptive noise cancellation of 
   residual gradient artifacts.
 
   Author: Rami K. Niazy
   Copyright (c) 2004, University of Oxford
*/

#include "math.h"
#include "mex.h"

void fastranc(double *refs, double *d, int N, double mu,
    double *out, double *y, int veclength) {
double *W;
int i, j;
W = (double *) mxMalloc((N+1) * sizeof(double));
    
    for (i = 0; i <= N; i++) {
        W[i]=0;
    }
    
    for (i = 0; i < veclength; i++) {
        if (i < N) { 
            out[i]=0; y[i]=0;
        } else {
            y[i]=0;
            
            for(j = 0; j <= N; j++) {
                y[i] += W[j]*refs[i-N+j]; 
            }
            
            out[i]=d[i]-y[i];
            
            for (j = 0; j <= N; j++) {
            W[j] += 2*mu*out[i]*refs[i-N+j];
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {
double *refs, *d, *out, *y, mu;
int    N, veclength;

/*  Check for proper number of arguments. */
if (nrhs != 4) 
    mexErrMsgTxt("Four inputs required.");
if (nlhs != 2) 
    mexErrMsgTxt("Two outputs required.");

/* Check to make sure the N input argument is a scalar. */
if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
  mxGetN(prhs[2])*mxGetM(prhs[2]) != 1) 
    mexErrMsgTxt("Input N must be a scalar of type double.");

/* Check to make sure the mu input argument is a scalar. */
if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
  mxGetN(prhs[3])*mxGetM(prhs[3]) != 1) 
    mexErrMsgTxt("Input mu must be a scalar of type double.");

/* Check that refs and d are doulbe */
if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
    mexErrMsgTxt("Inputs must be of type double.");

/* Check to make sure length refs and d is same and are col vectors. */
if (mxGetN(prhs[0]) != 1 ||  mxGetN(prhs[1]) != 1)
    mexErrMsgTxt("Reference and Input data must be column vectors.");

if (mxGetM(prhs[0]) != mxGetM(prhs[1])) 
    mexErrMsgTxt("Reference and Input must be of the same length.");


/* Get inputs */
refs = mxGetPr(prhs[0]);
d    = mxGetPr(prhs[1]);
N = (int) mxGetScalar(prhs[2]);
mu = mxGetScalar(prhs[3]);


/* Get the length of the inputs. */
veclength = mxGetM(prhs[0]);

/* Create an mxArray for the output */
plhs[0] = mxCreateDoubleMatrix(veclength,1, mxREAL);
plhs[1] = mxCreateDoubleMatrix(veclength,1, mxREAL);

/* Set the output pointer to the output matrix. */
out = (double *) mxMalloc(veclength * sizeof(double));
y = (double *) mxMalloc(veclength * sizeof(double));

/* Call the C subroutine. */
fastranc(refs,d,N,mu,out,y,veclength);

/* Assign the data array to the output array */
mxSetPr(plhs[0], out);
mxSetPr(plhs[1], y);
return;
}

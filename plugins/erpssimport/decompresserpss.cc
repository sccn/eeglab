// MEX function to decompress ERPSS data
// compile under matlab by typing
// mex decompresserpss.cc
// -------------------------------------

/*
 *	Jonathan C. Hansen	UCSD ERPL	Mar, 1990
 *  Distributed with the permission of the UCSD ERPSS project team 
 *  (Steve Hillyard and Matt Marlow (mmarlow@ucsd.edu))
 *  Contact persons above for copyright and license on this piece of text
 */

/*
 *  adapted from src/lib/libd/rdf/raw_expand.c
 *	Raw data expansion algorithm - takes data produced by raw_compress.
 *	Input and output areas CANNOT overlap...
 *	Return # shorts of encoded data processed (for error check).
 *
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
  char *x, *i;
  short *y;
  unsigned int m, n;
  char t;
  short *                 inptr;
  short *                outptr;
  int                     nch;
  int                     npts;
  int                     lendian;
  double *			pnch;
  double *			pnpts;
  double *			plendian;
  
  register int            code;           /* code nibble */
  register int            nba;        /* # bits available */
  register unsigned short            is;             /* input short */
  register signed short            os;             /* output short */
  register short *        rip;    /* in ptr */
  register double *        rop;   /* out ptr */
  register double *        rlp;  /* last out ptr */
  register double *        limp;  /* out limit ptr */

/* Check for proper number of arguments.  */
if (nrhs !=4)  {
	mexErrMsgTxt("Must have four input arguments.");}
else if (nlhs !=1) {
	mexErrMsgTxt("Only one output argument allowed.");
	}

/*Output matrix must be a column vector */
n = mxGetNumberOfDimensions(prhs[0]);
if (n > 2) {
	mexErrMsgTxt("Output matrix must be a column vector");}
n = mxGetN(prhs[0]); 
m = mxGetM(prhs[0]); 
if (n > 1) {
        mexErrMsgTxt("Output matrix must be a column vector");}
if (!mxIsUint8(prhs[0]) || mxIsComplex(prhs[0])) {
         mexErrMsgTxt("Input must be real uint8 format.");
        } 
if ((m%2) != 0) {
	 mexErrMsgTxt("Input array must have an even number of elements.");
	}

/* Assign pointers to each input.  */
x = (char *) mxGetPr(prhs[0]);
pnch = mxGetPr(prhs[1]);
pnpts = mxGetPr(prhs[2]);
plendian = mxGetPr(prhs[3]);
nch = *pnch;
npts = *pnpts;
lendian = *plendian;

if ((lendian != 0) && (lendian != 1)) {
	mexErrMsgTxt("lendian must equal 0 (big-endian) or 1 (little-endian)");}

/* Create a matrix for the return argument.  */
plhs[0] = mxCreateDoubleMatrix(nch*npts, 1, mxREAL);
if (plhs[0] == NULL) {
            mexErrMsgTxt("Could not create mxArray.\n");
	} 
/* Assign pointers to output.  */
y = (short *)mxGetPr(plhs[0]);  /* 05/27/05 Petr Janata added cast */

/* Perform byteswap on input variables. */
m = mxGetM(prhs[0]);
if (lendian == 1) {
	for(i=x; i<(x+m); i+=2) {
    	    t = *i;
        	*i = *(i+1);
        	*(i+1) = t;
    }
}

  	inptr = (short *) x;
	outptr = y;
	nba = 0;
	rip = inptr;
	rop = (double *) outptr;  /* 05/27/05 Petr Janata added cast */
	rlp = rop-nch;
	limp = rop + nch*npts;

	do {
		/*
		 *	Get 4 bits of code
		 */
		if ( nba == 0 ) {
			is = *rip++;
			nba = 16;
		}
		nba -= 4;
		code = ( is >> nba ) & 0x000f;

		/*
		 *	ensure >= 4 bits available
		 */
		if ( nba == 0 ) {
			is = *rip++;
			nba = 16;
		}
		/*
		 *	See if 3 bits
		 */
		if ( ( code & 0x0008 ) == 0 ) {
			os = code;
			if ( code & 0x0004 )
				os |= ~0x0003;
		} else
		/*
		 *	See if 6 bits
		 */
		if ( ( code & 0x0004 ) == 0 ) {
			/*
			 *	Get 4 bits - at least 4 available
			 */
			nba -= 4;
			os = (( is >> nba ) & 0x000f) | (( code & 01 ) << 4 );
			if ( code & 0x0002 )
				os |= ~0x001f;
		} else
		/*
		 *	See if 9 bits
		 */
		if ( ( code & 0x0002 ) == 0 ) {
			/*
			 *	Get 8 bits
			 */
			if ( nba >= 8 ) {
				nba -= 8;
				os = ( is >> nba ) & 0x00ff;
			} else {
				os = ( is << 4 ) & 0x00f0;
				is = *rip++;
				nba = 12;
				os |= ( is >> 12 ) & 0x000f;
			}
			if ( code & 0x0001 )
				os |= ~0x00ff;
		} else {
			/*
			 *	Get 12 and set
			 */
			if ( nba >= 12 ) {
				nba -= 12;
				os = ( is >> nba ) & 0x0fff;
			} else
			if ( nba == 8 ) {
				os = ( is << 4 ) & 0x0ff0;
				is = *rip++;
				nba = 12;
				os |= ( is >> 12 ) & 0x000f;
			} else {
				os = ( is << 8 ) & 0x0f00;
				is = *rip++;
				nba = 8;
				os |= ( is >> 8 ) & 0x00ff;
			}
			if ( os & 0x0800 )
				os |= 0xf000;
			*rop++ = os;
			rlp++;
			continue;
		}
		os += *rlp++;
		*rop++ = os;
	} while ( rop < limp );
	return;
}

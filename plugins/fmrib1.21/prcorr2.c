/* MEX-file (Plugin) for MATLAB 
   Created with: ELEMENT ITERATOR Version 1.0 Beta (C) 2001-2002 Peter Rydesäter
   Implements function: prcorr2
*/
#include "mex.h"
#include <math.h>
#include <string.h>

#define MATLAB_COMMENT(X) /* Dummy macro for matlab_comments */

#define _deftype_ double

#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))
#define LIMIT(X,A,B) MIN(MAX(X,B),A)
#define double(X) (double)(X)
#define single(X) (float)(X)
#define char(X) uint16(X)
#define int(X)  int32(X)
#define int8(X) (char)(X)
#define int16(X) (short int)(X)
#define int32(X) (long int)(X)
#define int64(X) (long long)(X)
#define uint8(X) (unsigned char)(X)
#define uint16(X) (unsigned short int)(X)
#define uint32(X) (unsigned long int)(X)
#define uint64(X) (unsigned long long)(X)

#define pointer(X) (&(X))

#define error(X) ErrStr=(X)
#define disp(X) fprintf(stderr,"At element %d:%s = %.15g\n",_iter_i_+1, #X , (double)(X))
#define dispstr(X) fprintf(stderr,"At element %d: %s\n",_iter_i_+1,(X))
#define mod(X,Y) ((X)%(Y))
#define round(X) (floor((X)-0.5))
#define power(X,Y) pow(X,Y)
#define prod_size(X) ( _ ## X ## _elements_ ) 
#define _iter_elements_  _iter_end_  
#define size(X,D) (((int)(D)> ( _ ## X ## _no_dims_ ) || (int)(D)<1 )?1:_ ## X ## _dims_[(int)(D)-1])
#define length(X) MAX(  MAX(size(X,1),size(X,2))  ,  size(X,3)  )
#define ndims(X) ( _ ## X ## _no_dims_ ) 
#define INDEXERR(X,mi,ma) ((X)<(mi)?((mi)+0*(int)(ErrStr="Index out of range")):((X)>=(ma)? ((ma)+0*((int)(ErrStr="Index out of range"))):(X))) 
#define _ERROR_FF_CALL(ERRSTR) \
mexErrMsgTxt(ERRSTR)

volatile const char *ErrStr=NULL;

const int ones_dims[]={1,1};


/* FUNCTION that combine two arrays of dimention sizes and allocate an array from it*/
mxArray *mk_iterator_NumericArray(int ndims1, const int *dims1,int ndims2,const int *dims2,mxClassID class)
{
    int ndims=0;
    int dims[512];
    if(ndims1+ndims2>512) mexErrMsgTxt("To many dimentions in array to create");
    memcpy(&dims[ndims],dims1,ndims1*sizeof(int));
    ndims+=ndims1;
    memcpy(&dims[ndims],dims2,ndims2*sizeof(int));
    ndims+=ndims2;
    return mxCreateNumericArray(ndims,dims,class,mxREAL);
}
int myGetNumberOfDimensions(const mxArray *ptr)
{
    const int *dims;
    int ndims=mxGetNumberOfDimensions(ptr);
    if(ndims>2) return ndims;
    dims=mxGetDimensions(ptr);
    if(dims[1]>1) return 2;
    if(dims[0]>1) return 1;
    return 0;
}

const double pi=3.141592653589793115997963468544185161590576171875;

void mexFunction(int nlhs_IN,mxArray *plhs_IN[],int nrhs_IN, const mxArray *prhs_IN[])
{
    /* Variables */
     mxArray *plhs[51],*prhs[51];
     int nlhs=0,nrhs=0;
     int _iter_start_=0, _iter_stop_=0, _iter_end_=0,no_arg_in=0,no_arg_ret=0;
     int _iter_no_dims_=2;
     const int *_iter_dims_=NULL;
     int _iter_make_dims[6]={1,1,1,1,1,1};
    

    /* Input and Return Variables ToFrom MATLAB WS */

    const double _def_A_[]={0};
    const double *_A_=_def_A_;
    int    _A_no_dims_ =0;
    int const   *_A_dims_    =ones_dims;
    int    _A_elements_=1;

    const double _def_B_[]={0};
    const double *_B_=_def_B_;
    int    _B_no_dims_ =0;
    int const   *_B_dims_    =ones_dims;
    int    _B_elements_=1;

    double _def_corr_[]={0};
    double *_corr_=_def_corr_;
    int    _corr_no_dims_ =0;
    int const   *_corr_dims_    =ones_dims;
    int    _corr_elements_=1;



    /*----Argument init section---*/ 
    {
        int idx=0;
             
        ErrStr=NULL;
        if(nrhs_IN==0 && nlhs_IN==1){plhs_IN[0]=mxCreateString(
"ELEMENT ITERATOR Version 1.0 Beta (C) 2001-2002 Peter Rydesäter -pipeline -force -fastmath -typecast [double corr]=prcorr2( A, B) init: elements=prod_size(iter) iter: suma=suma+A sumb=sumb+B post: ma=suma/elements mb=sumb/elements itersec: a=(A)-ma b=(B)-mb sab+=a*b saa+=a*a sbb+=b*b postsec: corr=sab/(sqrt(saa*sbb)) ");return;}

        do{        
            /* Start of breakable section */

            if( nrhs_IN>0)
                prhs[nrhs++]=(mxArray *)prhs_IN[0];
            else
               break;
            if( nrhs_IN>1)
                prhs[nrhs++]=(mxArray *)prhs_IN[1];
            else
               break;
        }while(0); /* End of breakable section */
        nlhs=MAX(nlhs_IN,0);

        if(nrhs_IN>2)
            mexWarnMsgTxt( "Too many arguments to function prcorr2. Argument(s) ignored.");

        /*==============GET DIMENSIONS SIZE FOR ITERATOR============*/
        /* If iterator dimensions not know, get it from first arguments shape/size */
        if(nrhs>idx+0 && _iter_dims_==NULL){
            _iter_no_dims_=myGetNumberOfDimensions(prhs[idx+0])-(0);
            _iter_dims_=mxGetDimensions(prhs[idx+0]);
            _iter_end_= mxGetNumberOfElements(prhs[idx+0])/(1);
        }
        /* If iterator dimensions not know, Compose it....*/
        if(_iter_dims_==NULL){
           _iter_dims_=_iter_make_dims;
           _iter_end_=1;
           _iter_no_dims_=0;
        }

        if(_iter_no_dims_<=0){_iter_no_dims_=1; _iter_dims_=ones_dims;}

        /*==============END GET DIMS========================*/
        do{
            /*========GET ARGUMENTS======================*/
            const mxArray *ptr=NULL;
              
            /* Get argument:  A    */
            if(nrhs<=idx) break; else no_arg_in=1;
            ptr=prhs[idx++];
            if(mxGetClassID(ptr)!=mxDOUBLE_CLASS){
                mxArray *cast_plhs[1]={ NULL };
                mexCallMATLAB(1,cast_plhs,1,(mxArray **)&ptr,"double");
                ptr=cast_plhs[0];
            }
            if(mxGetClassID(ptr)!=mxDOUBLE_CLASS)
                _ERROR_FF_CALL("Datatype double expected as input argument \"A\" (1)");
            if(mxIsComplex(ptr))
                _ERROR_FF_CALL("Expected REAL number as input argument \"A\" (1)");
            _A_ =(double *)mxGetPr(ptr);
            _A_no_dims_=myGetNumberOfDimensions(ptr);
            _A_elements_=mxGetNumberOfElements(ptr);
            _A_dims_=mxGetDimensions(ptr);
            /* Remove leading iterator dimensions */
            {   int n=0,e=1;
                for(n=0; ;n++){
                    if(e==_iter_end_ && _A_no_dims_-n==0)
                        break;
                    if( (n>=_A_no_dims_) || (e * _A_dims_[n]>_iter_end_) )
                        _ERROR_FF_CALL("Can not match size of dimensions for A to iterator");
                    e *=_A_dims_[n];
                }
                _A_elements_/=e;
                _A_no_dims_-=n;
                _A_dims_=&_A_dims_[n];
            }
            if( _A_elements_!=1 )
                _ERROR_FF_CALL("Expected 1 elements as iterator window dims as input argument \"A\" (1)");
              
            /* Get argument:  B    */
            if(nrhs<=idx) break; else no_arg_in=2;
            ptr=prhs[idx++];
            if(mxGetClassID(ptr)!=mxDOUBLE_CLASS){
                mxArray *cast_plhs[1]={ NULL };
                mexCallMATLAB(1,cast_plhs,1,(mxArray **)&ptr,"double");
                ptr=cast_plhs[0];
            }
            if(mxGetClassID(ptr)!=mxDOUBLE_CLASS)
                _ERROR_FF_CALL("Datatype double expected as input argument \"B\" (2)");
            if(mxIsComplex(ptr))
                _ERROR_FF_CALL("Expected REAL number as input argument \"B\" (2)");
            _B_ =(double *)mxGetPr(ptr);
            _B_no_dims_=myGetNumberOfDimensions(ptr);
            _B_elements_=mxGetNumberOfElements(ptr);
            _B_dims_=mxGetDimensions(ptr);
            /* Remove leading iterator dimensions */
            {   int n=0,e=1;
                for(n=0; ;n++){
                    if(e==_iter_end_ && _B_no_dims_-n==0)
                        break;
                    if( (n>=_B_no_dims_) || (e * _B_dims_[n]>_iter_end_) )
                        _ERROR_FF_CALL("Can not match size of dimensions for B to iterator");
                    e *=_B_dims_[n];
                }
                _B_elements_/=e;
                _B_no_dims_-=n;
                _B_dims_=&_B_dims_[n];
            }
            if( _B_elements_!=1 )
                _ERROR_FF_CALL("Expected 1 elements as iterator window dims as input argument \"B\" (2)");
        }while(0);

        do{
            /*========CREATE RETURN ARGUMENTS============*/
              
            /* Create return argument:  corr */
            no_arg_ret=1;
            {int di[]={1,1};plhs[0]=mxCreateNumericArray(2,di,mxDOUBLE_CLASS,mxREAL);}
            _corr_=(double *)mxGetPr(plhs[0]);
            _corr_no_dims_=myGetNumberOfDimensions(plhs[0]);
            _corr_elements_=mxGetNumberOfElements(plhs[0]);
            _corr_dims_=mxGetDimensions(plhs[0]);
        }while(0);


        if(no_arg_in!=2 ) mexErrMsgTxt("Wrong number of input variables");
        if(no_arg_ret!=1) mexErrMsgTxt("Wrong number of return arguments");
          _iter_stop_=_iter_end_;
    } /* End argument Init section  */
    /* THE ELEMENT ITERATOR LOOP  2 1 */
    {
        register int _iter_i_=0; 
#define A _A_[_iter_i_]
#define B _B_[_iter_i_]
#define corr _corr_[0]
#define _continue_this_loop_  _continue_this_loop_2_1_
        double elements=0;
        double suma=0;
        double sumb=0;
        double ma=0;
        double mb=0;
        double a=0;
        double b=0;
        double sab=0;
        double saa=0;
        double sbb=0;
        elements = prod_size(iter);
        for( _iter_i_=_iter_start_; _iter_i_<_iter_stop_; _iter_i_++){
            if(_iter_i_+16>=_iter_stop_){
                suma = suma+A;
                sumb = sumb+B;
            }else{
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                _iter_i_++;
                suma = suma+A;
                sumb = sumb+B;
                }
        }
        _iter_i_=_iter_stop_-1; /*Back to last element for post processing  */
        ma = suma/elements;
        mb = sumb/elements;
        for( _iter_i_=_iter_start_; _iter_i_<_iter_stop_; _iter_i_++){
            if(_iter_i_+16>=_iter_stop_){
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
            }else{
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                _iter_i_++;
                a = (A)-ma;
                b = (B)-mb;
                sab += a*b;
                saa += a*a;
                sbb += b*b;
                }
        }
        _iter_i_=_iter_stop_-1; /*Back to last element for post processing  */
        corr = sab/(sqrt(saa*sbb));
#undef _continue_this_loop_
#undef A
#undef B
#undef corr
    }
    do{        /* Start of breakable section */
    if(1)
        plhs_IN[0]=plhs[0];
    }while(0); /* End of breakable section */
    if(ErrStr) mexErrMsgTxt((const char *)ErrStr);
}

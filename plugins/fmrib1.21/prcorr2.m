% PRCORR2  - Compute correlation coefficient over all dimentions
%
%    This is a 10 time faster implementation of MATLABs corr2.
%    Implemented as a mex-file generated with the ITERATOR
%    tool that can be downloaded from www.mathworks.com
%    Re compile the source with MEX -O prcorr2.c
%    
%
%    C = PRCORR2(A,B) computes the correlation between A and B.
%    A and B are arrays witn any number of dims but with same 
%    number of elements
% 
%    Class Support
%    -------------
%    A and B can be any numeric type but not complex value.
%    C is a scalar double.
% 
%    See also CORR2, CORRCOEF, STD2, ITERATOR, MEX
%
%
%    (C) 2003 Peter Rydes?ter,  http://www.rydesater.com
%
%    Edits by Rami Niazy: If binary not available, use corrcoef.
%    NOTE: FOR VECTORS ONLY

function R=prcorr2(A,B)
tmp=corrcoef(A,B);
R=tmp(2);
return;

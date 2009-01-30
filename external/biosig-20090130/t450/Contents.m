% BIOSIG/T450 contains functions for multiple test statistics
%
% === MAIN FUNCTIONS === 
%  FDR.M	false discovery rate
%  FDP.M	false discovery proportion 
%  gFWE.M	generalized familiy-wise error 
%  GLOB.M	global hypothesis test
%
% === UTILITY FUNCTIONS === 
%    --- do not use them directly if you not know what to do. At least you are warned --- 

%
% REFERENCES: 
%  [1] Hemmelmann C, Horn M, Suesse T, Vollandt R, Weiss S.
%	New concepts of multiple tests and their use for evaluating high-dimensional EEG data.
%	J Neurosci Methods. 2005 Mar 30;142(2):209-17.
%  [2] Hemmelmann C, Horn M, Reiterer S, Schack B, Suesse T, Weiss S.
%	Multivariate tests for the evaluation of high-dimensional EEG data.
%	J Neurosci Methods. 2004 Oct 15;139(1):111-20. 
%

%	$Id: Contents.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%	Copyright (C) 2006,2007 by Alois Schloegl <a.schloegl@ieee.org>	
%	This is part of the BIOSIG project http://biosig.sf.net/


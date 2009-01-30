% BIOSIG/T300 contains signal processing functions for Matlab/Octave
%     It requires also the TSA-toolbox [1].
%
% Currently the following methods are supported:  
%   Multivariate Autoregressive Analysis   
%       baccala2001
%       mvar
%       mvfilter
%	mvfreqz
%   Time-varying Autoregressive spectral estimation  
%       aar, amarma
%       tfar
%	tvaar	wrapper for AAR estimators
%   Time-varying Multivariate Autoregressive Analysis
%       mvar
%       mvaar
%       tfmvar
%   EEG 
%       lumped model
%       kemp's feedback loop model
%	evoked potential 
%	arspectrum
%   EMG analysis
%       Paynter
%   ECG analysis
%	qrsdetect
%	berger
%       ecgbcorr
%       qrscorr
%       tvaar
%	heartratevariability
%   SaO2 - Oxygen saturation
%	desatur
%   Others:
%	Brainrate (including SEF90, SEF95)
%       Hjorth
%       Bandpower
%	Wackermann
% 
%
% REFERENCES: 
%  [1]  A. Schlögl, Time Series Analysis toolbox for Matlab and Octave. 1996-2004
%     available online: http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/download.html
%
%

%	$Id: Contents.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 1997-2004,2007 by Alois Schloegl <a.schloegl@ieee.org>	
%	This is part of the BIOSIG project http://biosig.sf.net/


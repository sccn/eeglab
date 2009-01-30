% Demo for Writing WAV files %
% DEMO5 is part of the biosig-toolbox
%     it demonstrates generating WAV files 
%     and contains a few tests 
% 

%	$Revision: 1.1 $
%	$Id: demo5.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 2003 by Alois Schloegl <a.schloegl@ieee.org>	

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


clear
F{1}='test1.wav';
s0 = randn(1000,6);	% Generate Test data
s0 = [-10:10]'*[1:size(s0,2)];

% File type, format specification
    HDR.TYPE='WAV';		% Define file format
% Filename
    HDR.FileName = F{1};	% Assign Filename
% Sampling frequency
    HDR.SampleRate = 8000;	% Sampling rate
% Scaling information
    HDR.PhysMax = max(abs(s0(:)));	% Physical maximum 
% normalize data
    s = (s0 / HDR.PhysMax);	% Wav does not include scaling
    HDR.FLAG.UCAL = 0;
% number of bits
    HDR.bits = 16;
% number of channels
    HDR.NS = size(s,2);     	
% number samples (per record)
    HDR.SPR = size(s,1);     	
    
%%%%%%% 1st way to generate WAV-file
HDR.FileName = F{1};	% Assign Filename
HDR = sopen(HDR,'w'); 	% OPEN BKR FILE
swrite(HDR,s); 
HDR = sclose(HDR);            % CLOSE BKR FILE

% test file 
[s1,H1]=sload(F{1});

plot([s-s1]);
 
% Demo for Writing BKR files %
% DEMO4 is part of the biosig-toolbox
%     it demonstrates generating BKR files 
%     and contains a few tests 
% 

%	$Revision: 1.1 $
%	$Id: demo4.m,v 1.1 2009-01-30 06:04:39 arno Exp $
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
F{1}='test1.bkr';
F{2}='test2.bkr';
s = randn(1000,5);	% Generate Test data

% File type, format specification
    HDR.TYPE='BKR';		% Define file format
% Filename
    HDR.FileName = F{1};	% Assign Filename
% Sampling frequency
    HDR.SampleRate = 100;	% Sampling rate
% Scaling information
    HDR.PhysMax = max(abs(s(:)));	% Physical maximum 
    HDR.DigMax  = max(2^15-1);	% Digital  maximum
% filters [Hz]
    HDR.Filter.LowPass  = 30;	% upper cutoff frequency
    HDR.Filter.HighPass = .5;	% lower cutoff frequency 

% number of records, HDR.NRec must be fixed before EEGCLOSE, in case of triggered data (FLAG.TRIGGERED==1)
    HDR.NRec = -1;		% number of trials (1 for continous data, >1 triggered data, <0 unknown),  
% FLAG.TRIGGERED indicates if triggered data is stored.
    HDR.FLAG.TRIGGERED = 0;		% 0: continous data, 1 Triggered data
% HDR.NRec must be known or FLAG.TRIGGERED must be set
% number of channels [must be fixed before calling EEGCLOSE]
    HDR.NS = size(s,2);	


%%%%%%% 1st way to generate BKR-file
HDR.FileName = F{1};	% Assign Filename
HDR = sopen(HDR,'w'); 	% OPEN BKR FILE
dig_values = s'*HDR.DigMax/HDR.PhysMax; 	% digital values without scaling, each row is a channel.
fwrite(HDR.FILE.FID,dig_values,'int16'); 
% number of blocks [must be fixed before calling SCLOSE]
HDR.NRec= 1;
HDR = sclose(HDR);            % CLOSE BKR FILE

%%%%%%% 2nd way to generate BKR-file
% number of records, HDR.NRec must be fixed before EEGCLOSE, in case of triggered data (FLAG.TRIGGERED==1)
    HDR.NRec = -1;		% number of trials (1 for continous data, >1 triggered data, <0 unknown),  
% FLAG.TRIGGERED indicates if triggered data is stored.
    HDR.FLAG.TRIGGERED = 1;		% 0: continous data, 1 Triggered data
% if HDR.NRec is not set, and FLAG.TRIGGERED = 0, HDR.NRec is set to 1; 

HDR.FileName = F{2};	% Assign Filename
HDR.Classlabel=[1,2,3,0,0];
HDR = sopen(HDR,'w'); 	% OPEN BKR FILE
HDR = swrite(HDR,s);  	% WRITE BKR FILE
HDR.NRec = 1; 
HDR = swrite(HDR,s);  	% WRITE BKR FILE
% number of blocks [must be fixed before calling EEGCLOSE]
HDR.NRec = HDR.NRec+1; 
HDR = sclose(HDR);            % CLOSE BKR FILE


% read file 
HDR = sopen(HDR,'r'); 	% OPEN BKR FILE
[s0,HDR] = sread(HDR);  	% WRITE BKR FILE
HDR = sclose(HDR);


[s1,H1]=sload('test1.bkr');
[s2,H2]=sload('test2.bkr');



 
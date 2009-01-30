% Demostration for generating EDF/BDF/GDF-files
% DEMO3 is part of the biosig-toolbox
%    and it tests also Matlab/Octave for its correctness. 
% 

%	$Id: demo3.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 2000-2005,2006,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.

x = randn(10000,6)+(1:1e4)'*ones(1,6); % test data
x = (1:1e4)'*ones(1,6)/1000; % test data
x = reshape(mod(1:6e4,100),6,1e4)'; x(:,6)=NaN;

clear HDR;

VER   = version;
cname = computer;

% select file format 
HDR.TYPE='GDF';
%HDR.TYPE='EDF';
%HDR.TYPE='BDF'; 
%HDR.TYPE='CFWB';
%HDR.TYPE='CNT';

% set Filename
HDR.FileName = ['TEST_',VER([1,3]),cname(1:3),'_e1.',HDR.TYPE];

% person identification, max 80 char
HDR.Patient.ID = 'P0000';	
HDR.Patient.Sex = 'F';
HDR.Patient.Birthday = [1951 05 13 0 0 0];
HDR.Patient.Name = 'X';		% for privacy protection  
HDR.Patient.Handedness = 0; 	% unknown, 1:left, 2:right, 3: equal

% description of recording device 
HDR.Manufacturer.Name = 'BioSig'; 
HDR.Manufacturer.Model = 'demo3.m'; 
HDR.Manufacturer.Version = '$Revision'; 
HDR.Manufacturer.SerialNumber = '00000000';

% recording identification, max 80 char.
HDR.RID = 'TestFile 001'; %StudyID/Investigation [consecutive number];
HDR.REC.Hospital   = 'BioSig Test Lab'; 
HDR.REC.Techician  = 'Mister Muster';
HDR.REC.Equipment  = 'biosig';
HDR.REC.IPaddr	   = [127,0,0,1];	% IP address of recording system 	
HDR.Patient.Name   = 'anonymous';  
HDR.Patient.Id     = '007';
HDR.Patient.Weight = 0; 	% undefined 
HDR.Patient.Height = 0; 	% undefined 
HDR.Patient.Sex    = 0; 	% 0: undefined,	1: male, 2: female 
HDR.Patient.Birthday = zeros(1,6); %    undefined 
HDR.Patient.Impairment.Heart = 0;  %	0: unknown 1: NO 2: YES 3: pacemaker 
HDR.Patient.Impairment.Visual = 0; %	0: unknown 1: NO 2: YES 3: corrected (with visual aid) 
HDR.Patient.Smoking = 0;           %	0: unknown 1: NO 2: YES 
HDR.Patient.AlcoholAbuse = 0; 	   %	0: unknown 1: NO 2: YES 
HDR.Patient.DrugAbuse = 0; 	   %	0: unknown 1: NO 2: YES 
HDR.Patient.Handedness = 0; 	   % 	unknown, 1:left, 2:right, 3: equal

% recording time [YYYY MM DD hh mm ss.ccc]
HDR.T0 = clock;	

% number of channels
HDR.NS = size(x,2);

% Duration of one block in seconds
HDR.SampleRate = 1024.6;
HDR.SPR = 1000;   
HDR.Dur = HDR.SPR/HDR.SampleRate;

% Samples within 1 block
HDR.AS.SPR = [1000;100;200;100;20;1];	% samples per block;
%HDR.AS.SampleRate = [1000;100;200;100;20;0];	% samplerate of each channel

% channel identification, max 80 char. per channel
HDR.Label=['chan 1  ';'chan 2  ';'chan 3  ';'chan 4  ';'chan 5  ';'NEQS    '];

% Transducer, mx 80 char per channel
HDR.Transducer = {'Ag-AgCl ';'Airflow ';'xyz     ';'        ';'        ';'Thermome'};

% define datatypes (GDF only, see GDFDATATYPE.M for more details)
HDR.GDFTYP = 3*ones(1,HDR.NS);

% define scaling factors 
HDR.PhysMax = [100;100;100;100;100;100];
HDR.PhysMin = [0;0;0;0;0;0];
HDR.DigMax  = [100;100;100;100;100;1000];
HDR.DigMin  = [0;0;0;0;0;0];
HDR.Filter.Lowpass = [0,0,0,NaN,NaN,NaN];
HDR.Filter.Highpass = [100,100,100,NaN,NaN,NaN];
HDR.Filter.Notch = [0,0,0,0,0,0];


% define physical dimension
HDR.PhysDim = {'uV';'mV';'%';'-';'-';'°C'};

t = [100:100:size(x,1)]';
%HDR.NRec = 100;
HDR.VERSION = 2.11; 
HDR = sopen(HDR,'w');
%HDR.SIE.RAW = 0; % [default] channel data mode, one column is one channel 
%HDR.SIE.RAW = 1; % switch to raw data mode, i.e. one column for one EDF-record

HDR = swrite(HDR,x);

HDR.EVENT.POS = t;
HDR.EVENT.TYP = t/100;
if 1, 
HDR.EVENT.CHN = repmat(0,size(t));
HDR.EVENT.DUR = repmat(1,size(t));
HDR.EVENT.VAL = repmat(NaN,size(t));
ix = 6:5:60; 
HDR.EVENT.CHN(ix) = 6; 
HDR.EVENT.VAL(ix) = 373+round(100*rand(size(ix))); % HDR.EVENT.TYP(ix) becomes 0x7fff
ix = 8; 
HDR.EVENT.CHN(ix) = 5; % not valid because #5 is not sparse sampleing
HDR.EVENT.VAL(ix) = 374; 
end; 

HDR = sclose(HDR);


%
[s0,HDR0] = sload(HDR.FileName);	% test file 

HDR0=sopen(HDR0.FileName,'r');
[s0,HDR0]=sread(HDR0);
HDR0=sclose(HDR0); 

%plot(s0-x)



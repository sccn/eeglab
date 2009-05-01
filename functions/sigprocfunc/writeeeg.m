% writeeeg - Generating CNT/EDF/BDF/GDF-files using BIOSIG toolbox. Note
%            that the CNT file format is not fully functional. See also the 
%            writecnt.m Matlab function (there is no fully working
%            Neuroscan writing function to our knowledge).
%
% writeeeg( filename, data, srate, 'key', 'val')
%
% Inputs:
%  filename   - [string] filename
%  data       - [float] continuous or epoched data (channels x times x
%               epochs)
%  srate      - [float] sampling rate in Hz
%
% Optional keys:
%  'TYPE'         - ['GDF'|'EDF'|'BDF'|'CFWB'|'CNT'] file format for writing
%                   default is 'EDF'.
%  'EVENT'        - event structure (BIOSIG or EEGLAB format)
%  'Label'        - cell array of channel labels. Warning, this input is
%                   case sensitive.
%  'SPR'          - [integer] sample per block, default is 1000
%
% Other optional keys:
%  'Patient.id'       - [string] person identification, max 80 char
%  'Patient.Sex'      - [0|1|2] 0: undefined (default),	1: male, 2: female
%  'Patient.Birthday' - Default [1951 05 13 0 0 0];
%  'Patient.Name'     - [string] default 'anonymous' for privacy protection  
%  'Patient.Handedness' - [0|1|2|3] 0: unknown (default), 1:left, 2:right, 
%                         3: equal
%  'Patient.Weight'            - [integer] default 0, undefined
%  'Patient.Height'            - [integer] default 0, undefined
%  'Patient.Impairment.Heart'  - [integer] 0: unknown 1: NO 2: YES 3: pacemaker 
%  'Patient.Impairment.Visual' - [integer] 0: unknown 1: NO 2: YES 
%                                3: corrected (with visual aid) 
%  'Patient.Smoking'           - [integer] 0: unknown 1: NO 2: YES 
%  'Patient.AlcoholAbuse'      - [integer] 0: unknown 1: NO 2: YES 
%  'Patient.DrugAbuse'         - [integer] 0: unknown 1: NO 2: YES 
%
%  'Manufacturer.Name          - [string] default is 'BioSig/EEGLAB' 
%  'Manufacturer.Model         - [string] default is 'writeeeg.m'
%  'Manufacturer.Version       - [string] default is current version in CVS repository
%  'Manufacturer.SerialNumber  - [string] default is '00000000'
%  
%  'T0'            - recording time [YYYY MM DD hh mm ss.ccc]
%  'RID'           - [string] StudyID/Investigation, default is empty
%  'REC.Hospital   - [string] default is empty
%  'REC.Techician  - [string] default is empty
%  'REC.Equipment  - [string] default is empty
%  'REC.IPaddr	   - [integer array] IP address of recording system. Default is
%                    [127,0,0,1]
%
% Author: Arnaud Delorme, SCCN, UCSD/CERCO, 2009
%         Based on BIOSIG, sopen and swrite

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 22 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $

function HDR = writeeeg(filename, x, srate, varargin)

% decode structures
% -----------------
HDR = [];
for i=1:2:length( varargin )
    str = varargin{i};
    val = varargin{i+1};
    ind = findstr('.', str);
    if length(ind) == 2
        HDR = setfield(HDR, varargin{i}(1:ind(1)-1), {1}, varargin{i}(ind(1)+1:ind(2)-1), ...
                                                     {1}, varargin{i}(ind(2)+1:end), val);
    elseif length(ind) == 1
        HDR = setfield(HDR, varargin{i}(1:ind-1), {1}, varargin{i}(ind+1:end), val);
    else
        HDR = setfield(HDR, varargin{i}, varargin{i+1});
    end;
end;

HDR.FileName = filename;

% select file format 
if ~isfield(HDR, 'TYPE'),                      HDR.TYPE ='GDF'; end;
if ~isfield(HDR, 'Patient')                    HDR.Patient = []; end;
if ~isfield(HDR.Patient, 'ID')                 HDR.Patient.ID = 'P0000'; end;
if ~isfield(HDR.Patient, 'Sex')                HDR.Patient.Sex = 0; end;
if ~isfield(HDR.Patient, 'Name')               HDR.Patient.Name = 'Anonymous'; end;
if ~isfield(HDR.Patient, 'Handedness')         HDR.Patient.Handedness = 0; end;% 	unknown, 1:left, 2:right, 3: equal

% description of recording device 
if ~isfield(HDR,'Manufacturer')                HDR.Manufacturer = []; end;
if ~isfield(HDR.Manufacturer, 'Name'),         HDR.Manufacturer.Name = 'BioSig/EEGLAB'; end;
if ~isfield(HDR.Manufacturer, 'Model'),        HDR.Manufacturer.Model = 'writeeeg.m'; end;
if ~isfield(HDR.Manufacturer, 'Version'),      HDR.Manufacturer.Version = '$Revision'; end;
if ~isfield(HDR.Manufacturer, 'SerialNumber'), HDR.Manufacturer.SerialNumber = '00000000'; end;

% recording identification, max 80 char.
if ~isfield(HDR,'RID')                         HDR.RID = ''; end; %StudyID/Investigation [consecutive number];
if ~isfield(HDR,'REC')                         HDR.REC = []; end;
if ~isfield(HDR.REC, 'Hospital')               HDR.REC.Hospital = ''; end; 
if ~isfield(HDR.REC, 'Techician')              HDR.REC.Techician = ''; end;
if ~isfield(HDR.REC, 'Equipment')              HDR.REC.Equipment = ''; end;
if ~isfield(HDR.REC, 'IPaddr')             	   HDR.REC.IPaddr = [127,0,0,1]; end; % IP address of recording system 	
if ~isfield(HDR.Patient, 'Weight')             HDR.Patient.Weight = 0; end; 	% undefined 
if ~isfield(HDR.Patient, 'Height')             HDR.Patient.Height = 0; end;	% undefined 
if ~isfield(HDR.Patient, 'Birthday')           HDR.Patient.Birthday = [1951 05 13 0 0 0]; end; %    undefined 
if ~isfield(HDR.Patient, 'Impairment')         HDR.Patient.Impairment = []; end; 
if ~isfield(HDR.Patient.Impairment, 'Heart')   HDR.Patient.Impairment.Heart = 0; end;  %	0: unknown 1: NO 2: YES 3: pacemaker 
if ~isfield(HDR.Patient.Impairment, 'Visual')  HDR.Patient.Impairment.Visual = 0; end; %	0: unknown 1: NO 2: YES 3: corrected (with visual aid) 
if ~isfield(HDR.Patient,'Smoking')             HDR.Patient.Smoking = 0; end;           %	0: unknown 1: NO 2: YES 
if ~isfield(HDR.Patient,'AlcoholAbuse')        HDR.Patient.AlcoholAbuse = 0; end; 	   %	0: unknown 1: NO 2: YES 
if ~isfield(HDR.Patient,'DrugAbuse')           HDR.Patient.DrugAbuse = 0; end;  	   %	0: unknown 1: NO 2: YES 
    
% recording time [YYYY MM DD hh mm ss.ccc]
if ~isfield(HDR,'T0')                          HDR.T0 = clock; end;
if ~isfield(HDR,'SPR')                         HDR.SPR = size(x,2); end;
if  isfield(HDR,'label')                       HDR.Label = HDR.label; HDR = rmfield(HDR, 'label'); end;
if HDR.SPR > 1000, HDR.SPR = round(srate); end;
   
% channel identification, max 80 char. per channel
if ~isfield(HDR,'Label')
    for i = 1:size(x,1)
        chans{i} = sprintf('Chan %d', i);
    end;
    HDR.Label = chans;
end;
if iscell(HDR.Label), HDR.Label = strvcat(HDR.Label{:}); end;

if ~isempty(HDR.EVENT) 
    if ~isfield(HDR.EVENT, 'TYP')
        EVENT = [];
        EVENT.CHN = zeros(length(HDR.EVENT),1);
        EVENT.TYP = zeros(length(HDR.EVENT),1);
        EVENT.POS = zeros(length(HDR.EVENT),1);
        EVENT.DUR = ones( length(HDR.EVENT),1);
        EVENT.VAL = zeros(length(HDR.EVENT),1)*NaN;
        if isfield(HDR.EVENT, 'type')
            alltypes = unique( { HDR.EVENT.type } );
        end;
        for i = 1:length(HDR.EVENT)
            if isfield(HDR.EVENT, 'type')
                ind = strmatch(HDR.EVENT(i).type, alltypes);
                EVENT.TYP(i) = ind;
            end;
            if isfield(HDR.EVENT, 'latency')
                EVENT.POS(i) = HDR.EVENT(i).latency;
            end;
            if isfield(HDR.EVENT, 'duration')
                EVENT.TYP(i) = HDR.EVENT(i).duration;
            end;
        end;
        HDR.EVENT = EVENT;
        HDR.EVENT.SampleRate = srate;
    end;
end;

% reformat data
HDR.NRec = size(x,3);
x = double(x);
x = reshape(x, [size(x,1) size(x,2)*size(x,3) ]);

% recreate event channel
if isempty(HDR.EVENT)
    HDR.EVENT.POS = [];
    HDR.EVENT.TYP = [];
    HDR.EVENT.POS = [];
    HDR.EVENT.DUR = [];
    HDR.EVENT.VAL = [];
    HDR.EVENT.CHN = [];
elseif ~strcmpi(HDR.TYPE, 'GDF') % GDF can save events
    disp('Recreating event channel');
    x(end+1,:) = 0;
    x(end,round(HDR.EVENT.POS')) = HDR.EVENT.TYP';
    HDR.Label = strvcat(HDR.Label, 'Status');
    HDR.EVENT.POS = [];
end;

% number of channels
HDR.NS = size(x,1);

if ~isfield(HDR,'Transducer')
    HDR.Transducer = cell(1,HDR.NS);
    HDR.Transducer(:) = { '' };
end;
if ~isfield(HDR,'PhysDim')
    HDR.PhysDim = cell(1,HDR.NS);
    HDR.PhysDim(:) = { 'uV' };
    % HDR.PhysDim = {'uV';'mV';'%';'-';'-';'°C'};
end;

% Duration of one block in seconds
HDR.SampleRate = srate;
HDR.Dur = HDR.SPR/HDR.SampleRate;

% Samples within 1 block
%HDR.AS.SPR = [1000;100;200;100;20;1];	% samples per block;
%HDR.AS.SampleRate = [1000;100;200;100;20;0];	% samplerate of each channel


% define datatypes (GDF only, see GDFDATATYPE.M for more details)
HDR.GDFTYP = 3*ones(1,HDR.NS);

% define scaling factors 
HDR.PhysMax = ones(HDR.NS,1)*100;
HDR.PhysMin = zeros(HDR.NS,1);
HDR.DigMax  = ones(HDR.NS,1)*100;
HDR.DigMin  = zeros(HDR.NS,1);
HDR.Filter.Lowpass  = zeros(1,HDR.NS)*NaN;
HDR.Filter.Highpass = zeros(1,HDR.NS)*NaN;
HDR.Filter.Notch    = zeros(1,HDR.NS)*NaN;

%HDR.NRec = 100;
HDR.VERSION = 2.11; 
%HDR.SIE.RAW = 0; % [default] channel data mode, one column is one channel 
%HDR.SIE.RAW = 1; % switch to raw data mode, i.e. one column for one
%EDF-record
    
    
HDR = sopen(HDR,'w');
HDR = swrite(HDR, x');
HDR = sclose(HDR);
return;

%
[s0,HDR0] = sload(HDR.FileName);	% test file 

HDR0=sopen(HDR0.FileName,'r');
[s0,HDR0]=sread(HDR0);
HDR0=sclose(HDR0); 

%plot(s0-x)



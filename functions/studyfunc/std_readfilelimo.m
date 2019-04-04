% std_readlimofile()- Read limo file into eeglab structure
%
% Usage:   
%   >> [MeasureData, parameters, MeasureRange1, MeasureRange2] = std_readfilelimo(filestruct,'key', val)
%
% Inputs:
% filestrcuct  - Strcuture with the following fields: 
%                datasets : Index of he datasets in datasetinfo (real)
%                trials   : index of the trials to pull out     (cell array)
%                value    : Cell array with the values of the independet variables of the set.                 
%                case     : name of the subject /case
%                filebase : Full Path and filename of the 
%
% Optional inputs:
%   'measure'     - ['itcbeta1' | 'itcbeta2' |'itcr2r' |'itcr2f' |'itcr2p' |
%                    'erpbeta1' | 'erpbeta2' |'erpr2r' |'erpr2f' |'erpr2p' |
%                    'erspbeta1'| 'erspbeta2'|'erspr2r'|'erspr2f'|'erspr2p'|
%                    'specbeta1'| 'specbeta2'|'specr2r'|'specr2f'|'specr2p ] 
%   'timelimits'  - [min max] ERSP time (latency in ms) range of interest
%   'freqlimits'  - [min max] ERSP frequency range (in Hz) of interest
%   'channels'    - [cell or integer] channel labels - for instance 
%                   { 'cz' 'pz' }
%
% Outputs:
%   MeasureData    - the multi-channel or multi-component data. The size of this
%                    output depends on the number of dimension (for ERSP or ERP
%                    for instance). The last dimension is channels or components.
%   parameters     - structure containing parameters saved with the data file.
%   MeasureRange1  - time points (ERP, ERSP) or frequency points (spectrum)
%   MeasureRange2  - frequency points (ERSP, ITCs)
%
%
%
% See also:
%
%
% Authors: Arnaud Delorme, Ramon Martinez-Cancino
%
% Copyright (C) 2014  Ramon Martinez-Cancino,INC, SCCN
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.
function [MeasureData, parameters, MeasureRange1, MeasureRange2] = std_readfilelimo(filestruct, varargin)

parameters = [];
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g = []; end
catch
    error('std_readfilelimo() error: calling convention {''key'', value, ... } error');
end

try, g.channels;                catch, g.channels          = [];        end; %
try, g.components;              catch, g.components        = [];        end; %
try, g.measure;                 catch, g.measure           = [];        end; %
try, g.timelimits;              catch, g.timelimits        = [];        end; %
try, g.triallimits;             catch, g.triallimits       = [];        end; %
try, g.freqlimits;              catch, g.freqlimits        = [];        end; %
try, g.getparamonly;            catch, g.getparamonly      = 'off';     end; %  work on this         

%  Checking entries
%..........................................................................



% Indetify if components or channels
%..........................................................................
if isempty(g.components)
    chanoric = 'Channels';
else
    chanoric = 'Components';
end

% Identify type of analysis for limo stuff path
%..........................................................................
if strcmp(g.measure(1:3),'erp') && strcmp(chanoric,'Channels')
    type = 'Channels_Time';
elseif strcmp(g.measure(1:3),'erp') && strcmp(chanoric,'Components')
    type = 'Components_Time';
elseif strcmp(g.measure(1:3),'spe') && strcmp(chanoric,'Channels')
    type = 'Channels_Frequency';
elseif strcmp(g.measure(1:3),'spe') && strcmp(chanoric,'Components')
    type = 'Components_Frequency';
elseif strcmp(g.measure(1:3),'ers') && strcmp(chanoric,'Channels')
    type = 'Channels_Time-Frequency';
elseif strcmp(g.measure(1:3),'ers') && strcmp(chanoric,'Components')
    type = 'Components_Time-Frequency';
end

% Study pathname
%..........................................................................
filebasetemp = filestruct.filebase;
[path, name, tmp] = fileparts(filebasetemp);
[path2study,tmp,tmp]= fileparts(path);

% Getting design index in STUDY
%..........................................................................
designindx = name(7:9);
if isempty(str2num(designindx))
    designindx = designindx(1:2);
    if isempty(str2num(designindx))
        designindx = designindx(1);
        if isempty(str2num(designindx))
            error('std_readfilelimo(): Fail to detect index of the design in STUDY');
        end
    end
end

% Creating path to LIMO.mat and loading it
%..........................................................................
filelist        = dir(path2study);
fileindxtmp     = strmatch('LIMO_', {filelist.name});
limofoldername  = filelist(fileindxtmp(1)).name;
path2limofiles  = [path2study filesep limofoldername filesep filestruct.case filesep 'GLM' num2str(designindx) '_' type];
structmp        = load ([path2limofiles filesep 'LIMO.mat']);
limostruct      = structmp.LIMO;

% Identifying condition
%..........................................................................
dims = []; cond_indxs = [];
condtemp = {limostruct.data.studyinfo.value};

for i =  1 : length(condtemp)
    dims(i)       = length(condtemp{i});
    tmp = condtemp{i};
    for j = 1: dims(i)
        tmp{j} = num2str(tmp{j});
    end
    condtemp{i} = tmp;
    cond_indxs(i) = find(strcmp(condtemp{i},num2str(filestruct.value{i})));
end
% Sub2ind stuff
command = 'sub2ind(dims';
for i = 1: length(dims)
    command = [command ',cond_indxs(' num2str(i) ')'];
end
command = [command ')'];
CondIndx = eval(command);

% Checking time and set TimeIndxVec if type = TIME
%..........................................................................
if all([~isempty(g.timelimits),isfield(limostruct.data,'sampling_rate')])
    if any([g.timelimits(1) < limostruct.data.start, g.timelimits(2) > limostruct.data.end])
        error('std_readfilelimo(): Time indices provided are out of bonds');
    else
        TimeVec = limostruct.data.start: 1000/limostruct.data.sampling_rate:limostruct.data.end
        [tmp,TimeIndx(1)] = min(abs(TimeVec - g.timelimits(1)));
        [tmp,TimeIndx(2)] = min(abs(TimeVec - g.timelimits(2)));
        TimeIndxVec = TimeIndx(1):TimeIndx(2);
    end
elseif all([isempty(g.timelimits),isfield(limostruct.data,'sampling_rate')]) % will use all the timepoints
    TimeIndxVec = 1:length(limostruct.data.start: 1000/limostruct.data.sampling_rate:limostruct.data.end);
end

% Checking frequencies and set FreqIndxVec if type = Frequency
%..........................................................................
 if all([~isempty(g.freqlimits),isfield(limostruct.data,'freqlist')])
     if any([g.freqlimits(1) < limostruct.data.lowf, g.freqlimits(2) > limostruct.data.highf])
         error('std_readfilelimo(): Frequency indices provided are out of bonds');
     else
         [tmp,FreqIndx(1)] = min(abs(limostruct.data.freqlist - g.freqlimits(1)));
         [tmp,FreqIndx(2)] = min(abs(limostruct.data.freqlist - g.freqlimits(2)));
         FreqIndxVec = FreqIndx(1):FreqIndx(2);
     end
 elseif all([isempty(g.freqlimits),isfield(limostruct.data,'freqlist')]) % will use all frequencies
     FreqIndxVec = 1 : length(limostruct.data.freqlist);
 end
 
 % Checking frequencies and time and set FreqIndxVec and TimeIndxVec if type = Time-Frequency
 %..........................................................................
 if all([isfield(limostruct.data,'tf_times'), isfield(limostruct.data,'tf_freq')])
     
     %1- Time
     if ~isempty(g.timelimits)
         if any([g.timelimits(1) < limostruct.data.start, g.timelimits(2) > limostruct.data.end])
             error('std_readfilelimo(): Time indices provided are out of bonds');
         else
             [tmp,TimeIndx(1)] = min(abs(limostruct.data.tf_times - g.timelimits(1)));
             [tmp,TimeIndx(2)] = min(abs(limostruct.data.tf_times - g.timelimits(2)));
             TimeIndxVec = TimeIndx(1):TimeIndx(2);
         end
     else % will use all timepoints
         TimeIndxVec = 1:length(limostruct.data.tf_times);
     end
     
     %2- Frequency
     if ~isempty(g.freqlimits)
         if any([g.freqlimits(1) < limostruct.data.lowf, g.freqlimits(2) > limostruct.data.highf])
             error('std_readfilelimo(): Frequency indices provided are out of bonds');
         else
             [tmp,FreqIndx(1)] = min(abs(limostruct.data.tf_freqs - g.freqlimits(1)));
             [tmp,FreqIndx(2)] = min(abs(limostruct.data.tf_freqs - g.freqlimits(2)));
             FreqIndxVec = FreqIndx(1):FreqIndx(2);
         end
     else % will use all frequencies
         FreqIndxVec = 1:length(limostruct.data.tf_freqs);
     end
     
 end

% Checking channels or ic to load
%..........................................................................
% 1- Channels
if isempty(g.channels) && strcmp(chanoric, 'Channels')
    ChanIcIndx = 1:length({limostruct.data.chanlocs.labels});
elseif ~isempty(g.channels) && strcmp(chanoric, 'Channels')
     for i = 1 : length(g.channels)
         if isempty(find(strcmp(lower(g.channels(i)),lower({limostruct.data.chanlocs.labels}))))
             error('std_readfilelimo(): Invalid channel label');
         else
             ChanIcIndx = find(strcmp(lower(g.channels(i)),lower({limostruct.data.chanlocs.labels})));    
         end
     end
end

% 2- Component
if isempty(g.components) && strcmp(chanoric, 'Components')
    ChanIcIndx = 1:length({limostruct.data.chanlocs.labels});
elseif ~isempty(g.channels) && strcmp(chanoric, 'Components')
     for i = 1 : length(g.channels)
         if isempty(find(strcmp(g.channels(i),{limostruct.data.chanlocs.labels})))
             error('std_readfilelimo(): Invalid channel label');
         else
             ChanIcIndx = find(strcmp(g.channels(i),{limostruct.data.chanlocs.labels}));    
         end
     end
end

% Retreiving data
%..........................................................................
if isfield(limostruct.data,'sampling_rate')
    TimeVec = limostruct.data.start: 1000/limostruct.data.sampling_rate:limostruct.data.end
end

switch g.measure
    % beta1
    %..................................................................
    case 'erpbeta1'
        datatmp = load([path2limofiles filesep 'Betas.mat']); % PUT THIS OUT OF THE LOOP
        data = datatmp.Betas;
        MeasureData   = data(ChanIcIndx,TimeIndxVec,CondIndx);   % Chan/IC x Time x Condition

    case 'specbeta1'
        datatmp = load([path2limofiles filesep 'Betas.mat']);
        data = datatmp.Betas;
        MeasureData   = data(ChanIcIndx,FreqIndxVec,CondIndx);   % Chan/IC x Freq x Condition        
        
    case {'erspbeta1','itcbeta1'}
        datatmp = load([path2limofiles filesep 'Betas.mat']);
        data = datatmp.Betas;
        MeasureData   = data(FreqIndxVec, TimeIndxVec,CondIndx); % Freq x Time x Condition


    % beta2
    %..................................................................
    case {'erpbeta2','erspbeta2','specbeta2','itcbeta2'}
        error('std_readfilelimo(): Invalid optin in this realease');
    
    % erpr2
    %..................................................................
    case 'erpr2r'
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(ChanIcIndx,TimeIndxVec,1);

     case 'erpr2f'
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(ChanIcIndx,TimeIndxVec,2);       

    case 'erpr2p'
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(ChanIcIndx,TimeIndxVec,3);        
        
    % erspr2 & itcr2
    %..................................................................
    case {'erspr2r','itcr2r'}
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(FreqIndxVec,TimeIndxVec,1);

     case {'erspr2f','itcr2f'}
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(FreqIndxVec,TimeIndxVec,2);

    case {'erspr2p','itcr2p'}
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(FreqIndxVec,TimeIndxVec,3); 
    
    %spec
    %..................................................................
    case {'specr2r'}
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(ChanIcIndx,FreqIndxVec,1);

     case {'specr2f'}
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(ChanIcIndx,FreqIndxVec,2);

    case {'specr2p'}
        datatmp = load([path2limofiles filesep 'R2.mat']);
        data = datatmp.R2;
        MeasureData   = data(ChanIcIndx,FreqIndxVec,3);
end

% Measures ranges
%..........................................................................
switch g.measure(1:3)
    case 'erp'
        MeasureRange1 = TimeVec(TimeIndxVec);
        MeasureRange2 = [];
    case {'ers','itc'}
        MeasureRange1 = limostruct.data.tf_times(TimeIndxVec);
        MeasureRange2 = limostruct.data.tf_freqs(FreqIndxVec);
    case 'spe'
        MeasureRange1 = limostruct.data.freqlist(FreqIndxVec);
        MeasureRange2 = [];
end

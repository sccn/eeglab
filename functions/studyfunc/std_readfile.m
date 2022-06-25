% std_readfile() - Read data file containing STUDY measures.
%
% Usage:
%   >>  [data param range1 range2] = std_readfile(filename, 'key', val);
%
% Inputs:
%   filename   - [string] read specific file, for instance 's1.daterp'
%                containing ERP data for dataset "s1.set". It is also
%                possible to provide only the "base" file name "s1" and
%                function will load the appropriate file based on the
%                selected input measure (measure input).
%
% Optional inputs:
%   'channels'       - [cell or integer] channel labels - for instance
%                      { 'cz' 'pz' } of channels to load from the data file.
%   'components'     - [integer] component index in the selected EEG dataset
%                      for which to read the data
%   'timelimits'     - [min max] ERSP/ERP time (latency in ms) range of interest
%   'freqlimits'     - [min max] ERSP/Spectrum frequency range (in Hz) of interest
%   'measure'        - ['erp'|'spec'|'ersp'|'itc'|'timef'|'erspbase'|'erspboot'
%                      'itcboot'|'erpim'] data measure to read. If a full file name
%                      is provided, the data measure is selected automatically.
%   'getparamsonly'  - ['on'|'off'] only read file parameters not data.
%   'trialselect'    - Cell of cells defining the field and the values
%                      of the field in the trialinfo structure to be used to pull out the trials.
%                      i.e. if values are string:  {'field1',{val1_fromfield1 val2_fromfield1}, 'field2'...}
%                           if values are numeric: {'field1',[val1 val2], 'field2'...}
%   'designvar'      - Structure of independent variables with fields 'label'
%                      and 'value'. Each independent variable should have these
%                      two fields so the function can use these values to pull
%                      out the trials. If empty, all the trials are read.
%   'singletrials'   - { 'on','off' } Extract single trials. Default 'off'
%   'getparamonly'   - { 'on','off'} Get only parameters. Default 'off'
%   'concatenate'    - { 'on','off'} In case of ERP images this function
%                      set 'singletrials' to 'on'. Default 'off'
%   'dataindices'    - obsolete input
%
% Outputs:
%   data                - the multi-channel or multi-component data. The size of this
%                         output depends on the number of dimension (for ERSP or ERP
%                         for instance). The last dimension is channels or components.
%   params              - structure containing parameters saved with the data file.
%   range1              - time points (ERP, ERSP) or frequency points (spectrum)
%   range2              - frequency points (ERSP, ITCs)
%   events              - Event read from the data structure
%
% Examples:
%   % the examples below read all data channels for the selected files
%   [ersp params times freqs] = std_readfiles('s1.datersp');
%   [erp params times] = std_readfiles('s1.daterp');
%   [erp params times] = std_readfiles('s1.daterp', timerange', [-100 500]);
%   [erp params times] = std_readfiles('s1.daterp', timerange', [-100 500],'trialselect','load',[1 2],'type', {'Y'});
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, May 2010

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2010, arno@sccn.ucsd.edu
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

% dimensions
% time x freqs x channel_comps x subjects_trials

function [measureData, parameters, measureRange1, measureRange2, events] = std_readfile(fileBaseName, varargin);

if nargin < 1
    help std_readfile;
    return;
end

if iscell(fileBaseName)
    % this is when the data of a subject is split across sessions
    % then we need to call std_readfile several times and combine
    % the results
    % ------------------------------------------------------------
    if length(fileBaseName) > 1
        for iFile = 1:length(fileBaseName)
            [measureDataTmp{iFile}, parameters, measureRange1, measureRange2, eventsTmp{iFile}] = std_readfile(fileBaseName{iFile}, varargin{:});
        end
    
        % get the size of each session cond x group and sum 3rd dim (trials) for each cond and group
        % ignore empty cells; these are always trials, even when processing subjects
        measureData = cell(size(measureDataTmp{1})); % cond x group
        sz          = cell(size(measureDataTmp{1})); % cond x group
        for iSess = 1:length(measureDataTmp) % scan session
            szTmp = cellfun(@size, measureDataTmp{iSess}, 'uniformoutput', false); 
            for iCond = 1:length(sz(:))
                if szTmp{iCond}(1) ~= 0
                    if isempty(sz{iCond})
                        sz{iCond} = szTmp{iCond};
                    else
                        sz{iCond}(2) = sz{iCond}(2)+szTmp{iCond}(2);
                    end
                end
            end
        end
        for iCond = 1:length(sz(:))
            measureData{iCond} = zeros(sz{iCond});
        end

        % combine arrays
        events = eventsTmp{1};
        for iCond = 1:length(measureDataTmp{1}(:))
            pointer = 1; 
            for iSess = 1:length(measureDataTmp)
                if ~isempty(measureDataTmp{iSess}{iCond})
                    % when trials concatenate them
                    measureData{iCond}(:, pointer:pointer+size(measureDataTmp{iSess}{iCond},2)-1,:,:) = measureDataTmp{iSess}{iCond};
                    pointer = pointer + size(measureDataTmp{iSess}{iCond},2);
                end
            end
        end
        
        return
    else
        fileBaseName = fileBaseName{1};
    end
end

limomeasures = {'itcbeta1' , 'itcbeta2' ,'itcr2r' ,'itcr2f' ,'itcr2p' ,...
    'erpbeta1' , 'erpbeta2' ,'erpr2r' ,'erpr2f' ,'erpr2p' ,...
    'erspbeta1', 'erspbeta2','erspr2r','erspr2f','erspr2p',...
    'specbeta1', 'specbeta2','specr2r','specr2f','specr2p'};
opt = finputcheck(varargin, { 'components'       'integer'  []    [];
    'getparamonly'     'string'   { 'on','off' }  'off';
    'trialselect'      'cell'     {}                 {};
    'trialinfo'        'struct'   {}                 struct([]);
    'designvar'        'struct'   []                 struct([]);
    'singletrials'     'string'   { 'on','off' }  'off';
    'concatenate'      'string'   { 'on','off' }  'off'; % ERPimage only
    'channels'         'cell'     []    {};
    'cache'            'struct'   []    struct([]);
    'function'         { 'function_handle' 'integer' } []  [];
    'measure'          'string'   {limomeasures{:} 'erp' 'spec' 'timef' 'topo'} 'erp';
    'timelimits'       'real'     []    []; % ERPimage, ERP, ERSP, ITC
    'triallimits'      'real'     []    []; % ERPimage only
    'freqlimits'       'real'     []    []; % SPEC, ERSP, ITC
    'dataindices'      'integer'  []    [] }, 'std_readfile');
if ischar(opt), error(opt); end

if ~isempty(opt.triallimits), opt.freqlimits = opt.triallimits; end
if strcmpi(opt.concatenate, 'on'), opt.singletrials = 'on'; end

% get file extension
% ------------------
if ~isempty(opt.channels) || (~isempty(opt.dataindices) && opt.dataindices(1) < 0) , dataType = 'chan';
else                                                                                 dataType = 'comp';
end
[~, ~, currentFileExt] = fileparts(fileBaseName);
if length(currentFileExt) > 3 && (strcmpi(currentFileExt(2:4), 'dat') || strcmpi(currentFileExt(2:4), 'ica'))
    opt.measure = currentFileExt(5:end);
    if strcmpi(currentFileExt(2:4), 'dat'), dataType = 'chan';
    else                                    dataType = 'comp';
    end
    fileExt = '';
else
    if strcmpi(dataType, 'chan'), fileExt = [ '.dat' opt.measure ];
    else                          fileExt = [ '.ica' opt.measure ];
    end
end
if nargin > 5 && strcmpi(opt.singletrials, 'on')
    indFlag = true;
else indFlag = false;
end

% get fields to read
% ------------------
v6Flag = testv6([ fileBaseName fileExt ]);
v6Flag = 1;
if v6Flag || strcmpi(opt.measure, 'topo')
    if ~isempty(opt.channels)
        fileData = load('-mat', [ fileBaseName fileExt ], 'labels');
    end
else
    fileData = matfile([ fileBaseName fileExt ]);
end

% get channel or component indices
% --------------------------------
if ~isempty(opt.channels) && isnumeric(opt.channels)
    opt.dataindices = opt.channels;
elseif ~isempty(opt.channels)
    chan.chanlocs = struct('labels', fileData.labels);
    opt.dataindices = std_chaninds(chan, opt.channels);
elseif ~isempty(opt.components)
    opt.dataindices = opt.components;
else opt.dataindices = abs(opt.dataindices);
end

if v6Flag
    for iChan = 1:length(opt.dataindices)
        chanList{iChan} = [ dataType int2str(opt.dataindices(iChan)) fastif(strcmpi(opt.measure, 'topo'), '_grid', '') ];
    end
    warning('off', 'MATLAB:load:variableNotFound');
    if ~isempty(opt.dataindices)
        fileData = load('-mat', [ fileBaseName fileExt ], chanList{:}, 'trialinfo', 'times', 'freqs', 'parameters', 'events', 'chanlocsforinterp', 'chanall');
    else fileData = load('-mat', [ fileBaseName fileExt ], 'trialinfo', 'times', 'freqs', 'parameters', 'events', 'chanlocsforinterp', 'chanall');
    end
    % read all channels ERPimage only
    if length(opt.channels) > 2 && isfield(fileData, 'chanall')
        chanlocsforinterp = fileData.chanlocsforinterp;
        tmpX = eval( fileData.chanall );
        for iChan = 1:size(tmpX,1)
            fileData.(sprintf('chan%d', iChan)) = squeeze(tmpX(iChan,:,:));
        end
    end
    if strcmpi(opt.measure, 'topo') && ~isfield(fileData, 'trialinfo'), error('Compatibilty issue. Recompute ICA topographic maps'); end
    warning('on', 'MATLAB:load:variableNotFound');
end

% scan datasets
% -------------
if strcmpi(opt.getparamonly, 'on'), opt.dataindices = 1; end
measureRange1  = [];
measureRange2  = [];
measureData    = [];
parameters     = [];
events         = {};

% get output for parameters and measure ranges
% --------------------------------------------
fileFields = fieldnames(fileData);
if ~isempty( strmatch('parameters', fileFields) )
    parameters = removedup(fileData.parameters);
    for index = 1:length(parameters), if iscell(parameters{index}), parameters{index} = { parameters{index} }; end; end
    parameters = struct(parameters{:});
end
if ~isempty(strmatch('times', fileFields, 'exact')),  measureRange1 = fileData.times; end
if ~isempty(strmatch('freqs', fileFields, 'exact')),  measureRange2 = fileData.freqs; end
if isempty(measureRange1) && isfield(fileData, 'chan1'), measureRange1 = 1:size(fileData.chan1,1); end
    
% if the function is only called to get parameters
% ------------------------------------------------
if strcmpi(opt.measure, 'spec'), measureRange1 = measureRange2; opt.timelimits = opt.freqlimits; end
if ~isempty(opt.timelimits)
    [measureRange1, indBegin1, indEnd1 ] = indicesselect(measureRange1, opt.timelimits);
else
    indBegin1 = 1;
    indEnd1   = length(measureRange1);
end
if ~isempty(opt.freqlimits)
    [measureRange2, indBegin2, indEnd2 ] = indicesselect(measureRange2, opt.freqlimits);
else
    indBegin2 = 1;
    indEnd2   = length(measureRange2);
end

if strcmpi(opt.getparamonly, 'on')
    return;
end

options = { opt.dataindices, opt.function, dataType, indBegin1, indEnd1, indBegin2, indEnd2 };
if isempty(opt.designvar)
    [ measureData, events ] = getfiledata(fileData, NaN, v6Flag, options{:}); % read all data
    measureData = { measureData };
    events      = { events };
else
    [ measureData, events ] = globalgetfiledata(fileData, opt.designvar, options, {}, v6Flag);
end

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
[tmp, indices] = unique_bc(cella(1:2:end));
if length(tmp) ~= length(cella)/2
    %fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
end
cella = cella(sort(union(indices*2-1, indices*2)));

% find indices for selection of measure
% -------------------------------------
function [measureRange, indBegin, indEnd] = indicesselect(measureRange, measureLimits)
indBegin = 1;
indEnd   = length(measureRange);
if ~isempty(measureRange) && ~isempty(measureLimits) && (measureLimits(1) > measureRange(1) || measureLimits(end) < measureRange(end))
    indBegin   = min(find(measureRange >= measureLimits(1)));
    indEnd     = max(find(measureRange <= measureLimits(end)));
    measureRange = measureRange(indBegin:indEnd);
end

% recursive function to load data
% -------------------------------
function [ measureData, eventVals ] = globalgetfiledata(fileData, designvar, options, trialselect, v6Flag)

if length(designvar) == 0
    [ measureData, eventVals ] = getfiledata(fileData, trialselect, v6Flag, options{:});
    measureData = { measureData };
    eventVals   = { eventVals   };
else
    % scan independent variable values
    if isfield(designvar(1), 'vartype') && strcmpi('continuous', designvar(1).vartype)
        if ~ischar(designvar(1).value), designvar(1).value = ''; end
        trialselect = { trialselect{:} designvar(1).label designvar(1).value };
        [ tmpMeasureData tmpEvents ] = globalgetfiledata(fileData, designvar(2:end), options, trialselect, v6Flag);
        measureData(1,:,:,:) = reshape(tmpMeasureData, [ 1 size(tmpMeasureData) ]);
        eventVals(  1,:,:,:) = reshape(tmpEvents     , [ 1 size(tmpEvents     ) ]);
    else
        trialselectOri = trialselect;
        for iField = 1:length(designvar(1).value)
            trialselect = { trialselectOri{:} designvar(1).label designvar(1).value{iField} };
            [ tmpMeasureData, tmpEvents ] = globalgetfiledata(fileData, designvar(2:end), options, trialselect, v6Flag);
            measureData(iField,:,:,:) = reshape(tmpMeasureData, [ 1 size(tmpMeasureData) ]);
            eventVals(  iField,:,:,:) = reshape(tmpEvents     , [ 1 size(tmpEvents     ) ]);
        end
    end
end

% load data from structure or file
% --------------------------------
function [ fieldData, events ] = getfiledata(fileData, trialselect, v6Flag, chan, func, dataType, indBegin1, indEnd1, indBegin2, indEnd2)

persistent tmpcache;
persistent hashcode;

if length(chan) > 1
    %    error('This function can only read one channel at a time');
end

% get trial indices
fieldData = [];
events    = [];
subTrials = [];
trials    = [];
if ~isempty(trialselect)
    if isnumeric(trialselect) && isnan(trialselect(1))
        trials = [1:length(fileData.trialinfo)]; % read all trials if NaN
    else
        [trials, events] = std_gettrialsind(fileData.trialinfo, trialselect{:});
        if length(unique(diff(trials))) > 1 && ~v6Flag
            temptrials = [trials(1):trials(end)];
            subTrials  = trials-trials(1)+1;
            trials     = temptrials;
        end
    end
end
% if ~isempty(subTrials)
%     trials = subTrials(trials(
    
for index = 1:length(chan)
    allfields   = fieldnames(fileData);
    topoFlag    = ~isempty(findstr(allfields{1}, '_grid'));
    fieldToRead = [ dataType int2str(chan(index)) fastif(topoFlag, '_grid', '') ];
    
    % find trials
    if isempty(trials)
        return;
        % trials = size(fileData.(fieldToRead), ndims(fileData.(fieldToRead))); % not sure what this does
    end
    
    % load data
    warning('off', 'MATLAB:MatFile:OlderFormat');
    if ischar(fileData.(fieldToRead)) % special ERP-image
        try
            fileData.chanlocsforinterp; % isfield does not work because fileData is a MatFile
        catch, error('Missing field in ERPimage STUDY file, try recomputing them');
        end
        chanlocsforinterp = fileData.chanlocsforinterp;
        
        % caching for ERPimage only
        if isequal(hashcode, fileData.(fieldToRead))
            tmpFieldData = tmpcache;
        else
            tmpFieldData = eval( fileData.(fieldToRead) );
            tmpcache = tmpFieldData;
            hashcode = fileData.(fieldToRead);
        end
        tmpFieldData = tmpFieldData(indBegin1:indEnd1,trials);
        if ~isempty(subTrials), tmpFieldData = tmpFieldData(:, subTrials); end
    else
        if topoFlag
            tmpFieldData = fileData.(fieldToRead);
            if isempty(trials), tmpFieldData = []; end
        elseif ndims(fileData.(fieldToRead)) == 2
            tmpFieldData = fileData.(fieldToRead)(indBegin1:indEnd1,trials);
            if ~isempty(subTrials), tmpFieldData = tmpFieldData(:, subTrials); end
        else
            tmpFieldData = fileData.(fieldToRead)(indBegin2:indEnd2,indBegin1:indEnd1,trials); % frequencies first here
            if ~isempty(subTrials), tmpFieldData = tmpFieldData(:, :, subTrials); end
        end
    end
    if isfield(fileData, 'events') && ~isempty(fileData.events)
        events = fileData.events(trials);
        if ~isempty(subTrials), events = events(subTrials); end
    end
    warning('on', 'MATLAB:MatFile:OlderFormat');
    
    % average single trials if necessary
    if ~isempty(func)
        tmpFieldData = func(tmpFieldData);
    end
    
    % store data
    if index == 1 && length(chan) == 1
        fieldData = tmpFieldData;
    else
        if index == 1
            fieldData = zeros([ size(tmpFieldData) length(chan) ]);
        end
        if ndims(tmpFieldData) == 2
            fieldData(:,:,index) = tmpFieldData;
        else
            fieldData(:,:,:,index) = tmpFieldData;
        end
    end
end

% see if a file is v6 of v7.3
function v6 = testv6(x)

fid=fopen(x);
if fid == -1
    error('File %s not found - this could be because your STUDY contains data files with relative path, try changing your MATLAB path to the STUDY folder', x);
end
txt=char(fread(fid,20,'uchar')');
tmp = fclose(fid);
txt=[txt,char(0)];
txt=txt(1:find(txt==0,1,'first')-1);
if ~isempty(strfind(txt, 'MATLAB 5.0')), v6 = true; else v6 = false; end


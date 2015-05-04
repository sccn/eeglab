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
%                      { 'cz' 'pz' }
%                      of channels to load from the data file.
%   'components'     - [integer] component index in the selected EEG dataset for which 
%                      to read the ERSP
%   'timelimits'     - [min max] ERSP time (latency in ms) range of interest
%   'freqlimits'     - [min max] ERSP frequency range (in Hz) of interest
%   'measure'        - ['erp'|'spec'|'ersp'|'itc'|'timef'|'erspbase'|'erspboot'
%                      'itcboot'|'erpim'] data measure to read. If a full file name
%                      is provided, the data measure is selected automatically.
%   'getparamsonly'  - ['on'|'off'] only read file parameters not data.
%   'trialselect'    - Cell of cells defining the field and the values
%                      of the field in the trialinfo structure to be used to pull out the trials. 
%                      i.e. if values are string:  {'field1',{val1_fromfield1 val2_fromfield1}, 'field2'...} 
%                           if values are numeric: {'field1',[val1 val2], 'field2'...}
%   'designvar'      - Structure of Independent Variables(IV) with fields 'label'
%                      and 'value'. Each IV should have these two fields so
%                      the function can use these values to pull out the trials. 
%   'singletrials'   - { 'on','off' } Extract single trials. Default 'off'
%   'getparamonly'   - { 'on','off'} Get only parameters. Default 'off'
%   'concatenate'    - { 'on','off'} In case of ERP images this function
%                      set 'singletrials' to 'on'. Default 'off'
%   'dataindices'    - (To be updated)
%
% Outputs:
%   data                - the multi-channel or multi-component data. The size of this
%                         output depends on the number of dimension (for ERSP or ERP
%                         for instance). The last dimension is channels or components.
%   params              - structure containing parameters saved with the data file.
%   range1              - time points (ERP, ERSP) or frequency points (spectrum)
%   range2              - frequency points (ERSP, ITCs)
%   events              - Event readed from the data structure
%   setinfoTrialIndices - Deprecated
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

% dimensions
% time x freqs x channel_comps x subjects_trials

function [measureData, parameters, measureRange1, measureRange2, events] = std_readfile(fileBaseName, varargin);

if nargin < 1
    help std_readfile;
    return;
end;

limomeasures = {'itcbeta1' , 'itcbeta2' ,'itcr2r' ,'itcr2f' ,'itcr2p' ,...
                'erpbeta1' , 'erpbeta2' ,'erpr2r' ,'erpr2f' ,'erpr2p' ,...
                'erspbeta1', 'erspbeta2','erspr2r','erspr2f','erspr2p',...
                'specbeta1', 'specbeta2','specr2r','specr2f','specr2p'};
opt = finputcheck(varargin, { 'components'       'integer'  []    [];
                              'getparamonly'     'string'   { 'on','off' }  'off';
                              'trialselect'      'cell'     {}                 {};
                              'designvar'        'struct'   []                 struct([]);
                              'singletrials'     'string'   { 'on','off' }  'off';
                              'concatenate'      'string'   { 'on','off' }  'off'; % ERPimage only
                              'channels'         'cell'     []    {};
                              'function'         { 'function_handle' 'integer' } []  [];
                              'measure'          'string'   {limomeasures{:} 'erp' 'spec' 'timef'} 'erp';                                                 
                              'timelimits'       'real'     []    []; % ERPimage, ERP, ERSP, ITC
                              'triallimits'      'real'     []    []; % ERPimage only
                              'freqlimits'       'real'     []    []; % SPEC, ERSP, ITC
                              'dataindices'      'integer'  []    [] }, 'std_readdatafile');
if isstr(opt), error(opt); end;

if ~isempty(opt.triallimits), opt.freqlimits = opt.triallimits; end;
if strcmpi(opt.concatenate, 'on'), opt.singletrials = 'on'; end;

% get file extension
% ------------------
if ~isempty(opt.channels) || (~isempty(opt.dataindices) && opt.dataindices(1) < 0) , dataType = 'chan';
else                                                                                 dataType = 'comp';
end;
[tmp1 tmp2 currentFileExt] = fileparts(fileBaseName);
if length(currentFileExt) > 3 && (strcmpi(currentFileExt(2:4), 'dat') || strcmpi(currentFileExt(2:4), 'ica'))
    opt.measure = currentFileExt(5:end);
    if strcmpi(currentFileExt(2:4), 'dat'), dataType = 'chan';
    else                                    dataType = 'comp';
    end;
    fileExt = '';
else
    if strcmpi(dataType, 'chan'), fileExt = [ '.dat' opt.measure ];
    else                          fileExt = [ '.ica' opt.measure ];
    end;
end;
if nargin > 5 && strcmpi(opt.singletrials, 'on')
     indFlag = true;
else indFlag = false;
end;

% get fields to read
% ------------------
fileData = matfile([ fileBaseName fileExt ]);

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
end;

% scan datasets
% -------------
if strcmpi(opt.getparamonly, 'on'), opt.dataindices = 1; end;
measureRange1  = [];
measureRange2  = [];
measureData    = [];
parameters     = [];
events         = {};

% get output for parameters and measure ranges
% --------------------------------------------
fileFields = fieldnames(fileData);
if any(strncmp('parameters', fileFields, 100))
    parameters = removedup(fileData.parameters);
    for index = 1:length(parameters), if iscell(parameters{index}), parameters{index} = { parameters{index} }; end; end;
    parameters = struct(parameters{:});
end;
if any(strncmp('times', fileFields, 100)),  measureRange1 = fileData.times; end;
if any(strncmp('freqs', fileFields, 100)),  measureRange2 = fileData.freqs; end;

% if the function is only called to get parameters
% ------------------------------------------------
if ~isempty(opt.timelimits)
    [measureRange1 indBegin1 indEnd1 ] = indicesselect(measureRange1, opt.timelimits);
else
    indBegin1 = 1;
    indEnd1   = length(measureRange1);
end;
if ~isempty(opt.freqlimits)
    [measureRange2 indBegin2 indEnd2 ] = indicesselect(measureRange2, opt.freqlimits);
else
    indBegin2 = 1;
    indEnd2   = length(measureRange2);
end;
if strcmpi(opt.measure, 'spec'), measureRange1 = measureRange2; indBegin2 = indBegin1; indEnd2 = indEnd1; end;
if strcmpi(opt.getparamonly, 'on'),
    return;
end;

options = { opt.dataindices, opt.function, dataType, indBegin1, indEnd1, indBegin2, indEnd2 };
[ measureData events ] = globalgetfiledata(fileData, opt.designvar, options, {});

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique_bc(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        %fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));
    
% find indices for selection of measure
% -------------------------------------
function [measureRange indBegin indEnd] = indicesselect(measureRange, measureLimits);
    indBegin = 1;
    indEnd   = length(measureRange);
    if ~isempty(measureRange) && ~isempty(measureLimits) && (measureLimits(1) > measureRange(1) || measureLimits(end) < measureRange(end))
        indBegin   = min(find(measureRange >= measureLimits(1)));
        indEnd     = max(find(measureRange <= measureLimits(end)));
        measureRange = measureRange(indBegin:indEnd);
    end;

% recursive function to load data
% -------------------------------
function [ measureData eventVals ] = globalgetfiledata(fileData, designvar, options, trialselect);

    if length(designvar) == 0
        [ measureData eventVals ] = getfiledata(fileData, trialselect, options{:});
        measureData = { measureData };
        eventVals   = { eventVals   };
    else
        % scan independent variable values
        if strcmpi('continuous', designvar(1).vartype)
            if ~isstr(designvar(1).value), designvar(1).value = ''; end;
            trialselect = { trialselect{:} designvar(1).label designvar(1).value };
            [ tmpMeasureData tmpEvents ] = globalgetfiledata(fileData, designvar(2:end), options, trialselect);
            measureData(1,:,:,:) = reshape(tmpMeasureData, [ 1 size(tmpMeasureData) ]);
            eventVals(  1,:,:,:) = reshape(tmpEvents     , [ 1 size(tmpEvents     ) ]);
        else
            for iField = 1:length(designvar(1).value)
                trialselect = { trialselect{:} designvar(1).label designvar(1).value{iField} };
                [ tmpMeasureData tmpEvents ] = globalgetfiledata(fileData, designvar(2:end), options, trialselect);
                measureData(iField,:,:,:) = reshape(tmpMeasureData, [ 1 size(tmpMeasureData) ]);
                eventVals(  iField,:,:,:) = reshape(tmpEvents     , [ 1 size(tmpEvents     ) ]);
            end;
        end;
    end;
    
% load data from structure or file
% --------------------------------
function [ fieldData events ] = getfiledata(fileData, trialselect, chan, func, dataType, indBegin1, indEnd1, indBegin2, indEnd2)

if length(chan) > 1
    error('This function can only read one channel at a time');
end;

% get trial indices
fieldData = [];
events    = [];
subTrials = [];
trials    = [];
if ~isempty(trialselect)
    [trials events] = std_gettrialsind(fileData.trialinfo, trialselect{:});
    if length(unique(diff(trials))) > 1
        temptrials = [trials(1):trials(end)];
        subTrials  = trials-trials(1)+1;
        trials     = temptrials;
    end;
end;

fieldToRead = [ dataType int2str(chan) ];

% find trials
if isempty(trials), 
    return;
    % trials = size(fileData.(fieldToRead), ndims(fileData.(fieldToRead))); % not sure what this does
end;

% load data
warning('off', 'MATLAB:MatFile:OlderFormat');
if ndims(fileData.(fieldToRead)) == 2
    fieldData = fileData.(fieldToRead)(indBegin1:indEnd1,trials);
    if ~isempty(subTrials), fieldData = fieldData(:, subTrials); end;
else fieldData = fileData.(fieldToRead)(indBegin1:indEnd1,indBegin2:indEnd2,trials);
    if ~isempty(subTrials), fieldData = fieldData(:, :, subTrials); end;
end;
warning('on', 'MATLAB:MatFile:OlderFormat');

% average single trials if necessary
if ~isempty(func)
    fieldData = func(fieldData);
end;


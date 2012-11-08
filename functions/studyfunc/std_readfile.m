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
%   'channels'   - [cell or integer] channel labels - for instance 
%                  { 'cz' 'pz' } - or indices - for instance [1 2 3]
%                  of channels to load from the data file.
%   'components' - [integer] component index in the selected EEG dataset for which 
%                  to read the ERSP
%   'timelimits' - [min max] ERSP time (latency in ms) range of interest
%   'freqlimits' - [min max] ERSP frequency range (in Hz) of interest
%   'measure'    - ['erp'|'spec'|'ersp'|'itc'|'timef'|'erspbase'|'erspboot'
%                  'itcboot'|'erpim'] data measure to read. If a full file name
%                  is provided, the data measure is selected automatically.
%   'getparamsonly' - ['on'|'off'] only read file parameters not data.
%
% Outputs:
%   data    - the multi-channel or multi-component data. The size of this
%             output depends on the number of dimension (for ERSP or ERP
%             for instance). The last dimension is channels or components.
%   params  - structure containing parameters saved with the data file.
%   range1  - time points (ERP, ERSP) or frequency points (spectrum)
%   range2  - frequency points (ERSP, ITCs)
%
% Examples:
%   % the examples below read all data channels for the selected files
%   [ersp params times freqs] = std_readfiles('s1.datersp');
%   [erp params times] = std_readfiles('s1.daterp');
%   [erp params times] = std_readfiles('s1.daterp', timerange', [-100 500]);
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

opt = finputcheck(varargin, { 'components'       'integer'  []    [];
                              'getparamonly'     'string'   { 'on','off' }  'off';
                              'singletrials'     'string'   { 'on','off' }  'off';
                              'concatenate'      'string'   { 'on','off' }  'off'; % ERPimage only
                              'channels'         'cell'     []    {};
                              'measure'          'string'   { 'erpim','ersp','erspboot','erspbase','itc','itcboot','spec','erp','timef' }  'erp';
                              'timelimits'       'real'     []    []; % ERPimage, ERP, ERSP, ITC
                              'triallimits'      'real'     []    []; % ERPimage only
                              'freqlimits'       'real'     []    []; % SPEC, ERSP, ITC
                              'dataindices'      'integer'  []    [] }, 'std_readdatafile');
if isstr(opt), error(opt); end;
if ~isempty(opt.triallimits), opt.freqlimits = opt.triallimits; end;
if strcmpi(opt.concatenate, 'on'), opt.singletrials = 'on'; end;
if isstruct(fileBaseName), fileBaseName = { fileBaseName.filebase }; 
else                       fileBaseName = { fileBaseName };
end;
if ~isempty(opt.channels) && length(opt.channels) < length(fileBaseName)
    opt.channels(2:length(fileBaseName)) = opt.channels(1);
end;

% get file extension
% ------------------
if ~isempty(opt.channels) || (~isempty(opt.dataindices) && opt.dataindices(1) < 0) , dataType = 'chan';
else                                                                                 dataType = 'comp';
end;
[tmp1 tmp2 currentFileExt] = fileparts(fileBaseName{1});
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

% get fields to read
% ------------------
erspFreqOnly = 0;
switch opt.measure
    case 'erpim'   , fieldExt = '';
    case 'erp'     , fieldExt = '';
    case 'spec'    , fieldExt = '';
    case 'ersp'    , fieldExt = '_ersp';
    case 'itc'     , fieldExt = '_itc';
    case 'timef'   , fieldExt = '_timef';
    case 'erspbase', fieldExt = '_erspbase'; fileExt = fileExt(1:end-4); erspFreqOnly = 1;
    case 'erspboot', fieldExt = '_erspboot'; fileExt = fileExt(1:end-4); erspFreqOnly = 1;
    case 'itcboot' , fieldExt = '_itcboot';  fileExt = fileExt(1:end-4); erspFreqOnly = 1;
end;

% get channel or component indices
% --------------------------------
if ~isempty(opt.channels) && isnumeric(opt.channels)
    opt.dataindices = opt.channels;
elseif ~isempty(opt.channels)
    %if length(fileBaseName) > 1, error('Cannot read channel labels when reading more than 1 input file'); end;
    for iFile = 1:length(fileBaseName)
        filename = [ fileBaseName{iFile} fileExt ];
        try, 
            warning('off', 'MATLAB:load:variableNotFound');
            fileData = load( '-mat', filename, 'labels', 'chanlabels' );
            warning('on', 'MATLAB:load:variableNotFound');
        catch, fileData = [];
        end;
        if isfield(fileData, 'labels'),        chan.chanlocs = struct('labels', fileData.labels);
        elseif isfield(fileData, 'chanlabels') chan.chanlocs = struct('labels', fileData.chanlabels);
        else error('Cannot use file to lookup channel names, the file needs to be recomputed');
        end;
        opt.dataindices(iFile) = std_chaninds(chan, opt.channels{iFile});
    end;
elseif ~isempty(opt.components)
     opt.dataindices = opt.components;
else opt.dataindices = abs(opt.dataindices);
end;

% file names and indices must have the same number of values
% ----------------------------------------------------------
if strcmpi(opt.getparamonly, 'on'), opt.dataindices = 1; end;
if length(fileBaseName   ) == 1, fileBaseName(   1:length(opt.dataindices)) = fileBaseName;    end;
if length(opt.dataindices) == 1, opt.dataindices(1:length(fileBaseName   )) = opt.dataindices; end;
if length(opt.dataindices) ~= length(fileBaseName) && ~isempty(opt.dataindices), error('Number of files and number of indices must be the same'); end;

% scan datasets
% -------------
measureRange1 = [];
measureRange2 = [];
measureData   = [];
parameters    = [];
events        = {};
    
% read only specific fields
% -------------------------
counttrial = 0;
for fInd = 1:length(opt.dataindices) % usually only one value

    fieldsToRead = [ dataType int2str(opt.dataindices(fInd)) fieldExt ]; 
    try,
        warning('off', 'MATLAB:load:variableNotFound');
        fileData = load( '-mat', [ fileBaseName{fInd} fileExt ], 'parameters', 'freqs', 'times', 'events', 'chanlocsforinterp', fieldsToRead );
        warning('on', 'MATLAB:load:variableNotFound');
    catch
        error( [ 'Cannot read file ''' fileBaseName{fInd} fileExt '''' ]);
    end;

    % get output for parameters and measure ranges
    % --------------------------------------------
    if isfield(fileData, 'chanlocsforinterp'), chanlocsforinterp = fileData.chanlocsforinterp; end;
    if isfield(fileData, 'parameters')
        parameters = removedup(fileData.parameters);
        for index = 1:length(parameters), if iscell(parameters{index}), parameters{index} = { parameters{index} }; end; end;
        parameters = struct(parameters{:});
    end;
    if isfield(fileData, 'times'),  measureRange1 = fileData.times; end;
    if isfield(fileData, 'freqs'),  measureRange2 = fileData.freqs; end;
    if isfield(fileData, 'events'), events{fInd}  = fileData.events; end;

    % if the function is only called to get parameters
    % ------------------------------------------------
    if strcmpi(opt.getparamonly, 'on'), 
        measureRange1 = indicesselect(measureRange1, opt.timelimits);
        measureRange2 = indicesselect(measureRange2, opt.freqlimits);
        if strcmpi(opt.measure, 'spec'), measureRange1 = measureRange2; end;
        
        parameters.singletrials = 'off';
        if strcmpi(opt.measure, 'timef')
            parameters.singletrials = 'on';
        elseif strcmpi(opt.measure, 'erp') || strcmpi(opt.measure, 'spec')
            if strcmpi(dataType, 'chan') 
                if size(fileData.chan1,1) > 1 && size(fileData.chan1,2) > 1
                    parameters.singletrials = 'on';
                end;
            else 
                if size(fileData.comp1,1) > 1 && size(fileData.comp1,2) > 1
                    parameters.singletrials = 'on';
                end;
            end;
        end;
        return; 
    end;

    % copy fields to output variables
    % -------------------------------
    if isfield(fileData, fieldsToRead)
        fieldData = getfield(fileData, fieldsToRead);
        if isstr(fieldData), eval( [ 'fieldData = ' fieldData ] ); end;
        
        % average single trials if necessary
        % ----------------------------------
        if strcmpi(opt.measure, 'erp') || strcmpi(opt.measure, 'spec')
            if strcmpi(opt.singletrials, 'off') && size(fieldData,1) > 1 && size(fieldData,2) > 1
                fieldData = mean(fieldData,2);
            end;
        end;
              
        % array reservation
        % -----------------
        if fInd == 1
            sizeFD = size(fieldData);
            if length(sizeFD) == 2 && (sizeFD(1) == 1 || sizeFD(2) == 1), sizeFD = sizeFD(1)*sizeFD(2); end;
            if strcmpi(opt.singletrials, 'off'), measureData = zeros([ sizeFD length(opt.dataindices) ], 'single');
            else                                 measureData = zeros([ sizeFD ], 'single');
            end;
            nDimData = length(sizeFD);
        end;
        
        % copy data to output variable
        % ----------------------------
        if nDimData == 1,     measureData(:,fInd)     = fieldData;
        else
            if strcmpi(opt.singletrials, 'off')
                if nDimData == 2, measureData(:,:,fInd) = fieldData;
                else            measureData(:,:,:,fInd) = fieldData;
                end;
            else
                if nDimData == 2, measureData(:,counttrial+1:counttrial+size(fieldData,2)) = fieldData;
                else              measureData(:,:,counttrial+1:counttrial+size(fieldData,2)) = fieldData;
                end;
                counttrial = counttrial+size(fieldData,2);
            end;
        end;
    elseif ~isempty(findstr('comp', fieldsToRead))
        error( sprintf([ 'Field "%s" not found in file %s' 10 'Try recomputing measure.' ], fieldsToRead, [ fileBaseName{fInd} fileExt ]));
    else
        % the case below is for the rare case where all the channels are read and the end of the array needs to be trimmed
        error('There is a problem with your data, please enter a bug report and upload your data at http://sccn.ucsd.edu/eeglab/bugzilla');
        if nDimData == 1,     measureData(:,1:(fInd-1))     = [];
        elseif nDimData == 2, measureData(:,:,1:(fInd-1))   = [];
        else                  measureData(:,:,:,1:(fInd-1)) = [];
        end;
        break;
    end;
end;

% special ERP image
% -----------------
if ~isempty(events)
    len    = length(events{1});
    events = [ events{:} ];
    if strcmpi(opt.singletrials, 'off') events = reshape(events, len, length(events)/len); end;
end;
    
% select plotting or clustering time/freq range
% ---------------------------------------------
if ~isempty(measureRange1) && ~erspFreqOnly
    [measureRange1 indBegin indEnd] = indicesselect(measureRange1, opt.timelimits);
    if ~isempty(measureData)
        if strcmpi(opt.measure, 'erp') || ( strcmpi(opt.measure, 'erpim') && strcmpi(opt.singletrials, 'on') )
             measureData = measureData(indBegin:indEnd,:,:);
        else measureData = measureData(:,indBegin:indEnd,:);
        end;
    end;
end;
if isempty(measureRange2) && size(measureData,1) > 1 && size(measureData,2) > 1 % for ERPimage
    measureRange2 = [1:size(measureData,1)];
end;
if ~isempty(measureRange2)
    [measureRange2 indBegin indEnd] = indicesselect(measureRange2, opt.freqlimits);
    if ~isempty(measureData)
        measureData = measureData(indBegin:indEnd,:,:);
    end;
    if strcmpi(opt.measure, 'spec'), measureRange1 = measureRange2; end;
end;

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique(cella(1:2:end));
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


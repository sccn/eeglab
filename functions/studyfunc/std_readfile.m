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
%   'channels'   - [cell or integet] channel labels - for instance 
%                  { 'cz' 'pz' } - or indices - for instance [1 2 3]
%                  of channels to load from the data file.
%   'components' - [integer] component index in the selected EEG dataset for which 
%                  to read the ERSP
%   'timelimits' - [min max] ERSP time (latency in ms) range of interest
%   'freqlimits' - [min max] ERSP frequency range (in Hz) of interest
%   'measure'    - ['erp'|'spec'|'ersp'|'itc'|'timef'|'erspbase'|'erspboot'
%                  'itcboot'] data measure to read. If a full file name
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
%   [ersp params times freqs] = std_readdatafiles('s1.datersp');
%   [erp params times] = std_readdatafiles('s1.daterp');
%   [erp params times] = std_readdatafiles('s1.daterp', timerange', [-100 500]);
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, May 2010

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: std_readdatafile.m,v $

% dimensions
% time x freqs x channel_comps x subjects_trials

function [measureData, parameters, measureRange1, measureRange2] = std_readfile(fileBaseName, varargin);

if nargin < 1
    help std_readfile;
    return;
end;

opt = finputcheck(varargin, { 'dataindices'      'integer'  []    [];
                              'components'       'integer'  []    [];
                              'getparamonly'     'string'   { 'on' 'off' }  'off';
                              'singletrials'     'string'   { 'on' 'off' }  'off';
                              'channels'         'cell'     []    {};
                              'measure'          'string'   { 'ersp' 'erspboot' 'erspbase' 'itc' 'itcboot' 'spec' 'erp' 'timef' }  'erp';
                              'timelimits'       'real'     []    [];
                              'freqlimits'       'real'     []    [] }, 'std_readdatafile');
if isstr(opt), error(opt); end;
if isstruct(fileBaseName), fileBaseName = { fileBaseName.filebase }; 
else                       fileBaseName = { fileBaseName };
end;

% get file extension
% ------------------
if ~isempty(opt.channels) || (~isempty(opt.dataindices) && opt.dataindices(1) < 0) , dataType = 'chan';
else                                                                                 dataType = 'comp';
end;
[tmp1 tmp2 currentFileExt] = fileparts(fileBaseName{1});
if ~isempty(currentFileExt)
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
switch opt.measure
    case 'erp'     , fieldExt = '';
    case 'spec'    , fieldExt = '';
    case 'ersp'    , fieldExt = '_ersp';
    case 'itc'     , fieldExt = '_itc';
    case 'timef'   , fieldExt = '_timef';
    case 'erspbase', fieldExt = '_erspbase'; fileExt = fileExt(1:end-4);
    case 'erspboot', fieldExt = '_erspboot'; fileExt = fileExt(1:end-4);
    case 'itcboot' , fieldExt = '_itcboot';  fileExt = fileExt(1:end-4);
end;

% get channel or component indices
% --------------------------------
if ~isempty(opt.channels) && isnumeric(opt.channels)
    opt.dataindices = opt.channels;
elseif ~isempty(opt.channels)
    if length(fileBaseName) > 1, error('Cannot read channel labels when reading more than 1 input file'); end;
    filename = [ fileBaseName{1} fileExt ];
    try, 
        warning('off', 'MATLAB:load:variableNotFound');
        fileData = load( '-mat', filename, 'labels', 'chanlabels' );
        warning('on', 'MATLAB:load:variableNotFound');
    catch, fileData = [];
    end;
    if ~isempty(fileData)
        if isfield(fileData, 'labels'), chan.chanlocs = struct('labels',     fileData.labels);
        else                            chan.chanlocs = struct('labels', fileData.chanlabels);
        end;
        opt.dataindices = std_chaninds(chan, opt.channels);
    else
        warning('Cannot use file to lookup channel names, the file needs to be recomputed')
        return;
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
if length(opt.dataindices) ~= length(fileBaseName), error('Number of files and number of indices must be the same'); end;

% scan datasets
% -------------
measureRange1 = [];
measureRange2 = [];
measureData   = [];
parameters    = [];
    
% read only specific fields
% -------------------------
for fInd = 1:length(opt.dataindices) % usually only one value

    fieldsToRead = [ dataType int2str(opt.dataindices(fInd)) fieldExt ]; 
    try,
        warning('off', 'MATLAB:load:variableNotFound');
        fileData = load( '-mat', [ fileBaseName{fInd} fileExt ], 'parameters', 'freqs', 'times', fieldsToRead );
        warning('on', 'MATLAB:load:variableNotFound');
    catch
        error( [ 'Cannot read file ''' fileBaseName{fInd} fileExt '''' ]);
    end;

    % get output for parameters and measure ranges
    % --------------------------------------------
    if isfield(fileData, 'parameters')
        parameters = removedup(fileData.parameters);
        parameters = struct(parameters{:});
    end;
    if isfield(fileData, 'times'),  measureRange1 = fileData.times; end;
    if isfield(fileData, 'freqs'),  measureRange2 = fileData.freqs; end;

    % if the function is only called to get parameters
    % ------------------------------------------------
    if strcmpi(opt.getparamonly, 'on'), 
        if isempty(measureRange1), measureRange1 = measureRange2; end;
        
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
            if nDimData == 2
                if strcmpi(opt.singletrials, 'off'), measureData(:,:,fInd)                      = fieldData;
                else                                 measureData(:,end+1:end+size(fieldData,2)) = fieldData;
                end;
            else
                if strcmpi(opt.singletrials, 'off'), measureData(:,:,:,fInd)                      = fieldData;
                else                                 measureData(:,:,end+1:end+size(fieldData,2)) = fieldData;
                end;
            end;
        end;
    else % the case below is for the rare case where all the channels are read
        if nDimData == 1,     measureData(:,1:(fInd-1))     = [];
        elseif nDimData == 2, measureData(:,:,1:(fInd-1))   = [];
        else                  measureData(:,:,:,1:(fInd-1)) = [];
        end;
        break;
    end;
end;

% select plotting or clustering time/freq range
% ---------------------------------------------
tminind = 1;
tmaxind = length(measureRange1);
if ~isempty(opt.timelimits) && (opt.timelimits(1) > measureRange1(1) || opt.timelimits(end) < measureRange1(end))
    tmaxind = max(find(measureRange1 <= opt.timelimits(end)));
    tminind = min(find(measureRange1 >= opt.timelimits(1)));
    if isempty(measureRange2), 
         measureData = measureData(:,tminind:tmaxind,:,:);
    else measureData = measureData(:,:,tminind:tmaxind,:);
    end;
end
fminind = 1;
fmaxind = length(measureRange2);
if ~isempty(opt.freqlimits) && (opt.freqlimits(1) > measureRange2(1) || opt.freqlimits(end) < measureRange2(end))
    fmaxind = max(find(measureRange2 <= opt.freqlimits(end)));
    fminind = min(find(measureRange2 >= opt.freqlimits(1)));
    measureData = measureData(:,fminind:fmaxind,:,:);
end
if isempty(measureRange2), measureRange1 = measureRange2; end;

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        %fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));

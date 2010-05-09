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

function [measureData, parameters, measureRange1, measureRange2] = std_readfile(fileBaseName, varargin);

if nargin < 1
    help std_readfile;
    return;
end;

opt = finputcheck(varargin, { 'dataIndices'      'integer'  []    [];
                              'components'       'integer'  []    [];
                              'getparamonly'     'string'   { 'on' 'off' }  'off';
                              'channels'         'cell'     []    {};
                              'measure'          'string'   { 'ersp' 'erspboot' 'erspbase' 'itc' 'itcboot' 'spec' 'erp' 'timef' 'auto' }  'auto';
                              'timeLimits'       'real'     []    [];
                              'freqLimits'       'real'     []    [] }, 'std_readdatafile');
if isstr(opt), return; end;

% get file extension
% ------------------
if ~isempty(opt.channels) || (~isempty(opt.dataIndices) && opt.dataIndices(1) < 0) , dataType = 'chan';
else                                                                                 dataType = 'comp';
end;
[tmp1 tmp2 currentFileExt] = fileparts(fileBaseName);
if ~isempty(currentFileExt)
    opt.measure = currentFileExt(5:end);
    if strcmpi(currentFileExt(2:4), 'dat'), dataType = 'chan';
    else                                    dataType = 'comp';
    end;
    fileExt = '';
else
    if strcmpi(dataType, 'chan'), fileExt    = [ '.dat' opt.measure ];
    else                          fileExt    = [ '.ica' opt.measure ];
    end;
end;

% get channel or component indices
% --------------------------------
if ~isempty(opt.channels) && isnumeric(opt.channels)
    opt.dataIndices = opt.channels;
elseif ~isempty(opt.channels)
    filename = [ fileBaseName fileExt ];
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
        opt.dataIndices = std_chaninds(chan, opt.channels);
    else
        warning('Cannot use file to lookup channel names, the file needs to be recomputed')
        return;
    end;
elseif ~isempty(opt.components)
     opt.dataIndices = opt.components;
else opt.dataIndices = abs(opt.dataIndices);
end;

% get fields to read
% ------------------
switch opt.measure
    case 'erp'     , fieldExt = '';
    case 'spec'    , fieldExt = '';
    case 'ersp'    , fieldExt = '_ersp';
    case 'itc'     , fieldExt = '_itc';
    case 'erspboot', fieldExt = '_erspboot';
    case 'itcboot' , fieldExt = '_itcboot';
    case 'timef'   , fieldExt = '_timef';
end;

% scan datasets
% -------------
measureRange1 = [];
measureRange2 = [];
measureData   = [];
parameters    = [];
    
% read whole file
% ---------------
filename = [ fileBaseName fileExt];
if isempty(opt.dataIndices) && strcmpi(opt.getparamonly, 'off')
    fileData = load( '-mat', filename);
    opt.dataIndices = 1:(length(fieldnames(fileData))-2);
else
    fileData = [];
end;

% read only specific fields
% -------------------------
fieldsToRead = {};
for fInd1 = 1:length(opt.dataIndices)
    fieldsToRead{fInd1} = [ dataType int2str(opt.dataIndices(fInd1)) fieldExt ]; 
end;
if isempty(fileData)
    try,
        warning('off', 'MATLAB:load:variableNotFound');
        fileData = load( '-mat', filename, 'parameters', 'freqs', 'times', fieldsToRead{:} );
        warning('on', 'MATLAB:load:variableNotFound');
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;
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
    if isempty(measureRange2), measureRange1 = measureRange2; end;
    return; 
end;

% copy fields to output variables
% -------------------------------
for fInd = 1:length(opt.dataIndices) % usually only one value
    if isfield(fileData, fieldsToRead{fInd})
        fieldData = getfield(fileData, fieldsToRead{fInd});
        if fInd == 1
            measureData = zeros([ size(fieldData) length(opt.dataIndices) ], 'single');
        end;
        measureData(:,:,fInd) = fieldData;
    else
        measureData = measureData(:,:,1:(fInd-1));
        break;
    end;
end;

% select plotting or clustering time/freq range
% ---------------------------------------------
tminind = 1;
tmaxind = length(measureRange1);
if ~isempty(opt.timeLimits) && (opt.timeLimits(1) > measureRange1(1) || opt.timeLimits(end) < measureRange1(end))
    tmaxind = max(find(measureRange1 <= opt.timeLimits(end)));
    tminind = min(find(measureRange1 >= opt.timeLimits(1)));
    if isempty(measureRange2), 
         measureData = measureData(:,tminind:tmaxind,:,:);
    else measureData = measureData(:,:,tminind:tmaxind,:);
    end;
end
fminind = 1;
fmaxind = length(measureRange2);
if ~isempty(opt.freqLimits) && (opt.freqLimits(1) > measureRange2(1) || opt.freqLimits(end) < measureRange2(end))
    fmaxind = max(find(measureRange2 <= opt.freqLimits(end)));
    fminind = min(find(measureRange2 >= opt.freqLimits(1)));
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

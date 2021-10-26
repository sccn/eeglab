% std_gettrialsind() - return the index of the trials that comply with the 
%                      defined values from trialinfo
%
% Usage:
%   >>  trialindsx = std_gettrialsind('testnewstudyfile.mat', 'type', { 'C' 'X' 'Y' }, ...
%                                     'load', [0 1 2 3 4 5 6], 'uncertainty1',[1 2])
%   >>  trialindsx = std_gettrialsind('testnewstudyfile.mat', 'load', [0 1 2 3 4 5 6])
%   >>  [trialindsx cont] = std_gettrialsind('testnewstudyfile.mat', 'load', [0 1 2 3 4 5 6], ...
%                                      'duration', '0<100', 'rt', '')
%
% Inputs:
% filename - Filename of the study file, the structure contained in that
% file or, the trialinfo structure for that specific subject
% Beware of the value type: cell or string is categorical string variables, 
% numerical array is categorical numerical variables, empty string or string
% containing "<" character is for continuous variables. For continuous variables, 
% a range may be entered. Enter empty for the whole value range.
%
% varargin - Pairs of inputs ,varname(string),varvalues . (see Usage)
%
% Optional inputs:
%   ['var']    - [cell or array] ['var'] is the name of one independent
%                variable followed by a value or a list of value for this
%                variable. To obtain values for a continuous variable,
%                use '' for all the values or a range '0<100' to get values
%                within a given range (see example below).
%
% Outputs:
%   trialind   - trial indices 
%   cont       - optional continuous variable value for these indices 
%   
% Example:
%   >>  trialindsx = std_gettrialsind('testnewstudyfile.mat', 'type', { 'C' 'X' 'Y' }, ...
%                                     'load', [0 1 2 3 4 5 6], 'uncertainty1',[1 2])
%   >>  trialindsx = std_gettrialsind('testnewstudyfile.mat', 'load', [0 1 2 3 4 5 6])
%   >>  [trialindsx cont] = std_gettrialsind('testnewstudyfile.mat', 'load', [0 1 2 3 4 5 6], ...
%                                      'duration', '0<100', 'rt', '')
% Authors: Arnaud Delorme
%          Ramon Martinez-Cancino
%          
% Copyright (C) 2015  Ramon Martinez-Cancino, UCSD, INC, SCCN
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

function [trialindsx, eventvals] = std_gettrialsind(filename,varargin)

trialindsx = [];
eventvals = [];
if nargin < 1
    help std_gettrialsind;
    return; 
end

% Input stuff
try
    if length(varargin) == 1
        varargin = varargin{:}; % Call from std_readfile or other func
    end
    varsin = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(varsin)
            queryvars.(varsin{i}) = varsin{i+1};
        end
    else queryvars = []; end
catch
    disp('std_gettrialsind() error: calling convention {''key'', value, ... } error'); return;
end

% Checking if first entry is string or struct
if isstruct(filename)
    if isfield(filename,'trialinfo')
        trialinfo = filename.trialinfo;
    else
        trialinfo = filename;
        %warning('std_gettrialsind: First input assumed as trialinfo structure');
    end
else
    % Loading file
    trialinfo = load(filename,'-mat', 'trialinfo');
    trialinfo = trialinfo.trialinfo;
end

% --- Query ---
varnames = fieldnames(queryvars);
hits = zeros(length(trialinfo),length(varnames));

for iVar = 1 :  length(varnames)
    dattrials = {trialinfo.(varnames{iVar})};
    indvarvals = queryvars.(varnames{iVar});
    
    % case of continuous variable
    if ischar(indvarvals) && (isempty(indvarvals) || any(indvarvals == '<'))
        if ischar(dattrials{1})
            error(sprintf('Type error - excepting string numerical for field %s in data file', varnames{iVar}));
        end
        if isempty(indvarvals)
            hits(:,iVar) = ~cellfun(@isempty, dattrials);
        else % test for range
            indGreater = find(indvarvals == '<');
            if length(indGreater) ~= 1, error('Undefined range for continuous var'); end
            if any(indvarvals == '='), error('<= or =< not allowed for continuous var'); end
            lowerBound = str2num(indvarvals(1:indGreater-1));
            upperBound = str2num(indvarvals(indGreater+1:end));
            hits(:,iVar) = cellfun(@(x)(~isempty(x)&&x>lowerBound&&x<upperBound), dattrials);
        end
        dattrials(~hits(:,iVar)) = { NaN };
        eventvals(end+1,:) = [ dattrials{:} ];
    else
        % case of categorical variable
        if ischar(dattrials{1})
            if ischar(indvarvals) indvarvals = { indvarvals }; end
            if ~iscell(indvarvals)
                error(sprintf('Type error - excepting numerical values for field %s', varnames{iVar}));
            end
            for iVal = 1:length(indvarvals) % programmed for speed - AD
                hits(strmatch(indvarvals{iVal}, dattrials, 'exact'),iVar) = 1;
                % older version not compatible with old Matlab
                % hits(:,iVar) = hits(:,iVar) || strncmp(indvarvals{iVal}, dattrials, 100)';
            end
        else
            if iscell(indvarvals)
                error(sprintf('Type error - excepting string values for field %s', varnames{iVar}));
            end
            excl = find(cellfun(@isempty, dattrials));
            if ~(isempty(excl))
                [dattrials{excl}] = deal(nan);
            end
            dattrials = [ dattrials{:} ];
            for iVal = 1:length(indvarvals) % programmed for speed - AD
                hits(:,iVar) = hits(:,iVar) | [ dattrials == indvarvals(iVal) ]';
            end
        end
    end
end

% Retrieving overlapped values
hits = sum(hits,2);
trialindsx = find(hits == length(varnames));
if ~isempty(eventvals)
    eventvals  = eventvals(:,trialindsx);
end

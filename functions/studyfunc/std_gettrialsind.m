% std_gettrialsind() - return the index of the trials that comply with the defined values from trialinfo
%
% Usage:
%   >>  trialindsx = std_gettrialsind('testnewstudyfile.mat', 'type', { 'C' 'X' 'Y' }, 'load', [0 1 2 3 4 5 6], 'uncertainty1',[1 2])
%   >>  trialindsx = std_gettrialsind('testnewstudyfile.mat', 'load', [0 1 2 3 4 5 6])
%
% Inputs:
% filename - Filename of the study file, the structure contained in that
% file or, the trialinfo strcuture for that specific subject
%
% varargin - Pairs of inputs ,varname(string),varvalues . (see Usage)
%
% Optional inputs:
%
% Outputs:
%   trialinfo   - 
%
% See also:
%   
% Authors: Arnaud Delorme
%          Ramon Martinez-Cancino
%          
% Copyright (C) 2015  Ramon Martinez-Cancino, UCSD, INC, SCCN
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

function trialindsx = std_gettrialsind(filename,varargin)

% Input stuff
try
    varsin = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(varsin)
            queryvars.(varsin{i}) = varsin{i+1};
        end
    else queryvars = []; end;
catch
    disp('std_gettrialsind() error: calling convention {''key'', value, ... } error'); return;
end;

% Checking if fisrt entry is string or struct
if isstruct(filename)
    if isfield(filename,'trialinfo')
        trialinfo = filename.trialinfo;
    else
        trialinfo = filename;
        warning('std_gettrialsind: First input assumed as trialinfo structure');
    end
else
    % Loading file
    datastruct  = load(filename,'-mat');
    trialinfo = datastruct.trialinfo;
end

% --- Query ---
varnames = fieldnames(queryvars);
hits = zeros(length(trialinfo),length(varnames));

for i = 1 :  length(varnames)
    
    % Checking if numeric
    tmpval = {trialinfo.(varnames{i})};
    for n = 1: length(tmpval)
        checknum(n) = isnumeric(tmpval{n});
    end
     % if nuemric then convert to string
    if sum(checknum)~=0
        tmpindx = find(checknum == 1);
        
        for m = 1: length(tmpindx)
            tmpval{tmpindx(m)} = num2str(tmpval{tmpindx(m)});
        end
    end
            
    % Getting indices of values that comply individually each value    
    for j = 1 : length(queryvars.(varnames{i}))
        query = queryvars.(varnames{i})(j);
        if isnumeric(query)
            query = num2str(query);
        end
        indxtmp = find(strcmp(tmpval,query));
        if ~isempty(indxtmp) 
            hits(indxtmp,i) = 1;
        end
    end
end

% Retreiving overlapped values
hits = sum(hits,2);
trialindsx = find(hits == length(varnames));
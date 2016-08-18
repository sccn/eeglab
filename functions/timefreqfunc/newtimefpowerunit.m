% newtimefpowerunit() - Find power unit for y-axis based on input structure.
%
% Usage:    
%   >>  str = newtimefpowerunit(paramstruct); 
%
% Inputs:
%   paramstruct - [structure] structure with fields .
%
% Outputs:
%   str         - string containing the 
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, August 2016

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2016, arno@sccn.ucsd.edu
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

function unitPower = newtimefpowerunit(tmpparams)

if nargin < 1
    help newtimefunit;
    return;
end;

if ~isfield(tmpparams, 'baseline'), tmpparams.baseline = 0;     end;
if ~isfield(tmpparams, 'scale'   ), tmpparams.scale    = 'log'; end;
if ~isfield(tmpparams, 'basenorm'), tmpparams.basenorm = 'off'; end;
if strcmpi(tmpparams.scale, 'log')
    if strcmpi(tmpparams.basenorm, 'on')
        unitPower = '10*log(std.)'; % impossible
    elseif isnan(tmpparams.baseline)
        unitPower = '10*log10(\muV^{2}/Hz)';
    else
        unitPower = 'dB';
    end;
else
    if strcmpi(tmpparams.basenorm, 'on')
        unitPower = 'std.';
    elseif isnan(tmpparams.baseline)
        unitPower = '\muV^{2}/Hz';
    else
        unitPower = '% of baseline';
    end;
end;
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

function unitPower = newtimefpowerunit(tmpparams)

if nargin < 1
    help newtimefunit;
    return;
end

if ~isfield(tmpparams, 'baseline'), tmpparams.baseline = 0;     end
if ~isfield(tmpparams, 'scale'   ), tmpparams.scale    = 'log'; end
if ~isfield(tmpparams, 'basenorm'), tmpparams.basenorm = 'off'; end
if strcmpi(tmpparams.scale, 'log')
    if strcmpi(tmpparams.basenorm, 'on')
        unitPower = '10*log(std.)'; % impossible
    elseif isnan(tmpparams.baseline)
        unitPower = '10*log10(\muV^{2}/Hz)';
    else
        unitPower = 'dB';
    end
else
    if strcmpi(tmpparams.basenorm, 'on')
        unitPower = 'std.';
    elseif isnan(tmpparams.baseline)
        unitPower = '\muV^{2}/Hz';
    else
        unitPower = '% of baseline';
    end
end

% std_indvarmatch - match independent variable value in a list of values
%            
% Usage:
%   indices = std_indvarmatch(value, valuelist);
%
% Input:
%   value     - [string|real|cell] value to be matched
%   valuelist - [cell array] cell array of string, numerical values or
%               cell array
%
% Output:
%   indices   - [integer] numerical indices
%
% Example:
%   std_indvarmatch( 3, { 3 4 [2 3] });
%   std_indvarmatch( [2 3], { 3 4 [2 3] });
%   std_indvarmatch( [2 3], { 3 4 2 4 3 });
%   std_indvarmatch( 'test1', { 'test1' 'test2' { 'test1' 'test2' } });
%   std_indvarmatch( { 'test1' 'test2' }, { 'test1' 'test2' { 'test1' 'test2' } });
%   std_indvarmatch( { 'test1' 'test2' }, { 'test1' 'test2' 'test3' 'test1' });
%
% Author: Arnaud Delorme, CERCO/CNRS, UCSD, 2009-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function indices = std_indvarmatch(val, allvals);

indices = [];
if nargin < 1
    help std_indvarmatch;
    return;
end

% match values
% ------------
if all(cellfun(@isstr, allvals)) % string
    if ~iscell(val)
        indices = strmatch( val, allvals, 'exact')';
    else
        for indcell = 1:length(val)
            indices = [ indices std_indvarmatch(val{indcell}, allvals) ];
        end
    end
elseif all(cellfun(@length, allvals) == 1) % numerical
    if length(val) == 1
        indices = find( val == [allvals{:}]);
    else
        for ind = 1:length(val)
            indices = [ indices std_indvarmatch(val(ind), allvals) ];
        end
    end
else
    % mixed with cell array
    indices = [];
    for index = 1:length(allvals)
        if isequal(val, allvals{index})
            indices = [indices index]; 
        end
    end
end
return;

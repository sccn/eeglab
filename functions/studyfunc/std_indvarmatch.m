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

function indices = std_indvarmatch(val, allvals);

indices = [];
if nargin < 1
    help std_indvarmatch;
    return;
end;

% match values
% ------------
if all(cellfun(@isstr, allvals)) % string
    if ~iscell(val)
        indices = strmatch( val, allvals, 'exact')';
    else
        for indcell = 1:length(val)
            indices = [ indices std_indvarmatch(val{indcell}, allvals) ];
        end;
    end;
elseif all(cellfun(@length, allvals) == 1) % numerical
    if length(val) == 1
        indices = find( val == [allvals{:}]);
    else
        for ind = 1:length(val)
            indices = [ indices std_indvarmatch(val(ind), allvals) ];
        end;
    end;
else
    % mixed with cell array
    indices = [];
    for index = 1:length(allvals)
        if isequal(val, allvals{index})
            indices = [indices index]; 
        end;
    end;
end;
return;

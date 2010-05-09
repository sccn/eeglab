% std_buildsetcomps - rebuild sets and comps array for component clusters
%                     this is original building for 
%            
% Usage:
%   STUDY = std_buildsetcomps(STUDY, ALLEEG);
%
% Author: Arnaud Delorme, CERCO/CNRS, UCSD, 2010-

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

% $Log: std_cell2setcomps.m,v $
% Revision 1.1  2009/10/20 02:28:35  arno
% Updated conversion between sets and indices formats
%

function [ STUDY ] = std_buildsetcomps(STUDY, ALLEEG, cind)

if nargin < 2
    help std_buildsetcomps;
    return;
end;
if nargin < 3
    cind = 1:length(STUDY.cluster);
end;

% find datasets with common ICA decompositions
icaclusters = std_findsameica(ALLEEG);

for c = 2 %1:length(STUDY.cluster)
   for comp = 1:length(STUDY.cluster(c).comps)
       ind = find(cellfun(@(x)(~isempty(find(x == STUDY.cluster(c).sets(1,comp)))), icaclusters));
       if length(ind) ~= 1, error('Dataset not found'); end;
       tmpsets = icaclusters{ind};
       newsets(1:length(tmpsets),comp) = tmpsets';
   end;
   newsets(find(newsets == 0)) = NaN;
   STUDY.cluster(c).sets = newsets;
end;

    
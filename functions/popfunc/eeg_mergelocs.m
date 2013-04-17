% eeg_mergelocs() - merge channel structure while preserving channel
%                   order
%
%      >> mergedlocs = eeg_mergelocs(loc1, loc2, loc3, ...);
%
% Inputs: 
%     loc1     - EEGLAB channel location structure
%     loc2     - second EEGLAB channel location structure
%
% Output: 
%     mergedlocs - merged channel location structure
%     warning    - [0|1] dissimilar structures found (0=false, 1=true)
%
% Author: Arnaud Delorme, August 2006

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
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

function [alllocs warn] = eeg_mergelocs(varargin)

persistent warning_shown;
warn = 0;

try
    % sort by length
    % --------------
    len = cellfun(@length, varargin);
    [tmp so] = sort(len, 2, 'descend');
    varargin = varargin(so);

    alllocs = varargin{1};
    for index = 2:length(varargin)

        % fuse while preserving order (assumes the same channel order)
        % ------------------------------------------------------------
        tmplocs = varargin{index};
        newlocs = myunion(alllocs, tmplocs);

        if length(newlocs) > length(union({ alllocs.labels }, { tmplocs.labels }))

            warn = 1;
            if isempty(warning_shown)
                disp('Warning: different channel montage or electrode order for the different datasets');
                warning_shown = 1;
            end;

            % trying to preserve order of the longest array
            %----------------------------------------------
            if length(alllocs) < length(tmplocs)
                tmp     = alllocs;
                alllocs = tmplocs;
                tmplocs = tmp;
            end;
            allchans = { alllocs.labels tmplocs.labels };
            [uniquechan ord1 ord2 ]  = unique_bc( allchans );

            [tmp rminds] = intersect_bc( uniquechan, { alllocs.labels });
            ord1(rminds) = [];
            tmplocsind = ord1-length(alllocs);

            newlocs = [ alllocs tmplocs(tmplocsind) ];

        end;
        alllocs = newlocs;
    end;
catch,
    % temporary fix for dissimilar structures
    % should check channel structure consistency instead
    % using checkchan function
    disp('Channel merging warning: dissimilar fields in the two structures');
    [alllocs warn ] = eeg_mergelocs_diffstruct(varargin{:});
end;

% union of two channel location structure
% without loosing the order information
% ---------------------------------------
function alllocs = myunion(locs1, locs2)

   labs1 = { locs1.labels };
   labs2 = { locs2.labels };
   
   count1 = 1;
   count2 = 1;
   count3 = 1;
   alllocs = locs1; alllocs(:) = [];
   while count1 <= length(locs1) || count2 <= length(locs2)
       
       if count1 > length(locs1)
           alllocs(count3) = locs2(count2);
           count2 = count2 + 1;
           count3 = count3 + 1;
       elseif count2 > length(locs2)
           alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count3 = count3 + 1;
       elseif strcmpi(labs1{count1}, labs2{count2})
           alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count2 = count2 + 1;
           count3 = count3 + 1;
       elseif isempty(strmatch(labs1{count1}, labs2, 'exact'))
           alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count3 = count3 + 1;
       else 
           alllocs(count3) = locs2(count2);
           count2 = count2 + 1;
           count3 = count3 + 1;
       end;
       
   end;

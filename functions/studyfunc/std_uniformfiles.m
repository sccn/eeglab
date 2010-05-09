% std_uniformfiles() - Check uniform channel distribution accross data files
%
% Usage:    
%     >> boolval = std_uniformfiles(STUDY, ALLEEG);   
%
% Inputs:
%   STUDY      - EEGLAB STUDY
%
% Outputs:
%   boolval    - [-1|0|1] 0 if non uniform, 1 if uniform, -1 if error
%
% Authors: Arnaud Delorme, SCCN/UCSD, CERCO/CNRS, 2010-

% Copyright (C) Arnaud Delorme
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

function uniformlabels = std_uniformfiles( STUDY, ALLEEG );

if nargin < 2
    help std_checkconsist;
    return;
end;

filetypes = { 'daterp' 'datspec' 'datersp' 'datitc' 'dattimef' };

% set channel interpolation mode
% ------------------------------
uniformlabels = [];

count = 1;
selectedtypes = {};
for index = 1:length(filetypes)
    
    % scan datasets
    % -------------
    [ tmpstruct compinds filepresent ] = std_fileinfo(ALLEEG, filetypes{index});
    
    % check if the structures are equal
    % ---------------------------------
    if ~isempty(tmpstruct)
        if any(filepresent == 0)
            uniformlabels = -1;
            return;
        else
            firstval = tmpstruct(1).labels;
            uniformlabels(count) = length(firstval);
            for dat = 2:length(tmpstruct)
                tmpval = tmpstruct(dat).labels;
                if ~isequal(firstval, tmpval)
                     uniformlabels(count) = 0;
                end;
            end;
            selectedtypes{count} = filetypes{index};
            count = count+1;
        end;
    end;
end;

% check output
% ------------
if any(uniformlabels == 0)
    if any(uniformlabels) == 1
        disp('Warning: conflict between datafiles - all should have the same channel distribution');
        for index = 1:length(selectedtypes)
            if uniformlabels(index)
                 fprintf('  Data files with extention "%s" have interpolated channels\n',  selectedtypes{index});
            else fprintf('  Data files with extention "%s" have non-interpolated channels\n', selectedtypes{index});
            end;
        end;
        uniformlabels = -1;
    else
        uniformlabels = 0;
    end;
else
    if ~(length(unique(uniformlabels)) == 1)
        fprintf('Warning: Data files have different number of channel in each file\n');
        for index = 1:length(selectedtypes)
            fprintf('  Data files with extention "%s" have %d channels\n',  selectedtypes{index}, uniformlabels(index));
        end;
    end;
    uniformlabels = 1;
end;

% std_renamestudyfiles() - rename files for design 1 if necessary. In design
%                          1, for backward compatibility, files could have
%                          legacy names. For consistency these files now
%                          need to be renamed. Note that the STUDY is
%                          automatically resave on disk to avoid any potential
%                          inconsistency.
%
% Usage:
%   >> STUDY = std_renamestudyfiles(STUDY, ALLEEG)
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of the EEG datasets included in the STUDY structure
%
% Output:
%      STUDY - The input STUDY with new design files.
%
% Author: Arnaud Delorme, Institute for Neural Computation UCSD, 2013-

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

function STUDY = std_renamestudyfiles(STUDY, ALLEEG)

if nargin < 2
    help std_renamestudyfiles;
    return;
end;

STUDY2 = std_makedesign(STUDY, ALLEEG, 1, STUDY.design(1), 'defaultdesign', 'forceoff', 'verbose', 'off');
allCell1 = {  STUDY.design(1).cell.filebase };
allCell2 = { STUDY2.design(1).cell.filebase };

fileExtensions = { 'daterp' 'datspec' 'datersp' 'daterpim' 'dattimef' 'datitc' 'daterpim' ...
    'icaerp' 'icaspec' 'icaersp' 'icaerpim' 'icatimef' 'icaitc' 'icaerpim' };

if ~isequal(allCell1, allCell2)
    
    thereIsAFileNotDesign = false;
    for index = 1:length(allCell1), if length(allCell1{index}) < 6 || all(allCell1{index}(1:6) == 'design'), thereIsAFileNotDesign = true; end; end;
    for index = 1:length(allCell1), if length(allCell1{index}) < 6 || all(allCell1{index}(1:6) == 'design'), thereIsAFileNotDesign = true; end; end;
    
    if thereIsAFileNotDesign
        res = questdlg2(['Old STUDY design data files have been detected.' 10 ...
            'EEGLAB wants to rename these files to improve consistency' 10 ...
            'and stability. No dataset will be renamed, only preprocessed' 10 ...
            'STUDY data files.' ], 'Rename STUDY data files', 'Cancel', ...
            'Rename', 'Rename');
        
        if strcmpi(res, 'rename')
            STUDY = pop_savestudy(STUDY, ALLEEG, 'savemode', 'resave');
            for iCell = 1:length(allCell1)
                
                % scan file extensions
                for iExt = 1:length(fileExtensions)
                    files = dir( [ allCell1{iCell} '.' fileExtensions{iExt} ]);
                    if ~isempty(files) && ~strcmpi(allCell1{iCell}, allCell2{iCell})
                        movefile( [ allCell1{iCell} '.' fileExtensions{iExt} ], ...
                            [ allCell2{iCell} '.' fileExtensions{iExt} ]);
                        disp([ 'Moving ' [ allCell1{iCell} '.' fileExtensions{iExt} ] ' to ' [ allCell2{iCell} '.' fileExtensions{iExt} ] ]);
                    end;
                end;
            end;
            STUDY = STUDY2;
            STUDY = pop_savestudy(STUDY, ALLEEG, 'savemode', 'resave');
        end;
    end;
end;

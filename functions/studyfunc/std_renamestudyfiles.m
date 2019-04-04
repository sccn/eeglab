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

function STUDY = std_renamestudyfiles(STUDY, ALLEEG)

if nargin < 2
    help std_renamestudyfiles;
    return;
end

if ~isfield(STUDY.design(1), 'cell') || isempty(STUDY.design(1).cell), return; end
STUDY2 = std_makedesign(STUDY, ALLEEG, 1, STUDY.design(1), 'defaultdesign', 'forceoff', 'verbose', 'off');
allCell1 = {  STUDY.design(1).cell.filebase };
allCell2 = { STUDY2.design(1).cell.filebase };

fileExtensions = { 'daterp' 'datspec' 'datersp' 'daterpim' 'dattimef' 'datitc' 'daterpim' ...
    'icaerp' 'icaspec' 'icaersp' 'icaerpim' 'icatimef' 'icaitc' 'icaerpim' };

if ~isequal(allCell1, allCell2)
    
    thereIsAFileNotDesign = false;
    for index = 1:length(allCell1), if length(allCell1{index}) < 6 || all(allCell1{index}(1:6) == 'design'), thereIsAFileNotDesign = true; end; end
    for index = 1:length(allCell1), if length(allCell1{index}) < 6 || all(allCell1{index}(1:6) == 'design'), thereIsAFileNotDesign = true; end; end
    
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
                    end
                end
            end
            STUDY = STUDY2;
            STUDY = pop_savestudy(STUDY, ALLEEG, 'savemode', 'resave');
        end
    end
end

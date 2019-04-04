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

function uniformlabels = std_uniformfiles( STUDY, ALLEEG );

if nargin < 2
    help std_checkconsist;
    return;
end

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
                end
            end
            selectedtypes{count} = filetypes{index};
            count = count+1;
        end
    end
end

% check output
% ------------
if any(uniformlabels == 0)
    if any(uniformlabels) == 1
        disp('Warning: conflict between datafiles - all should have the same channel distribution');
        for index = 1:length(selectedtypes)
            if uniformlabels(index)
                 fprintf('  Data files with extention "%s" have interpolated channels\n',  selectedtypes{index});
            else fprintf('  Data files with extention "%s" have non-interpolated channels\n', selectedtypes{index});
            end
        end
        uniformlabels = -1;
    else
        uniformlabels = 0;
    end
else
    if ~(length(unique(uniformlabels)) == 1)
        fprintf('Warning: Data files have different number of channel in each file\n');
        for index = 1:length(selectedtypes)
            fprintf('  Data files with extention "%s" have %d channels\n',  selectedtypes{index}, uniformlabels(index));
        end
    end
    uniformlabels = 1;
end

% std_fileinfo() - Check uniform channel distribution across datasets
%
% Usage:    
%   >> [struct filepresent] = std_fileinfo(ALLEEG);   
% Inputs:
%   ALLEEG    - EEGLAB ALLEEG structure
%
% Outputs:
%   struct      - structure of the same length as the ALLEEG variable
%                 containing all the fields in the datafiles
%   filepresent - array of 0 and 1 indicating if the file is present for
%                 each dataset
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

function [tmpstructout compinds filepresent] = std_fileinfo( ALLEEG, filetype );

    firstpass   = 1;
    notequal = 0;
    compinds = {};
    tmpstructout = [];
    filepresent = zeros(1,length(ALLEEG));
    for dat = 1:length(ALLEEG)
        
        filename = fullfile(ALLEEG(dat).filepath, [ ALLEEG(dat).filename(1:end-3) filetype ]);
        thisfilepresent = exist(filename);
        if thisfilepresent && firstpass == 1
            %fprintf('Files of type "%s" detected, checking...',  filetype);
        elseif firstpass == 1
            notequal = 1;
        end
        firstpass = 0;
        filepresent(dat) = thisfilepresent;
        
        if filepresent(dat)
            try
                warning('off', 'MATLAB:load:variableNotFound');
                tmptmpstruct = load( '-mat', filename, 'times', 'freqs', 'parameters', 'labels', 'chanlabels' );
                warning('on', 'MATLAB:load:variableNotFound');
            catch
                passall = 0;
                fprintf(' Error\n');
                break;
            end
            
            % rename chanlabels and store structure
            if isfield(tmptmpstruct, 'chanlabels')
                tmptmpstruct.labels = tmptmpstruct.chanlabels;
                tmptmpstruct = rmfield(tmptmpstruct, 'chanlabels');
            end
            
            if filetype(1) ~= 'd' % ICA components
                allvars = whos('-file', filename);
                tmpinds = [];
                for cind = 1:length(allvars)
                    str = allvars(cind).name(5:end);
                    ind_ = find(str == '_');
                    if ~isempty(ind_), str(ind_:end) = []; end
                    tmpinds = [ tmpinds str2num(str) ];
                end
                compinds(dat) = { unique(tmpinds) };
            elseif ~isfield(tmptmpstruct, 'labels')
                allvars = whos('-file', filename);
                allvarnames = { allvars.name };
                tmptmpstruct.labels = allvarnames(strmatch('chan', allvarnames));
            end
            try,
                tmpstruct(dat) = tmptmpstruct;
            catch,
                passall = 0;
                break;
            end
        end
    end
    if exist('tmpstruct') == 1
        tmpstructout = tmpstruct;
    end

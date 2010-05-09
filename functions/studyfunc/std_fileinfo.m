% std_fileinfo() - Check uniform channel distribution accross datasets
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
        end;
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
            end;
            
            % rename chanlabels and store structure
            if isfield(tmptmpstruct, 'chanlabels')
                tmptmpstruct.labels = tmptmpstruct.chanlabels;
                tmptmpstruct = rmfield(tmptmpstruct, 'chanlabels');
            end;
            
            if filetype(1) ~= 'd' % ICA components
                allvars = whos('-file', filename);
                tmpinds = [];
                for cind = 1:length(allvars)
                    str = allvars(cind).name(5:end);
                    ind_ = find(str == '_');
                    if ~isempty(ind_), str(ind_:end) = []; end;
                    tmpinds = [ tmpinds str2num(str) ];
                end;
                compinds(dat) = { unique(tmpinds) };
            elseif ~isfield(tmptmpstruct, 'labels')
                allvars = whos('-file', filename);
                allvarnames = { allvars.name };
                tmptmpstruct.labels = allvarnames(strmatch('chan', allvarnames));
            end;
            try,
                tmpstruct(dat) = tmptmpstruct;
            catch,
                passall = 0;
                break;
            end;
        end;
    end;
    if exist('tmpstruct') == 1
        tmpstructout = tmpstruct;
    end;

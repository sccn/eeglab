% readelp() - read electrode locations from an .elp (electrode positions)
%             file as generated, for example, by a Polhemus tracking device 
% Usage:
%   >> [eloc, elocnames, X, Y, Z] = readelp(filename);
%
% Inputs:
%   filename - name of the .elp file containing cartesian (XYZ) electrode
%              locations
% Outputs:
%   eloc      - structure containing the names and locations of the channels
%   elocnames - cell array containing the names of the electrodes
%   X,Y,Z     - vectors containing the xyz locations of the electrodes
%
% Author: Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
%
% Note: ignores comments and the sensor type field
% Note: convert output XYZ locations to polar coordinates using cart2pol()
%
% See also: readlocs(), snapread(), floatread(), cart2pol()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
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

function [eloc, names, x, y, z] = readelp( filename ); 

if nargin < 1
	help readelp;
	return;
end

% open file
% ---------
fid = fopen(filename, 'r');
if fid == -1
  disp('Cannot open file'); return;
end

index = 1;
countfid  = 1;
fidlabels = { 'Nz' 'LPA' 'RPA' };

needToReadNumbers = 0;
tmpstr = fgetl(fid);

while 1
    if needToReadNumbers==1 % Has sparsed $N.
        if (~isempty(str2num(tmpstr)))
            tmp = sscanf(tmpstr, '%f');
            
            eloc(index).X  = tmp(1); x(index) = tmp(1);
            eloc(index).Y  = tmp(2); y(index) = tmp(2);
            eloc(index).Z  = tmp(3); z(index) = tmp(3);
            eloc(index).type = 'EEG';
            
            index = index + 1;
            needToReadNumbers = 0;
        end
    elseif tmpstr(1) == '%'
        if tmpstr(2) == 'F' % fiducial
            tmp = sscanf(tmpstr(3:end), '%f');
                    
            eloc(index).labels = fidlabels{countfid};
            eloc(index).X  = tmp(1); x(index) = tmp(1);
            eloc(index).Y  = tmp(2); y(index) = tmp(2);
            eloc(index).Z  = tmp(3); z(index) = tmp(3);
            eloc(index).type = 'FID';
            
            index     = index    + 1;
            countfid  = countfid + 1;
                    
        elseif tmpstr(2) == 'N' % regular channel
            nm = strtok(tmpstr(3:end));
            if ~(strcmp(nm, 'Name') || strcmp(nm, 'Status'))
                eloc(index).labels = nm;
                needToReadNumbers = 1;
            end
        end
    end

    % Get the next line
    tmpstr = fgetl(fid);
    while isempty(tmpstr)
        tmpstr = fgetl(fid);
        if ~ischar(tmpstr) && tmpstr == -1
            break;
        end; 
    end
    
    if ~ischar(tmpstr) && tmpstr == -1
        break;
    end;  
end;  

names = { eloc.labels };

fclose(fid);    % close the file

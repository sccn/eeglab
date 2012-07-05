% pop_loadascinstep() - import an INStep ASC (ASCII) file
%
% Usage:
%   >> OUTEEG = pop_loadascinstep( filename );
%
% Inputs:
%   filename       - file name
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, July 2006

% Copyright (C) 2006 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] = pop_loadascinstep( filename );

if nargin < 1
    [filename, filepath] = uigetfile('*.*', 'Choose an INStep ASC file -- pop_loadascinstep'); 
    drawnow;
    if filename == 0 return; end;
    filename = fullfile(filepath, filename);
end;

EEG = eeg_emptyset;
fid = fopen(filename, 'r');
if fid == -1, error('Cannot open file'); end;
num1 = fscanf(fid, '%f', 1);tmp   = fgetl(fid);
num2 = fscanf(fid, '%f', 1);tmp   = fgetl(fid);
num3 = fscanf(fid, '%f', 1);tmp   = fgetl(fid);

if ~isempty(num3)
    EEG.nbchan = num1;
    EEG.pnts   = num2;
    EEG.srate  = 1000/num3;
    tline = fgetl(fid);
    tline = fgetl(fid);
else 
    EEG.nbchan = num1;
    EEG.srate  = 1000/num2;
    tline = tmp;
end;

allf  = parsetxt(tline);

EEG.data = fscanf(fid, '%f', [EEG.nbchan+5 Inf]);
fclose(fid);

EEG.xmin = EEG.data(3,1)/1000;
EEG.data(4,:) = EEG.data(4,:).*EEG.data(2,:); % stimulus times category
EEG.data([1:3],:) = []; % trial column
EEG.chanlocs = struct('labels', allf(4:end));
EEG = eeg_checkset(EEG);

EEG = pop_chanevent(EEG, 1, 'edge', 'leading', 'delchan', 'on');
EEG = pop_chanevent(EEG, 1, 'edge', 'leading', 'delchan', 'on', 'delevent', 'off', 'nbtype', 1, 'typename', 'resp' );
EEG = eeg_checkset(EEG, 'eventconsistency');

com = sprintf('EEG = pop_loadascinstep(''%s'');',filename);

% scanfold() - scan folder content
%
% Usage:    
%    >> [cellres textres] = scanfold(foldname);
%    >> [cellres textres] = scanfold(foldname, ignorelist, maxdepth);
%
% Inputs:
%   foldname   - [string] name of the folder
%   ignorelist - [cell] list of folders to ignore
%   maxdepth   - [integer] maximum folder depth
%
% Outputs:
%   cellres   - cell array containing all the files
%   textres   - string array containing all the names preceeded by "-a"
% 
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2009

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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

function [ cellres, textres ] = scanfold(foldname, ignorelist, maxdepth)

if nargin < 2, ignorelist = {}; end;
if nargin < 3, maxdepth = 100; end;
foldcontent = dir(foldname);
textres = '';
cellres = {};
if maxdepth == 0, return; end;
for i = 1:length(foldcontent)
    if (exist(foldcontent(i).name) == 7 || strcmpi(foldcontent(i).name, 'functions')) && ~ismember(foldcontent(i).name, ignorelist)
        if ~strcmpi(foldcontent(i).name, '..') && ~strcmpi(foldcontent(i).name, '.')
            disp(fullfile(foldname, foldcontent(i).name));
            [tmpcellres tmpres] = scanfold(fullfile(foldname, foldcontent(i).name), ignorelist, maxdepth-1);
            textres = [ textres tmpres ];
            cellres = { cellres{:} tmpcellres{:} };
        end;
    elseif length(foldcontent(i).name) > 2
        if strcmpi(foldcontent(i).name(end-1:end), '.m')
            textres = [ textres ' -a ' foldcontent(i).name ];
            cellres = { cellres{:} foldcontent(i).name };
        end;
    else 
        disp( [ 'Skipping ' fullfile(foldname, foldcontent(i).name) ]);
    end;
end;
return;

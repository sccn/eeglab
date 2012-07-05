% eeg_getversion() - obtain EEGLAB version number
%
% Usage:
%     >> vers = eeg_getversion;
%     >> [vers vnum] = eeg_getversion;
%
% Outputs:
%    vers = [string] EEGLAB version number
%    vnum = [float] numerical value for the version. For example 11.3.2.4b
%           is converted to 11.324
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, 2010

% Copyright (C) 2010  Arnaud Delorme, SCCN/INC/UCSD
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

function [vers versnum releaseDate] = eeg_getversion;

vers = '';
filepath = fileparts(which('eeglab.m'));
filename = dir(fullfile(filepath, 'Contents.m'));

releaseDate = filename.date;

if isempty(filename), return; end;

fid = fopen(fullfile(filepath, filename.name), 'r');
fgetl(fid);
versionline = fgetl(fid);
vers = versionline(11:end);
fclose(fid);

tmpvers = vers;
if isempty(str2num(tmpvers(end))), tmpvers(end) = []; end;
indsDot = find(tmpvers == '.' );
tmpvers(indsDot(2:end)) = [];
versnum = str2num(tmpvers);
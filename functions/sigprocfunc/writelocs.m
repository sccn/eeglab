% writelocs() - write electrode positions (expanded from the ICA 
%              toolbox function)
%             
% Usage:
%   >> writelocs( chanstruct, filename );
%   >> writelocs( chanstruct, filename, 'key', 'val' );
%
% Inputs:
%   chanstruct - channel structure retuned by readlocs()
%   filename   - file name for saving electrode locations
%
% Optional inputs:
%   'filetype'  - ['loc'|'sph'|'xyz'|'besa'|'chanedit'] type of 
%                 the file to be read. By default the file type is
%                 determined using the file extension.
%                  'loc' is the eeglab() format (see below)
%                  'sph' is a matlab spherical coordinate file (spherical
%                        coordinates in Matlab are different from Besa 
%                        spherical coordinates).
%                  'xyz' contain carhtesian coordinates of the electrodes
%                  'besa' besa '.elp' coordinate file (contact author for
%                        reading other Besa files). '.elp' Besa files
%                        have only 3 columns of 'theta' angle, 'phi' angle
%                        and 'weight' (which is ignored). Also note that 
%                        'theta' and 'phi' meaning is also the inverse
%                        as the one they have in Matlab (in Besa 'theta' =
%                        azimut and 'phi' = horiz).
%                  'chanedit' files created by pop_chanedit().
%                  'custom' specify your own format using the 'format' entry
%                 The Default file type is 'loc'.
%   'format'    - [cell array] for 'custom' file types only (or if no
%                 file type is defined). The input cell array contains one
%                 string entry per column. The string name can be 
%                     'channum'   for channel number. 
%                     'labels'    for channel name. 
%                     'theta'     for channel angle polar coordinate. 
%                     'radius'    for channel radius polar coordinate. 
%                     'X'         for channel cathesian coordinate X. 
%                     'Y'         for channel cathesian coordinate Y. 
%                     'Z'         for channel cathesian coordinate Z. 
%                     '-X'        for minus channel cathesian coordinate X. 
%                     '-X'        for minus channel cathesian coordinate X. 
%                     '-Z'        for minus channel cathesian coordinate Z.  
%                     'sph_theta' for Matlab channel spherical theta angle
%                                 (=horizontal). 
%                     'sph_phi'   for Matlab channel spherical phi angle
%                                 (=azimut). 
%                     'sph_radius' for channel spherical radius. 
%                     'sph_theta_besa' for BESA channel spherical theta
%                                 angle (=horizontal).  
%                     'sph_phi_besa' for BESA channel spherical phi 
%                                 angle (=azimut). 
%                     'calib'     for the channel calibration value. 
%                     'gain'      for the channel gain. 
%                     'custom1'   custom field 1.
%                     'custom2', 'custom3', 'custom4' other custom fields.
%   'header'    - ['on'|'off'] add a header with the name of each column.
%                 Default is 'off'.
%   'customheader' - [string] add a custom header at the beginning of the
%                 file. If used with 'header' set to 'on', the column names
%                 will be insterted after the custom header.
%   'elecind'   - [integer array] indices of electrode to export. Default 
%                 is all ellectrodes.
%
% Note: for file formats, see readlocs() help (>> help readlocs)
%
% Author: Arnaud Delorme, Salk Institute, 16 Dec 2002
%
% See also: readlocs(), readelp()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
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

% $Log: not supported by cvs2svn $

function writelocs( chans, filename, varargin ); 

if nargin < 2
	help writelocs;
	return;
end;

chans = convertelocs(chans, 'auto');

% get infos from readlocs
% -----------------------
[listtype formatinfo listcolformat formatskip] = readlocs('getinfos');

g = finputcheck( varargin, ...
   { 'filetype'	'string'	 listtype 			'loc';
     'header'     'string' { 'on' 'off' } 	'off';
     'customheader' 'string' [] 					'';
     'elecind'    'integer' [1 Inf]				[];
     'format'		'cell'	 []					{} }, 'writelocs');
if isstr(g), error(g); end;  

% select channels
% ---------------
if ~isempty(g.elecind)
	chans = chans(g.elecind);
end;

% finding types of input
% ----------------------
if isempty(g.format)
   indexformat = strmatch(lower(g.filetype), listtype, 'exact');
   g.format = formatinfo{indexformat};
   g.skipline = formatskip(indexformat);
else 
   g.skipline = 0;   
end;

% creating file
% -------------
fid = fopen(filename, 'w');

% exporting header
% ----------------
if ~isempty(g.customheader)
   fprintf(fid, '%s\n', g.customheader);
end;
if  strcmpi(g.header, 'on') | g.skipline == 2
   for index=1:length(g.format)
      fprintf(fid, '%s', g.format{index});
      if index == length(g.format)
         fprintf(fid, '\t');
      end;         
   end;
   for index=1:length(g.format)
      fprintf(fid, '%s', char(ones(1,length(g.format{index}))*45));
      if index == length(g.format)
         fprintf(fid, '\t');
      end;         
   end;
end;
if g.skipline == 1
   fprintf(fid, '%d\n', length(chans));
end;         

% writing infos
% -------------
for indexchan = 1:length(chans)
   for index=1:length(g.format)
      [str, mult] = checkformat(g.format{index});
      if strcmpi(str, 'channum')
         fprintf(fid, '%d', indexchan);
      else
         if ~isfield(chans, str)
            error([ 'Non-existant field: ''' str '''' ]);
         end;
         eval( [ 'chanval = chans(indexchan).' str ';' ] );
         if   isstr(chanval), fprintf(fid, '%s', chanval);
         else   					fprintf(fid, '%s', num2str(mult*chanval));
         end;
      end;
      if index ~= length(g.format)
         fprintf(fid, '\t');
      end;         
   end;
   fprintf(fid, '\n');
end;
fclose(fid);

return;

% check field format
% ------------------
function [str, mult] = checkformat(str)
	mult = 1;
	if strcmpi(str, 'labels'), str = lower(str); return; end;
	if strcmpi(str, 'channum'), str = lower(str); return; end;
	if strcmpi(str, 'theta'), str = lower(str); return; end;
	if strcmpi(str, 'radius'), str = lower(str); return; end;
	if strcmpi(str, 'sph_theta'), str = lower(str); return; end;
	if strcmpi(str, 'sph_phi'), str = lower(str); return; end;
	if strcmpi(str, 'sph_radius'), str = lower(str); return; end;
	if strcmpi(str, 'sph_theta_besa'), str = lower(str); return; end;
	if strcmpi(str, 'sph_phi_besa'), str = lower(str); return; end;
	if strcmpi(str, 'gain'), str = lower(str); return; end;
	if strcmpi(str, 'calib'), str = lower(str); return; end;
	if strcmpi(str, 'X'), str = upper(str); return; end;
	if strcmpi(str, 'Y'), str = upper(str); return; end;
	if strcmpi(str, 'Z'), str = upper(str); return; end;
	if strcmpi(str, '-X'), str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Y'), str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Z'), str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, 'custum1'), return; end;
	if strcmpi(str, 'custum2'), return; end;
	if strcmpi(str, 'custum3'), return; end;
	if strcmpi(str, 'custum4'), return; end;
   error(['writelocs: undefined field ''' str '''']);
   

% readlocs() - read electrode positions (expanded from the ICA 
%              toolbox function)
%             
% Usage:
%   >> [eloc, labels, theta, radius] = readlocs( filename );
%   >> [eloc, labels, theta, radius] = readlocs( filename, ...
%                                          'key', 'val' );
%
% Inputs:
%   filename   - name of the file containing the electrode locations
%                and convert in polar coordinates
%
% Optional inputs:
%   'filetype'  - ['loc'|'sph'|'xyz'|'polhemus'|'besa'|'chanedit'|'custom'] 
%                 type of the file to be read. By default the file type is
%                 determined using the file extension.
%                  'loc' is the eeglab() format (see below)
%                  'sph' is a matlab spherical coordinate file (spherical
%                        coordinates in Matlab are different from Besa 
%                        spherical coordinates).
%                  'xyz' contain carhtesian coordinates of the electrodes
%                  'polhemus' or 'polhemusx' polhemus file with electrodes 
%                        recorded using X as the main axis (see readelp() )
%                  'polhemusy' polhemus file with electrodes recorded
%                        using Y as the main axis (see readelp() )
%                  'besa' besa '.elp' coordinate file (contact author for
%                        reading other Besa files). '.elp' Besa files
%                        have only 3 columns of 'theta' angle, 'phi' angle
%                        and 'weight' (which is ignored). Also note that 
%                        'theta' and 'phi' meaning is also the inverse
%                        as the one they have in Matlab (in Besa 'theta' =
%                        azimut and 'phi' = horiz).
%                  'chanedit' files created by pop_chanedit().
%                  'custom' specify your own format using the 'format' entry
%   'format'    - [cell array] for 'custom' file types (or if no
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
%   'skipline'  - [integer] number of header lines to skip for custom file
%                 types only.
%   'elecind'   - [integer array] indices of electrode to read 
%   'center'    - [integer array or 'auto'] for X, Y, Z coordinates, specify
%                the center of the sphere. 'auto' uses the best matching
%                that passes through all the electrodes. Default is [0 0 0].
%
% File formats:
%   The extension of the file determines its type
%   if the 'filetype' parameter is not specified
%   '.loc' or '.locs':
%               polar format. Example
%               1    -18    .352       Fp1
%               2     18    .352       Fp2
%               3    -90    .181       C3
%               4     90    .181       C4
%                 more lines ...
%               Note that in previous releases, the channel labels
%               had to contain 4 characters (spaces replaced by '.').
%               This format still works but the dots are no longer
%               required.
%   '.sph':
%               spherical coordinate file. Example
%               1    -63.36    -72      Fp1
%               2     63.36    72       Fp2
%               3     32.58    0        C3
%               4     32.58    0        C4
%                 more lines ...
%   '.xyz': 
%               cartesian coordinate file. Example
%               1   -0.8355   -0.2192   -0.5039      Fp1
%               2   -0.8355    0.2192    0.5039      Fp2
%               3    0.3956         0   -0.9184      C3
%               4    0.3956         0    0.9184      C4
%                 more lines ...
%   '.txt':   
%               read ascii files saved using pop_chanedit()
%   '.elp':     
%               Polhemus coordinate file (uses readelp())
%
% Outputs:
%   eloc      - structure containing the channel names and locations.
%               It has three fields 'labels', 'theta' and 'radius'.
%   labels    - cell array of strings giving the  names of the electrodes
%   theta     - vector of polar angles of the electrodes in degrees
%   radius    - vector of polar norms of the electrodes
%
% Note: If the channel label are not given, they are given the name
%       'E1', 'E2', etc.
%       ************************** TO CHANGE HERE
%       the function cart2topo() is used to convert cartesian to polar
%       coordinates. By default the current function uses cart2topo()
%       options to recompute the best center coordinates.
%
% Author: Arnaud Delorme, Salk Institute, 8 Dec 2002
%
% See also: readelp(), writelocs(), topo2sph(), sph2topo(), cart2topo(), sph2cart()

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

% To DO: remove the cart2topo, use cart2sph and sph2topo instead
% use chancenter to recenter data
%

function [eloc, labels, theta, radius] = readlocs( filename, varargin ); 

if nargin < 1
	help readlocs;
	return;
end;

% to add a new channel format
% ---------------------------
% 1) add the format name at the end of the listtype variable list
% 2) enter the column type in a new list at the end of the listimportformat variable
% 3) enter the number of lines to skip at the end of the listskipline array
% Note: these infos are also used by writelocs() and pop_readlocs() but
% you do not have to edit these functions.

listtype = { ...
   		'polhemus' ...
         'polhemusX' ...
         'polhemusY' ...
         'besa' ...
         'xyz' ...
         'loc' ...
         'sph' ...
         'chanedit' ...
         'custom' };
   
listimportformat = { ...
   	{ } ... % polhemus (specific non-columnar implementation)	
      { } ... % polhemus (specific non-columnar implementation)
      { } ... % polhemus (specific non-columnar implementation)
      { 'sph_theta_besa' 'sph_phi_besa' } ... % BESA/EGI format
      { 'channum' 'X' 'Y' 'Z' } ... % xyz format
      { 'channum' 'theta' 'radius' 'labels' } ... % loc format
      { 'channum' 'sph_theta' 'sph_radius' 'labels' } ... % sph format
      { 'labels'  'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' 'sph_radius' } }; %chanedit format

listcolformat = { 'labels' 'channum' 'theta' 'radius' 'sph_theta' 'sph_phi' ...
      'sph_radius' 'sph_theta_besa' 'sph_phi_besa' 'gain' 'calib' ...
      'X' 'Y' 'Z' '-X' '-Y' '-Z' 'custom1' 'custom2' 'custom3' 'custom4' 'not def' };

listskipline = [ ...
   0 ... % polhemus, not applicable
   0 ... % polhemus, not applicable
   0 ... % polhemus, not applicable
   1 ...  % besa
   0 ...
   0 ...
   0 ...
   2 ]; % skip the 2 lines header for the chanedit format

% ------------------------------------------------------
% special mode for getting the infos
% ----------------------------------
if isstr(filename) & strcmp(filename, 'getinfos')
   eloc = listtype;
   labels = listimportformat;
   theta = listcolformat;
   radius = listskipline;
   return;
end;

g = finputcheck( varargin, ...
   { 'filetype'	'string'	 listtype '';
     'skipline'   'integer' [0 Inf] 			[];
     'elecind'    'integer' [1 Inf]				[];
     'format'		'cell'	 []					{} }, 'readlocs');
if isstr(g), error(g); end;  

if isstr(filename)
   
   % format auto detection
	% ---------------------
   periods = find(filename == '.');
   fileextension = filename(periods(end)+1:end);
	if isempty(g.filetype)
		switch lower(fileextension),
		 case {'loc' 'locs' }, g.filetype = 'loc';
		 case 'xyz', g.filetype = 'xyz';
		 case 'sph', g.filetype = 'sph';
		 case 'txt', g.filetype = 'chanedit';
		 case 'elp', g.filetype = 'polhemus';
		 case 'eps', g.filetype = 'besa';
		 otherwise, g.filetype =  ''; 
      end;
      fprintf('Readlocs: ''%s'' format detected from file extension\n', g.filetype); 
	end;
   
   % assign format from filetype
   % ---------------------------
   if isempty(g.format) & ~isempty(g.filetype)
      indexformat = strmatch(lower(g.filetype), listtype, 'exact');
      g.format = listimportformat{indexformat};
      if isempty(g.skipline)
         g.skipline = listskipline(indexformat);
      end;
      if isempty(g.filetype) 
         error( ['Readlocs error: filetype can not be detected from' ...
               'file extension and custom format not specified']);
      end;
   end;
   
   % import file
   % -----------
   if strcmp(lower(g.filetype(1:end-1)), 'polhemus') | ...
         strcmp(lower(g.filetype), 'polhemus')
	      [eloctmp labels X Y Z]= readelp( filename );
 			if strcmp(lower(g.filetype), 'polhemusy')
             [theta radius X Y Z] = cart2topo( Y', X', Z','optim',1);  
         else
             [theta radius X Y Z] = cart2topo( X', Y', Z','optim',1);  
         end;
         for index = 1:length( eloctmp )
			  eloc(index).labels = labels{index};
			  eloc(index).theta  = theta(index);
			  eloc(index).radius = radius(index);
			  if strcmp(lower(g.filetype), 'polhemusy')
				  eloc(index).X = Y(index);
				  eloc(index).Y = X(index);	
				  eloc(index).Z = Z(index);	
			  else
				  eloc(index).X = X(index);
				  eloc(index).Y = Y(index);	
				  eloc(index).Z = Z(index);	
			  end;
         end;
	else      
      % importing file
      % --------------
      array = load_file_or_array( filename, g.skipline);
      if size(array,2) < length(g.format)
      	fprintf('Readlocs warning: # of columns in file inferior to # format entries');
   	elseif size(array,2) > length(g.format)
      	fprintf('Readlocs warning: # of columns in file superior to # format entries');
   	end;
      
      % removing comments and empty lines
      % ---------------------------------
      indexbeg = 1;
      while isempty(array{indexbeg,1}) | ...
            (isstr(array{indexbeg,1}) & array{indexbeg,1}(1) == '#' )
				indexbeg = indexbeg+1;
      end;
      array = array(indexbeg:end,:);
      
      % converting file
      % ---------------
   	for indexcol = 1:min(size(array,2), length(g.format))
      	for indexrow = 1:size( array, 1)
         	eval ( [ 'eloc(indexrow).' g.format{indexcol} '= array{indexrow, indexcol};' ]);   
      	end;
   	end;
   end;
   
   % handling BESA coordinates
   % -------------------------
   if isfield(eloc, 'sph_phi_besa') & isfield(eloc, 'sph_theta_besa')
		for index = 1:length(eloc)
			[ tmp eloc(index).theta eloc(index).radius] ...
				= sph2topo( [ 1 eloc(index).sph_phi_besa eloc(index).sph_theta_besa] , 1, 1);
      end;
      fprintf('Readlocs: BESA spherical coords. converted, now deleting BESA fields\n');   
      fprintf('          to avoid confusion (these field can though be exported)\n');   
		eloc = rmfield(eloc, 'sph_phi_besa');
		eloc = rmfield(eloc, 'sph_theta_besa');
   end;
   
   % converting XYZ coordinates to polar
   % -----------------------------------
   if isfield(eloc, 'sph_phi_besa') & isfield(eloc, 'sph_theta_besa')
      	[theta radius X Y Z] = cart2topo( Y', X', Z','optim',1);  
         for index = 1:length( eloc )
			  eloc(index).labels = labels{index};
			  eloc(index).theta  = theta(index);
			  eloc(index).radius = radius(index);
			  eloc(index).X = X(index);
			  eloc(index).Y = Y(index);	
			  eloc(index).Z = Z(index);	
			end;
   end;
   
   % inserting labels if no labels
   % -----------------------------
   if ~isfield(eloc, 'labels')
      fprintf('Readlocs: Automatically inserting electrode labels\n');
      for index = 1:length(eloc)
         eloc(index).labels = [ 'E' int2str(index) ];
      end;
   end;
   
   % resorting electrodes if number not-sorted
   % -----------------------------------------
   if isfield(eloc, 'channum')
      allchannum = cell2mat( { eloc.channum } );
      if any( sort(allchannum) ~= allchannum )
         fprintf('Readlocs: Resorting channel number based on ''channum'' column indices\n');
         [tmp newindices] = sort(allchannum);
         eloc = eloc(newindices);
      end;
		eloc = rmfield(eloc, 'channum');      
   end;
else
    if isstruct(filename)
        eloc = filename;
    else
        disp('Readlocs: input variable must be a string or a structure');
    end;        
end;
if ~isempty(g.elecind)
	eloc = eloc(g.elecind);
end;
theta = cell2mat({ eloc.theta });
radius  = cell2mat({ eloc.radius });
if isnumeric(eloc(1).labels)
    for index = 1:length(eloc)
        eloc(index).labels = int2str(eloc(index).labels);
    end;
end;
labels = { eloc.labels };

return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skipline );
	 if isempty(skipline),
       skipline = 0;
    end;
    if exist( varname ) == 2
        array = loadtxt(varname,'verbose','off','skipline',skipline);
    else % variable in the global workspace
         % --------------------------
         try, array = evalin('base', varname);
	     catch, error('readlocs: cannot find file or variable, check syntax');
		 end;
    end;     
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
	if strcmpi(str, '-X'), str = upper(str); mult = -1; return; end;
	if strcmpi(str, '-Y'), str = upper(str); mult = -1; return; end;
	if strcmpi(str, '-Z'), str = upper(str); mult = -1; return; end;
	if strcmpi(str, 'custum1'), return; end;
	if strcmpi(str, 'custum2'), return; end;
	if strcmpi(str, 'custum3'), return; end;
	if strcmpi(str, 'custum4'), return; end;
   error(['Readlocs: undefined field ''' str '''']);
   
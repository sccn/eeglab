% readlocs() - read polar electrode positions (expanded from the ICA 
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
%   'filetype'  - ['loc'|'sph'|'xyz'|'polhemus'|'besa'|'chanedit'] type of the 
%                 file to be read. 
%                  'loc' is the eeglab() format (see below)
%                  'sph' is a matlab spherical coordinate file (spherical
%                        coordinates in Matlab are different from Besa 
%                        spherical coordinates).
%                  'xyz' contain carhtesian coordinates of the electrodes
%                  'polhemus' polhemus file (see readelp() )
%                  'besa' besa '.elp' coordinate file (contact author for
%                         reading other Besa files). '.elp' Besa files
%                         have only 3 columns of 'theta' angle, 'phi' angle
%                         and 'weight' (which is ignored). Also note that 
%                         'theta' and 'phi' meaning is also the inverse
%                         as the one they have in Matlab (in Besa 'theta' =
%                         azimut and 'phi' = horiz).
%                  'chanedit' files created by pop_chanedit()
%   'skipline'  - [integer] number of lines to skip 
%   'elecind'   - [integer array] indices of electrode to read 
%   'maindir'   - ['X'|'Y'] Direction pointing toward the subject in 
%                the Polhemus .elp file. Default is 'X'.  Used to 
%                convert locations from cartesian to polar coordinates.
%
% File formats:
%   The extension of the file determines its type
%   if filetype is not specified
%   '.loc' or '.locs'   - polar format. Example:
%               1    -18    .352       Fp1
%               2     18    .352       Fp2
%               3    -90    .181       C3
%               4     90    .181       C4
%                 more lines ...
%               Note that in previous releases, the channel labels
%               had to contain 4 characters (spaces replaced by '.').
%               This format still works but the dots are no longer
%               required.
%   '.sph' - spherical coordinate file. Example:
%               1    -63.36    -72      Fp1
%               2     63.36    72       Fp2
%               3     32.58    0        C3
%               4     32.58    0        C4
%                 more lines ...
%   '.xyz' - cartesian coordinate file. Example:
%               1   -0.8355   -0.2192   -0.5039      Fp1
%               2   -0.8355    0.2192    0.5039      Fp2
%               3    0.3956         0   -0.9184      C3
%               4    0.3956         0    0.9184      C4
%                 more lines ...
%   '.txt' - read ascii files saved using pop_chanedit()
%   '.elp' - Polhemus coordinate file (uses readelp())
%
% Outputs:
%   eloc      - structure containing the channel names and locations.
%               It has three fields 'labels', 'theta' and 'radius'.
%   labels    - cell array of strings giving the  names of the electrodes
%   theta     - vector of polar angles of the electrodes in degrees
%   radius    - vector of polar norms of the electrodes
%
% Note: the function cart2topo() is used to convert cartesian to polar
%       coordinates. By default the current function uses cart2topo()
%       options to recompute the best center coordinates.
%
% Author: Arnaud Delorme & Scott Makeig CNL / Salk Institute, 28 Feb 2002
%
% See also: readelp(), topo2sph(), sph2topo(), cart2topo(), sph2cart()

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
% Revision 1.16  2002/11/12 23:39:57  arno
% convert numerical labels to string
%
% Revision 1.15  2002/10/23 15:36:43  arno
% header comment
%
% Revision 1.14  2002/08/17 18:20:22  scott
% help msg update
%
% Revision 1.13  2002/07/02 00:30:23  arno
% updating header
%
% Revision 1.12  2002/07/02 00:26:28  arno
% updating for reading besa .elp files
%
% Revision 1.11  2002/05/19 13:40:07  scott
% turned off loadtxt 'verbose' -sm
%
% Revision 1.10  2002/05/18 19:03:58  scott
% typo -sm
%
% Revision 1.9  2002/05/02 01:34:56  arno
% getting XYZ back
%
% Revision 1.8  2002/05/02 00:33:43  arno
% remove minus
%
% Revision 1.7  2002/05/02 00:22:57  arno
% same
%
% Revision 1.6  2002/05/02 00:06:19  arno
% correcting polhemus reading
%
% Revision 1.5  2002/05/01 03:27:40  arno
% removing spherical
%
% Revision 1.4  2002/05/01 02:11:56  arno
% new .txt format
%
% Revision 1.3  2002/05/01 01:15:31  arno
% removing topo optimization
%
% Revision 1.2  2002/04/18 15:17:25  scott
% editted help msg -sm
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%

function [eloc, labels, theta, radius] = readlocs( filename, varargin ); 

if nargin < 1
	help readlocs;
	return;
end;

g = [];
if ~isempty(varargin)
	for index=1:length(varargin)
		if iscell(varargin{index})
			varargin{index} = {varargin{index}};
		end;
	end;
	try
		g = struct(varargin{:});
	catch
		error('readlocs: error in ''key'', ''val'' sequence');
	end;
end;

try, g.filetype = g.filetype;    catch, g.filetype = []; end;
try, g.skipline;                 catch, g.skipline = 0; end;
try, g.maindir;                  catch, g.maindir = 'X'; end;
try, g.elecind;                  catch, g.elecind = []; end;
g.filetype = lower( g.filetype );

allfields = fieldnames(g);
for index=1:length(allfields)
	switch allfields{index}
	 case { 'filetype' 'skipline' 'maindir' 'elecind' };
	 otherwise, error(['readlocs: unknow option ''' allfields{index} '''']);
	end;
end;

if isstr(filename)
	% open file
	% ---------
	array = load_file_or_array( filename, g.skipline);

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
		 otherwise, g.filetype = 'loc';
		end;
	end;
	
	% scan file
	% ---------
    switch g.filetype,
        case 'xyz', 
			for index = 1:size( array, 1)
			  eloc(index).X = array{index, 2};
			  eloc(index).Y = array{index, 3};
			  eloc(index).Z = array{index, 4};
			  eloc(index).labels  = array{index, 5};
			  [ eloc(index).theta eloc(index).radius] ...
			         = cart2topo( eloc(index).X, eloc(index).Y, eloc(index).Z);
			end;  
        case 'sph', 
			for index = 1:size( array, 1)
			  eloc(index).sph_theta    = array{index, 2};
			  eloc(index).sph_phi = array{index, 3};
			  eloc(index).sph_radius = 1;
			  eloc(index).labels  = array{index, 4};
			end;  
        case 'loc', 
			for index = 1:size( array, 1)
			  eloc(index).theta = array{index, 2};
			  eloc(index).radius  = array{index, 3};
			  eloc(index).labels  = array{index, 4};
			  eloc(index).labels( find(eloc(index).labels == '.' )) = ' ';
			end;
        case 'polhemus' 
            [eloctmp labels X Y Z]= readelp( filename );
            if exist('elpmaindir') ~= 1, elpmaindir = 'X'; end;
 			if strcmp(lower(elpmaindir), 'x')
                [theta radius X Y Z] = cart2topo( X', Y', Z','optim',1);  
            else
                [theta radius X Y Z] = cart2topo( Y', X', Z','optim',1);  
            end;
			for index = 1:length( eloctmp )
			  eloc(index).labels = labels{index};
			  eloc(index).theta  = theta(index);
			  eloc(index).radius = radius(index);
			  if strcmp(lower(elpmaindir), 'x')
				  eloc(index).X = X(index);
				  eloc(index).Y = Y(index);	
				  eloc(index).Z = Z(index);	
			  else
				  eloc(index).X = Y(index);
				  eloc(index).Y = X(index);	
				  eloc(index).Z = Z(index);	
			  end;
            end;
     	 case 'chanedit', 
		    if isempty(array(end,1)), totlines = size( array, 1)-1; else totlines = size( array, 1); end;
			for index = 2:totlines
				if ~isempty(array{index,2}) eloc(index-1).labels  = array{index, 2}; end;
				if ~isempty(array{index,3}) eloc(index-1).theta = array{index, 3}; end;
				if ~isempty(array{index,4}) eloc(index-1).radius  = array{index, 4}; end;
				if ~isempty(array{index,5}) eloc(index-1).X = array{index, 5}; end;
				if ~isempty(array{index,6}) eloc(index-1).Y = array{index, 6}; end;
				if ~isempty(array{index,7}) eloc(index-1).Z = array{index, 7}; end;
				if ~isempty(array{index,8}) eloc(index-1).sph_theta = array{index, 8}; end;
				if ~isempty(array{index,9}) eloc(index-1).sph_phi   = array{index, 9}; end;
				if ~isempty(array{index,10}) eloc(index-1).sph_radius   = array{index, 10}; end;
			end;
	     case 'besa',
	        indexbeg = 1;
		    while isempty(array{indexbeg,1}) | isstr(array{indexbeg,1})
				indexbeg = indexbeg+1;
			end;
			for index = indexbeg:size( array, 1)
				shifted_i = index-indexbeg+1;
				eloc(shifted_i).labels  = num2str(shifted_i);
				eloc(shifted_i).besathloc  = array{index, 1};
				eloc(shifted_i).besaphloc = array{index, 2};
				[ tmp eloc(shifted_i).theta eloc(shifted_i).radius] ...
					= sph2topo( [ 1 eloc(shifted_i).besathloc eloc(shifted_i).besaphloc] , 1, 1);
			end;
			eloc = rmfield(eloc, 'besathloc');
			eloc = rmfield(eloc, 'besaphloc');
        otherwise, error('Readlocs(): unrecognized file extension');
    end;
    for index = 1:length( eloc )
        if ~isstr(eloc(index).labels)
            eloc(index).labels = num2str( eloc(index).labels );
            if ~isempty(findstr( '0.', eloc(index).labels ))
                eloc(index).labels = eloc(index).labels(3:end);
            end;    
        else
            eloc(index).labels = deblank(num2str( eloc(index).labels ));
        end;
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

    if exist( varname ) == 2
        array = loadtxt(varname,'verbose','off');
    else % variable in the global workspace
         % --------------------------
         try, array = evalin('base', varname);
	     catch, error('readlocs: cannot find file or variable, check syntax');
		 end;
    end;     
return;

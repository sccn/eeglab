% readlocs() - read polar electrode positions from ICA toolbox
%             
% Usage:
%   >> [eloc, labels, theta, radius] = readlocs( filename, elpmaindir );
%
% Inputs:
%   filename   - name of the file containing the electrode locations
%                and convert in polar coordinates
%   elpmaindir - direction pointing toward the subject in the Polhemus 
%                elp file. Can be 'X' or 'Y'. Default is 'X'. This 
%                information is used to convert carthesian to polar
%                coordinates.
%
% Note on suported formats:
% The extension of the file determine its type
%   '.loc' or '.locs'   - polar format, for example
%               1    -18    .352       Fp1
%               2     18    .352       Fp2
%               3    -90    .181       C3
%               4     90    .181       C4
%               more lines ...
%   '.sph' - spherical coordinate file, for example
%               1    -63.36    -72      Fp1
%               2     63.36    72       Fp2
%               3     32.58    0        C3
%               4     32.58    0        C4
%               more lines ...
%   '.xyz' - carthesian coordinate file, for example
%               1   -0.8355   -0.2192   -0.5039      Fp1
%               2   -0.8355    0.2192    0.5039      Fp2
%               3    0.3956         0   -0.9184      C3
%               4    0.3956         0    0.9184      C4
%               more lines ...
%   '.elp' - Polhemus coordinate file (use the readelp() function)
%
% Outputs:
%   eloc      - structure containing the channel names and locations.
%               It has three fields 'labels', 'theta' and 'radius'.
%   labels    - names of the electrodes
%   theta     - vector of polar angles of the electrodes in degree
%   radius    - vector of polar norms of the electrodes
%
% Note: the function cart2topo() is used to convert carthesian to polar
%       coordinates. By default the current function uses cart2topo()
%       options to recompute the best center and optimize the squeezing
%       parameter.
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

function [eloc, labels, theta, radius] = readlocs( filename, elpmaindir ); 

if nargin < 1
	help readlocs;
	return;
end;

if isstr(filename)
	% open file
	% ---------
	array = load_file_or_array( filename, 0);

    periods = find(filename == '.');
    fileextension = filename(periods(end)+1:end);

	% scan file
	% ---------
    switch lower(fileextension),
        case { '' 'chanlocs'}, eloc = array; 
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
			  eloc(index).sperical_az    = array{index, 2};
			  eloc(index).sperical_horiz = array{index, 3};
			  eloc(index).labels  = array{index, 4};
			  [ tmp eloc(index).theta eloc(index).radius] ...
			         = sph2topo( index, eloc(index).sperical_az, eloc(index).sperical_horiz);
			  [eloc.X eloc.Y eloc.Z] = sph2cart(eloc(index).sperical_horiz'/180*pi, eloc(index).sperical_az'/180*pi, 1);
			end;  
        case { 'loc' 'locs' }, 
			for index = 1:size( array, 1)
			  eloc(index).theta = array{index, 2};
			  eloc(index).radius  = array{index, 3};
			  eloc(index).labels  = array{index, 4};
			  eloc(index).labels( find(eloc(index).labels == '.' )) = ' ';
              [eloc(index).sperical_az eloc(index).sperical_horiz] = topo2sph( [eloc(index).theta eloc(index).radius]);
			  [eloc(index).X eloc(index).Y eloc(index).Z] = sph2cart(eloc(index).sperical_horiz'/180*pi, eloc(index).sperical_az'/180*pi, 1);
			end;
        case 'elp', 
            [eloc labels X Y Z]= readelp( filename );
            if exist('elpmaindir') ~= 1, elpmaindir = 'X'; end;
 			if strcmp(lower(elpmaindir), 'x')
                [theta radius] = cart2topo( -X', -Y', Z',[],-1);  
            else
                [theta radius] = cart2topo( -Y', -X', Z',[],-1);  
            end;                   
			for index = 1:length( eloc )
			  eloc(index).theta  = theta(index);
			  eloc(index).radius = radius(index);
			  eloc(index).labels = labels{index};
            end;
        otherwise, error('Readlocs: unrecognised file extension');
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
theta = cell2mat({ eloc.theta });
radius  = cell2mat({ eloc.radius });
labels = { eloc.labels };

return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skipline );

    if exist( varname ) == 2
        array = loadtxt(varname);
    else % variable in the global workspace
         % --------------------------
         array = evalin('base', varname);
    end;     
return;

% readlocs() - read electrode location coordinates and other information from a file. 
%              Several standard file formats are supported. Users may also specify 
%              a custom column format. Defined format examples are given below 
%              (see File Formats).
% Usage:
%   >>  eloc = readlocs( filename );
%   >>  EEG.chanlocs = readlocs( filename, 'key', 'val', ... ); 
%   >>  [eloc, labels, theta, radius, indices] = ...
%                                               readlocs( filename, 'key', 'val', ... );
% Inputs:
%   filename   - Name of the file containing the electrode locations
%                {default: 2-D polar coordinates} (see >> help topoplot )
%
% Optional inputs:
%   'filetype'  - ['loc'|'sph'|'sfp'|'xyz'|'asc'|'polhemus'|'besa'|'chanedit'|'custom'] 
%                 Type of the file to read. By default the file type is determined 
%                 using the file extension (see below under File Formats),
%                  'loc'   an EEGLAB 2-D polar coordinates channel locations file 
%                          Coordinates are theta and radius (see definitions below).
%                  'sph'   Matlab spherical coordinates (Note: spherical
%                          coordinates used by Matlab functions are different 
%                          from spherical coordinates used by BESA - see below).
%                  'sfp'   EGI Cartesian coordinates (NOT Matlab Cartesian - see below).
%                  'xyz'   Matlab/EEGLAB Cartesian coordinates (NOT EGI Cartesian).
%                          z is toward nose; y is toward left ear; z is toward vertex
%                  'asc'   Neuroscan polar coordinates.
%                  'polhemus' or 'polhemusx' - Polhemus electrode location file recorded 
%                          with 'X' on sensor pointing to subject (see below and readelp()).
%                  'polhemusy' - Polhemus electrode location file recorded with 
%                          'Y' on sensor pointing to subject (see below and readelp()).
%                  'besa' BESA-'.elp' spherical coordinates. (Not MATLAB spherical -
%                           see below).
%                  'chanedit' - EEGLAB channel location file created by pop_chanedit().
%                  'custom' - Ascii file with columns in user-defined 'format' (see below).
%   'importmode' - ['eeglab'|'native'] for location files containing 3-D cartesian electrode
%                  coordinates, import either in EEGLAB format (nose pointing toward +X). 
%                  This may not always be possible since EEGLAB might not be able to 
%                  determine the nose direction for scanned electrode files. 'native' import
%                  original carthesian coordinates (user can then specify the position of
%                  the nose when calling the topoplot() function; in EEGLAB the position
%                  of the nose is stored in the EEG.chaninfo structure). {default 'eeglab'}
%   'format'    -  [cell array] Format of a 'custom' channel location file (see above).
%                  {default: if no file type is defined. The cell array contains
%                  labels defining the meaning of each column of the input file.
%                           'channum'   [positive integer] channel number.
%                           'labels'    [string] channel name (no spaces).
%                           'theta'     [real degrees] 2-D angle in polar coordinates.
%                                       positive => rotating from nose (0) toward left ear
%                           'radius'    [real] radius for 2-D polar coords; 0.5 is the head
%                                       disk radius and limit for topoplot() plotting).
%                           'X'         [real] Matlab-Cartesian X coordinate (to nose).
%                           'Y'         [real] Matlab-Cartesian Y coordinate (to left ear).
%                           'Z'         [real] Matlab-Cartesian Z coordinate (to vertex).
%                           '-X','-Y','-Z' Matlab-Cartesian coordinates pointing opposite
%                                       to the above.
%                           'sph_theta' [real degrees] Matlab spherical horizontal angle.
%                                       positive => rotating from nose (0) toward left ear.
%                           'sph_phi'   [real degrees] Matlab spherical elevation angle.
%                                       positive => rotating from horizontal (0) upwards.
%                           'sph_radius' [real] distance from head center (unused).
%                           'sph_phi_besa' [real degrees] BESA phi angle from vertical.
%                                       positive => rotating from vertex (0) towards right ear.
%                           'sph_theta_besa' [real degrees] BESA theta horiz/azimuthal angle.
%                                       positive => rotating from right ear (0) toward nose.
%                           'ignore'    ignore column}.
%     The input file may also contain other channel information fields.
%                           'type'      channel type: 'EEG', 'MEG', 'EMG', 'ECG', others ...
%                           'calib'     [real near 1.0] channel calibration value.
%                           'gain'      [real > 1] channel gain.
%                           'custom1'   custom field #1.
%                           'custom2', 'custom3', 'custom4', etc.    more custom fields
%   'skiplines' - [integer] Number of header lines to skip (in 'custom' file types only).
%                 Note: Characters on a line following '%' will be treated as comments.
%   'readchans' - [integer array] indices of electrodes to read. {default: all}
%   'center'    - [(1,3) real array or 'auto'] center of xyz coordinates for conversion 
%                 to spherical or polar, Specify the center of the sphere here, or 'auto'. 
%                 This uses the center of the sphere that best fits all the electrode 
%                 locations read. {default: [0 0 0]}
% Outputs:
%   eloc        - structure containing the channel names and locations (if present).
%                 It has three fields: 'eloc.labels', 'eloc.theta' and 'eloc.radius' 
%                 identical in meaning to the EEGLAB struct 'EEG.chanlocs'.
%   labels      - cell array of strings giving the names of the electrodes. NOTE: Unlike the
%                 three outputs below, includes labels of channels *without* location info.
%   theta       - vector (in degrees) of polar angles of the electrode locations.
%   radius      - vector of polar-coordinate radii (arc_lengths) of the electrode locations 
%   indices     - indices, k, of channels with non-empty 'locs(k).theta' coordinate
%
% File formats:
%   If 'filetype' is unspecified, the file extension determines its type.
%
%   '.loc' or '.locs' or '.eloc': 
%               polar coordinates. Notes: angles in degrees: 
%               right ear is 90; left ear -90; head disk radius is 0.5. 
%               Fields:   N    angle  radius    label
%               Sample:   1    -18    .511       Fp1   
%                         2     18    .511       Fp2  
%                         3    -90    .256       C3
%                         4     90    .256       C4
%                           ...
%               Note: In previous releases, channel labels had to contain exactly 
%               four characters (spaces replaced by '.'). This format still works, 
%               though dots are no longer required.
%   '.sph':
%               Matlab spherical coordinates. Notes: theta is the azimuthal/horizontal angle
%               in deg.: 0 is toward nose, 90 rotated to left ear. Following this, performs
%               the elevation (phi). Angles in degrees.
%               Fields:   N    theta    phi    label
%               Sample:   1      18     -2      Fp1
%                         2     -18     -2      Fp2
%                         3      90     44      C3
%                         4     -90     44      C4
%                           ...
%   '.elc':
%               Cartesian 3-D electrode coordinates scanned using the EETrak software. 
%               See readeetraklocs().
%   '.elp':     
%               Polhemus-.'elp' Cartesian coordinates. By default, an .elp extension is read
%               as PolhemusX-elp in which 'X' on the Polhemus sensor is pointed toward the 
%               subject. Polhemus files are not in columnar format (see readelp()).
%   '.elp':
%               BESA-'.elp' spherical coordinates: Need to specify 'filetype','besa'.
%               The elevation angle (phi) is measured from the vertical axis. Positive 
%               rotation is toward right ear. Next, perform azimuthal/horizontal rotation 
%               (theta): 0 is toward right ear; 90 is toward nose, -90 toward occiput. 
%               Angles are in degrees.  If labels are absent or weights are given in 
%               a last column, readlocs() adjusts for this. Default labels are E1, E2, ...
%               Fields:   label      phi  theta   
%               Sample:   Fp1        -92   -72    
%                         Fp2         92    72   
%                         C3         -46    0  
%                         C4          46    0 
%                           ...
%   '.xyz': 
%               Matlab/EEGLAB Cartesian coordinates. Here. x is towards the nose, 
%               y is towards the left ear, and z towards the vertex. Note that the first
%               column (x) is -Y in a Matlab 3-D plot, the second column (y) is X in a 
%               matlab 3-D plot, and the third column (z) is Z.
%               Fields:   channum   x           y         z     label
%               Sample:   1       .950        .308     -.035     Fp1
%                         2       .950       -.308     -.035     Fp2
%                         3        0           .719      .695    C3
%                         4        0          -.719      .695    C4
%                           ...
%   '.asc', '.dat':     
%               Neuroscan-.'asc' or '.dat' Cartesian polar coordinates text file.
%   '.sfp': 
%               BESA/EGI-xyz Cartesian coordinates. Notes: For EGI, x is toward right ear, 
%               y is toward the nose, z is toward the vertex. EEGLAB converts EGI 
%               Cartesian coordinates to Matlab/EEGLAB xyz coordinates. 
%               Fields:   label   x           y          z
%               Sample:   Fp1    -.308        .950      -.035    
%                         Fp2     .308        .950      -.035  
%                         C3     -.719        0          .695  
%                         C4      .719        0          .695  
%                           ...
%   '.ced':   
%               ASCII file saved by pop_chanedit(). Contains multiple MATLAB/EEGLAB formats.
%               Cartesian coordinates are as in the 'xyz' format (above).
%               Fields:   channum  label  theta  radius   x      y      z    sph_theta   sph_phi  ...
%               Sample:   1        Fp1     -18    .511   .950   .308  -.035   18         -2       ...
%                         2        Fp2      18    .511   .950  -.308  -.035  -18         -2       ...
%                         3        C3      -90    .256   0      .719   .695   90         44       ...
%                         4        C4       90    .256   0     -.719   .695  -90         44       ...
%                           ...
%               The last columns of the file may contain any other defined fields (gain,
%               calib, type, custom).
%
%    Fieldtrip structure: 
%               If a Fieltrip structure is given as input, an EEGLAB
%               chanlocs structure is returned
%
% Author: Arnaud Delorme, Salk Institute, 8 Dec 2002
%
% See also: readelp(), writelocs(), topo2sph(), sph2topo(), sph2cart()

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


function [eloc, labels, theta, radius, indices] = readlocs( filename, varargin ); 

if nargin < 1
	help readlocs;
	return;
end;

% NOTE: To add a new channel format:
% ----------------------------------
% 1) Add a new element to the structure 'chanformat' (see 'ADD NEW FORMATS HERE' below):
% 2)  Enter a format 'type' for the new file format, 
% 3)  Enter a (short) 'typestring' description of the format
% 4)  Enter a longer format 'description' (possibly multiline, see ex. (1) below)
% 5)  Enter format file column labels in the 'importformat' field (see ex. (2) below)
% 6)  Enter the number of header lines to skip (if any) in the 'skipline' field
% 7)  Document the new channel format in the help message above.
% 8)  After testing, please send the new version of readloca.m to us
%       at eeglab@sccn.ucsd.edu with a sample locs file.
% The 'chanformat' structure is also used (automatically) by the writelocs() 
% and pop_readlocs() functions. You do not need to edit these functions.

chanformat(1).type         = 'polhemus';
chanformat(1).typestring   = 'Polhemus native .elp file';
chanformat(1).description  = [ 'Polhemus native coordinate file containing scanned electrode positions. ' ...
                               'User must select the direction ' ...
                               'for the nose after importing the data file.' ];
chanformat(1).importformat = 'readelp() function';
% ---------------------------------------------------------------------------------------------------
chanformat(2).type         = 'besa';
chanformat(2).typestring   = 'BESA spherical .elp file';
chanformat(2).description  = [ 'BESA spherical coordinate file. Note that BESA spherical coordinates ' ...
                               'are different from Matlab spherical coordinates' ];
chanformat(2).skipline     = 0; % some BESA files do not have headers
chanformat(2).importformat = { 'type' 'labels' 'sph_theta_besa' 'sph_phi_besa' 'sph_radius' };
% ---------------------------------------------------------------------------------------------------
chanformat(3).type         = 'xyz';
chanformat(3).typestring   = 'Matlab .xyz file';
chanformat(3).description  = [ 'Standard 3-D cartesian coordinate files with electrode labels in ' ...
                               'the first column and X, Y, and Z coordinates in columns 2, 3, and 4' ];
chanformat(3).importformat = { 'channum' '-Y' 'X' 'Z' 'labels'};
% ---------------------------------------------------------------------------------------------------
chanformat(4).type         = 'sfp';
chanformat(4).typestring   = 'BESA or EGI 3-D cartesian .sfp file';
chanformat(4).description  = [ 'Standard BESA 3-D cartesian coordinate files with electrode labels in ' ...
                               'the first column and X, Y, and Z coordinates in columns 2, 3, and 4.' ...
                               'Coordinates are re-oriented to fit the EEGLAB standard of having the ' ...
                               'nose along the +X axis.' ];
chanformat(4).importformat = { 'labels' '-Y' 'X' 'Z' };
chanformat(4).skipline     = 0;
% ---------------------------------------------------------------------------------------------------
chanformat(5).type         = 'loc';
chanformat(5).typestring   = 'EEGLAB polar .loc file';
chanformat(5).description  = [ 'EEGLAB polar .loc file' ];
chanformat(5).importformat = { 'channum' 'theta' 'radius' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(6).type         = 'sph';
chanformat(6).typestring   = 'Matlab .sph spherical file';
chanformat(6).description  = [ 'Standard 3-D spherical coordinate files in Matlab format' ];
chanformat(6).importformat = { 'channum' 'sph_theta' 'sph_phi' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(7).type         = 'asc';
chanformat(7).typestring   = 'Neuroscan polar .asc file';
chanformat(7).description  = [ 'Neuroscan polar .asc file, automatically recentered to fit EEGLAB standard' ...
                               'of having ''Cz'' at (0,0).' ];
chanformat(7).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(8).type         = 'dat';
chanformat(8).typestring   = 'Neuroscan 3-D .dat file';
chanformat(8).description  = [ 'Neuroscan 3-D cartesian .dat file. Coordinates are re-oriented to fit ' ...
                               'the EEGLAB standard of having the nose along the +X axis.' ];
chanformat(8).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(9).type         = 'elc';
chanformat(9).typestring   = 'ASA .elc 3-D file';
chanformat(9).description  = [ 'ASA .elc 3-D coordinate file containing scanned electrode positions. ' ...
                               'User must select the direction ' ...
                               'for the nose after importing the data file.' ];
chanformat(9).importformat = 'readeetraklocs';
% ---------------------------------------------------------------------------------------------------
chanformat(10).type         = 'chanedit';
chanformat(10).typestring   = 'EEGLAB complete 3-D file';
chanformat(10).description  = [ 'EEGLAB file containing polar, cartesian 3-D, and spherical 3-D ' ...
                               'electrode locations.' ];
chanformat(10).importformat = { 'channum' 'labels'  'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' ...
                               'sph_radius' 'type' };
chanformat(10).skipline     = 1;
% ---------------------------------------------------------------------------------------------------
chanformat(11).type         = 'custom';
chanformat(11).typestring   = 'Custom file format';
chanformat(11).description  = 'Custom ASCII file format where user can define content for each file columns.';
chanformat(11).importformat = '';
% ---------------------------------------------------------------------------------------------------
% ----- ADD MORE FORMATS HERE -----------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------

listcolformat = { 'labels' 'channum' 'theta' 'radius' 'sph_theta' 'sph_phi' ...
      'sph_radius' 'sph_theta_besa' 'sph_phi_besa' 'gain' 'calib' 'type' ...
      'X' 'Y' 'Z' '-X' '-Y' '-Z' 'custom1' 'custom2' 'custom3' 'custom4' 'ignore' 'not def' };

% ----------------------------------
% special mode for getting the info
% ----------------------------------
if isstr(filename) & strcmp(filename, 'getinfos')
   eloc = chanformat;
   labels = listcolformat;
   return;
end;

g = finputcheck( varargin, ...
   { 'filetype'	   'string'  {}                 '';
     'importmode'  'string'  { 'eeglab','native' } 'eeglab';
     'defaultelp'  'string'  { 'besa','polhemus' } 'polhemus';
     'skiplines'   'integer' [0 Inf] 			[];
     'elecind'     'integer' [1 Inf]	    	[];
     'format'	   'cell'	 []					{} }, 'readlocs');
if isstr(g), error(g); end;  
if ~isempty(g.format), g.filetype = 'custom'; end;

if isstr(filename)
   
   % format auto detection
	% --------------------
   if strcmpi(g.filetype, 'autodetect'), g.filetype = ''; end;
   g.filetype = strtok(g.filetype);
   periods = find(filename == '.');
   fileextension = filename(periods(end)+1:end);
   g.filetype = lower(g.filetype);
   if isempty(g.filetype)
       switch lower(fileextension),
        case {'loc' 'locs' }, g.filetype = 'loc';
        case 'xyz', g.filetype = 'xyz'; 
          fprintf( [ 'WARNING: Matlab Cartesian coord. file extension (".xyz") detected.\n' ... 
                  'If importing EGI Cartesian coords, force type "sfp" instead.\n'] );
        case 'sph', g.filetype = 'sph';
        case 'ced', g.filetype = 'chanedit';
        case 'elp', g.filetype = g.defaultelp;
        case 'asc', g.filetype = 'asc';
        case 'dat', g.filetype = 'dat';
        case 'elc', g.filetype = 'elc';
        case 'eps', g.filetype = 'besa';
        case 'sfp', g.filetype = 'sfp';
        otherwise, g.filetype =  ''; 
       end;
       fprintf('readlocs(): ''%s'' format assumed from file extension\n', g.filetype); 
   else 
       if strcmpi(g.filetype, 'locs'),  g.filetype = 'loc'; end
       if strcmpi(g.filetype, 'eloc'),  g.filetype = 'loc'; end
   end;
   
   % assign format from filetype
   % ---------------------------
   if ~isempty(g.filetype) & ~strcmpi(g.filetype, 'custom') ...
           & ~strcmpi(g.filetype, 'asc') & ~strcmpi(g.filetype, 'elc') & ~strcmpi(g.filetype, 'dat')
      indexformat = strmatch(lower(g.filetype), { chanformat.type }, 'exact');
      g.format = chanformat(indexformat).importformat;
      if isempty(g.skiplines)
         g.skiplines = chanformat(indexformat).skipline;
      end;
      if isempty(g.filetype) 
         error( ['readlocs() error: The filetype cannot be detected from the \n' ...
                 '                  file extension, and custom format not specified']);
      end;
   end;
   
   % import file
   % -----------
   if strcmp(g.filetype, 'asc') | strcmp(g.filetype, 'dat')
       eloc = readneurolocs( filename );
       eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
       eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
       if isfield(eloc, 'type')
           for index = 1:length(eloc)
               type = eloc(index).type;
               if type == 69,     eloc(index).type = 'EEG';
               elseif type == 88, eloc(index).type = 'REF';
               elseif type >= 76 & type <= 82, eloc(index).type = 'FID';
               else eloc(index).type = num2str(eloc(index).type);
               end;
           end;
       end;
   elseif strcmp(g.filetype, 'elc')
       eloc = readeetraklocs( filename );
       %eloc = read_asa_elc( filename ); % from fieldtrip
       %eloc = struct('labels', eloc.label, 'X', mattocell(eloc.pnt(:,1)'), 'Y', ...
       %                        mattocell(eloc.pnt(:,2)'), 'Z', mattocell(eloc.pnt(:,3)'));
       eloc = convertlocs(eloc, 'cart2all');
       eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
       eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
   elseif strcmp(lower(g.filetype(1:end-1)), 'polhemus') | ...
           strcmp(g.filetype, 'polhemus')
       try, 
           [eloc labels X Y Z]= readelp( filename );
           if strcmp(g.filetype, 'polhemusy')
               tmp = X; X = Y; Y = tmp;
           end;
           for index = 1:length( eloc )
               eloc(index).X = X(index);
               eloc(index).Y = Y(index);	
               eloc(index).Z = Z(index);	
           end;
       catch, 
           disp('readlocs(): Could not read Polhemus coords. Trying to read BESA .elp file.');
           [eloc, labels, theta, radius, indices] = readlocs( filename, 'defaultelp', 'besa', varargin{:} );
       end;
   else      
       % importing file
       % --------------
       if isempty(g.skiplines), g.skiplines = 0; end;
       if strcmpi(g.filetype, 'chanedit')
           array = loadtxt( filename, 'delim', 9, 'skipline', g.skiplines, 'blankcell', 'off');
       else
           array = load_file_or_array( filename, g.skiplines);
       end;
       if size(array,2) < length(g.format)
           fprintf(['readlocs() warning: Fewer columns in the input than expected.\n' ...
                    '                    See >> help readlocs\n']);
       elseif size(array,2) > length(g.format)
           fprintf(['readlocs() warning: More columns in the input than expected.\n' ...
                    '                    See >> help readlocs\n']);
       end;
       
       % removing lines BESA
       % -------------------
       if isempty(array{1,2})
           disp('BESA header detected, skipping three lines...');
           array = load_file_or_array( filename, g.skiplines-1);
           if isempty(array{1,2})
               array = load_file_or_array( filename, g.skiplines-1);
           end;
       end;

       % xyz format, is the first col absent
       % -----------------------------------
       if strcmp(g.filetype, 'xyz')
           if size(array, 2) == 4
               array(:, 2:5) = array(:, 1:4);
           end;
       end;
       
       % removing comments and empty lines
       % ---------------------------------
       indexbeg = 1;
       while isempty(array{indexbeg,1}) | ...
               (isstr(array{indexbeg,1}) & array{indexbeg,1}(1) == '%' )
           indexbeg = indexbeg+1;
       end;
       array = array(indexbeg:end,:);
       
       % converting file
       % ---------------
       for indexcol = 1:min(size(array,2), length(g.format))
           [str mult] = checkformat(g.format{indexcol});
           for indexrow = 1:size( array, 1)
               if mult ~= 1
                   eval ( [ 'eloc(indexrow).'  str '= -array{indexrow, indexcol};' ]);
               else
                   eval ( [ 'eloc(indexrow).'  str '= array{indexrow, indexcol};' ]);
               end;
           end;
       end;
   end;
   
   % handling BESA coordinates
   % -------------------------
   if isfield(eloc, 'sph_theta_besa')
       if isfield(eloc, 'type')
           if isnumeric(eloc(1).type)
               disp('BESA format detected ( Theta | Phi )');
               for index = 1:length(eloc)
                   eloc(index).sph_phi_besa   = eloc(index).labels;
                   eloc(index).sph_theta_besa = eloc(index).type;
                   eloc(index).labels         = '';
                   eloc(index).type           = '';
               end;
               eloc = rmfield(eloc, 'labels');
           end;
       end;
       if isfield(eloc, 'labels')       
           if isnumeric(eloc(1).labels)
               disp('BESA format detected ( Elec | Theta | Phi )');
               for index = 1:length(eloc)
                   eloc(index).sph_phi_besa   = eloc(index).sph_theta_besa;
                   eloc(index).sph_theta_besa = eloc(index).labels;
                   eloc(index).labels         = eloc(index).type;
                   eloc(index).type           = '';
                   eloc(index).radius         = 1;
               end;           
           end;
       end;
       
       try
           eloc = convertlocs(eloc, 'sphbesa2all');
           eloc = convertlocs(eloc, 'topo2all'); % problem with some EGI files (not BESA files)
       catch, disp('Warning: coordinate conversion failed'); end;
       fprintf('Readlocs: BESA spherical coords. converted, now deleting BESA fields\n');   
       fprintf('          to avoid confusion (these fields can be exported, though)\n');   
       eloc = rmfield(eloc, 'sph_phi_besa');
       eloc = rmfield(eloc, 'sph_theta_besa');

       % converting XYZ coordinates to polar
       % -----------------------------------
   elseif isfield(eloc, 'sph_theta')
       try
           eloc = convertlocs(eloc, 'sph2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   elseif isfield(eloc, 'X')
       try
           eloc = convertlocs(eloc, 'cart2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   else 
       try
           eloc = convertlocs(eloc, 'topo2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   end;
   
   % inserting labels if no labels
   % -----------------------------
   if ~isfield(eloc, 'labels')
       fprintf('readlocs(): Inserting electrode labels automatically.\n');
       for index = 1:length(eloc)
           eloc(index).labels = [ 'E' int2str(index) ];
       end;
   else 
       % remove trailing '.'
       for index = 1:length(eloc)
           if isstr(eloc(index).labels)
               tmpdots = find( eloc(index).labels == '.' );
               eloc(index).labels(tmpdots) = [];
           end;
       end;
   end;
   
   % resorting electrodes if number not-sorted
   % -----------------------------------------
   if isfield(eloc, 'channum')
       if ~isnumeric(eloc(1).channum)
           error('Channel numbers must be numeric');
       end;
       allchannum = [ eloc.channum ];
       if any( sort(allchannum) ~= allchannum )
           fprintf('readlocs(): Re-sorting channel numbers based on ''channum'' column indices\n');
           [tmp newindices] = sort(allchannum);
           eloc = eloc(newindices);
       end;
       eloc = rmfield(eloc, 'channum');      
   end;
else
    if isstruct(filename)
        % detect Fieldtrip structure and convert it
        % -----------------------------------------
        if isfield(filename, 'pnt')
            neweloc = [];
            for index = 1:length(filename.label)
                neweloc(index).labels = filename.label{index};
                neweloc(index).X      = filename.pnt(index,1);
                neweloc(index).Y      = filename.pnt(index,2);
                neweloc(index).Z      = filename.pnt(index,3);
            end;
            eloc = neweloc;
            eloc = convertlocs(eloc, 'cart2all');
        else
            eloc = filename;
        end;
    else
        disp('readlocs(): input variable must be a string or a structure');
    end;        
end;
if ~isempty(g.elecind)
	eloc = eloc(g.elecind);
end;
if nargout > 2
    if isfield(eloc, 'theta')
         tmptheta = { eloc.theta }; % check which channels have (polar) coordinates set
    else tmptheta = cell(1,length(eloc));
    end;
    if isfield(eloc, 'theta')
         tmpx = { eloc.X }; % check which channels have (polar) coordinates set
    else tmpx = cell(1,length(eloc));
    end;
    
    indices           = find(~cellfun('isempty', tmptheta));
    indices           = intersect(find(~cellfun('isempty', tmpx)), indices);
    indices           = sort(indices);
    
    indbad            = setdiff(1:length(eloc), indices);
    tmptheta(indbad)  = { NaN };
    theta             = [ tmptheta{:} ];
end;
if nargout > 3
    if isfield(eloc, 'theta')
         tmprad = { eloc.radius }; % check which channels have (polar) coordinates set
    else tmprad = cell(1,length(eloc));
    end;
    tmprad(indbad)    = { NaN };
    radius            = [ tmprad{:} ];
end;

%tmpnum = find(~cellfun('isclass', { eloc.labels }, 'char'));
%disp('Converting channel labels to string');
for index = 1:length(eloc)
    if ~isstr(eloc(index).labels)
        eloc(index).labels = int2str(eloc(index).labels);
    end;
end;
labels = { eloc.labels };
if isfield(eloc, 'ignore')
    eloc = rmfield(eloc, 'ignore');
end;

% process fiducials if any
% ------------------------
fidnames = { 'nz' 'lpa' 'rpa' 'nasion' 'left' 'right' 'nazion' 'fidnz' 'fidt9' 'fidt10' 'cms' 'drl' };
for index = 1:length(fidnames)
    ind = strmatch(fidnames{index}, lower(labels), 'exact');
    if ~isempty(ind), eloc(ind).type = 'FID'; end;
end;

return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skiplines );
	 if isempty(skiplines),
       skiplines = 0;
    end;
    if exist( varname ) == 2
        array = loadtxt(varname,'verbose','off','skipline',skiplines,'blankcell','off');
    else % variable in the global workspace
         % --------------------------
         try, array = evalin('base', varname);
	     catch, error('readlocs(): cannot find the named file or variable, check syntax');
		 end;
    end;     
return;

% check field format
% ------------------
function [str, mult] = checkformat(str)
	mult = 1;
	if strcmpi(str, 'labels'),         str = lower(str); return; end;
	if strcmpi(str, 'channum'),        str = lower(str); return; end;
	if strcmpi(str, 'theta'),          str = lower(str); return; end;
	if strcmpi(str, 'radius'),         str = lower(str); return; end;
	if strcmpi(str, 'ignore'),         str = lower(str); return; end;
	if strcmpi(str, 'sph_theta'),      str = lower(str); return; end;
	if strcmpi(str, 'sph_phi'),        str = lower(str); return; end;
	if strcmpi(str, 'sph_radius'),     str = lower(str); return; end;
	if strcmpi(str, 'sph_theta_besa'), str = lower(str); return; end;
	if strcmpi(str, 'sph_phi_besa'),   str = lower(str); return; end;
	if strcmpi(str, 'gain'),           str = lower(str); return; end;
	if strcmpi(str, 'calib'),          str = lower(str); return; end;
	if strcmpi(str, 'type') ,          str = lower(str); return; end;
	if strcmpi(str, 'X'),              str = upper(str); return; end;
	if strcmpi(str, 'Y'),              str = upper(str); return; end;
	if strcmpi(str, 'Z'),              str = upper(str); return; end;
	if strcmpi(str, '-X'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Y'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Z'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, 'custom1'), return; end;
	if strcmpi(str, 'custom2'), return; end;
	if strcmpi(str, 'custom3'), return; end;
	if strcmpi(str, 'custom4'), return; end;
    error(['readlocs(): undefined field ''' str '''']);
   

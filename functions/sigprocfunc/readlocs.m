% readlocs() - read electrode location information from a file. Several standard file
%              formats are supported. Users may also specify a custom column format.
%              Examples of the defined formats are given below (File Formats).
% Usage:
%   >>  eloc = readlocs( filename );
%   >>  EEG.chanlocs = readlocs( filename, 'key', 'val', ... ); 
%   >>  [eloc, labels, theta, radius] = readlocs( filename, 'key', 'val', ... );
%
% Inputs:
%   filename   - Name of the file containing the electrode locations
%                Default is 2-D polar coordinates (see >> help topoplot )
%
% Optional inputs:
%   'filetype'  - ['loc'|'sph'|'sfp'|'xyz'|'asc'|'polhemus'|'besa'|'chanedit'|'custom'] 
%                 Type of the file to read. By default the file type is determined 
%                 using the file extension (see below under File Formats):
%                  'loc' - an EEGLAB 2-D polar coordinates channel locations file 
%                          Coordinates are theta and radius (see definitions below).
%                  'sph' - a Matlab spherical coordinates file (Note: spherical
%                          coordinates used by Matlab functions are different 
%                          from spherical coordinates used in BESA - see below).
%                  'sfp' - EGI cartesian coordinates (not Matlab cartesian - see below).
%                  'xyz' - MATLAB/EEGLAB cartesian coordinates (Not EGI cartesian; 
%                          z is toward nose; y is toward left ear; z is toward vertex).
%                  'asc' - Neuroscan polar coordinate file.
%                  'polhemus' or 'polhemusx' - Polhemus electrode location file recorded with 
%                          'X' on sensor pointing to subject (see below and readelp()).
%                  'polhemusy' - Polhemus electrode location file recorded with 
%                          'Y' on sensor pointing to subject (see below and readelp()).
%                  'besa' - BESA-'.elp' spherical coordinate file. (Not MATLAB spherical
%                           - see below).
%                  'chanedit' - EEGLAB channel location files created by pop_chanedit().
%                  'custom' - Ascii files with columns in user-defined 'format' (see below).
%   'format'    - [cell array] Format of a 'custom' channel location file (see above).
%                          Default if no file type is defined. The cell array contains
%                          labels defining the meaning of each column of the input file:
%                           'channum'   [positive integer] channel number 
%                           'labels'    [string] channel name (no spaces)
%                           'theta'     [real degrees] 2-D angle in polar coordinates; 
%                                       positive = rotating from nose (0) toward left ear 
%                           'radius'    [real] radius in 2-D polar coords (0.5 is disk limits)
%                           'X'         [real] Matlab-cartesian X coordinate (to nose)
%                           'Y'         [real] Matlab-cartesian Y coordinate (to left ear)
%                           'Z'         [real] Matlab-cartesian Z coordinate (to vertex)
%                           '-X','-Y','-Z' Matlab-cartesian coordinates pointing away from above
%                           'sph_theta' [real degrees] Matlab spherical horizontal angle; 
%                                       positive = rotating from nose (0) toward left ear.
%                           'sph_phi'   [real degrees] Matlab spherical elevation angle;
%                                       positive = rotating from horizontal (0) upwards.
%                           'sph_radius' [real] distance from head center (unused) 
%                           'sph_phi_besa' [real degrees] BESA phi angle from vertical; 
%                                       positive = rotating from vertex (0) towards right ear.
%                           'sph_theta_besa' [real degrees] BESA theta horiz/azimuthal angle; 
%                                       positive = rotating from right ear (0) toward nose.
%                           'ignore'    ignore column
%     The input file may also contain other channel information fields
%                           'type'      channel type: 'EEG', 'MEG', 'EMG', 'ECG', others ...
%                           'calib'     [real near 1.0] channel calibration value.
%                           'gain'      [real > 1] channel gain. 
%                           'custom1'   custom field #1.
%                           'custom2', 'custom3', 'custom4' more custom fields.
%   'skiplines' - [integer] Number of header lines to skip (in 'custom' file types only).
%                 All characters following '%' will be treated as comments and ignored.
%   'readchans' - [integer array] indices of electrodes to read. Default is all.
%   'center'    - [(1,3) real array or 'auto'] of xyz coordinates for conversion to 
%                 spherical or polar, Specify the center of the sphere here, or 
%                'auto'. This uses the center of the sphere that best fits all 
%                 the electrode locations read. Default is [0 0 0].
% Outputs:
%   eloc      - structure containing the channel names and locations.
%               It has three fields: 'labels', 'theta' and 'radius', identical
%               to the EEGLAB struct EEG.chanlocs.
%   labels    - cell array of strings giving the  names of the electrodes
%   theta     - vector of polar angles for the electrodes (in degrees).
%   radius    - vector of polar coordinate radius values for the electrodes
%
% File formats:
%   The extension of the file determines its type if 'filetype' is unspecified:
%
%   '.loc' or '.locs': 
%               polar coordinates. Notes: angle in degrees: right ear is 90, 
%               left ear -90; head disk radius is 0.5. 
%               Fields:   N    angle  radius    label
%               Sample:   1    -18    .511       Fp1   
%                         2     18    .511       Fp2  
%                         3    -90    .256       C3
%                         4     90    .256       C4
%                           ...
%               Note: In previous releases, channel labels had to contain exactly 
%               four characters (spaces replaced by '.'). This format still works 
%               but dots are no longer required.
%   '.sph':
%               Matlab spherical coordinates. Notes: theta is the azimuthal/horizontal angle
%               in deg.: 0 is toward nose, 90 rotated to left ear. Following this, perform
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
%               Polhemus-.'elp' cartesian coordinates. By default, an .elp extension is read
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
%               Matlab/EEGLAB cartesian coordinates. Here. x is towards the nose, 
%               y is towards the left ear, and z towards the vertex.
%               Fields:   channum   x           y         z     label
%               Sample:   1       .950        .308     -.035     Fp1
%                         2       .950       -.308     -.035     Fp2
%                         3        0           .719      .695    C3
%                         4        0          -.719      .695    C4
%                           ...
%   '.asc', '.dat':     
%               Neuroscan-.'asc' or '.dat' cartesian polar coordinates text file.
%   '.sfp': 
%               BESA/EGI-xyz cartesian coordinates. Notes: For EGI, x is toward right ear, 
%               y is toward the nose, z is toward the vertex. EEGLAB converts EGI 
%               cartesian coordinates to Matlab/EEGLAB xyz coordinates. 
%               Fields:   label   x           y          z
%               Sample:   Fp1    -.308        .950      -.035    
%                         Fp2     .308        .950      -.035  
%                         C3     -.719        0          .695  
%                         C4      .719        0          .695  
%                           ...
%   '.ced':   
%               ASCII file saved by pop_chanedit(). Contains multiple MATLAB/EEGLAB formats.
%               Carthesian coordinates are the same as the ones for the 'xyz' format.
%               Fields:   channum  label  theta  radius   x      y      z    sph_theta   sph_phi  ...
%               Sample:   1        Fp1     -18    .511   .950   .308  -.035   18         -2       ...
%                         2        Fp2      18    .511   .950  -.308  -.035  -18         -2       ...
%                         3        C3      -90    .256   0      .719   .695   90         44       ...
%                         4        C4       90    .256   0     -.719   .695  -90         44       ...
%                           ...
%               The last columns of the file may contain any other defined field (gain,
%               calib, type, custom).
%
% Author: Arnaud Delorme, Salk Institute, 8 Dec 2002 (expanded from the ICA 
%         toolbox function)
%
% See also: readelp(), writelocs(), topo2sph(), sph2topo(), sph2cart()

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
% Revision 1.55  2003/12/02 03:21:39  arno
% neuroscan format
%
% Revision 1.54  2003/11/27 00:38:13  arno
% conversion elc
%
% Revision 1.53  2003/11/27 00:31:30  arno
% debuging elc format
%
% Revision 1.52  2003/11/27 00:25:51  arno
% automatically detecting elc files
%
% Revision 1.51  2003/11/05 17:20:23  arno
% first convert spherical instead of carthesian
%
% Revision 1.50  2003/09/18 00:07:05  arno
% further checks for neuroscan
%
% Revision 1.49  2003/07/16 18:52:21  arno
% allowing file type locs
%
% Revision 1.48  2003/06/30 15:00:43  arno
% fixing inputcheck problem
%
% Revision 1.47  2003/05/13 23:31:25  arno
% number of lines to skip in chanedit format
%
% Revision 1.46  2003/05/13 22:09:01  arno
% updating sph format
%
% Revision 1.45  2003/05/13 22:07:07  arno
% removing labels in sfp format
%
% Revision 1.44  2003/05/13 21:14:11  arno
% only write a subset of file format
%
% Revision 1.43  2003/03/10 16:28:12  arno
% removing help for elc
%
% Revision 1.42  2003/03/10 16:26:59  arno
% adding then removing .elc format
%
% Revision 1.41  2003/03/08 17:36:13  arno
% import spherical EGI files correctly
%
% Revision 1.40  2003/03/05 15:38:15  arno
% fixing '.' bug
%
% Revision 1.39  2003/03/04 20:04:44  arno
% adding neuroscan .asc format
%
% Revision 1.38  2003/01/30 16:45:12  arno
% debugging ced format
%
% Revision 1.37  2003/01/10 17:40:11  arno
% removing trailing dots
%
% Revision 1.36  2003/01/03 22:47:00  arno
% typo in warning messages
%
% Revision 1.35  2003/01/03 22:45:48  arno
% adding another warning message
%
% Revision 1.34  2003/01/03 22:41:38  arno
% autodetect format .sfp
%
% Revision 1.33  2003/01/03 22:38:39  arno
% adding warning message
%
% Revision 1.32  2002/12/29 23:04:00  scott
% header
%
% Revision 1.31  2002/12/29 22:37:15  arno
% txt -> ced
%
% Revision 1.30  2002/12/29 22:35:35  arno
% adding coords. info for file format in header, programming .sph, ...
%
% Revision 1.29  2002/12/29 22:00:10  arno
% skipline -> skiplines
%
% Revision 1.28  2002/12/28 23:46:45  scott
% header
%
% Revision 1.27  2002/12/28 02:02:35  scott
% header details
%
% Revision 1.26  2002/12/28 01:32:41  scott
% worked on header information - axis details etcetc. -sm & ad
%
% Revision 1.25  2002/12/27 23:23:35  scott
% edit header msg - NEEDS MORE DETAILS -sm
%
% Revision 1.24  2002/12/27 22:57:23  arno
% debugging polhemus
%
% Revision 1.23  2002/12/27 17:47:32  arno
% compatible with more BESA formats
%
% Revision 1.22  2002/12/26 16:41:23  arno
% new release
%
% Revision 1.21  2002/12/24 02:51:22  arno
% new version of readlocs
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
% 3) enter the number of lines to skip at the end of the listskiplines array
% Note: these infos are also used by writelocs() and pop_readlocs() but
% you do not have to edit these functions.

listtype = { ...
   		 'polhemus' ...
         'polhemusX' ...
         'polhemusY' ...
         'besa' ...
         'xyz' ...
         'sfp' ...
         'loc' ...
         'sph' ...
         'asc' ...
         'dat' ...
         'elc' ...
         'chanedit' ...
         'custom' };
   
listimportformat = { ...
   	{ } ... % polhemus (specific non-columnar implementation)	
      { } ... % polhemus (specific non-columnar implementation)
      { } ... % polhemus (specific non-columnar implementation)
      { 'labels' 'sph_theta_besa' 'sph_phi_besa' 'sph_radius' } ... % BESA/EGI format
      { 'channum' 'X' 'Y' 'Z' 'labels' } ... % xyz format
      { 'labels' '-Y' 'X' 'Z' } ... % sfp format
      { 'channum' 'theta' 'radius' 'labels' } ... % loc format
      { 'channum' 'sph_theta' 'sph_phi' 'labels' } ... % sph format
      { } ... % ascii Neuroscan format
      { } ... % ascii Neuroscan format
      { } ... % eetrak format
      { 'channum' 'labels'  'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' 'sph_radius' } }; %chanedit format

listcolformat = { 'labels' 'channum' 'theta' 'radius' 'sph_theta' 'sph_phi' ...
      'sph_radius' 'sph_theta_besa' 'sph_phi_besa' 'gain' 'calib' 'type' ...
      'X' 'Y' 'Z' '-X' '-Y' '-Z' 'custom1' 'custom2' 'custom3' 'custom4' 'ignore' 'not def' };

listskipline = [ ...
   0 ... % polhemus, not applicable
   0 ... % polhemus, not applicable
   0 ... % polhemus, not applicable
   -1 ...  % besa
   0 ...
   0 ...
   0 ...
   0 ...
   0 ...
   0 ...
   0 ...
   1 ]; % skip the 2 lines header for the chanedit format

% ------------------------------------------------------
% special mode for getting the infos
% ----------------------------------
if isstr(filename) & strcmp(filename, 'getinfos')
   eloc = listtype;
   labels = listimportformat;
   theta = listcolformat;
   radius = listskipline;
   return;
elseif isstr(filename) & strcmp(filename, 'getinfoswrite')
   eloc = listtype([4:8,10]);
   labels = listimportformat([4:8,10]);
   theta = listcolformat([4:8,10]);
   radius = listskipline([4:8,10]);
   return;
end;

g = finputcheck( varargin, ...
   { 'filetype'	   'string'  {}                 '';
     'skiplines'   'integer' [0 Inf] 			[];
     'elecind'     'integer' [1 Inf]	    	[];
     'format'	   'cell'	 []					{} }, 'readlocs');
if isstr(g), error(g); end;  

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
        case 'xyz', g.filetype = 'xyz'; disp( [ 'WARNING: Matlab carthesian coords "xyz" file extension' ... 
                                       ' detected; if importing EGI cartesian coords, force to type "sfp" instead'] );
        case 'sph', g.filetype = 'sph';
        case 'ced', g.filetype = 'chanedit';
        case 'elp', g.filetype = 'polhemus';disp( [ 'WARNING: Polhemus carthesian coords "elp" file extension' ... 
                                       ' detected; if importing BESA spherical coords. force to type "besa" instead'] );
        case 'asc', g.filetype = 'asc';
        case 'dat', g.filetype = 'dat';
        case 'elc', g.filetype = 'elc';
        case 'eps', g.filetype = 'besa';
        case 'sfp', g.filetype = 'sfp';
        otherwise, g.filetype =  ''; 
       end;
       fprintf('Readlocs: ''%s'' format detected from file extension\n', g.filetype); 
   else 
       if strcmpi(g.filetype, 'locs'),  g.filetype = 'loc'; end;
   end;
   
   % assign format from filetype
   % ---------------------------
   if ~isempty(g.filetype) & ~strcmpi(g.filetype, 'custom') ...
           & ~strcmpi(g.filetype, 'asc') & ~strcmpi(g.filetype, 'elc') & ~strcmpi(g.filetype, 'dat')
      indexformat = strmatch(lower(g.filetype), lower(listtype), 'exact');
      g.format = listimportformat{indexformat};
      if isempty(g.skiplines)
         g.skiplines = listskipline(indexformat);
      end;
      if isempty(g.filetype) 
         error( ['Readlocs error: filetype can not be detected from' ...
               'file extension and custom format not specified']);
      end;
   end;
   
   % import file
   % -----------
   if strcmp(g.filetype, 'asc') | strcmp(g.filetype, 'dat')
       eloc = readneurolocs( filename );
       eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
       eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
   elseif strcmp(g.filetype, 'elc')
       eloc = readeetraklocs( filename );
       eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
       eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
   elseif strcmp(lower(g.filetype(1:end-1)), 'polhemus') | ...
           strcmp(g.filetype, 'polhemus')
       try, 
           [eloc labels X Y Z]= readelp( filename );
       catch, error('Error while reading Polhemus (for BESA .elp file force file type to BESA)'); end;
       if strcmp(g.filetype, 'polhemusy')
           tmp = X; X = Y; Y = tmp;
       end;
       for index = 1:length( eloc )
           eloc(index).X = X(index);
           eloc(index).Y = Y(index);	
           eloc(index).Z = Z(index);	
       end;
   else      
       % importing file
       % --------------
       array = load_file_or_array( filename, max(g.skiplines,0));
       if size(array,2) < length(g.format)
           fprintf('Readlocs warning: # of columns in file inferior to # format entries\n');
       elseif size(array,2) > length(g.format)
           fprintf('Readlocs warning: # of columns in file superior to # format entries\n');
       end;
       
       % removing lines BESA
       % -------------------
       if g.skiplines == -1
           if isempty(array{1,2})
               disp('BESA header detected, skipping 3 lines');
               array = load_file_or_array( filename, -2);
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
       if isnumeric(eloc(1).labels)
           disp('Alternate BESA format detected ( Theta | Phi )');
           for index = 1:length(eloc)
               eloc(index).sph_phi_besa   = eloc(index).sph_theta_besa;
               eloc(index).sph_theta_besa = eloc(index).labels;
           end;
           eloc = rmfield(eloc, 'labels');
       elseif isstr(eloc(1).sph_theta_besa)
           disp('Alternate BESA format detected ( E_type| Elec | Theta | Phi )');
           for index = 1:length(eloc)
               eloc(index).labels         = eloc(index).sph_theta_besa;
               eloc(index).sph_theta_besa = eloc(index).sph_phi_besa;
               eloc(index).sph_phi_besa   = eloc(index).sph_radius;
               eloc(index).radius         = 1;
           end;           
       end;
       try
           eloc = convertlocs(eloc, 'sphbesa2all');
           eloc = convertlocs(eloc, 'topo2all'); % problem with some EGI files (not BESA files)
       catch, disp('Warning: coordinate conversion failed'); end;
       fprintf('Readlocs: BESA spherical coords. converted, now deleting BESA fields\n');   
       fprintf('          to avoid confusion (these field can be exported though)\n');   
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
       fprintf('Readlocs: Automatically inserting electrode labels\n');
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
if nargout > 2
    theta = cell2mat({ eloc.theta });
end;
if nargout > 3
    radius  = cell2mat({ eloc.radius });
end;
%tmpnum = find(~cellfun('isclass', { eloc.labels }, 'char'));
%disp('Converting channel labels to string');
for index = 1:length(eloc)
    if ~isstr(eloc(index).labels)
        eloc(index).labels = int2str(eloc(index).labels);
    end;
end;
if nargout > 1
    labels = { eloc.labels };
end;
if isfield(eloc, 'ignore')
    eloc = rmfield(eloc, 'ignore');
end;
return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skiplines );
	 if isempty(skiplines),
       skiplines = 0;
    end;
    if exist( varname ) == 2
        array = loadtxt(varname,'verbose','off','skipline',skiplines);
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
	if strcmpi(str, 'custum1'), return; end;
	if strcmpi(str, 'custum2'), return; end;
	if strcmpi(str, 'custum3'), return; end;
	if strcmpi(str, 'custum4'), return; end;
    error(['Readlocs: undefined field ''' str '''']);
   

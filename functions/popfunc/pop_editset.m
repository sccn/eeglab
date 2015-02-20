% pop_editset() - Edit EEG dataset structure fields.
%
% Usage:
%   >> EEGOUT = pop_editset( EEG ); % pops-up a data entry window
%   >> EEGOUT = pop_editset( EEG, 'key', val,...); % no pop-up window
%
% Graphic interface:
%   "Dataset name" - [Edit box] Name for the new dataset. 
%                  In the right column of the graphic interface, the "EEG.setname"
%                  text indicates which field of the EEG structure this parameter
%                  corresponds to (in this case, .'setname').
%                  Command line equivalent: 'setname'. 
%   "Data sampling rate" - [Edit box] In Hz. Command line equivalent: 'srate'
%   "Time points per epoch" - [Edit box] Number of data frames (points) per epoch.
%                  Changing this value will change the number of data epochs.
%                  Command line equivalent: 'pnts'. 
%   "Start time" - [Edit box]  This edit box is only present for epoched data
%                  and specifies the epoch start time in ms. Epoch end time
%                  is automatically calculated. Command line equivalent: 'xmin'
%   "Number of channels" - [Edit box] Number of data channels. Command line 
%                  equivalent: 'nbchan'. This edit box cannot be edited.
%   "Ref. channel indices or mode" - [edit box] current reference. This edit box
%                  cannot be edited. To change the data reference, use menu item,
%                  'Tools > Re-reference', calling function pop_reref(). The 
%                  reference can be either a string (channel name), 'common', 
%                  indicating an unknown common reference, 'averef' indicating 
%                  average reference, or an array of integers containing indices 
%                  of the reference channel(s).
%   "Subject code" - [Edit box] subject code. For example, 'S01'. The command 
%                  line equivalent is 'subject'.
%   "Task Condition" - [Edit box] task condition. For example, 'Targets'. 
%                    The command line equivalent 'condition'.
%   "Session number" - [Edit box] session number (from the same subject). 
%                   All datasets from the same subject and session will be 
%                   assumed to use the same ICA decomposition. The command 
%                   line equivalent 'session'.
%   "Subject group" - [Edit box] subject group. For example 'Patients' or 
%                   'Control'. The command line equivalent is 'group'.
%   "About this dataset" - [Edit box] Comments about the dataset. Command line 
%                   equivalent is 'comments'.
%   "Channel locations file or array" - [Edit box] For channel data formats, see 
%                  >> readlocs help     Command line equivalent: 'chanlocs'
%   "ICA weights array or text/binary file" - [edit box] Import ICA weights from 
%                  other decompositions (e.g., same session, different conditions). 
%                  To use the ICA weights from another loaded dataset (n), enter 
%                  ALLEEG(n).icaweights. Command line equivalent: 'icaweights'
%   "ICA sphere array or text/binary file" - [edit box] Import ICA sphere matrix. 
%                  In EEGLAB, ICA decompositions require a sphere matrix and
%                  an unmixing weight matrix (see above). To use the sphere 
%                  matrix from a loaded dataset (n), enter ALLEEG(n).icasphere 
%                  Command line equivalent: 'icasphere'.
%   "From other dataset" - [push button] Press this button to enter the index
%                  of another dataset. This will update the channel locations or 
%                  the ICA edit box.
% Inputs:
%   EEG          - EEG dataset structure
%
% Optional inputs:
%   'setname'    - Name of the EEG dataset
%   'data'       - ['varname'|'filename'] Import data from a Matlab variable or 
%                  mat file into an EEG data structure 
%   'dataformat' - ['array|matlab|ascii|float32le|float32be'] Input data format.
%                  'array' is a Matlab array in the global workspace.
%                  'matlab' is a Matlab file (containing a single variable).
%                  'ascii' is an ascii file. 
%                  'float32le' and 'float32be' are 32-bit float data files 
%                  with little-endian and big-endian byte order, respectively.
%                  Data must be organised as 2-D (channels, timepoints), i.e. 
%                  channels = rows, timepoints = columns; else as 3-D (channels, 
%                  timepoints, epochs). For convenience, the data file is 
%                  transposed if the number of rows is larger than the number 
%                  of columns, as the program assumes that there are more 
%                  channels than data points. 
%   'subject'    - [string] subject code. For example, 'S01'.
%                   {default: none -> each dataset from a different subject}
%   'condition'  - [string] task condition. For example, 'Targets'
%                   {default: none -> all datasets from one condition}
%   'group'      - [string] subject group. For example 'Patients' or 'Control'.
%                   {default: none -> all subjects in one group}
%   'session'    - [integer] session number (from the same subject). All datasets
%                   from the same subject and session will be assumed to use the
%                   same ICA decomposition {default: none -> each dataset from
%                   a different session}
%   'chanlocs'   - ['varname'|'filename'] Import a channel location file.
%                   For file formats, see >> help readlocs
%   'nbchan'     - [int] Number of data channels. 
%   'xmin'       - [real] Data epoch start time (in seconds).
%                   {default: 0}
%   'pnts'       - [int] Number of data points per data epoch. The number of 
%                  data trials is automatically calculated.
%                   {default: length of the data -> continuous data assumed}
%   'srate'      - [real] Data sampling rate in Hz {default: 1Hz}
%   'ref'        - [string or integer] reference channel indices; 'averef' 
%                  indicates average reference. Note that this does not perform 
%                  referencing but only sets the initial reference when the data 
%                  are imported.
%   'icaweight'  - [matrix] ICA weight matrix. 
%   'icasphere'  - [matrix] ICA sphere matrix. By default, the sphere matrix 
%                  is initialized to the identity matrix if it is left empty.
%   'comments'   - [string] Comments on the dataset, accessible through the 
%                  EEGLAB main menu using ('Edit > About This Dataset'). 
%                  Use this to attach background information about the 
%                  experiment or the data to the dataset.
% Outputs:
%   EEGOUT       - Modified EEG dataset structure
%
% Note:
%   To create a new dataset:
%   >> EEG = pop_editset( eeg_emptyset ); % eeg_emptyset() returns an empty dataset
%
%   To erase a variable, use '[]'. The following suppresses channel locations:
%   >> EEG = pop_editset( EEG, 'chanlocs', '[]');
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_importdata(), pop_select(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 03-16-02 text interface editing -sm & ad 
% 03-16-02 remove EEG.xmax et EEG.xmin (for continuous) -ad & sm
% 03-31-02 changed interface, reprogrammed all function -ad
% 04-02-02 recompute event latencies when modifying xmin -ad

function [EEGOUT, com] = pop_editset(EEG, varargin);
   
com = '';
if nargin < 1
   help pop_editset;
   return;
end;   

EEGOUT = EEG;
if nargin < 2                 % if several arguments, assign values 
   % popup window parameters	
   % -----------------------
   % popup window parameters	
   % -----------------------
    geometry    = { [2 3.38] [1] [2.5 1 1.5 1.5] [2.5 1 1.5 1.5] [2.5 1 1.5 1.5] [2.5 1 1.5 1.5] [2.5 1 1.5 1.5] ...
                    [1] [1.4 0.7 .8 0.5] [1] [1.4 0.7 .8 0.5] [1.4 0.7 .8 0.5] [1.4 0.7 .8 0.5] };
    editcomments = [ 'tmp = pop_comments(get(gcbf, ''userdata''), ''Edit comments of current dataset'');' ...
                     'if ~isempty(tmp), set(gcf, ''userdata'', tmp); end; clear tmp;' ];
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename(1) ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    '   if strcmpi(tagtest,  ''weightfile''),' ...
                    '       set(findobj( ''parent'', gcbf, ''tag'', ''sphfile''), ''string'', ''eye(' num2str(EEG.nbchan) ')'');' ...
                    '   end;' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    commandselica = [ 'res = inputdlg2({ ''Enter dataset number'' }, ''Select ICA weights and sphere from other dataset'', 1, { ''1'' });' ...
                      'if ~isempty(res),' ...
                      '   set(findobj( ''parent'', gcbf, ''tag'', ''weightfile''), ''string'', sprintf(''ALLEEG(%s).icaweights'', res{1}));' ...
                      '   set(findobj( ''parent'', gcbf, ''tag'', ''sphfile'')   , ''string'', sprintf(''ALLEEG(%s).icasphere'' , res{1}));' ...
                      '   set(findobj( ''parent'', gcbf, ''tag'', ''icainds'')   , ''string'', sprintf(''ALLEEG(%s).icachansind'' , res{1}));' ...
                      'end;' ];
    commandselchan = [ 'res = inputdlg2({ ''Enter dataset number'' }, ''Select channel information from other dataset'', 1, { ''1'' });' ...
                      'if ~isempty(res),' ...
                      '   set(findobj( ''parent'', gcbf, ''tag'', ''chanfile''), ' ...
                      '                ''string'', sprintf(''{ ALLEEG(%s).chanlocs ALLEEG(%s).chaninfo ALLEEG(%s).urchanlocs }'', res{1}, res{1}, res{1}));' ...
                      'end;' ];
    if isstr(EEGOUT.ref)
        curref = EEGOUT.ref;
    else
        if length(EEGOUT.ref) > 1
            curref = [ int2str(abs(EEGOUT.ref)) ];
        else
            curref = [ int2str(abs(EEGOUT.ref)) ];
        end;
    end;
                        
    uilist = { ...
         { 'Style', 'text', 'string', 'Dataset name', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, { 'Style', 'edit', 'string', EEG.setname }, { } ...
         ...
         { 'Style', 'text', 'string', 'Data sampling rate (Hz)', 'horizontalalignment', 'right', 'fontweight', ...
		   'bold' }, { 'Style', 'edit', 'string', num2str(EEGOUT.srate) }, ...
         { 'Style', 'text', 'string', 'Subject code', 'horizontalalignment', 'right', ...
		    },   { 'Style', 'edit', 'string', EEG.subject }, ...
         { 'Style', 'text', 'string', 'Time points per epoch (0->continuous)', 'horizontalalignment', 'right', ...
		   },  { 'Style', 'edit', 'string', num2str(EEGOUT.pnts) }, ...
         { 'Style', 'text', 'string', 'Task condition', 'horizontalalignment', 'right', ...
		   },   { 'Style', 'edit', 'string', EEG.condition }, ...
         { 'Style', 'text', 'string', 'Start time (sec) (only for data epochs)', 'horizontalalignment', 'right', ...
		   }, { 'Style', 'edit', 'string', num2str(EEGOUT.xmin) }, ...
         { 'Style', 'text', 'string', 'Session number', 'horizontalalignment', 'right', ...
		   },   { 'Style', 'edit', 'string', EEG.session }, ...
         { 'Style', 'text', 'string', 'Number of channels (0->set from data)', 'horizontalalignment', 'right', ...
		    },   { 'Style', 'edit', 'string', EEG.nbchan 'enable' 'off' }, ...
         { 'Style', 'text', 'string', 'Subject group', 'horizontalalignment', 'right', ...
		   },   { 'Style', 'edit', 'string', EEG.group }, ...
         { 'Style', 'text', 'string', 'Ref. channel indices or mode (see help)', 'horizontalalignment', 'right', ...
		   }, { 'Style', 'edit', 'string', curref 'enable' 'off' }, ...
         { 'Style', 'text', 'string', 'About this dataset', 'horizontalalignment', 'right', ...
		   },   { 'Style', 'pushbutton', 'string', 'Enter comments' 'callback' editcomments }, ...
         { } ...
         { 'Style', 'text', 'string', 'Channel location file or info', 'horizontalalignment', 'right', 'fontweight', ...
		   'bold' }, {'Style', 'pushbutton', 'string', 'From other dataset', 'callback', commandselchan }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'chanfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''chanfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', ...
           ' Note: The file format may be auto-detected from its file extension. See menu "Edit > Channel locations" for other options.', ...
           'horizontalalignment', 'right' }, ...
         ...
         { 'Style', 'text', 'string', 'ICA weights array or text/binary file (if any):', 'horizontalalignment', 'right' }, ...
         { 'Style', 'pushbutton' 'string' 'from other dataset' 'callback' commandselica }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'weightfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''weightfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA sphere array:', 'horizontalalignment', 'right' },  ...
         { 'Style', 'pushbutton' 'string' 'from other dataset' 'callback' commandselica }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'sphfile' } ...
         { } ...
         ...
         { 'Style', 'text', 'string', 'ICA channel indices (by default all):', 'horizontalalignment', 'right' },  ...
         { 'Style', 'pushbutton' 'string' 'from other dataset' 'callback' commandselica }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'icainds' } ...
         { } };

    [ results newcomments ] = inputgui( geometry, uilist, 'pophelp(''pop_editset'');', 'Edit dataset information - pop_editset()', ...
                                         EEG.comments);
    if length(results) == 0, return; end;
	args = {};

    i = 1;
	if ~strcmp( results{i  },         EEG.setname   ) , args = { args{:}, 'setname',           results{i  }  }; end;    
	if ~strcmp( results{i+1}, num2str(EEG.srate)    ) , args = { args{:}, 'srate',     str2num(results{i+1}) }; end;
	if ~strcmp( results{i+2},         EEG.subject   ) , args = { args{:}, 'subject',           results{i+2}  }; end;
	if ~strcmp( results{i+3}, num2str(EEG.pnts)     ) , args = { args{:}, 'pnts',      str2num(results{i+3}) }; end;
	if ~strcmp( results{i+4},         EEG.condition ) , args = { args{:}, 'condition',         results{i+4}  }; end;
    if ~strcmp( results{i+5}, num2str(EEG.xmin)     ) , args = { args{:}, 'xmin',      str2num(results{i+5}) }; end;
    if ~strcmp( results{i+6}, num2str(EEG.session)  ) , args = { args{:}, 'session',   str2num(results{i+6}) }; end;
	if ~strcmp( results{i+7}, num2str(EEG.nbchan)   ) , args = { args{:}, 'nbchan',    str2num(results{i+7}) }; end;
    if ~strcmp( results{i+8},        EEG.group      ) , args = { args{:}, 'group',             results{i+8}  }; end;
    if ~strcmp( results{i+9}, num2str(EEG.ref)      ) , args = { args{:}, 'ref',               results{i+9}  }; end;
    if ~strcmp(EEG.comments, newcomments)             , args = { args{:}, 'comments' , newcomments }; end;
    
    if abs(str2num(results{i+5})) > 10,
        fprintf('WARNING: are you sure the epoch start time (%3.2f) is in seconds\n', str2num(results{i+5}));
    end;
    
	if ~isempty( results{i+12} ) , args = { args{:}, 'icachansind',  results{i+13} }; end;
	if ~isempty( results{i+10} ) , args = { args{:}, 'chanlocs' ,  results{i+10} }; end;
	if ~isempty( results{i+11} ),  args = { args{:}, 'icaweights', results{i+11} }; end;
	if ~isempty( results{i+12} ) , args = { args{:}, 'icasphere',  results{i+12} }; end;
    
else % no interactive inputs
    args = varargin;
    % Do not copy varargin
    % --------------------
    %for index=1:2:length(args)
    %    if ~isempty(inputname(index+2)) & ~isstr(args{index+1}) & length(args{index+1})>1, 
	%		args{index+1} = inputname(index+1); 
	%	end;
    %end;                
end;

% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('Setevent: wrong syntax in function arguments'); return; end;
else
    g = [];
end;

% test the presence of variables
% ------------------------------
try, g.dataformat;	 	  catch, g.dataformat = 'ascii'; end;

% assigning values
% ----------------
tmpfields = fieldnames(g);
for curfield = tmpfields'
    switch lower(curfield{1})
        case {'dataformat' }, ; % do nothing now
        case 'setname'   , EEGOUT.setname   = getfield(g, {1}, curfield{1});
        case 'subject'   , EEGOUT.subject   = getfield(g, {1}, curfield{1});
        case 'condition' , EEGOUT.condition = getfield(g, {1}, curfield{1});
        case 'group'     , EEGOUT.group     = getfield(g, {1}, curfield{1});
        case 'session'   , EEGOUT.session   = getfield(g, {1}, curfield{1});
        case 'setname'   , EEGOUT.setname   = getfield(g, {1}, curfield{1});
        case 'setname'   , EEGOUT.setname   = getfield(g, {1}, curfield{1});
        case 'pnts'      , EEGOUT.pnts      = getfield(g, {1}, curfield{1});
        case 'comments'  , EEGOUT.comments  = getfield(g, {1}, curfield{1});
        case 'nbchan'    , tmp              = getfield(g, {1}, curfield{1});
                           if tmp ~=0, EEGOUT.nbchan = tmp; end;
	    case 'averef'    , disp('The ''averef'' argument is obsolete; use function pop_reref() instead');
        case 'ref'       , EEGOUT.ref       = getfield(g, {1}, curfield{1});
                           disp('WARNING: CHANGING REFERENCE DOES NOT RE-REFERENCE THE DATA, use menu Tools > Rereference instead');
                           if ~isempty(str2num( EEGOUT.ref )), EEG,ref = str2num(EEG.ref); end;
        case 'xmin'    , oldxmin = EEG.xmin;
                         EEGOUT.xmin = getfield(g, {1}, curfield{1});
                         if oldxmin ~= EEGOUT.xmin
                             if ~isempty(EEG.event)
                                 if nargin < 2
                                     if ~popask( ['Warning: changing the starting point of epochs will' 10 'lead to recomputing epoch event latencies, Continue?'] )
                                         com = ''; warndlg2('pop_editset(): transformation cancelled by user'); return; 
                                     end;
                                 end;
                                 if isfield(EEG.event, 'latency')
                                     for index = 1:length(EEG.event)
                                         EEG.event(index).latency = EEG.event(index).latency - (EEG.xmin-oldxmin)*EEG.srate;
                                     end;
                                 end;       
                             end;    
                         end;
        case 'srate'   , EEGOUT.srate = getfield(g, {1}, curfield{1});
        case 'chanlocs', varname = getfield(g, {1}, curfield{1});
                         if isempty(varname)
                             EEGOUT.chanlocs = [];
                         elseif isstr(varname) & exist( varname ) == 2
                            fprintf('pop_editset(): channel locations file ''%s'' found\n', varname); 
                            [ EEGOUT.chanlocs lab theta rad ind ] = readlocs(varname);
                         elseif isstr(varname)
                            EEGOUT.chanlocs = evalin('base', varname, 'fprintf(''pop_editset() warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
                            if iscell(EEGOUT.chanlocs)
                                if length(EEGOUT.chanlocs) > 1, EEGOUT.chaninfo   = EEGOUT.chanlocs{2}; end;
                                if length(EEGOUT.chanlocs) > 2, EEGOUT.urchanlocs = EEGOUT.chanlocs{3}; end;
                                EEGOUT.chanlocs = EEGOUT.chanlocs{1};
                            end;
                         else
                             EEGOUT.chanlocs = varname;
                         end;
        case 'icaweights', varname = getfield(g, {1}, curfield{1});
                         if isstr(varname) & exist( varname ) == 2
                            fprintf('pop_editset(): ICA weight matrix file ''%s'' found\n', varname);
                            if ~isempty(EEGOUT.icachansind), nbcol = length(EEGOUT.icachansind); else nbcol = EEG.nbchan; end;                            
							try, EEGOUT.icaweights = load(varname, '-ascii');
								EEGOUT.icawinv = [];
                            catch, 
                                try
                                    EEGOUT.icaweights = floatread(varname, [1 Inf]);
                                    EEGOUT.icaweights = reshape( EEGOUT.icaweights, [length(EEGOUT.icaweights)/nbcol nbcol]);
                                catch
                                    fprintf('pop_editset() warning: error while reading filename ''%s'' for ICA weight matrix\n', varname); 
                                end;
                            end;
                         else
							 if isempty(varname) 
								 EEGOUT.icaweights = [];
							 elseif isstr(varname)
								 EEGOUT.icaweights = evalin('base', varname, 'fprintf(''pop_editset() warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
								 EEGOUT.icawinv = [];
                             else
								 EEGOUT.icaweights = varname;
								 EEGOUT.icawinv = [];                                 
							 end;
						 end;
                         if ~isempty(EEGOUT.icaweights) & isempty(EEGOUT.icasphere)
                            EEGOUT.icasphere = eye(size(EEGOUT.icaweights,2));
                         end;
        case 'icachansind', varname = getfield(g, {1}, curfield{1});
							 if isempty(varname) 
								 EEGOUT.icachansind = [];
							 elseif isstr(varname)
								 EEGOUT.icachansind = evalin('base', varname, 'fprintf(''pop_editset() warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
                             else
  								 EEGOUT.icachansind = varname;
							 end;
        case 'icasphere', varname = getfield(g, {1}, curfield{1});
                         if isstr(varname) & exist( varname ) == 2
                            fprintf('pop_editset(): ICA sphere matrix file ''%s'' found\n', varname); 
                            if ~isempty(EEGOUT.icachansind), nbcol = length(EEGOUT.icachansind); else nbcol = EEG.nbchan; end;
                            try, EEGOUT.icasphere = load(varname, '-ascii');
								EEGOUT.icawinv = [];
                            catch,
                                try
                                    EEGOUT.icasphere = floatread(varname, [1 Inf]);
                                    EEGOUT.icasphere = reshape( EEGOUT.icasphere, [length(EEGOUT.icasphere)/nbcol nbcol]);
                                catch
                                    fprintf('pop_editset() warning: erro while reading filename ''%s'' for ICA weight matrix\n', varname); 
                                end;
                            end;
                         else
							 if isempty(varname) 
								 EEGOUT.icasphere = [];
							 elseif isstr(varname)
								 EEGOUT.icasphere = evalin('base', varname, 'fprintf(''pop_editset() warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
								 EEGOUT.icawinv = [];
                             else
  								 EEGOUT.icaweights = varname;
								 EEGOUT.icawinv = [];                                 
							 end;
                         end;
                         if ~isempty(EEGOUT.icaweights) & isempty(EEGOUT.icasphere)
                            EEGOUT.icasphere = eye(size(EEGOUT.icaweights,2));
                         end;
	    case 'data'    , varname = getfield(g, {1}, curfield{1});
                         if isnumeric(varname)
                             EEGOUT.data = varname;
                         elseif exist( varname ) == 2 & ~strcmp(lower(g.dataformat), 'array');
                            fprintf('pop_editset(): raw data file ''%s'' found\n', varname); 
                            switch lower(g.dataformat)
							 case 'ascii' , 
							  try, EEGOUT.data = load(varname, '-ascii');
							  catch, disp(lasterr); error(['pop_editset() error: cannot read ascii file ''' varname ''' ']); 
							  end;
							  if ndims(EEGOUT.data)<3 & size(EEGOUT.data,1) > size(EEGOUT.data,2), EEGOUT.data = transpose(EEGOUT.data); end;
							 case 'matlab', 
							  try,
								  x = whos('-file', varname);
								  if length(x) > 1, 
									  error('pop_editset() error: .mat file must contain a single variable'); 
								  end;
								  tmpdata = load(varname, '-mat');									  
								  EEGOUT.data = getfield(tmpdata,{1},x(1).name);
								  clear tmpdata;
							  catch, error(['pop_editset() error: cannot read .mat file ''' varname ''' ']); 
							  end;
							  if ndims(EEGOUT.data)<3 & size(EEGOUT.data,1) > size(EEGOUT.data,2), EEGOUT.data = transpose(EEGOUT.data); end;
							 case {'float32le' 'float32be'}, 
							  if EEGOUT.nbchan == 0,
								  error(['pop_editset() error: to read float32 data you must first specify the number of channels']);
							  end;     
							  try, EEGOUT.data = floatread(varname, [EEGOUT.nbchan Inf], ...
														   fastif(strcmpi(g.dataformat, 'float32le'), 'ieee-le', 'ieee-be'));
							  catch, error(['pop_editset() error: cannot read float32 data file ''' varname ''' ']); 
							  end;
							 otherwise, error('pop_editset() error: unrecognized file format');
                            end;
                         elseif isstr(varname)
                             % restoration command
                             %--------------------
                             try 
                                 res = evalin('base', ['exist(''' varname ''') == 1']);
                             catch
                                 disp('pop_editset() warning: cannot find specified variable in global workspace!');
                             end;
                             if ~res, 
                                 error('pop_editset(): cannot find specified variable.'); 
                             end;
                             warning off;
                             try,
                                 testval = evalin('base', ['isglobal(' varname ')']);
                             catch, testval = 0; end;
                             if ~testval
                                 commandrestore = [ ' tmpp = '  varname '; clear global ' varname ';'   varname '=tmpp;clear tmpp;' ]; 
                             else
                                 commandrestore = [];
                             end;		  
                             % make global, must make these variable global, if you try to evaluate them direclty in the base
                             % workspace, with a large array the computation becomes incredibly slow.  
                             %--------------------------------------------------------------------
                             comglobal = sprintf('global %s;', varname);
                             evalin('base', comglobal);
                             eval(comglobal);
                             eval( ['EEGOUT.data = ' varname ';' ]);
                             try, evalin('base', commandrestore); catch, end;
                             warning on;
                         else 
                             EEGOUT.data = varname;
                             if ndims(EEGOUT.data)<3 & size(EEGOUT.data,1) > size(EEGOUT.data,2), EEGOUT.data = transpose(EEGOUT.data); end;
                         end;
         otherwise, error(['pop_editset() error: unrecognized field ''' curfield{1} '''']); 
    end;
end;

EEGOUT = eeg_checkset(EEGOUT);

% generate the output command
% ---------------------------
if nargout > 1
    com = sprintf( '%s = pop_editset(%s', inputname(1), inputname(1) );
    for i=1:2:length(args)
        if ~isempty( args{i+1} )
            if isstr( args{i+1} ) com = sprintf('%s, ''%s'', %s', com, args{i}, vararg2str(args{i+1}) );
            else                  com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
            end;
        else
            com = sprintf('%s, ''%s'', []', com, args{i} );
        end;
    end;
    com = [com ');'];
end;
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;

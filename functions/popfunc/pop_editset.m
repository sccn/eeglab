% pop_editset() - Edit EEG dataset structure fields.
%
% Usage:
%   >> EEGOUT = pop_editset( EEG ); % pops-up a data entry window
%   >> EEGOUT = pop_editset( EEG, 'key', val,...); % no pop-up window
%
% Graphic interface:
%   "EEGLAB dataset name" - [Edit box] Name for the new dataset. 
%                  In the last column of the graphic interface, the "EEG.setname"
%                  text indicates which field of the EEG structure this parameter
%                  is corresponding to (in this case 'setname').
%                  Command line equivalent: 'setname'. 
%   "Time points per epoch" - [Edit box] Number of data frames (points) per epoch.
%                  Command line equivalent: 'pnts'
%   "Data sampling rate" - [Edit box] In Hz. Command line equivalent: 'srate'
%   "Optional epoch start time" - [Edit box]  This edit box is only present for 
%                  data epoch and specify the epochs start time in ms. Epoch upper
%                  time limit is automatically calculated. 
%                  Command line equivalent: 'xmin'
%   "Channel locations file or array" - [Edit box] For channel data formats, see 
%                  >> readlocs help     Command line equivalent: 'chanlocs'
%   "ICA weights array or text file" - [edit box] Import ICA weights from other 
%                  decompositions (e.g., same data, different conditions). 
%                  To use the ICA weights from another loaded dataset (n), enter 
%                  ALLEEG(n).icaweights. Command line equivalent: 'icaweights'
%   "ICA sphere array or text file" - [edit box] Import ICA sphere matrix. 
%                  Infomax ICA decompositions may be defined by a sphere matrix 
%                  and an unmixing weight matrix (see above).  To use the sphere 
%                  matrix from another loaded dataset (n), enter ALLEEG(n).icasphere 
%                  Command line equivalent: 'icasphere'.
%   "Averaged referenced data" - [checkbox] Re-reference data to average reference 
%                  by checking the checkbox. Transform back to common reference 
%                  by unchecking the checkbox. Command line equivalent: 'averef' 
%                  See also pop_reref().
% Inputs:
%   EEG          - EEG dataset structure
%
% Optional inputs:
%   'setname'    - Name of the EEG dataset
%   'data'       - ['varname'|'filename'] Import data from a Matlab variable or file
%                  into an EEG data structure 
%   'dataformat' - ['array|matlab|ascii|float32le|float32be'] Input data format.
%                  'array' is a Matlab array in the global workspace.
%                  'matlab' is a Matlab file (which must contain a single variable).
%                  'ascii' is an ascii file. 'float32le' and 'float32be' are 32-bits
%                  float data files (little endian or big endian byte ordering).
%                  Data must be organised as (channels, timepoints) i.e. 
%                  channels = rows and timepoints = columns or (channels, timepoints, 
%                  epochs). For convenience, The data file is transposed if the number
%                  of rows is larger than the number of columns.
%   'chanlocs'   - ['varname'|'filename'] Import a channel location file.
%                  For file formats, see >> help readlocs
%   'nbchan'     - [int] Number of data channels. 
%   'xmin'       - [real] Data start time (in seconds).
%   'averef'     - ['Yes'|'No'] 'Yes' if data are average-reference. 
%   'pnts'       - [int] Number of data points per epoch (epoched data only)
%   'srate'      - [real] Data sampling rate in Hz. 
%   'icaweight'  - [matrix] ICA weight matrix. 
%   'icasphere'  - [matrix] ICA sphere matrix. By default, the sphere matrix 
%                  is initialized to the identity matrix if it is left empty.
%   'comments'   - [string] Comments on the dataset accessible through the EEGLAB
%                  main menu (Edit > About This Dataset). Use this to attach 
%                  background information about the data to the dataset.
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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.34  2003/02/25 01:00:53  scott
% header edit -sm
%
% Revision 1.33  2003/02/24 16:26:38  arno
% resolving ???
%
% Revision 1.32  2003/02/22 17:11:31  scott
% header edit -sm
%
% Revision 1.31  2003/02/21 22:55:09  arno
% adding gui info
%
% Revision 1.30  2003/01/23 21:34:59  scott
% header edits -sm
%
% Revision 1.29  2003/01/22 22:43:53  scott
% edit header -sm
%
% Revision 1.28  2002/11/14 18:29:03  arno
% new average reference
%
% Revision 1.27  2002/11/14 17:38:31  arno
% gui text
%
% Revision 1.26  2002/11/08 19:05:05  arno
% test if array=string before importing it
%
% Revision 1.25  2002/09/25 23:49:58  arno
% correcting float-le problem
%
% Revision 1.24  2002/09/04 18:28:25  luca
% debug command line big variable passed as text - arno
%
% Revision 1.23  2002/08/29 23:17:01  arno
% debugging icaweights and sphere
%
% Revision 1.22  2002/08/26 22:03:01  arno
% average reference more message
%
% Revision 1.21  2002/08/12 18:31:31  arno
% questdlg2
%
% Revision 1.20  2002/07/31 18:12:35  arno
% reading floats le and be
%
% Revision 1.19  2002/05/21 20:47:02  scott
% removed ; from evalin() commands -sm
%
% Revision 1.18  2002/05/01 01:23:44  luca
% same
%
% Revision 1.17  2002/05/01 01:22:56  luca
% same
%
% Revision 1.16  2002/05/01 01:21:18  luca
% transpose bug
%
% Revision 1.15  2002/04/30 18:38:21  arno
% adding about button
%
% Revision 1.14  2002/04/18 16:19:38  scott
% EEG.averef -sm
%
% Revision 1.13  2002/04/18 16:13:20  scott
% working on EEG.averef -sm
%
% Revision 1.12  2002/04/18 14:43:02  scott
% edited error msgs -sm
%
% Revision 1.11  2002/04/18 02:52:05  scott
% [same] -sm
%
% Revision 1.10  2002/04/18 02:44:08  scott
% edited error messages -sm
%
% Revision 1.9  2002/04/18 02:29:24  arno
% further checks for matlab file import
%
% Revision 1.8  2002/04/11 18:28:58  arno
% adding average reference input
%
% Revision 1.7  2002/04/11 18:01:24  arno
% removing warning when removing ICA components
%
% Revision 1.6  2002/04/10 22:42:11  arno
% debuging variable name
%
% Revision 1.5  2002/04/08 02:29:33  scott
% *** empty log message ***
%
% Revision 1.4  2002/04/08 02:26:28  scott
% *** empty log message ***
%
% Revision 1.3  2002/04/08 02:25:00  scott
% *** empty log message ***
%
% Revision 1.2  2002/04/08 02:22:35  scott
% improved workding, moved 'EEG.data' placement -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

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
    geometry    = { [1.5 0.6 1 0.7] [1.1 1 1 0.7]   [2 0.1 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] ...
					[2 0.1 1 0.7] [1.5 0.6 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] };
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
	editcomments = 'set(gcf, ''userdata'', pop_comments(get(gcbf, ''userdata''), ''Edit comments of current dataset''));';
		
    uilist = { ...
         { 'Style', 'text', 'string', 'Dataset name (optional):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
		 { 'Style', 'pushbutton', 'string', 'About', 'callback', editcomments },  ...
		 { 'Style', 'edit', 'string', EEG.setname }, { 'Style', 'text', 'string', '(EEG.setname)' }...
         ...
         { 'Style', 'text', 'string', 'Data file or array', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { }, { }, { 'Style', 'text', 'string', '(EEG.data)' }, ...
         ...
         { 'Style', 'text', 'string', 'Channels in data:',  'horizontalalignment', 'right', 'fontweight', 'bold' }, {},  ...
         { 'Style', 'text', 'string', num2str(EEG.nbchan),  'horizontalalignment', 'center' }, ...
		 { 'Style', 'text', 'string', '(EEG.nbchan)' } ...
		 ...
         { 'Style', 'text', 'string', 'Time points per epoch:', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
		 { },  { 'Style', 'edit', 'string', num2str(EEG.pnts) }, {'Style', 'text', 'string', '(EEG.pnts)' } ...
		 ...
         { 'Style', 'text', 'string', 'Data sampling rate (Hz):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { }, ... 
		 { 'Style', 'edit', 'string', num2str(EEG.srate) }, {'Style', 'text', 'string', '(EEG.srate)' },...
         ...
		 { 'Style', 'text', 'string', 'Epoch start time (sec):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { }, ...
		 { 'Style', 'edit', 'string', num2str(EEG.xmin) }, {'Style', 'text', 'string', '(EEG.xmin)' },...
         ...
         { 'Style', 'text', 'string', 'Channel position file or array:', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
		 {'Style', 'pushbutton', 'string', 'Help', 'callback', 'pophelp(''readlocs'');' }, ...
		 { 'Style', 'edit', 'string', fastif(isempty(EEG.chanlocs), '', 'EEG.chanlocs'), 'horizontalalignment', 'left', 'tag',  'chanfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''chanfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA weight array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
         { 'Style', 'edit', 'string', fastif(isempty(EEG.icaweights), '', 'EEG.icaweights'), 'horizontalalignment', 'left', 'tag',  'weightfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''weightfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA sphere array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
         { 'Style', 'edit', 'string', fastif(isempty(EEG.icasphere), '', 'EEG.icasphere'), 'horizontalalignment', 'left', 'tag',  'sphfile' } ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''sphfile'';' commandload ] } ...
	     ...
		 { 'Style', 'text', 'string', 'Averaged referenced data ?'} { } ...
		 { 'Style', 'checkbox', 'string', '(set = Yes)', 'value', ~strcmpi(EEG.ref, 'common') }, { 'Style', 'text', 'string', '(EEG.ref)'} ...
			 };

    if EEG.trials == 1,  uilist(21:24) = []; geometry(6) = []; end;
    [results newcomments] = inputgui( geometry, uilist, 'pophelp(''pop_editset'');', ...
						fastif(isempty(EEG.data), 'Import dataset info -- pop_editset()', 'Edit dataset info -- pop_editset()'), EEG.comments);
    if length(results) == 0, return; end;

	args = {};
	if ~strcmp( results{1}, EEG.setname ), args = { args{:}, 'setname', results{1} }; end;
	if ~strcmp( results{2}, num2str(EEG.pnts) )  , args = { args{:}, 'pnts', str2num(results{2}) }; end;
	if ~strcmp( results{3}, num2str(EEG.srate) )  , args = { args{:}, 'srate', str2num(results{3}) }; end;
	if EEG.trials ~= 1,
	   if ~strcmp( results{4}, num2str(EEG.xmin) ), args = { args{:}, 'xmin', str2num(results{4}) }; end;
	   i = 5;
	else i = 4;
	end;
	if ~strcmp( results{i  }, fastif(isempty(EEG.chanlocs), '', 'EEG.chanlocs')  ) , args = { args{:}, 'chanlocs' , results{i} }; end;
	if ~strcmp( results{i+1}, fastif(isempty(EEG.icaweights), '', 'EEG.icaweights') ), args = { args{:}, 'icaweights', results{i+1} }; end;
	if ~strcmp( results{i+2}, fastif(isempty(EEG.icasphere), '', 'EEG.icasphere') ) , args = { args{:}, 'icasphere', results{i+2} }; end;
	if strcmpi(EEG.ref, 'averef') ~= results{i+3}, args = { args{:}, 'averef', fastif(results{i+3}, 'yes', 'no') }; end;
	if ~strcmp(EEG.comments, newcomments), args = { args{:}, 'comments' , newcomments }; end;
else % no interactive inputs
    args = varargin;
    for index=1:2:length(args)
        if ~isempty(inputname(index+2)) & ~isstr(args{index+1}) & length(args{index+1})>1, 
			args{index+1} = inputname(index+2); 
		end;
    end;                
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
        case 'setname' , EEGOUT.setname = getfield(g, {1}, curfield{1});
        case 'pnts'    , EEGOUT.pnts = getfield(g, {1}, curfield{1});
        case 'comments', EEGOUT.comments = getfield(g, {1}, curfield{1});
	    case 'averef'  , reref = getfield(g, {1}, curfield{1});
	                     if (strcmpi(EEGOUT.ref, 'common') & strcmpi(reref, 'yes')) | (~strcmpi(EEGOUT.ref, 'common') & strcmpi(reref, 'no'))
							 EEGOUT = pop_reref(EEG, []);
						 end;
        case 'nbchan'  , EEGOUT.nbchan = getfield(g, {1}, curfield{1});
        case 'xmin'    , oldxmin = EEG.xmin;
                         EEGOUT.xmin = getfield(g, {1}, curfield{1});
                         if ~isempty(EEG.event)
                             if nargin < 2
                                if ~popask( ['Warning: changing the starting point of epochs will' 10 'lead to recomputing epoch event latencies ?'] )
                                    error('Pop_editset: transformation cancelled by user');
                                end;
                             end;
                             if isfield(EEG.event, 'latency')
                                for index = 1:length(EEG.event)
                                    EEG.event(index).latency = EEG.event(index).latency - (EEG.xmin-oldxmin)*EEG.srate;
                                end;
                             end;       
                         end;    
        case 'srate'   , EEGOUT.srate = getfield(g, {1}, curfield{1});
        case 'chanlocs', varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2
                            fprintf('Pop_editset: channel locations file ''%s'' found\n', varname); 
                            EEGOUT.chanlocs = readlocs(varname);
                         else
                            EEGOUT.chanlocs = evalin('base', varname, 'fprintf(''Pop_editset warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
                         end;
        case 'icaweights', varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2
                            fprintf('Pop_editset: ICA weight matrix file ''%s'' found\n', varname); 
							try, EEGOUT.icaweights = load(varname, '-ascii');
								EEGOUT.icawinv = [];
                            catch, fprintf('Pop_editset warning: error while reading filename ''%s'' for ICA weight matrix\n', varname); 
                            end;
                         else
							 if isempty(varname) 
								 EEGOUT.icaweights = [];
							 else
								 EEGOUT.icaweights = evalin('base', varname, 'fprintf(''Pop_editset warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
								 EEGOUT.icawinv = [];
							 end;
						 end;
                         if ~isempty(EEGOUT.icaweights) & isempty(EEGOUT.icasphere)
                            EEGOUT.icasphere = eye(size(EEGOUT.icaweights,2));
                         end;
        case 'icasphere', varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2
                            fprintf('Pop_editset: ICA sphere matrix file ''%s'' found\n', varname); 
                            try, EEGOUT.icasphere = load(varname, '-ascii');
								EEGOUT.icawinv = [];
                            catch, fprintf('Pop_editset warning: erro while reading filename ''%s'' for ICA weight matrix\n', varname); 
                            end;
                         else
							 if isempty(varname) 
								 EEGOUT.icasphere = [];
							 else
								 EEGOUT.icasphere = evalin('base', varname, 'fprintf(''Pop_editset warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
								 EEGOUT.icawinv = [];
							 end;
                         end;
	    case 'data'    , varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2 & ~strcmp(lower(g.dataformat), 'array');
                            fprintf('Pop_editset: raw data file ''%s'' found\n', varname); 
                            switch lower(g.dataformat)
							 case 'ascii' , 
							  try, EEGOUT.data = load(varname, '-ascii');
							  catch, error(['Pop_editset error: cannot read ascii file ''' varname ''' ']); 
							  end;
							  if ndims(EEGOUT.data)<3 & size(EEGOUT.data,1) > size(EEGOUT.data,2), EEGOUT.data = transpose(EEGOUT.data); end;
							 case 'matlab', 
							  try,
								  x = whos('-file', varname);
								  if length(x) > 1, 
									  error('Pop_editset error: .mat file must contain a single variable'); 
								  end;
								  tmpdata = load(varname, '-mat');									  
								  EEGOUT.data = getfield(tmpdata,{1},x(1).name);
								  clear tmpdata;
							  catch, error(['Pop_editset error: cannot read .mat file ''' varname ''' ']); 
							  end;
							  if ndims(EEGOUT.data)<3 & size(EEGOUT.data,1) > size(EEGOUT.data,2), EEGOUT.data = transpose(EEGOUT.data); end;
							 case {'float32le' 'float32be'}, 
							  if EEGOUT.nbchan == 0,
								  error(['Pop_editset error: to read float32 data you must first specify the number of channels']);
							  end;     
							  try, EEGOUT.data = floatread(varname, [EEGOUT.nbchan Inf], ...
														   fastif(strcmpi(g.dataformat, 'float32le'), 'ieee-le', 'ieee-be'));
							  catch, error(['Pop_editset error: cannot read float32 data file ''' varname ''' ']); 
							  end;
							 otherwise, error('Pop_editset error: unrecognized file format');
                            end;
                         elseif isstr(varname)
                             % restoration command
                             %--------------------
                             try 
                                 res = evalin('base', ['exist(''' varname ''') == 1']);
                             catch
                                 disp('Pop_editset warning: cannot find specified variable in global workspace!');
                             end;
                             if ~res, 
                                 error('Pop_editset: cannot find specified variable.'); 
                             end;
                             testval = evalin('base', ['isglobal(' varname ')']);
                             warning off;
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
         otherwise, error(['Pop_editset error: unrecognized field ''' curfield{1} '''']); 
    end;
end;

% generate the output command
% ---------------------------
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
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;

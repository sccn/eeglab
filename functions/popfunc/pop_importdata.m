% pop_importdata() - Import data from a Matlab variable or from a file. 
%
% Usage:
%   >> EEGOUT = pop_importdata( EEG ); % pop-up window mode
%   >> EEGOUT = pop_importdata( 'key', val,...);
%
% Graphical interface:
%   "EEGLAB dataset name" - [Edit box] name for the new dataset. Command line
%                  equivalent: 'setname'
%   "Data file/array" - [Edit box] Data file or Matlab variable name to import
%                  to EEGLAB. Command line equivalent: 'data'
%   "Data file/array" - [list box] select data format from listbox. If you
%                  browse for a data file, the graphical interface might be
%                  able to detect the file format from the file extension and
%                  his list box accordingly. Note that you have to click on
%                  the option to make it active. Command line equivalent:
%                  'dataformat'
%   "Number of channels" - [Edit box] number of data channel. Command line
%                  equivalent: 'nbchan'
%   "Time points per epoch" - [Edit box] Number of points per data frame for 
%                  data epochs. Command line equivalent: 'pnts'
%   "Data sampling rate" - [Edit box] command line equivalent: 'srate'
%   "Optional epoch start time" - [Edit box] command line equivalent: 'xmin'
%   "Channel locations file or array" - [Edit box] see readlocs() help for
%                  data channel format. Command line equivalent: 'chanlocs'
%   "ICA weights array or text file" - [edit box] use this option to import
%                  ICA weight from other decompositions (for instance: same
%                  data, different conditions). To use the ICA weights from
%                  an other dataset (i.e. dataset 2) enter "ALLEEG(2).icaweights"
%                  in this edit box. Command line equivalent: 'icaweight'
%   "ICA sphere array or text file" - [edit box] import ICA sphere matrix. For
%                  computational reasons, an ICA decomposition is defined by
%                  a sphere matrix and an unmixing matrix (see previous option).
%                  To use the ICA weights from an other dataset (i.e. dataset 2)
%                  enter "ALLEEG(2).icasphere" in this edit box. Command line 
%                  equivalent: 'icasphere'.
%
% Optional inputs:
%   'setname'    - name of the new EEGLAB dataset
%   'data'       - ['varname'|'filename'] Data file or variable to import
%                  into EEGLAB.
%   'dataformat' - ['array|matlab|ascii|float32le|float32be'] input data file format.
%                  The ata file is transposed if the number of rows is greater
%                  than the number of columns. Note: for type 'float32le' and
%                  'float32be' (little endian and big endian byte ordering), data
%                  must be organised in the format (channels x timepoints x epochs).
%   'chanlocs'   - ['varname'|'filename'] Import a file containing
%                  electrode locations (see >> help readlocs for file format).
%   'nbchan'     - Number of data channels 
%   'xmin'       - Starting time in seconds
%   'pnts'       - Number of points per data frame (epoched data only)
%   'srate'      - Data sampling rate
%   'icaweight'  - ICA weight matrix. 
%   'icasphere'  - ICA sphering matrix (if [], eye(nchans)).
% 
% Outputs:
%   EEGOUT      - modified dataset structure
%
% Note: this function call pop_editset() to modify parameter values.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_editset(), pop_select(), eeglab()

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
% Revision 1.14  2002/12/18 22:25:46  arno
% Automatic file format detection debug
%
% Revision 1.13  2002/12/05 02:26:41  arno
% adding an additionnal warning for clicking on the selected option
%
% Revision 1.12  2002/09/25 23:50:01  arno
% correcting float-le problem
%
% Revision 1.11  2002/09/04 18:30:31  luca
% same
%
% Revision 1.10  2002/09/04 18:30:11  luca
% same
%
% Revision 1.9  2002/09/04 18:28:23  luca
% 'debug command line big variable passed as text - arno
%
% Revision 1.8  2002/09/04 18:24:09  luca
% debug for command line - arno
%
% Revision 1.7  2002/07/31 18:02:26  arno
% adding more options
%
% Revision 1.6  2002/05/02 23:55:22  arno
% auto file type selection
%
% Revision 1.5  2002/04/20 23:53:14  scott
% editted screen items -sm
%
% Revision 1.4  2002/04/18 02:35:24  arno
% put default Matlab file read
%
% Revision 1.3  2002/04/11 22:22:52  arno
% removing comment
%
% Revision 1.2  2002/04/11 21:18:34  arno
% *** empty log message ***
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-16-02 text interface editing -sm & ad 
% 03-16-02 remove EEG.xmax et EEG.xmin (for continuous) -ad & sm
% 04-02-02 debugging command line calls -ad & lf

function [EEGOUT, com] = pop_importdata( varargin);

com = '';
EEGOUT = eeg_emptyset;
if nargin < 1                 % if several arguments, assign values 
   % popup window parameters	
   % -----------------------
    geometry    = { [2 0.1 0.8 0.5] [1.3 0.8 .8 0.5] [2 0.1 0.8 0.5] [2 0.1 0.8 0.5] [2 0.1 0.8 0.5] ...
                    [2 0.1 0.8 0.5] [1.5 0.6 0.8 0.5] [2 0.1 0.8 0.5] [2 0.1 0.8 0.5] };
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
	commandsetfiletype = [ 'filename = get( findobj(''parent'', gcbf, ''tag'', ''globfile''), ''string'');' ...
					'tmpext = findstr(filename,''.'');' ...
					'tmpext = lower(filename(tmpext(end)+1:end));' ...
					'switch tmpext, ' ...
					'  case ''mat'', set(findobj(gcbf,''tag'', ''loclist''), ''value'',5);' ...
					'  case ''fdt'', set(findobj(gcbf,''tag'', ''loclist''), ''value'',3);' ...
					'  case ''txt'', set(findobj(gcbf,''tag'', ''loclist''), ''value'',2);' ...
					'end; clear tmpext filename;' ];
    uilist = { ...
         { 'Style', 'text', 'string', 'EEGLAB dataset name (optional):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', '' }, { }...
         ...
         { 'Style', 'text', 'string', 'Data file/array (click on selected option)', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { 'Style', 'listbox', 'string', 'Matlab variable|ASCII text file|float32 le file|float32 be file|Matlab .mat file', ...
		   'fontweight', 'bold', 'tag','loclist' } ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'globfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', ...
		   [ 'tagtest = ''globfile'';' commandload commandsetfiletype ] }, ...
         ...
         { 'Style', 'text', 'string', 'Number of channels (0->set from data):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, {},  { 'Style', 'edit', 'string', '0' }, { } ...
         { 'Style', 'text', 'string', 'Time points per epoch (0=continuous data):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', num2str(EEGOUT.pnts) }, { } ...
         { 'Style', 'text', 'string', 'Data sampling rate (Hz):', 'horizontalalignment', 'right', 'fontweight', ...
		   'bold' }, { }, { 'Style', 'edit', 'string', num2str(EEGOUT.srate) }, { },...
         { 'Style', 'text', 'string', 'Optional epoch start time for data epochs (sec):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, { }, { 'Style', 'edit', 'string', num2str(EEGOUT.xmin) }, { },...
         ...
         { 'Style', 'text', 'string', 'Channel locations file or array:', 'horizontalalignment', 'right', 'fontweight', ...
		   'bold' }, {'Style', 'pushbutton', 'string', 'Help', 'callback', 'pophelp(''readlocs'');' }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'chanfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''chanfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA weights array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'weightfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''weightfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA sphere array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'sphfile' } ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''sphfile'';' commandload ] } };

    results = inputgui( geometry, uilist, 'pophelp(''pop_importdata'');', 'Import dataset info -- pop_importdata()');
    if length(results) == 0, return; end;

	args = {};
	if ~isempty( results{1} ), args = { args{:}, 'setname', results{1} }; end;
	switch results{2}
	   case 1, args = { args{:}, 'dataformat', 'array' };
	   case 2, args = { args{:}, 'dataformat', 'ascii' };
	   case 3, args = { args{:}, 'dataformat', 'float32le' };
	   case 4, args = { args{:}, 'dataformat', 'float32be' };
	   case 5, args = { args{:}, 'dataformat', 'matlab' };
	end;
	if ~isempty( results{3} ) , args = { args{:}, 'data', results{3} }; end;
	if ~isempty( results{4} ) , args = { 'nbchan', str2num(results{4}) args{:} }; end;
	if ~isempty( results{5} ) , args = { args{:}, 'pnts', str2num(results{5}) }; end;
	if ~isempty( results{6} ) , args = { args{:}, 'srate', str2num(results{6}) }; end;
    if ~isempty( results{7} ) , args = { args{:}, 'xmin', str2num(results{7}) }; end;
	i = 8;
	if ~isempty( results{i  } ) , args = { args{:}, 'chanlocs' , results{i} }; end;
	if ~isempty( results{i+1} ),  args = { args{:}, 'icaweights', results{i+1} }; end;
	if ~isempty( results{i+2} ) , args = { args{:}, 'icasphere', results{i+2} }; end;
else % no interactive inputs
    args = varargin;
    for index=1:2:length(args)
        if ~isempty(inputname(index+1)) & ~isstr(args{index+1}) & length(args{index+1})>1, 
			args{index+1} = inputname(index+1); 
		end;
    end;                
end;

% generate the output command
% ---------------------------
com = '';
for i=1:2:length(args)
    if ~isempty( args{i+1} )
        if isstr( args{i+1} ) com = sprintf('%s, ''%s'', ''%s''', com, args{i}, char(args{i+1}) );
        else                  com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
        end;
    else
        com = sprintf('%s, ''%s'', []', com, args{i} );
    end;
end;

eval( [ 'EEGOUT = pop_editset(EEGOUT' com ');'] );
com = [ 'EEG = pop_importdata(' com(2:end) ');'];
return;

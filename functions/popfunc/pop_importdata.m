% pop_importdata() - Import data from a Matlab variable or from a file. 
%
% Usage:
%   >> EEGOUT = pop_importdata( 'key', val,...);
%
% Optional inputs:
%   'setname'    - string with the name of the dataset
%   'data'       - ['varname'|'filename'] import a file or variable
%                  into the EEG structure of EEGLAB.
%   'dataformat' - ['array|matlab|ascii|float32'] format of the input data file. The
%                  data file is transposed if the number of row is greater than
%                  the number of columns. Note that for 'float32', data must
%                  be organised in the format channels x ndata.
%   'chanlocs'   - ['varname'|'filename'] import a file containing
%                  electrodes locations (see >> help readlocs for file format).
%   'nbchan'     - number of channel in data
%   'xmin'       - starting time in second
%   'pnts'       - number of point per frame in the data (for epoched data only)
%   'srate'      - data sampling rate
%   'icaweight'  - ica weight matrix. By default, the sphering matrix is set to
%                  the identity matrix if it is empty.
%   'icasphere'  - ica sphering matrix
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
    geometry    = { [2 0.1 1 0.7] [1.1 1 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] [1.5 0.6 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] };
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    uilist = { ...
         { 'Style', 'text', 'string', 'Dataset name (optional):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', '' }, { }...
         ...
         { 'Style', 'text', 'string', 'Data file or array', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { 'Style', 'listbox', 'string', 'Matlab array|text ASCII file format|Float_32 file format|Matlab file format' } ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'globfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''globfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'Number of channel (0=set from data):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, {},  { 'Style', 'edit', 'string', '0' }, { } ...
         { 'Style', 'text', 'string', 'Time points per epoch (0=continuous data):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', num2str(EEGOUT.pnts) }, { } ...
         { 'Style', 'text', 'string', 'Data sampling rate (Hz):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', num2str(EEGOUT.srate) }, { },...
         { 'Style', 'text', 'string', 'Optional epoch start time for data epochs (sec):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { }, { 'Style', 'edit', 'string', num2str(EEGOUT.xmin) }, { },...
         ...
         { 'Style', 'text', 'string', 'Channel position file or array:', 'horizontalalignment', 'right', 'fontweight', 'bold' }, {'Style', 'pushbutton', 'string', 'Help', 'callback', 'pophelp(''readlocs'');' }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'chanfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''chanfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA weight array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'weightfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''weightfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'CA sphere array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
        { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'sphfile' } ...
        { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''sphfile'';' commandload ] } };

    results = inputgui( geometry, uilist, 'pophelp(''pop_importdata'');', 'Import dataset info -- pop_editset()');
    if length(results) == 0, return; end;

	args = {};
	if ~isempty( results{1} ), args = { args{:}, 'setname', results{1} }; end;
	switch results{2}
	   case 1, args = { args{:}, 'dataformat', 'array' };
	   case 2, args = { args{:}, 'dataformat', 'ascii' };
	   case 3, args = { args{:}, 'dataformat', 'float32' };
	   case 4, args = { args{:}, 'dataformat', 'matlab' };
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
        if ~isempty(inputname(index+1)) & ~isstr(args{index+1}) , args{index+1} = { inputname(index+1) }; end;
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

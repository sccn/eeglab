% pop_editset() - Edit dataset info. 
%
% Usage:
%   >> EEGOUT = pop_editset( EEG, 'key', val,...);
%
% Inputs:
%   EEG             - dataset structure
%
% Optional inputs:
%   'setname'    - Name of the dataset
%   'data'       - ['varname'|'filename'] Import a file or variable
%                   into the EEG structure of EEGLAB.
%   'dataformat' - ['array|matlab|ascii|float32'] Format of the input data file. The
%                   data file is transposed if the number of rows is larger than
%                   the number of columns. Note that for 'float32', data must
%                   be organised in the format channels x times.
%   'chanlocs'   - ['varname'|'filename'] Import a file containing electrode locations 
%                  (See >> help readlocs for file format).
%   'nbchan'     - Number of channels in data
%   'xmin'       - Starting time (in seconds)
%   'pnts'       - Number of points per epoch in the data (for epoched data only)
%   'srate'      - Data sampling rate
%   'icaweight'  - ICA weight matrix. By default, the sphering matrix is set to
%                   the identity matrix if it is empty.
%   'icasphere'  - ICA sphering matrix
% 
% Outputs:
%   EEGOUT       - modified dataset structure
%
% Note:
%   To create a new dataset, type 
%   >> EEG = pop_editset( eeg_emptyset );
%   To erase a variable, use '[]'. The folowing command suppresses channel locations.
%   Ex: EEG = pop_editset( EEG, 'chanlocs', '[]');
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
    geometry    = { [2 0.1 1 0.7] [1.1 1 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] [1.5 0.6 1 0.7] [2 0.1 1 0.7] [2 0.1 1 0.7] };
    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    uilist = { ...
         { 'Style', 'text', 'string', 'Dataset name (optional):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', EEG.setname }, { 'Style', 'text', 'string', '(EEG.setname)' }...
         ...
         { 'Style', 'text', 'string', 'Data file or array', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { } ...
         { }, ...
         { 'Style', 'text', 'string', 'EEG.data' }, ...
         ...
         { 'Style', 'text', 'string', 'Channels in data:', 'horizontalalignment', 'right', 'fontweight', 'bold' }, {},  { 'Style', 'text', 'string', num2str(EEG.nbchan) }, { 'Style', 'text', 'string', '(EEG.nbchan)' } ...
         { 'Style', 'text', 'string', 'Time points per epoch (0=continuous data):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', num2str(EEG.pnts) }, {'Style', 'text', 'string', '(EEG.pnts)' } ...
         { 'Style', 'text', 'string', 'Data sampling rate (Hz):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', num2str(EEG.srate) }, {'Style', 'text', 'string', '(EEG.srate)' },...
         { 'Style', 'text', 'string', 'Epoch start time (sec):', 'horizontalalignment', 'right', 'fontweight', 'bold' }, { }, { 'Style', 'edit', 'string', num2str(EEG.xmin) }, {'Style', 'text', 'string', '(EEG.xmin)' },...
         ...
         { 'Style', 'text', 'string', 'Channel position file or array:', 'horizontalalignment', 'right', 'fontweight', 'bold' }, {'Style', 'pushbutton', 'string', 'Help', 'callback', 'pophelp(''readlocs'');' }, ...
         { 'Style', 'edit', 'string', fastif(isempty(EEG.chanlocs), '', 'EEG.chanlocs'), 'horizontalalignment', 'left', 'tag',  'chanfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''chanfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA weight array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
         { 'Style', 'edit', 'string', fastif(isempty(EEG.icaweights), '', 'EEG.icaweights'), 'horizontalalignment', 'left', 'tag',  'weightfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''weightfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'CA sphere array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
        { 'Style', 'edit', 'string', fastif(isempty(EEG.icasphere), '', 'EEG.icasphere'), 'horizontalalignment', 'left', 'tag',  'sphfile' } ...
        { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''sphfile'';' commandload ] } };

    if EEG.trials == 1,  uilist(21:24) = []; geometry(6) = []; end;
    results = inputgui( geometry, uilist, 'pophelp(''pop_editset'');', fastif(isempty(EEG.data), 'Import dataset info -- pop_editset()', 'Edit dataset info -- pop_editset()'));
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
else % no interactive inputs
    args = varargin;
    for index=1:2:length(args)
        if ~isempty(inputname(index+2)), args{index+1} = { inputname(index+2) }; end;
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
                            fprintf('Pop_editset: filename ''%s'' found for channel location\n', varname); 
                            EEGOUT.chanlocs = readlocs(varname);
                         else
                            EEGOUT.chanlocs = evalin('base', varname, 'fprintf(''Pop_editset warning: variable name ''''%s'''' not found, ignoring\n'', varname);' );
                         end;
        case 'icaweights', varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2
                            fprintf('Pop_editset: filename ''%s'' found for ICA weight matrix\n', varname); 
                            try, EEGOUT.icaweights = load(varname, '-ascii');
                            catch, fprintf('Pop_editset warning: erro while reading filename ''%s'' for ICA weight matrix\n', varname); 
                            end;
                         else
                            EEGOUT.icaweights = evalin('base', varname, 'fprintf(''Pop_editset warning: variable name ''''%s'''' not found, ignoring\n'', varname);' );
                         end;
                         if ~isempty(EEGOUT.icaweight) & isempty(EEGOUT.icasphere)
                            EEGOUT.icasphere = eye(size(EEGOUT.icaweights,2));
                         end;
        case 'icasphere', varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2
                            fprintf('Pop_editset: filename ''%s'' found for ICA sphere matrix\n', varname); 
                            try, EEGOUT.icasphere = load(varname, '-ascii');
                            catch, fprintf('Pop_editset warning: erro while reading filename ''%s'' for ICA weight matrix\n', varname); 
                            end;
                         else
                            EEGOUT.icasphere = evalin('base', varname, 'fprintf(''Pop_editset warning: variable name ''''%s'''' not found, ignoring\n'', varname);' );
                         end;
        case 'data'    , varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2 & ~strcmp(lower(g.dataformat), 'array');
                            fprintf('Pop_editset: filename ''%s'' found for raw data\n', varname); 
                            switch lower(g.dataformat)
                                case 'ascii' , try, EEGOUT.data = load(varname, '-ascii');
                                              catch, error(['Pop_editset error: can not read ascii file ''' varname ''' for raw data']); 
                                              end;
                                              if size(EEGOUT.data,1) > size(EEGOUT.data,2), EEGOUT.data = transpose(EEGOUT.data); end;
                                case 'matlab', try, EEGOUT.data = load(varname, '-mat');
                                              catch, error(['Pop_editset error: can not read matlab file ''' varname ''' for raw data']); 
                                              end;
                                              if size(EEGOUT.data,1) > size(EEGOUT.data,2), EEGOUT.data = transpose(EEGOUT.data); end;
                                case 'float32', if EEGOUT.nbchan == 0,
                                                    error(['Pop_editset error: can not read file float_32 without the number of channel']);
                                                end;     
                                                try, EEGOUT.data = floatread(varname, [EEGOUT.nbchan Inf]);
                                                catch, error(['Pop_editset error: can not read float_32 file ''' varname ''' for raw data']); 
                                                end;
                                otherwise, error('Pop_editset error: unrecognized file format');
                            end;
                         else
							% restoration command
							%--------------------
							try 
							     res = evalin('base', ['exist(''' varname ''') == 1']);
							catch
							     error('Pop_editset: cannot evaluate variable, that may be a filename ?');
							end;
							if ~res, error('Pop_editset: cannot evaluate variable, that may be a filename ?'); end;
						    testval = evalin('base', ['isglobal(' varname ')']);
							warning off;
							if ~testval
								commandrestore = [ ' tmpp = '  varname '; clear global ' varname ';'   varname '=tmpp;clear tmpp;' ]; 
							else
                                commandrestore = [];
							end;		  
							% make global, must make these variable global, if you try to evaluate them direclty in the base
							% workspace, with large array the computation become incredibly slow.  
							%--------------------------------------------------------------------
				         	comglobal = sprintf('global %s;', varname);
				      		evalin('base', comglobal);
				      		eval(comglobal);
				      		eval( ['EEGOUT.data = ' varname ';' ]);
				      		try, evalin('base', commandrestore); catch, end;
				      		warning on;
                            if size(EEGOUT.data,1) > size(EEGOUT.data,2), EEGOUT.data = transpose(EEGOUT.data); end;
                         end;
         otherwise, error(['Pop_editset error: unrecognized field ''' curfield{1} '''']); 
    end;
end;

% generate the output command
% ---------------------------
com = sprintf( '%s = pop_editset(%s', inputname(1), inputname(1) );
for i=1:2:length(args)
    if ~isempty( args{i+1} )
        if isstr( args{i+1} ) com = sprintf('%s, ''%s'', ''%s''', com, args{i}, args{i+1} );
        else                  com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
        end;
    else
        com = sprintf('%s, ''%s'', []', com, args{i} );
    end;
end;
com = [com ');'];
return;

function num = popask( text )
	 ButtonName=questdlg( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;

% pop_editoptions() - Edit memory-saving eeglab() options. These are stored in 
%                     a file 'eeg_options.m'. With no argument, pop up a window 
%                     to allow the user to set/unset these options. Store
%                     user choices in a new 'eeg_options.m' file in the 
%                     working directory.
%
% Usage: >> pop_editoptions;
%        >> pop_editoptions( 'key1', value1, 'key2', value2, ...);
%
% Graphic interface inputs:
%   "If set, keep at most one dataset in memory ..." - [checkbox] If set, EEGLAB will only retain the current
%                   dataset in memory. All other datasets will be automatically
%                   read and writen to disk. All EEGLAB functionalities are preserved
%                   even for dataset stored on disk. 
%   "If set, write data in same file as dataset ..." - [checkbox] Set -> dataset data (EEG.data) are 
%                   saved in the EEG structure in the standard Matlab dataset (.set) file. 
%                   Unset -> The EEG.data are saved as a transposed stream of 32-bit 
%                   floats in a separate binary file. As of Matlab 4.51, the order 
%                   of the data in the binary file is as in the transpose of EEG.data 
%                   (i.e., as in EEG.data', frames by channels). This allows quick 
%                   reading of single channels from the data, e.g. when comparing 
%                   channels across datasets. The stored files have the extension 
%                   .dat instead of the pre-4.51, non-transposed .fdt. Both file types 
%                   are read by the dataset load function. Command line equivalent is
%                   option_savematlab.
%   "Precompute ICA activations" - [checkbox] If set, all the ICA activation
%                   time courses are precomputed (this requires more RAM). 
%                   Command line equivalent: option_computeica.
%   "If set, remember old folder when reading dataset" - [checkbox] this option
%                   is convinient if the file you are working on are not in the 
%                   current folder.
%
% Commandline keywords:
%   'option_computeica' - [0|1] If 1, compute the ICA component activitations and
%                   store them in a new variable. If 0, compute ICA activations
%                   only when needed (& only partially, if possible) and do not
%                   store the results).
%   NOTE: Turn OFF the options above when working with very large datasets or on 
%                   computers with limited memory.
%   'option_savematlab' - [0|1] If 1, datasets are saved as single Matlab .set files. 
%                   If 0, dataset data are saved in separate 32-bit binary float 
%                   .dat files.  See the corresponding GUI option above for details. 
% Outputs:
%   In the output workspace, variables 'option_computeica', 
%   and 'option_savematlab'  are updated, and a new 'eeg_options.m' file may be
%   written to the working directory. The copy of 'eeg_options.m' placed in your 
%   working directory overwrites system defaults whenever EEGLAB operates in this
%   directory (assuming your working directory is in your MATLABPATH - see path()).
%   To adjust these options system-wide, edit the master "eeg_options.m" file in the
%   EEGLAB directory heirarchy.
%
% Author: Arnaud Delorme, SCCN / INC / UCSD, March 2002
%
% See also: eeg_options(), eeg_readoptions()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 09 March 2002, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function com = pop_editoptions(varargin)

com = '';
argsoutput = {};

datasets_in_memory = 0;
if nargin > 0
    if ~ischar(varargin{1})
        datasets_in_memory = varargin{1};
        varargin = {};
    end
end

% parse the eeg_options file
% ----------------------------
eeglab_options;
if isdeployed || (exist('ismcc') && ismcc)
    filename = which('eeg_options.txt');
    eegoptionbackup = which('eeg_optionsbackup.txt');
else
    % folder for eeg_options file (also update the eeglab_options)
    if ~isempty(EEGOPTION_PATH)
         homefolder = EEGOPTION_PATH;
    elseif ispc
         if ~exist('evalc'), eval('evalc = @(x)(eval(x));'); end
         homefolder = deblank(evalc('!echo %USERPROFILE%'));
    else homefolder = '~';
    end
    filename = fullfile(homefolder, 'eeg_options.m');
    eegoptionbackup = which('eeg_optionsbackup.m');
end

fid = fopen( filename, 'r+'); % existing file
storelocal = 0;
if	fid == -1
    filepath = homefolder;
    filename = 'eeg_options.m';
    fid = fopen( fullfile(filepath, filename), 'w'); % new file possible?
    if fid == -1
        error([ 'Cannot write into HOME folder: ' homefolder 10 'You may specify another folder for the eeg_option.m' 10 'file by editing the icadefs.m file' ]);
    end
    fclose(fid);
    delete(fullfile(filepath, filename));

    % read variables values and description
    % --------------------------------------
    [ header, opt ] = eeg_readoptions( eegoptionbackup ); 
else 
    [filepath, filename, ext] = fileparts(filename);
    filename  = [ filename ext ];
    fprintf('Using option file in directory %s\n', filepath);
    
    % read variables values and description
    % --------------------------------------
    [ header, opt ] = eeg_readoptions( eegoptionbackup ); 
    [ header, opt ] = eeg_readoptions( fid, opt  ); % use opt from above as default
end

optionsToShow = {
    'option_storedisk' ...
    'option_savetwofiles'  ...
    'option_computeica'  ...
    'option_rememberfolder' ...
    'option_allmenus'  ...
    'option_checkversion' ...
    'option_showadvanced' ...
    'option_cachesize' };

% remove advanced options if necessary
if isempty(varargin)
    if ~option_showadvanced
        % remove options 
        for iOpt = length(opt):-1:1
            if ~isempty(opt(iOpt).varname) && ~ismember(opt(iOpt).varname, optionsToShow)
                opt(iOpt) = [];
            end
        end
        % remove header not serving any option
        for iOpt = length(opt)-1:-1:1
            if isempty(opt(iOpt).varname) && isempty(opt(iOpt+1).varname)
                opt(iOpt) = [];
            end
        end
    end
end

if nargin < 2
    geometry = { [6 1] };
    tmpfile = fullfile(filepath, filename);
    
    cb_file = [ '[filename, filepath] = uiputfile(''eeg_options.txt'', ''Pick a folder to save option file'');' ...
                'if filename(1) ~= 0,' ...
                '   filepath = fullfile(filepath, ''eeg_options.m'');' ...
                '   set(gcf, ''userdata'', filepath);' ...
                '   if length(filepath) > 100,' ...
                '        filepath =  [ ''...'' filepath(end-100:end) ];' ...
                '   end;' ...
                '   set(findobj(gcf, ''tag'', ''filename''), ''string'', filepath);' ...
                'end;' ...
                'clear filepath;' ];
            
    uilist = { ...
         { 'Style', 'text', 'string', '', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Set/Unset', 'fontweight', 'bold'   } };

    % add all fields to graphic interface
    % -----------------------------------
    for index = 1:length(opt)
        % format the description to fit a help box
        % ----------------------------------------
        descrip = { 'string', opt(index).description }; % strmultiline(description{ index }, 80, 10) };
           
        % create the gui for this variable
        % --------------------------------
        if strcmpi(opt(index).varname, 'option_storedisk') && datasets_in_memory
            cb_nomodif = [ 'set(gcbo, ''value'', ~get(gcbo, ''value''));' ...
                           'warndlg2(strvcat(''This option may only be modified when at most one dataset is stored in memory.''));' ];
            
        elseif strcmpi(opt(index).varname, 'option_memmapdata')
            cb_nomodif = [ 'if get(gcbo, ''value''), warndlg2(strvcat(''Matlab memory is beta, use at your own risk'')); end;' ];
        elseif strcmpi(opt(index).varname, 'option_donotusetoolboxes')
            cb_nomodif = [ 'if get(gcbo, ''value''), warndlg2([''You have selected the option to disable'' 10 ''Matlab toolboxes. Use with caution.'' 10 ''Matlab toolboxes will be removed from'' 10 ''your path. Unlicking this option later will not'' 10 ''add back the toolboxes. You will need'' 10 ''to add them back manually. If you are unsure'' 10 ''if you want to disable Matlab toolboxes'' 10 ''deselect the option now.'' ]); end;' ];
        else
            cb_nomodif = '';
        end
        
        if ~isempty(opt(index).value)
            if opt(index).value <= 1
                uilist   = { uilist{:}, { 'Style', 'text', descrip{:}, 'horizontalalignment', 'left' }, ...
                             { 'Style', 'checkbox', 'string', '    ', 'value', opt(index).value 'callback' cb_nomodif } { } }; 
                geometry = { geometry{:} [4 0.3 0.1] };
            else
                uilist   = { uilist{:}, { 'Style', 'text', descrip{:}, 'horizontalalignment', 'left' }, ...
                             { 'Style', 'edit', 'string', num2str(opt(index).value), 'callback' cb_nomodif } { } }; 
                geometry = { geometry{:} [3 0.5 0.1] };
            end
        else
            uilist   = { uilist{:}, { 'Style', 'text', descrip{:}, 'fontweight' 'bold', 'horizontalalignment', 'left' }, { } { } }; 
            geometry = { geometry{:} [4 0.3 0.1] };
        end
    end

    % change option file
    uilist = { uilist{:} {} ...
                 { 'Style', 'text', 'string', 'Edit the EEGOPTION_PATH variable of functions/sigprocfunc/icadefs.m to change where the option file is saved' } };
    geometry = { geometry{:} [1] [1] };
    [results, userdat ] = inputgui( geometry, uilist, 'pophelp(''pop_editoptions'');', 'Memory options - pop_editoptions()', ...
                        [], 'normal');
    if isempty(results), return; end
   
    % decode inputs
    % -------------
    args = {};
    count = 1;
    for index = 1:length(opt)
        if ~isempty(opt(index).varname)
            args = {  args{:}, opt(index).varname, results{count} }; 
            count = count+1;
        end
    end
else 
    % no interactive inputs
    % ---------------------
    args = varargin;
end

% change default folder option
% ----------------------------
W_MAIN = findobj('tag', 'EEGLAB');
if ~isempty(W_MAIN)
    tmpuserdata    = get(W_MAIN, 'userdata');
    tmpuserdata{3} = filepath;
    set(W_MAIN, 'userdata', tmpuserdata);
end

% decode inputs
% -------------
for index = 1:2:length(args)
    ind = strmatch(args{index}, { opt.varname }, 'exact');
    if isempty(ind)
        if strcmpi(args{index}, 'option_savematlab')
            disp('pop_editoptions: option_savematlab is obsolete, use option_savetwofiles instead');
            ind = strmatch('option_savetwofiles', { opt.varname }, 'exact');
        else
            error(['Variable name ''' args{index} ''' is invalid']);
        end
    end
    
    % case for 'option_cachesize'
    if strcmpi(args{index}, 'option_cachesize') && ischar(args{index+1})
        args{index+1}  = str2num(args{index+1});
    end
    
    % overwrite only if different
    if args{index+1} ~= opt(ind).value 
        opt(ind).value    = args{index+1};
        argsoutput{end+1} = args{index};   % for history
        argsoutput{end+1} = args{index+1}; % for history
    end  
end

% write to eeg_options file
% -------------------------
fid = fopen( fullfile(filepath, filename), 'w');
addpath(filepath);
if fid == -1
	error('File writing error, check writing permission');
end
fprintf(fid, '%s\n', header);
for index = 1:length(opt)
    if isempty(opt(index).varname)
        fprintf( fid, '%% %s\n', opt(index).description);
    else
        fprintf( fid, '%s = %d ; %% %s\n', opt(index).varname, opt(index).value, opt(index).description);
    end
end
fclose(fid);    
% clear it from the MATLAB function cache
clear(fullfile(filepath,filename));

% generate the output text command
% --------------------------------
if ~isempty(argsoutput)
    com = 'pop_editoptions(';
    for index = 1:2:length(argsoutput)
        com = sprintf( '%s ''%s'', %d,', com, argsoutput{index}, argsoutput{index+1});
    end
    com = [com(1:end-1) ');'];
else
    disp('pop_editoptions: Options were not modified');
end
wtmp = warning;
warning off;
clear functions
warning(wtmp);

% ---------------------------
function  chopedtext = choptext( tmptext )
    chopedtext = '';
    while length(tmptext) > 30
          blanks = findstr( tmptext, ' ');
          [tmp, I] = min( abs(blanks - 30) );
          chopedtext = [ chopedtext ''' 10 ''' tmptext(1:blanks(I)) ];
          tmptext  = tmptext(blanks(I)+1:end);
    end    
    chopedtext = [ chopedtext ''' 10 ''' tmptext];
    chopedtext = chopedtext(7:end);
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName)
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end

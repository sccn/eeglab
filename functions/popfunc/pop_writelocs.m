% pop_writelocs() - load a EGI EEG file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_writelocs(chanstruct);             % a window pops up
%   >> EEG = pop_writelocs(chanstruct, filename, 'key', val, ...);
%
% Inputs:
%   chanstruct     - channel structure. See readlocs()
%   filename       - Electrode location file name
%   'key',val      - same as writelocs()
% 
% Author: Arnaud Delorme, CNL / Salk Institute, 17 Dec 2002
%
% See also: writelocs()

% Copyright (C) Arnaud Delorme, Salk Institute, arno@salk.edu
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

function com = pop_writelocs(chans, filename, varargin); 
    
com = '';
if nargin < 1
   help pop_writelocs;
   return;
end

if isfield(chans, 'shrink')
    chans = rmfield(chans, 'shrink');
    disp('Warning: shrink factor ignored');
end

disp('WARNING: ELECTRODE COORDINATES MUST BE WITH NOSE ALONG THE +X DIMENSION TO BE EXPORTED')
disp('         IF NOT, THE EXPORTED FILE COORDINATES MAY BE INACURATE')

% get infos from readlocs
% -----------------------
[chanformat listcolformat] = readlocs('getinfos');
chanformat(end)    = [];
listcolformat(end) = []; % remove chanedit
chanformat(end)    = [];
listcolformat(end) = []; % remove chanedit
indformat  = [];
for index = 1:length(chanformat), 
    if ~ischar(chanformat(index).importformat)
        indformat = [ indformat index ];
    end
    if isempty(chanformat(index).skipline), chanformat(index).skipline = 0; end
end
listtype   = { chanformat(indformat).type };
formatinfo = { chanformat(indformat).importformat };
formatskip = [ chanformat(indformat).skipline ];
   
%[listtype formatinfo listcolformat formatskip] = readlocs('getinfoswrite');

% GUI support of `custom` filetype, removed until `custom` can pass readlocs() check
%listtype{end+1} = 'custom';
%formatinfo{end+1} = {};
%formatskip = [ formatskip 0];


if nargin < 2
   updatefields = [ 'tmpdata = get(gcf, ''userdata'');' ...
                    'tmpobj = findobj(gcf, ''tag'', ''list2'');' ...
                    'set(tmpobj, ''string'', strvcat(tmpdata{2}));' ...
                    'clear tmpobj tmpdata;' ];
   addfieldcom = [ 'tmpdata = get(gcbf, ''userdata'');' ...
                   'tmpobj = findobj(gcf, ''tag'', ''list1'');' ...
                   'tmpdata{2}{end+1} = tmpdata{1}{get(tmpobj, ''value'')};' ...
                   'set(gcbf, ''userdata'', tmpdata);' ...
                   updatefields ];
   rmfieldcom  = [ 'tmpdata = get(gcbf, ''userdata'');' ...
                   'tmpobj = findobj(gcbf, ''tag'', ''list2'');' ...
                   'try, tmpdata{2}(get(tmpobj, ''value'')) = [];' ...
                   '    set(tmpobj, ''value'', 1);' ...
                   'catch, end;' ...
                   'set(gcbf, ''userdata'', tmpdata);' ...
                   updatefields ];                      
   filetypecom = [ 'tmpdata = get(gcf, ''userdata'');' ...
                   'tmpobj = findobj(gcf, ''tag'', ''formatlist'');' ...
                   'tmpval = get(tmpobj, ''value'');' ...
                   'try, tmpdata{2} = tmpdata{3}{tmpval}; catch, end;' ... %try and catch for custom
                   'set(gcf, ''userdata'', tmpdata);' ...
                   updatefields ...
                   'tmpdata = get(gcf, ''userdata'');' ...
						 'tmpobj1 = findobj(gcf, ''tag'', ''insertcol'');' ... % the lines below
                   'tmpobj2 = findobj(gcf, ''tag'', ''inserttext'');' ... % update the checkbox
                   'try, ' ...                                     % and the edit text box
                   '  if tmpdata{4}(tmpval) == 2,' ...
                   '     set(tmpobj1, ''value'', 1);' ...
                   '  else,' ...
                   '     set(tmpobj1, ''value'', 0);' ...
                   '  end;' ...
                   '  if tmpval == 1,' ... % besa only
                   '     set(tmpobj2, ''string'', ''' int2str(length(chans)) ''');' ...
                   '  else,' ...
                   '     set(tmpobj2, ''string'', '''');' ...
                   '  end;' ...
                   'catch, end;' ... % catch for custom case
                   'tmpobj = findobj(gcf, ''userdata'', ''setfield'');' ...
                   'if tmpval == ' int2str(length(listtype)) ',' ... % disable if non-custom type
                   '   set(tmpobj, ''enable'', ''on'');' ...
                   'else,' ...
                   '   set(tmpobj, ''enable'', ''off'');' ...
                   'end; clear tmpobj tmpobj2 tmpdata tmpval;' ];
                
   geometry = { [1 1 1] [1 1] [1] [1 1 1] [1 1] [1 1 1] [1 0.3 0.7] [1] [1] };
   listui = { ...
         { 'style' 'text' 'string' 'Filename' } ...
         { 'style' 'edit' 'string' '' 'tag' 'filename' 'horizontalalignment' 'left' } ...
         { 'style' 'pushbutton' 'string' 'Browse' 'callback' ...
           [  '[tmpfile tmppath] = uiputfile(''*'', ''Exporting electrode location file -- pop_writelocs()'');' ... 
              'set(findobj(gcbf, ''tag'', ''filename''), ''string'', char([tmppath tmpfile ]));' ...
              'clear tmpfile tmppath;' ] } ...
         { 'style' 'text' 'string' strvcat('Select output file type', ' ', ' ') } ...
         { 'style' 'listbox' 'tag' 'formatlist' 'string' strvcat(listtype) ...
            'value' length(listtype) 'callback' filetypecom } ...
         { 'style' 'text' 'string' 'Select fields to export below' } ...
         { } { 'style' 'pushbutton' 'string' '-> ADD' 'callback' addfieldcom 'userdata' 'setfield' } { } ...
         { 'style' 'listbox' 'tag' 'list1' 'string' strvcat(fieldnames(chans)) 'userdata' 'setfield' } ...
         { 'style' 'listbox' 'tag' 'list2' 'string' '' 'userdata' 'setfield2' } ...
         { } { 'style' 'pushbutton' 'string' 'REMOVE <-' 'callback' rmfieldcom 'userdata' 'setfield' } { } ...
         { 'style' 'text' 'string' 'Insert column names' } ...
         { 'style' 'checkbox' 'tag' 'insertcol' 'value' 1 'userdata' 'setfield' } { } ...
         { 'style' 'text' 'string' 'Enter custom header below' } ...
         { 'style' 'edit' 'userdata' 'setfield' 'tag' 'inserttext' 'horizontalalignment' 'left' 'max' 2 } ...
	};
   
   inputgui(geometry, listui, 'pophelp(''writelocs'');', ...
      'Exporting electrode location file -- pop_writelocs()', { fieldnames(chans) {} formatinfo formatskip }, 'plot', [1 3 1 1 3 1 1 1 3 ]);
   fig = gcf;
   
   % set default format
   tmpobj = findobj(fig, 'tag', 'formatlist'); 
   set(tmpobj, 'value', 6);
   eval(get(tmpobj, 'callback')); 
   
   res = inputgui(geometry, listui, 'pophelp(''writelocs'');', ...
      'Exporting electrode location file -- pop_writelocs()', { listcolformat {} formatinfo formatskip }, fig, [1 3 1 1 3 1 1 1 3 ]);
   if gcf ~= fig, return; end
   exportfields = get(fig, 'userdata');
   exportfields = exportfields{2};
   close(fig);
   
   % decode the inputs
   filename = res{1};
   if isempty(filename), 
      errordlg2('Error: Empty file name', 'Error');
      return;
   end
   options = { 'filetype' listtype{res{2}} 'format' exportfields ...
         'header' fastif(res{5}, 'on', 'off') 'customheader' res{6} };
else
	options = varargin;   
end

% generate history
% ----------------
if isempty(inputname(1)) % not a variable name -> probably the structure from pop_chanedit
    writelocs(chans, filename, options{:});
   com = sprintf('pop_writelocs( EEG.chanlocs, ''%s'', %s);', filename, vararg2str(options));
else
    if strcmpi(inputname(1), 'chantmp')
        % do not write file (yet)
        com = sprintf('pop_writelocs( chans, ''%s'', %s);', filename, vararg2str(options));
    else
        writelocs(chans, filename, options{:});
        com = sprintf('pop_writelocs( %s, ''%s'', %s);', inputname(1), filename, vararg2str(options));
    end
end

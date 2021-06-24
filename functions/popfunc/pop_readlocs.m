% pop_readlocs() - load channel location file (pop up an interactive window 
%                  if no arguments).
%
% Usage:
%   >> chanlocs = pop_readlocs;                             % a window pops up
%   >> chanlocs = pop_readlocs( filename, 'key', val, ...); % no window
%
% Inputs:
%   filename       - Electrode location file name
%   'key',val      - Same options as readlocs() (see >> help readlocs)
% 
% Outputs: same as readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 17 Dec 2002
%
% See also: readlocs()

% Copyright (C) 17 Dec 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [tmplocs, command] = pop_readlocs(filename, varargin); 
    
tmplocs = [];
command = '';

% get infos from readlocs
% -----------------------
[chanformat listcolformat] = readlocs('getinfos');
listtype = { chanformat.type };

if nargin < 1 
   [filename, filepath] = uigetfile('*.*', 'Importing electrode location file -- pop_readlocs()'); 
   %filename = 'chan32.locs';
   %filepath = 'c:\eeglab\eeglab4.0\';
   drawnow;
   if filename == 0 return; end
   filename = [filepath filename];
   tmpfile = loadtxt(filename);
   nbcols = cellfun('isempty', tmpfile(end,:));
   nbcols = ones(1,length(find(~nbcols)));
   
   % decoding file type
   % ------------------
   periods = find(filename == '.');
   fileextension = filename(periods(end)+1:end);
   switch lower(fileextension),
		 case {'loc' 'locs' }, filetype = 'loc';
		 case 'xyz', filetype = 'xyz';
		 case 'sph', filetype = 'sph';
		 case 'txt', filetype = 'chanedit';
		 case 'elp', filetype = 'polhemus';
		 case 'eps', filetype = 'besa';
		 otherwise, filetype =  ''; 
   end
   indexfiletype = strmatch(filetype, listtype, 'exact'); 
   
   % convert format info
   % -------------------
   formatinfo     = { chanformat.importformat };
   formatskipcell = { chanformat.skipline };
   rmindex    = [];
   count      = 1;
   for index = 1:length(formatinfo)
       if ~ischar(formatinfo{index})
           for index2 = 1:length(formatinfo{index})
               indexformat = strmatch(formatinfo{index}{index2}, listcolformat, 'exact');
               indexlist(count, index2) = indexformat;
           end
           if isempty(formatskipcell{index}), formatskip(count) = 0; 
           else                               formatskip(count) = formatskipcell{index}; 
           end
           count = count+1;
       else
           rmindex = [ rmindex index ];
       end
   end
   listtype(rmindex) = [];
   listtype  (end+1) = { 'custom' };
   formatskip(end+1) = 0;
   indexlist(end+1,:) = -1;
   
   % ask user
   formatcom = [  ...
         'indexformat = get(findobj(gcf, ''tag'', ''format''), ''value'');' ...
         'tmpformat = get(findobj(gcf, ''tag'', ''format''), ''userdata'');' ...
         'for tmpindex=1:' int2str(length(nbcols)) ...
         '   tmpobj = findobj(gcf, ''tag'', [ ''col'' int2str(tmpindex) ]);' ...
         '   if tmpformat{1}(indexformat,tmpindex) == -1,' ...
         '       set(tmpobj, ''enable'', ''on'');' ...
         '   else ' ...
         '       set(tmpobj, ''enable'', ''off'');' ...
         '       set(tmpobj, ''value'', tmpformat{1}(indexformat,tmpindex));' ...
         '   end;' ...
         'end;' ...
         'strheaderline = fastif(tmpformat{2}(indexformat)<0, ''auto'', int2str(tmpformat{2}(indexformat)));' ...
         'set(findobj(gcf, ''tag'', ''headlines''), ''string'', strheaderline);' ...
         'eval(get(findobj(gcf, ''tag'', ''headlines''), ''callback''));' ...
         'clear tmpindex indexformat tmpformat tmpobj strheaderline;' ];
         
   headercom = [ ...
         'tmpheader = str2num(get(findobj(gcf, ''tag'', ''headlines''), ''string''));' ...
         'if ~isempty(tmpheader), ' ...
        '   for tmpindex=1:' int2str(length(nbcols)) ...
         '     tmpobj = findobj(gcf, ''tag'', [ ''text'' int2str(tmpindex) ]);' ...
         '     tmpstr = get(tmpobj, ''userdata'');' ...
         '     set(tmpobj, ''string'', tmpstr(tmpheader+1:min(tmpheader+10, size(tmpstr,1)),:));' ...
         '  end;' ...
         'end;' ...
         'clear tmpobj tmpstr tmpheader tmpindex;' ];
            
   geometry = { [1 1 1] [1 2] [nbcols] [nbcols] [nbcols] [1] [1 1 1]};
   listui = { ...
         { 'style' 'text' 'string' 'Predefined format' } ...
         { 'style' 'text' 'string' 'Header lines' } ...
         { 'style' 'edit' 'tag' 'headlines' 'string' '0' 'callback' headercom } ...
         { 'style' 'listbox' 'tag' 'format' 'string' strvcat(listtype{:}) ...
            'callback' formatcom 'userdata' { indexlist;formatskip } } ...
         { 'style' 'pushbutton' 'string' 'preview' ...
            'callback' 'set(findobj(gcbf, ''tag'', ''Import''), ''userdata'', ''preview'')'} };
   
   % custom columns
   % --------------
   for index = 1:length(nbcols)
      listui{end+1} = { 'style' 'text' 'string' [ 'column ' int2str(index) ] };
   end
   for index = 1:length(nbcols)
      listui{end+1} = { 'style' 'listbox' 'tag' [ 'col' int2str(index) ] 'string' ...
            strvcat(listcolformat{:}) };
   end
   for index = 1:length(nbcols)
      listui{end+1} = { 'style' 'text' 'string' formatstr(tmpfile(1:min(10, size(tmpfile,1)),index)) ...
            'tag' ['text' int2str(index) ] 'userdata' formatstr(tmpfile(:,index)) };      
   end
   listui = { listui{:} ...
         {} ...
         { 'style' 'pushbutton' 'string' 'Cancel', 'callback', 'close gcbf;' } ...
         { 'style' 'pushbutton' 'string' 'Help', 'callback', 'pophelp(''pop_readlocs'');' } ...
         { 'style' 'pushbutton' 'string' 'Import' 'tag' 'Import' ...
            'callback', 'set(findobj(gcbf, ''tag'', ''Import''), ''userdata'', ''import'')'} };
   
   fig = figure('name', 'Importing electrode location file -- pop_readlocs()', 'visible', 'off');
   supergui( fig, geometry, [1 2 1 3 7 1 1], listui{:});
   
   % update figure
   set(findobj(gcf, 'tag', 'format'), 'value', indexfiletype);
   eval(formatcom);
   
   % this loop is necessary for decoding Preview
   cont = 1;
   while cont
      waitfor( findobj('parent', fig, 'tag', 'Import'), 'userdata');
      if isempty( findobj('parent', fig, 'tag', 'Import') ), return; end
      
   	% decode inputs
   	% -------------
   	tmpheader = str2num(get(findobj(fig, 'tag', 'headlines'), 'string'));
   	tmpobj = findobj(fig, 'tag', 'format');
   	tmpformatstr = get(tmpobj, 'string');
   	tmpformatstr = tmpformatstr(get(tmpobj, 'value'),:);
   	for tmpindex=1:length(nbcols)
      	tmpobj = findobj(fig, 'tag', [ 'col' int2str(tmpindex) ]);
      	tmpcolstr{tmpindex} = get(tmpobj, 'string');
      	tmpcolstr{tmpindex} = deblank(tmpcolstr{tmpindex}(get(tmpobj,'value'),:));
   	end
   
   	% take appropriate measures
   	% -------------------------
   	res = get(findobj('parent', fig, 'tag', 'Import'), 'userdata');
   	set(findobj('parent', fig, 'tag', 'Import'), 'userdata', '');
      if strcmp(res, 'preview')
         try, 
            tmplocs = readlocs( filename, 'filetype', tmpformatstr, ...
             'skiplines', tmpheader, 'format', tmpcolstr);
         	figure; topoplot([],tmplocs, 'style', 'blank', 'electrodes', 'labelpoint');
         catch,
            errordlg2(strvcat('Error while importing locations:', lasterr), 'Error');
         end
      else 
      	cont = 0;   
      end
   end
   close(fig);
   
   % importing files
   % ---------------
   tmplocs = readlocs( filename, 'filetype', tmpformatstr, ...
             'skiplines', tmpheader, 'format', tmpcolstr);
   
else 
	tmplocs = readlocs( filename, varargin{:});	   
end

command = sprintf('EEG.chanlocs = pop_readlocs(''%s'', %s);', filename, vararg2str(varargin)); 
return;

% format string for text file
% ---------------------------
function alltxt = formatstr( celltxt )
for index1 = 1:size(celltxt,1)
   alltxt{index1} = '';
   for index2 = 1:size(celltxt,2)
      alltxt{index1} = sprintf('%s\t%s', alltxt{index1}, num2str(celltxt{index1,index2}));
   end
   alltxt{index1} = sprintf(alltxt{index1}, '%s\n', alltxt{index1});
end
alltxt = strvcat( alltxt{:});

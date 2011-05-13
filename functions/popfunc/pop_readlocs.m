% pop_readlocs() - load a EGI-format EEG file (pop up an interactive window if no arguments).
%
% Usage:
%   >> EEG = pop_readlocs;                             % a window pops up
%   >> EEG = pop_readlocs( filename, 'key', val, ...); % no window
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
   if filename == 0 return; end;
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
   end;
   indexfiletype = strmatch(filetype, listtype, 'exact'); 
   
   % convert format info
   % -------------------
   formatinfo     = { chanformat.importformat };
   formatskipcell = { chanformat.skipline };
   rmindex    = [];
   count      = 1;
   for index = 1:length(formatinfo)
       if ~isstr(formatinfo{index})
           for index2 = 1:length(formatinfo{index})
               indexformat = strmatch(formatinfo{index}{index2}, listcolformat, 'exact');
               indexlist(count, index2) = indexformat;
           end;
           if isempty(formatskipcell{index}), formatskip(count) = 0; 
           else                               formatskip(count) = formatskipcell{index}; 
           end;
           count = count+1;
       else
           rmindex = [ rmindex index ];
       end;
   end;
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
   end;
   for index = 1:length(nbcols)
      listui{end+1} = { 'style' 'listbox' 'tag' [ 'col' int2str(index) ] 'string' ...
            strvcat(listcolformat{:}) };
   end;
   for index = 1:length(nbcols)
      listui{end+1} = { 'style' 'text' 'string' formatstr(tmpfile(1:min(10, size(tmpfile,1)),index)) ...
            'tag' ['text' int2str(index) ] 'userdata' formatstr(tmpfile(:,index)) };      
   end;
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
      if isempty( findobj('parent', fig, 'tag', 'Import') ), return; end;
      
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
   	end;
   
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
         end;
      else 
      	cont = 0;   
      end;
   end;
   close(fig);
   
   % importing files
   % ---------------
   tmplocs = readlocs( filename, 'filetype', tmpformatstr, ...
             'skiplines', tmpheader, 'format', tmpcolstr);
   
else 
	tmplocs = readlocs( filename, varargin{:});	   
end;

command = sprintf('EEG.chanlocs = pop_readlocs(''%s'', %s);', filename, vararg2str(varargin)); 
return;

% format string for text file
% ---------------------------
function alltxt = formatstr( celltxt )
for index1 = 1:size(celltxt,1)
   alltxt{index1} = '';
   for index2 = 1:size(celltxt,2)
      alltxt{index1} = sprintf('%s\t%s', alltxt{index1}, num2str(celltxt{index1,index2}));
   end;
   alltxt{index1} = sprintf(alltxt{index1}, '%s\n', alltxt{index1});
end;
alltxt = strvcat( alltxt{:});

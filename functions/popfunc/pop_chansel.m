% pop_chansel() - select channel graphical interface
%
% Usage:
%   >> chanlist = pop_chansel(chanstruct); % a window pops up
%
% Inputs:
%   chanstruct     - channel structure. See readlocs()
% 
% Author: Arnaud Delorme, CNL / Salk Institute, 3 March 2003

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 3 March 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.3  2003/03/05 18:33:27  arno
% handling cancel
%
% Revision 1.2  2003/03/04 15:06:27  roberto
% no change
%
% Revision 1.1  2003/03/03 19:32:31  arno
% Initial revision
%

function chanlist = pop_chansel(chans); 
    
    if nargin < 1
        help pop_chansel;
        return;
    end;
    
    % get infos from readlocs
    % -----------------------
    updatefields = [ 'tmpdata = get(gcf, ''userdata'');' ...
                     'tmpobj = findobj(gcf, ''tag'', ''list2'');' ...
                     'set(tmpobj, ''string'', strvcat(tmpdata{2}));' ...
                     'clear tmpobj tmpdata;' ];
    addfieldcom = [ 'tmpdata = get(gcbf, ''userdata'');' ...
                    'tmpobj = findobj(gcf, ''tag'', ''list1'');' ...
                    'if strmatch(  tmpdata{1}{get(tmpobj, ''value'')}, tmpdata{2} ),' ...
                    '   clear tmpobj tmpdata; return;' ...
                    'end;' ...
                    'tmpdata{2}{end+1} = tmpdata{1}{get(tmpobj, ''value'')};' ...
                    'set(gcbf, ''userdata'', tmpdata);' ...
                    updatefields ];
    rmfieldcom  = [ 'tmpdata = get(gcbf, ''userdata'');' ...
                    'tmpobj = findobj(gcbf, ''tag'', ''list2'');' ...
                    'if get(tmpobj, ''value'') == length(tmpdata{2}) &' ...
                    '   length(tmpdata{2}) > 1,' ...
                    '   set(tmpobj, ''value'', get(tmpobj, ''value'')-1);' ...
                    'end;' ...
                    'try, tmpdata{2}(get(tmpobj, ''value'')) = [];' ...
                    'catch, end;' ...
                    'set(gcbf, ''userdata'', tmpdata);' ...
                    updatefields ];                      

   channelnames = strvcat({chans.labels});
   geometry = { [1] [1 4 1] [1 1] [1 4 1] };
   listui = { ...
         { 'style' 'text' 'string' 'Select channels to remove' } ...
         { } { 'style' 'pushbutton' 'string' '-> ADD TO LIST' 'callback' addfieldcom 'userdata' 'setfield' } { } ...
         { 'style' 'listbox' 'tag' 'list1' 'string' channelnames 'userdata' 'setfield' } ...
         { 'style' 'listbox' 'tag' 'list2' 'string' '' 'userdata' 'setfield' } ...
         { } { 'style' 'pushbutton' 'string' 'REMOVE FROM LIST <-' 'callback' rmfieldcom 'userdata' 'setfield' } { } ...
	};
   
   [outparams userdat] = inputgui(geometry, listui, 'pophelp(''pop_chansel'');', ...
                                  'Select channels -- pop_chansel()', { {chans.labels} {} }, 'normal', [ 1 1 10 1 1]);
   
   % output
   % ------
   if isempty(userdat), chanlist = []; return; end;
   chanliststr = userdat{2};
   chanlist = [];
   for index = 1:length(chanliststr)
       i = strmatch (chanliststr{index}, channelnames);
       chanlist  = [chanlist i];
   end;
   chanlist = sort(chanlist);
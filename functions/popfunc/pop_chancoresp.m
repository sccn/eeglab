% pop_chancoresp() - define correspondances between two channel locations structures
%                    (EEG.chanlocs) automatically (by matching channel labels) 
%                    else using a user input gui.
% Usage:
%   >> [chanlist1 chanlist2] = pop_chancoresp(chanstruct1, chanstruc2, 'key', 'val', ...); 
%
% Inputs:
%   chanstruct1     - first (new) channel locations structure (EEG.chanlocs). 
%                      For details, >> help readlocs
%   chanstruct2     - second (reference) chanlocs structure.
%
% Optional parameters:
%   'gui'           - ['on'|'off'] display gui or not ('on' -> yes)
%   'autoselect'    - ['none'|'fiducials'|'all'] automatically pair channels
%   'chaninfo1'     - EEG.chaninfo structure for first (new) EEG.chanlocs
%   'chaninfo2'     - EEG.chaninfo structure for second (reference) EEG.chanlocs
%   'chanlist1'     - [integer] selected channel to pair in the graphic interface
%                     for the first channel structure. This requires the input 
%                     'chanlist2' below.
%   'chanlist2'     - [integer] selected channel to pair in the graphic interface
%                     for the first channel structure. This requires the input 
%                     requires the input 'chanlist1' that must be of the same length.
% Output:
%   chanlist1       - [int vector] indices of selected channels from first (new) EEG.chanlocs
%   chanlist2       - [int vector] selected channels from second (reference) EEG.chanlocs
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2005

% Copyright (C) 2005 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [chanlistout1, chanlistout2, thirdout, outfourth] = pop_chancoresp(chans1, chans2, varargin); 
    
    if nargin < 2
        help pop_chancoresp;
        return;
    end;
    chanlistout1 = [];
    chanlistout2 = [];
    
    % process sub command
    % -------------------
    if isstr(chans1) 
       if strcmpi(chans1, 'pair')
           [chanlistout1, chanlistout2, thirdout, outfourth] = pair(chans2, varargin{:});
       elseif strcmpi(chans1, 'unpair')
           [chanlistout1, chanlistout2, thirdout, outfourth] = unpair(chans2, varargin{:});
       elseif strcmpi(chans1, 'clear')
           [chanlistout1, chanlistout2, thirdout, outfourth] = clearchans(chans2, varargin{:});
       elseif strcmpi(chans1, 'auto')
           [chanlistout1, chanlistout2, thirdout, outfourth] = autoselect(chans2, varargin{:});
       end;
       return;
    end;
    
    g = finputcheck(varargin, { 'autoselect'     'string'  {'none';'fiducials';'all'}   'all';
                                'chanlist1'      'integer' [1 Inf] []; 
                                'chanlist2'      'integer' [1 Inf] []; 
                                'chaninfo1'      ''        []      []; 
                                'chaninfo2'      ''        []      []; 
                                'gui'            'string'  { 'on';'off' }  'on' } );
    if isstr(g), error(g); end;
    g.chanstruct1 = chans1;
    g.chanstruct2 = chans2;
    if length(g.chanlist1) ~= length(g.chanlist2)
        error('input arguments ''chanlist1'' and ''chanlist2'' must have the same length');
    end;
    
    % decode different input formats
    % ------------------------------
    if isstruct(chans1)
        if isfield(chans1, 'label') % fieldtrip
            chanstr1 = chans1.label;
            chanstr2 = chans2.label;
        else % EEGLAB
            chanstr1 = { chans1.labels };
            chanstr2 = { chans2.labels };
        end;
    else % only channel labels
        chanstr1 = chans1;
        chanstr2 = chans2;
    end;
    
    % convert selection to integer
    % ----------------------------
    if isempty(g.chanlist1)
        if strcmpi(g.autoselect, 'fiducials')
            % find fiducials in both channel location strustures
            % --------------------------------------------------
            naz1 = strmatch('nz', lower( chanstr1 ), 'exact'); if isempty(naz1), naz1 = strmatch('nasion', lower( chanstr1 ), 'exact'); end; if isempty(naz1), naz1 = strmatch('fidnz', lower( chanstr1 ), 'exact'); end;
            naz2 = strmatch('nz', lower( chanstr2 ), 'exact'); if isempty(naz2), naz2 = strmatch('nasion', lower( chanstr2 ), 'exact'); end; if isempty(naz2), naz2 = strmatch('fidnz', lower( chanstr2 ), 'exact'); end;
            lpa1 = strmatch('lpa', lower( chanstr1 ), 'exact'); if isempty(lpa1), lpa1 = strmatch('left',  lower( chanstr1 ), 'exact'); end; if isempty(lpa1), lpa1 = strmatch('fidt10', lower( chanstr1 ), 'exact'); end;
            lpa2 = strmatch('lpa', lower( chanstr2 ), 'exact'); if isempty(lpa2), lpa2 = strmatch('left',  lower( chanstr2 ), 'exact'); end; if isempty(lpa2), lpa2 = strmatch('fidt10', lower( chanstr2 ), 'exact'); end;
            rpa1 = strmatch('rpa', lower( chanstr1 ), 'exact'); if isempty(rpa1), rpa1 = strmatch('right', lower( chanstr1 ), 'exact'); end; if isempty(rpa1), rpa1 = strmatch('fidt9', lower( chanstr1 ), 'exact'); end;
            rpa2 = strmatch('rpa', lower( chanstr2 ), 'exact'); if isempty(rpa2), rpa2 = strmatch('right', lower( chanstr2 ), 'exact'); end; if isempty(rpa2), rpa2 = strmatch('fidt9', lower( chanstr2 ), 'exact'); end;
            g.chanlist1 = [ naz1 lpa1 rpa1 ];
            g.chanlist2 = [ naz2 lpa2 rpa2 ];
            if length(g.chanlist1) ~= length(g.chanlist2) | length(g.chanlist1) == 0
                disp('Warning: could not find fiducials in at least one of the channel location structure');
                g.chanlist1 = [];
                g.chanlist2 = [];
            end;
        elseif strcmpi(g.autoselect, 'all')
            % find common channels in both channel location strustures
            % --------------------------------------------------------
            chanstr2low = lower( chanstr2 );
            chanstr1low = lower( chanstr1 );
            for index = 1:length( chanstr1 )
                ind = strmatch(chanstr1low{index}, chanstr2low, 'exact' );
                if ~isempty(ind)
                    g.chanlist1(end+1) = index;
                    g.chanlist2(end+1) = ind;
                end;
            end;
        end;
    end;
    
    % plot
    % ----
    if strcmpi(g.gui, 'off')
        chanlistout1 = g.chanlist1;
        chanlistout2 = g.chanlist2;
        return;
    end;
    try,  g.promptstring;  catch, g.promptstring = ''; end;
    try,  g.selectionmode; catch, g.selectionmode = 'multiple'; end;
    try,  g.listsize;      catch, g.listsize = []; end;
    try,  g.initialvalue;  catch, g.initialvalue = []; end;
    try,  g.name;          catch, g.name = ''; end;
	g.chanstr1 = chanstr1;
	g.chanstr2 = chanstr2;

    fig = figure('visible', 'off');
    set(fig, 'name', 'Select corresponding channels to pair');
    
    % make text for list
    % ------------------
    [ g.newchanstr1 g.newchanstr2 ] = makelisttext( chanstr1, chanstr2, g.chanlist1, g.chanlist2);
    
    % callback for paring and unpairing
    % ---------------------------------
    cb_pair = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                'tmpval1 = get(findobj(gcbf, ''tag'', ''list1''), ''value'');' ...
                'tmpval2 = get(findobj(gcbf, ''tag'', ''list2''), ''value'');' ...
                '[tmpdat.newchanstr1{tmpval1}, tmpdat.newchanstr2{tmpval2}, tmpdat.chanlist1, tmpdat.chanlist2] = pop_chancoresp(''pair'', tmpval1, tmpval2, tmpdat.chanstr1, tmpdat.chanstr2, tmpdat.chanlist1, tmpdat.chanlist2, tmpdat.newchanstr1{tmpval1}, tmpdat.newchanstr2{tmpval2});' ...
                'set(findobj(gcbf, ''tag'', ''list1''), ''string'', tmpdat.newchanstr1);' ...
                'set(findobj(gcbf, ''tag'', ''list2''), ''string'', tmpdat.newchanstr2);' ...
                'set(gcbf, ''userdata'', tmpdat);' ...
                'clear tmpdat tmpval1 tmpval2;' ];
    cb_unpair = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                'tmpval1 = get(findobj(gcbf, ''tag'', ''list1''), ''value'');' ...
                'tmpval2 = get(findobj(gcbf, ''tag'', ''list2''), ''value'');' ...
                '[tmpdat.newchanstr1{tmpval1}, tmpdat.newchanstr2{tmpval2}, tmpdat.chanlist1, tmpdat.chanlist2] = pop_chancoresp(''unpair'', tmpval1, tmpval2, tmpdat.chanstr1, tmpdat.chanstr2, tmpdat.chanlist1, tmpdat.chanlist2, tmpdat.newchanstr1{tmpval1}, tmpdat.newchanstr2{tmpval2});' ...
                'set(findobj(gcbf, ''tag'', ''list1''), ''string'', tmpdat.newchanstr1);' ...
                'set(findobj(gcbf, ''tag'', ''list2''), ''string'', tmpdat.newchanstr2);' ...
                'set(gcbf, ''userdata'', tmpdat);' ...
                'clear tmpdat tmpval1 tmpval2;' ];
    cb_plot1 = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                 'figure; topoplot([], tmpdat.chanstruct1, ''style'', ''blank'', ''drawaxis'', ''on'', ' ...
                 '''electrodes'', ''labelpoint'', ''chaninfo'', tmpdat.chaninfo1);' ];
    cb_plot2 = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                 'figure; topoplot([], tmpdat.chanstruct2, ''style'', ''blank'', ''drawaxis'', ''on'', ' ...
                 '''electrodes'', ''labelpoint'', ''chaninfo'', tmpdat.chaninfo2);' ];
    cb_list1 = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                'tmpval = get(findobj(gcbf, ''tag'', ''list1''), ''value'');' ...
                'tmppos = find(tmpdat.chanlist1 == tmpval);' ...
                'if ~isempty(tmppos), set(findobj(gcbf, ''tag'', ''list2''), ''value'', tmpdat.chanlist2(tmppos)); end;' ...
                'clear tmpdat tmpval tmppos;' ];
    cb_list2 = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                'tmpval = get(findobj(gcbf, ''tag'', ''list2''), ''value'');' ...
                'tmppos = find(tmpdat.chanlist2 == tmpval);' ...
                'if ~isempty(tmppos), set(findobj(gcbf, ''tag'', ''list1''), ''value'', tmpdat.chanlist1(tmppos)); end;' ...
                'clear tmpdat tmpval tmppos;' ];
    cb_clear = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                 '[tmpdat.newchanstr1, tmpdat.newchanstr2, tmpdat.chanlist1, tmpdat.chanlist2] = pop_chancoresp(''clear'', tmpdat.chanstr1, tmpdat.chanstr2);' ...
                 'set(gcbf, ''userdata'', tmpdat);' ...
                 'set(findobj(gcbf, ''tag'', ''list1''), ''string'', tmpdat.newchanstr1);' ...
                 'set(findobj(gcbf, ''tag'', ''list2''), ''string'', tmpdat.newchanstr2);' ...
                 'set(gcbf, ''userdata'', tmpdat);' ...
                 'clear tmpdat;' ];
    cb_auto  = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                 '[tmpdat.newchanstr1, tmpdat.newchanstr2, tmpdat.chanlist1, tmpdat.chanlist2] = pop_chancoresp(''auto'', tmpdat.chanstr1, tmpdat.chanstr2, tmpdat.newchanstr1, tmpdat.newchanstr2, tmpdat.chanlist1, tmpdat.chanlist2);' ...
                 'set(gcbf, ''userdata'', tmpdat);' ...
                 'set(findobj(gcbf, ''tag'', ''list1''), ''string'', tmpdat.newchanstr1);' ...
                 'set(findobj(gcbf, ''tag'', ''list2''), ''string'', tmpdat.newchanstr2);' ...
                 'set(gcbf, ''userdata'', tmpdat);' ...
                 'clear tmpdat;' ];
    geometry = {[1 1] [1 1] [1 1] [1 1] [1 1]};
    geomvert = [ 1 min(max(length(chanstr1),length(chanstr2)), 10) 1 1 1];
    listui = { ...
              { 'Style', 'pushbutton', 'string', 'Plot new montage', 'callback', cb_plot1 }  ...
              { 'Style', 'pushbutton', 'string', 'Plot ref montage', 'callback', cb_plot2 }  ...
              { 'Style', 'listbox', 'tag', 'list1', 'string', g.newchanstr1, 'value', 1, 'min', 1, 'max', 2, 'callback', cb_list1 } ...
              { 'Style', 'listbox', 'tag', 'list2', 'string', g.newchanstr2, 'value', 1, 'min', 1, 'max', 2, 'callback', cb_list2 } ...
              { 'Style', 'pushbutton', 'string', 'Pair channels'  , 'callback', cb_pair } ...     
              { 'Style', 'pushbutton', 'string', 'Clear this pair', 'callback', cb_unpair } ...     
              { 'Style', 'pushbutton', 'string', 'Clear all pairs' , 'callback', cb_clear } ...     
              { 'Style', 'pushbutton', 'string', 'Auto select', 'callback', cb_auto } ...     
              { 'Style', 'pushbutton', 'string', 'Cancel', 'callback', 'close(gcbf);' }  ...
              { 'Style', 'pushbutton', 'string', 'Ok'    , 'tag', 'ok', 'callback', ['set(gcbo, ''userdata'', ''ok'');'] } };

    [tmp tmp2 allobj] = supergui( fig, geometry, geomvert, listui{:} );
    set(fig, 'userdata', g);
    
    % decode output
    % -------------
    okbut = findobj( 'parent', fig, 'tag', 'ok');
    figure(fig);
    drawnow;
    waitfor( okbut, 'userdata');
    try,
        tmpdat = get(fig, 'userdata');
        chanlistout1 = tmpdat.chanlist1;
        chanlistout2 = tmpdat.chanlist2;
        close(fig);
        drawnow;
    end;

% unpair channels
% ---------------
function [ str1, str2, chanlist1, chanlist2 ] = unpair(ind1, ind2, chanstr1, chanstr2, chanlist1, chanlist2, str1, str2);
    if nargin > 4
        if isempty(find(chanlist2 == ind2)), disp('Channels not associated'); return; end;
        if isempty(find(chanlist1 == ind1)), disp('Channels not associated'); return; end;
    end;
    
    str1 = sprintf('%2d - %3s', ind1, chanstr1{ind1});
    str2 = sprintf('%2d - %3s', ind2, chanstr2{ind2});
    if nargout > 2
        tmppos = find( chanlist1 == ind1);
        chanlist1(tmppos) = [];
        chanlist2(tmppos) = [];
    end;
    
% pair channels
% -------------
function [ str1, str2, chanlist1, chanlist2 ] = pair(ind1, ind2, chanstr1, chanstr2, chanlist1, chanlist2, str1, str2);
    if nargin > 4
        if ~isempty(find(chanlist2 == ind2)), disp('Channel in second structure already associated'); return; end;
        if ~isempty(find(chanlist1 == ind1)), disp('Channel in first structure already associated'); return; end;
    end;
    
    str1 = sprintf('%2d - %3s   -> %2d - %3s', ind1, chanstr1{ind1}, ind2, chanstr2{ind2});
    str2 = sprintf('%2d - %3s   -> %2d - %3s', ind2, chanstr2{ind2}, ind1, chanstr1{ind1});
    if nargout > 2
        chanlist1 = [ chanlist1 ind1 ];
        chanlist2 = [ chanlist2 ind2 ];
    end;

% make full channel list    
% ----------------------
function [ newchanstr1, newchanstr2 ] = makelisttext( chanstr1, chanstr2, chanlist1, chanlist2);    
    for index = 1:length(chanstr1)
        if ismember(index, chanlist1)
            pos = find(chanlist1 == index);
            newchanstr1{index} = pair( chanlist1(pos), chanlist2(pos), chanstr1, chanstr2 );
        else
            newchanstr1{index} = sprintf('%2d - %3s', index, chanstr1{index});
        end;
    end;
    for index = 1:length(chanstr2)
        if ismember(index, chanlist2)
            pos = find(chanlist2 == index);
            [tmp newchanstr2{index}] = pair( chanlist1(pos), chanlist2(pos), chanstr1, chanstr2 );
        else
            newchanstr2{index} = sprintf('%2d - %3s', index, chanstr2{index});
        end;
    end;

% clear channel pairs
% -------------------
function [ newchanstr1, newchanstr2, chanlist1, chanlist2 ] = clearchans(chanstr1, chanstr2);
    
    chanlist1 = [];
    chanlist2 = [];
    [ newchanstr1, newchanstr2 ] = makelisttext( chanstr1, chanstr2, chanlist1, chanlist2);

% autoselect channel pairs
% ------------------------
function [ newchanstr1, newchanstr2, chanlist1, chanlist2 ] = autoselect(chanstr1, chanstr2, newchanstr1, newchanstr2, chanlist1, chanlist2);
    
    % GUI for selecting pairs
    % -----------------------
    listoptions = { 'All channels (same labels)' 'Fiducials (same labels)' 'BIOSEMI -> 10-20' ...
                    'EGI -> 10-20' }; 
    listui = { { 'style' 'text' 'string' [ 'How to pair channels (click to select)' 10 ] } ...
               { 'style' 'listbox' 'string' listoptions 'value' 1 } };
    results = inputgui({ [1 1] }, listui, '', 'Auto select channel pairs', [], 'normal', 2);
    
    % decode results
    % --------------
    if isempty(results), return; end;
    if results{1} == 1 % select all pairs
        [chanlist1 chanlist2] = pop_chancoresp(chanstr1, chanstr2, 'autoselect', 'all', 'gui', 'off');
    elseif results{1} == 2
        [chanlist1 chanlist2] = pop_chancoresp(chanstr1, chanstr2, 'autoselect', 'fiducials', 'gui', 'off');
    else
        disp('Not implemented yet'); return;
    end;
    [ newchanstr1, newchanstr2 ] = makelisttext( chanstr1, chanstr2, chanlist1, chanlist2);
        
    
    

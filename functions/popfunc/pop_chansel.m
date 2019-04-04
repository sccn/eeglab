% pop_chansel() - pop up a graphic interface to select channels
%
% Usage:
%   >> [chanlist] = pop_chansel(chanstruct); % a window pops up
%   >> [chanlist strchannames cellchannames] = ...
%                        pop_chansel(chanstruct, 'key', 'val', ...);
%
% Inputs:
%   chanstruct     - channel structure. See readlocs()
%
% Optional input:
%   'withindex'      - ['on'|'off'] add index to each entry. May also a be 
%                      an array of indices
%   'select'         - selection of channel. Can take as input all the
%                      outputs of this function.
%   'field'          - ['type'|'labels'] information to select. Default is 
%                      'labels'
%   'handle'         - [handle] update handle (GUI)
%   'selectionmode'  - selection mode 'multiple' or 'single'. See listdlg2().
%
% Output:
%   chanlist      - indices of selected channels
%   strchannames  - names of selected channel names in a concatenated string
%                   (channel names are separated by space characters)
%   cellchannames - names of selected channel names in a cell array
%
% Author: Arnaud Delorme, CNL / Salk Institute, 3 March 2003

% Copyright (C) 3 March 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [chanlist,chanliststr, allchanstr] = pop_chansel(chans, varargin); 
    
    if nargin < 1
        help pop_chansel;
        return;
    end
    if isempty(chans), disp('Empty input'); return; end
    if isnumeric(chans),
        for c = 1:length(chans)
            newchans{c} = num2str(chans(c));
        end
        chans = newchans;
    end
    chanlist    = [];
    chanliststr = {};
    allchanstr  = '';
    
    g = finputcheck(varargin, { 'withindex'     '' [] 'off';
                                'select'        '' [] [];
                                'handle'        '' [] [];
                                'field'         'string' [] 'labels';
                                'selectionmode' 'string' { 'single';'multiple' } 'multiple'});
    if ischar(g), error(g); end
    if ~ischar(g.withindex), chan_indices = g.withindex; g.withindex = 'on';
    else                    chan_indices = 1:length(chans);
    end
    if isstruct(chans), chans = { chans.(g.field) }; end
    if strcmpi(g.field, 'type'), chans = unique_bc(chans); end
        
    % convert selection to integer
    % ----------------------------
    if ischar(g.select) && ~isempty(g.select)
        g.select = parsetxt(g.select);
    end
    if iscell(g.select) && ~isempty(g.select)
        if ischar(g.select{1})
            tmplower = lower( chans );
            for index = 1:length(g.select)
                matchind = strmatch(lower(g.select{index}), tmplower, 'exact');
                if ~isempty(matchind), g.select{index} = matchind;
                else error( [ 'Cannot find ''' g.select{index} '''' ] );
                end
            end
        end
        g.select = [ g.select{:} ];
    end
    if ~isnumeric( g.select ), g.select = []; end
    
    % add index to channel name
    % -------------------------
	tmpstr = {chans};
    if isnumeric(chans{1})
        tmpstr = [ chans{:} ];
        tmpfieldnames = cell(1, length(tmpstr));
        for index=1:length(tmpstr), 
            if strcmpi(g.withindex, 'on')
                tmpfieldnames{index} = [ num2str(chan_indices(index)) '  -  ' num2str(tmpstr(index)) ]; 
            else
                tmpfieldnames{index} = num2str(tmpstr(index)); 
            end
        end
    else
        tmpfieldnames = chans;
        if strcmpi(g.withindex, 'on')
            for index=1:length(tmpfieldnames), 
                tmpfieldnames{index} = [ num2str(chan_indices(index)) '  -  ' tmpfieldnames{index} ]; 
            end
        end
    end
    [chanlist,tmp,chanliststr] = listdlg2('PromptString',strvcat('(use shift|Ctrl to', 'select several)'), ...
                'ListString', tmpfieldnames, 'initialvalue', g.select, 'selectionmode', g.selectionmode);   
    if tmp == 0
        chanlist = [];
        chanliststr = '';
        return;
    else
        allchanstr = chans(chanlist);
    end
    
    % test for spaces
    % ---------------
    spacepresent = 0;
    if ~isnumeric(chans{1})
        tmpstrs = [ allchanstr{:} ];
        if ~isempty( find(tmpstrs == ' ')) || ~isempty( find(tmpstrs == 9))
            spacepresent = 1;
        end
    end
    
    % get concatenated string (if index)
    % -----------------------
    if strcmpi(g.withindex, 'on') || spacepresent
        if isnumeric(chans{1})
            chanliststr = num2str(celltomat(allchanstr));
        else
            chanliststr = '';
            for index = 1:length(allchanstr)
                if spacepresent
                    chanliststr = [ chanliststr '''' allchanstr{index} ''' ' ];
                else
                    chanliststr = [ chanliststr allchanstr{index} ' ' ];
                end
            end
            chanliststr = chanliststr(1:end-1);
        end
    end
    
    if ~isempty(g.handle)
        set(g.handle, 'string', chanliststr);
    end
       
    return;

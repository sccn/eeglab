% forcelocs() - rotate location in 3-D so specified electrodes
%               match specified locations.
%
% Usage:
%   >> chanlocs = forcelocs( chanlocs ); % pop-up window mode
%   >> chanlocs = forcelocs( filename, loc1, loc2, ... );
%
% Inputs:
%   chanlocs  - EEGLAB channel structure. See help readlocs()
%
% Optional inputs:
%   loc1      - cell array. First element is spherical horizontal angle;
%               second element is spherical elevation angle 
%               (90 = vertical); other elements are channel indices or
%               name. If several channel names are given, the function
%               set the average of channel horizontal angle and elevation
%               to the new values given as input.
%   loc2      - same as loc1
%
% Outputs:
%   chanlocs  - updated EEGLAB channel structure.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 April 2003
%
% See also: readlocs()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.2  2003/04/16 01:56:47  arno
% adding hidden output
%
% Revision 1.1  2003/04/16 01:48:32  arno
% Initial revision
%

function [chanlocs,options] = forcelocs( chanlocs, varargin)
    
    NENTRY = 1; % number of lines in GUI
    FIELDS = { 'X' 'Y' };
    
    options = [];
    if nargin < 1
        help forcelocs;
        return;
    end;
    if nargin < 2
        geom = { [0.3 1 1 0.3] };
        uilist = { { 'style' 'text' 'string' 'Value' } ...
                   { 'style' 'text' 'string' 'Coordinate' } ...
                   { 'style' 'text' 'string' 'Electrode list' } ...
                   { } };
        for index = 1:NENTRY
            tag = [ 'c' int2str(index) ];
            geom = { geom{:}  [0.3 1 1 0.3] };
            uilist = { uilist{:} { 'style' 'edit' 'string' fastif(index==1, '0','') } ...
                       { 'style' 'listbox' 'string' 'X (rotate X-Z plane)|Y (rotate Y-Z plane)' } ...
                       { 'style' 'edit' 'string'  fastif(index==1, 'Cz','') 'tag' tag } ...
                       { 'style' 'pushbutton' 'string' 'Pick' ...
                         'callback', [ '[tmp1 tmp2 tmp3] = pop_chansel(EEG.chanlocs);' ...
                                       'if ~isempty(tmp3) set(findobj(gcbf, ''tag'', ''' tag '''), ''string'', tmp3); end;' ...
                                       'clear tmp1 tmp2 tmp3;' ] } };
        end;
        
        results = inputgui( geom, uilist, 'pophelp(''forcelocs'');', 'Force electrode location -- forcelocs()' );
        if length(results) == 0, return; end;
        
        options = {};
        for index = 1:NENTRY
            tmpi = 3*(index-1)+1;
            if ~isempty(results{tmpi})
                tmpchans = parsetxt(results{tmpi+2});
                options = { options{:} { str2num(results{tmpi}) FIELDS{results{tmpi+1}} tmpchans{:} }};
            end;
        end;    
    else 
        options = varargin;
    end;

    % scan all locations
    % ------------------
    channelnames = lower(strvcat({chanlocs.labels}));
    for index = 1:length(options)
        
        val   = options{index}{1};
        type  = options{index}{2};
        chans = getchans(options{index}(3:end), channelnames);

        % rotate X-Z plane 
        % ----------------
        if strcmpi(type, 'x')
            curx   = mean(cell2mat( { chanlocs(chans).X }));
            curz   = mean(cell2mat( { chanlocs(chans).Z }));
            newx = val;
            rotangle = solvesystem(curx, curz, newx);
            
            for chanind = 1:length(chanlocs)
                [chanlocs(chanind).X chanlocs(chanind).Z]= rotation(chanlocs(chanind).X, chanlocs(chanind).Z, rotangle);
            end;
            chanlocs = convertlocs(chanlocs, 'cart2all');
        end;
        
        % rotate Y-Z plane 
        % ----------------
        if strcmpi(type, 'y')
            cury   = mean(cell2mat( { chanlocs(chans).Y }));
            curz   = mean(cell2mat( { chanlocs(chans).Z }));
            newy = val;
            rotangle = solvesystem(cury, curz, newy);
            
            for chanind = 1:length(chanlocs)
                [chanlocs(chanind).Y chanlocs(chanind).Z]= rotation(chanlocs(chanind).Y, chanlocs(chanind).Z, rotangle);
            end;
            chanlocs = convertlocs(chanlocs, 'cart2all');
        end;
    
    end;
        

% get channel indices
% -------------------
function chanlist = getchans(chanliststr, channelnames);
    chanlist = [];
    for index = 1:length(chanliststr)
        i = strmatch (lower(chanliststr{index}), channelnames, 'exact');
        chanlist  = [chanlist i];
    end;

% function rotate coordinates
% ---------------------------
function [X,Y] = rotation(x,y,rotangle);
    X = real((x+j*y)*exp(j*rotangle));
    Y = imag((x+j*y)*exp(j*rotangle));
    
% function solvesyst
% ------------------
function theta = solvesystem(x,y,nx);
    eq(1,:) = [x -y]; res(1) = nx;
    eq(2,:) = [y x];  res(2) = sqrt(x^2+y^2-nx^2);
    sol = eq\res';
    theta = atan2(sol(2), sol(1));
    
    % simplier solution
    ny = sqrt(x^2+y^2-nx^2);
    ang1 = angle(x+j*y);
    ang2 = angle(nx+j*ny);
    theta = ang2-ang1;
    
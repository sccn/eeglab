% pop_chancenter() - recenter cartesian X,Y,Z channel coordinates
%
% Usage:  
%    >> chanlocs = pop_chancenter(chanlocs); % pop up interactive window
%    >> [chanlocs centerloc] = pop_chancenter(chanlocs, center, omitchan); 
%
% Inputs:
%    chanlocs  = eeglab channel location structure (see readlocs())
%    center    = [X Y Z] known center different from [0 0 0]
%                [] will optimize the center location according
%                to the best sphere. Default is [0 0 0].
%    omitchan  = indices of channel to omit when computing center
%
% Outputs:
%    chanlocs  = updated channel location structure
%    centerloc = 3-D location of the new center (in the old coordinate
%                frame).
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, Feb 2004
%
% See also: chancenter(), spherror(), cart2topo()

% Copyright (C) 2004, Arnaud Delorme, SCCN/INC/UCSD, arno@salk.edu
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

function [ chanlocs, newcenter, com] = pop_chancenter( chanlocs, center, omitchans)

optim = 0;

if nargin<1
    help pop_chancenter
    return;
end;

com = '';
newcenter = [];
if nargin < 3
    omitchans = [];
end;
if nargin < 2
    cb_browse = [ 'tmpchans = get(gcbf, ''userdata'');' ...
                  'set(findobj(gcbf, ''tag'', ''chans''), ''string'', ' ...
                  'int2str(pop_chansel( { tmpchans.labels } )));' ];
    cb_opt    = [ 'if get(gcbo, ''value''), ' ...
                  '    set(findobj(gcbf, ''tag'', ''center''), ''enable'', ''off'');' ...
                  'else,' ...
                  '    set(findobj(gcbf, ''tag'', ''center''), ''enable'', ''on'');' ...
                  'end;' ];
    geometry = { [1.3 0.28 1 1] [1] [1] [2 1] };
    uilist = { { 'Style', 'text',       'string', 'Optimize center location', 'fontweight', 'bold'   } ...
			   { 'Style', 'checkbox',   'value',  1  'callback' cb_opt } ... 
               { 'Style', 'text',       'string', 'or specify center', 'fontweight', 'bold'  } ...
               { 'Style', 'edit',       'string', '0 0 0', 'tag' 'center' 'enable' 'off' } ...
               { } ...
			   { 'Style', 'text',       'string', 'Channel indices to ignore for best-sphere matching'  } ...
			   { 'Style', 'edit',       'string', '', 'tag', 'chans' } ...
			   { 'Style', 'pushbutton', 'string', 'Browse', 'callback', cb_browse } };
    
    results = inputgui( geometry, uilist, 'pophelp(''pop_chancenter'');', ...
                        'Convert channel locations -- pop_chancenter()', chanlocs );
	if isempty(results), return; end;
	if results{1}
        center = [];
    else
        center  = eval( [ '[' results{2} ']' ] );
    end;
    if ~isempty(results{3})
        omitchans =  eval( [ '[' results{3} ']' ] );
    end;
end;

% remove channels
% ---------------
c = setdiff([1:length(chanlocs)], union(omitchans, find(cellfun('isempty', { chanlocs.theta }))));

% optimize center
% ---------------
[X Y Z newcenter]= chancenter( [ chanlocs(c).X ]', [ chanlocs(c).Y ]', [ chanlocs(c).Z ]', center);
for index = 1:length(c)
    chanlocs(c(index)).X  = X(index);
    chanlocs(c(index)).Y  = Y(index);
    chanlocs(c(index)).Z  = Z(index);
end;
disp('Note: automatically convert XYZ coordinates to spherical and polar');
chanlocs = convertlocs(chanlocs, 'cart2all');
if ~isempty(omitchans)
    disp('Important warning: the location of omitted channels has not been modified');
end;

com = sprintf('%s = pop_chancenter( %s, %s);', inputname(1), inputname(1), vararg2str({ center omitchans }));

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

function [ chanlocs, newcenter, com] = pop_chancenter( chanlocs, center, omitchans)

optim = 0;

if nargin<1
    help pop_chancenter
    return;
end

com = '';
newcenter = [];
if nargin < 3
    omitchans = [];
end
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
	if isempty(results), return; end
	if results{1}
        center = [];
    else
        center  = eval( [ '[' results{2} ']' ] );
    end
    if ~isempty(results{3})
        omitchans =  eval( [ '[' results{3} ']' ] );
    end
end

% remove channels
% ---------------
c = setdiff_bc([1:length(chanlocs)], union(omitchans, find(cellfun('isempty', { chanlocs.theta }))));

% optimize center
% ---------------
[X, Y, Z, newcenter]= chancenter( [ chanlocs(c).X ]', [ chanlocs(c).Y ]', [ chanlocs(c).Z ]', center);
for index = 1:length(c)
    chanlocs(c(index)).X  = X(index);
    chanlocs(c(index)).Y  = Y(index);
    chanlocs(c(index)).Z  = Z(index);
end
disp('Note: automatically convert XYZ coordinates to spherical and polar');
chanlocs = convertlocs(chanlocs, 'cart2all');
if ~isempty(omitchans)
    disp('Important warning: the location of omitted channels has not been modified');
end
if nargout > 2
    com = sprintf('%s = pop_chancenter( %s, %s);', inputname(1), inputname(1), vararg2str({ center omitchans }));
end

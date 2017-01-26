% pop_firma() - Filter data using moving average FIR filter
%
% Usage:
%   >> [EEG, com] = pop_firma(EEG); % pop-up window mode
%   >> [EEG, com] = pop_firma(EEG, 'forder', order);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   'forder'  - scalar filter order. Mandatory even
%
% Outputs:
%   EEG       - filtered EEGLAB EEG structure
%   com       - history string
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   firfilt, plotfresp

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function [EEG, com] = pop_firma(EEG, varargin)

    com = '';
    if nargin < 1
        help pop_firma;
        return;
    end
    if isempty(EEG.data)
        error('Cannot process empty dataset');
    end

    if nargin < 2
        drawnow;
        uigeom = {[1 1 1] [1] [1 1 1]};
        uilist = {{'style' 'text' 'string' 'Filter order (mandatory even):'} ...
                  {'style' 'edit' 'string' '' 'tag' 'forderedit'} {} ...
                  {} ...
                  {} {} {'Style' 'pushbutton' 'string' 'Plot filter responses' 'callback' {@complot, EEG.srate}}};
        result = inputgui(uigeom, uilist, 'pophelp(''pop_firma'')', 'Filter the data -- pop_firma()');
        if length(result) == 0, return; end

        if ~isempty(result{1})
            args = [{'forder'} {str2num(result{1})}];
        else
            error('Not enough input arguments');
        end
    else
        args = varargin;
    end

    % Convert args to structure
    args = struct(args{:});

    % Filter coefficients
    b = ones(1, args.forder + 1) / (args.forder + 1);

    % Filter
    disp('pop_firma() - filtering the data');
    EEG = firfilt(EEG, b);

    % History string
    com = sprintf('%s = pop_firma(%s', inputname(1), inputname(1));
    for c = fieldnames(args)'
        if ischar(args.(c{:}))
            com = [com sprintf(', ''%s'', ''%s''', c{:}, args.(c{:}))];
        else
            com = [com sprintf(', ''%s'', %s', c{:}, mat2str(args.(c{:})))];
        end
    end
    com = [com ');'];

% Callback plot filter properties
function complot(obj, evt, srate)
    args.forder = str2num(get(findobj(gcbf, 'tag', 'forderedit'), 'string'));
    if isempty(args.forder)
        error('Not enough input arguments');
    end
    b = ones(1, args.forder + 1) / (args.forder + 1);
    H = findobj('tag', 'filter responses', 'type', 'figure');
    if ~isempty(H)
        figure(H);
    else
        H = figure;
        set(H, 'color', [.93 .96 1], 'tag', 'filter responses');
    end
    plotfresp(b, 1, [], srate);

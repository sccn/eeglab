% pop_firpm() - Filter data using Parks-McClellan FIR filter
%
% Usage:
%   >> [EEG, com, b] = pop_firpm(EEG); % pop-up window mode
%   >> [EEG, com, b] = pop_firpm(EEG, 'key1', value1, 'key2', ...
%                                value2, 'keyn', valuen);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   'fcutoff' - vector or scalar of cutoff frequency/ies (~-6 dB; Hz)
%   'ftrans'  - scalar transition band width
%   'ftype'   - char array filter type. 'bandpass', 'highpass',
%               'lowpass', or 'bandstop'
%   'forder'  - scalar filter order. Mandatory even
%
% Optional inputs:
%   'wtpass'  - scalar passband weight
%   'wtstop'  - scalar stopband weight
%
% Outputs:
%   EEG       - filtered EEGLAB EEG structure
%   com       - history string
%   b         - filter coefficients
%
% Note:
%   Requires the signal processing toolbox.
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   firfilt, pop_firpmord, plotfresp, firpm, firpmord

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

function [EEG, com, b] = pop_firpm(EEG, varargin)

    if ~(exist('firpm', 'file') == 2 || exist('firpm', 'file') == 6)
       error('Requires the signal processing toolbox.');
    end

    com = '';
    if nargin < 1
        help pop_firpm;
        return;
    end
    if isempty(EEG.data)
        error('Cannot process empty dataset');
    end

    if nargin < 2
        drawnow;
        ftypes = {'bandpass' 'highpass' 'lowpass' 'bandstop'};
        uigeom = {[1 0.75 0.75] [1 0.75 0.75] [1 0.75 0.75] 1 [1 0.75 0.75] [1 0.75 0.75] [1 0.75 0.75] 1 [1 0.75 0.75]};
        uilist = {{'Style' 'text' 'String' 'Cutoff frequency(ies) [hp lp] (~-6 dB; Hz):'} ...
                  {'Style' 'edit' 'String' '' 'Tag' 'fcutoffedit'} {} ...
                  {'Style' 'text' 'String' 'Transition band width:'} ...
                  {'Style' 'edit' 'String' '' 'Tag' 'ftransedit'} {} ...
                  {'Style' 'text' 'String' 'Filter type:'} ...
                  {'Style' 'popupmenu' 'String' ftypes 'Tag' 'ftypepop'} {} ...
                  {} ...
                  {'Style' 'text' 'String' 'Passband weight:'} ...
                  {'Style' 'edit' 'String' '' 'Tag' 'wtpassedit'} {} ...
                  {'Style' 'text' 'String' 'Stopband weight:'} ...
                  {'Style' 'edit' 'String' '' 'Tag' 'wtstopedit'} {} ...
                  {'Style' 'text' 'String' 'Filter order (mandatory even):'} ...
                  {'Style' 'edit' 'String' '' 'Tag' 'forderedit'} ...
                  {'Style' 'pushbutton' 'String' 'Estimate' 'Tag' 'orderpush' 'Callback' {@comcb, ftypes, EEG.srate}} ...
                  {} ...
                  {} {} {'Style' 'pushbutton' 'String', 'Plot filter responses' 'Tag' 'plotpush' 'Callback' {@comcb, ftypes, EEG.srate}}};
        result = inputgui(uigeom, uilist, 'pophelp(''pop_firpm'')', 'Filter the data -- pop_firpm()');
        if isempty(result), return; end

        args = {};
        if ~isempty(result{1})
            args = [args {'fcutoff'} {str2num(result{1})}];
        end
        if ~isempty(result{2})
            args = [args {'ftrans'} {str2double(result{2})}];
        end
        args = [args {'ftype'} ftypes(result{3})];
        if ~isempty(result{4})
            args = [args {'wtpass'} {str2double(result{4})}];
        end
        if ~isempty(result{5})
            args = [args {'wtstop'} {str2double(result{5})}];
        end
        if ~isempty(result{6})
            args = [args {'forder'} {str2double(result{6})}];
        end
    else
        args = varargin;
    end

    % Convert args to structure
    args = struct(args{:});

    c = parseargs(args, EEG.srate);
    if ~isfield(args, 'forder') || isempty(args.forder)
        error('Not enough input arguments');
    end
    b = firpm(args.forder, c{:});

    % Filter
    disp('pop_firpm() - filtering the data');
    EEG = firfilt(EEG, b);

    % History string
    com = sprintf('%s = pop_firpm(%s', inputname(1), inputname(1));
    for c = fieldnames(args)'
        if ischar(args.(c{:}))
            com = [com sprintf(', ''%s'', ''%s''', c{:}, args.(c{:}))];
        else
            com = [com sprintf(', ''%s'', %s', c{:}, mat2str(args.(c{:})))];
        end
    end
    com = [com ');'];

% Convert structure args to cell array firpm parameters
function c = parseargs(args, srate)

    if ~isfield(args, 'fcutoff') || ~isfield(args, 'ftype') || ~isfield(args, 'ftrans') || isempty(args.fcutoff) || isempty(args.ftype) || isempty(args.ftrans)
        error('Not enough input arguments.');
    end

    % Cutoff frequencies
    args.fcutoff = [args.fcutoff - args.ftrans / 2 args.fcutoff + args.ftrans / 2];
    args.fcutoff = sort(args.fcutoff / (srate / 2)); % Sorting and normalization
    if any(args.fcutoff < 0)
        error('Cutoff frequencies - transition band width / 2 must not be < DC');
    elseif any(args.fcutoff > 1)
        error('Cutoff frequencies + transition band width / 2 must not be > Nyquist');
    end
    c = {[0 args.fcutoff 1]};

    % Filter type
    switch args.ftype
        case 'bandpass'
            c = [c {[0 0 1 1 0 0]}];
        case 'bandstop'
            c = [c {[1 1 0 0 1 1]}];
        case 'highpass'
            c = [c {[0 0 1 1]}];
        case 'lowpass'
            c = [c {[1 1 0 0]}];
    end

    %Filter weights
    if all(isfield(args, {'wtpass', 'wtstop'})) && ~isempty(args.wtpass) && ~isempty(args.wtstop)
        w = [args.wtstop args.wtpass];
        c{3} = w(c{2}(1:2:end) + 1);
    end

% Callback
function comcb(obj, evt, ftypes, srate)

    args.fcutoff = str2num(get(findobj(gcbf, 'Tag', 'fcutoffedit'), 'String'));
    args.ftype = ftypes{get(findobj(gcbf, 'Tag', 'ftypepop'), 'Value')};
    args.ftrans = str2double(get(findobj(gcbf, 'Tag', 'ftransedit'), 'String'));
    args.wtpass = str2double(get(findobj(gcbf, 'Tag', 'wtpassedit'), 'String'));
    args.wtstop = str2double(get(findobj(gcbf, 'Tag', 'wtstopedit'), 'String'));
    c = parseargs(args, srate);

    switch get(obj, 'Tag')
        case 'orderpush'
            [args.forder, args.wtpass, args.wtstop] = pop_firpmord(c{1}(2:end - 1), c{2}(1:2:end));
            if ~isempty(args.forder) || ~isempty(args.wtpass) || ~isempty(args.wtstop)
                set(findobj(gcbf, 'Tag', 'forderedit'), 'String', ceil(args.forder / 2) * 2);
                set(findobj(gcbf, 'Tag', 'wtpassedit'), 'String', args.wtpass);
                set(findobj(gcbf, 'Tag', 'wtstopedit'), 'String', args.wtstop);
            end

        case 'plotpush'
            args.forder = str2double(get(findobj(gcbf, 'Tag', 'forderedit'), 'String'));
            if isempty(args.forder)
                error('Not enough input arguments');
            end
            b = firpm(args.forder, c{:});
            H = findobj('Tag', 'filter responses', 'Type', 'figure');
            if ~isempty(H)
                figure(H);
            else
                H = figure;
                set(H, 'color', [.93 .96 1], 'Tag', 'filter responses');
            end
            plotfresp(b, 1, [], srate);
    end

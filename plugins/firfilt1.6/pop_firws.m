% pop_firws() - Filter data using windowed sinc FIR filter
%
% Usage:
%   >> [EEG, com, b] = pop_firws(EEG); % pop-up window mode
%   >> [EEG, com, b] = pop_firws(EEG, 'key1', value1, 'key2', ...
%                                value2, 'keyn', valuen);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   'fcutoff' - vector or scalar of cutoff frequency/ies (-6 dB; Hz)
%   'forder'  - scalar filter order. Mandatory even
%
% Optional inputs:
%   'ftype'       - char array filter type. 'bandpass', 'highpass',
%                   'lowpass', or 'bandstop' {default 'bandpass' or
%                   'lowpass', depending on number of cutoff frequencies}
%   'wtype'       - char array window type. 'rectangular', 'bartlett',
%                   'hann', 'hamming', 'blackman', or 'kaiser' {default
%                   'blackman'} 
%   'warg'        - scalar kaiser beta
%   'minphase'    - scalar boolean minimum-phase converted causal filter
%                   {default false}
%
% Outputs:
%   EEG       - filtered EEGLAB EEG structure
%   com       - history string
%   b         - filter coefficients
%
% Note:
%   Window based filters' transition band width is defined by filter
%   order and window type/parameters. Stopband attenuation equals
%   passband ripple and is defined by the window type/parameters. Refer
%   to table below for typical parameters. (Windowed sinc) symmetric FIR
%   filters have linear phase and can be made zero phase (non-causal) by
%   shifting the data by the filters group delay (what firfilt does by
%   default). Pi phase jumps noticable in the phase reponse reflect a
%   negative frequency response and only occur in the stopband. pop_firws
%   also allows causal filtering with minimum-phase (non-linear!) converted
%   filter coefficients with similar properties. Non-linear causal
%   filtering is NOT recommended for most use cases.
%
%               Beta    Max stopband    Max passband    Max passband    Transition width    Mainlobe width
%                       attenuation     deviation       ripple (dB)     (normalized freq)   (normalized rad freq)
%                       (dB)
%   Rectangular         -21             0.0891          1.552           0.9 / m*             4 * pi / m
%   Bartlett            -25             0.0562          0.977                                8 * pi / m
%   Hann                -44             0.0063          0.109           3.1 / m              8 * pi / m
%   Hamming             -53             0.0022          0.038           3.3 / m              8 * pi / m
%   Blackman            -74             0.0002          0.003           5.5 / m             12 * pi / m
%   Kaiser      5.653   -60             0.001           0.017           3.6 / m
%   Kaiser      7.857   -80             0.0001          0.002           5.0 / m
%   * m = filter order
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   firfilt, firws, pop_firwsord, pop_kaiserbeta, plotfresp, windows

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

function [EEG, com, b] = pop_firws(EEG, varargin)

    com = '';
    if nargin < 1
        help pop_firws;
        return;
    end
    if isempty(EEG.data)
        error('Cannot process empty dataset');
    end

    if nargin < 2
        drawnow;
        ftypes = {'bandpass', 'highpass', 'lowpass', 'bandstop'};
        ftypesStr = {'Bandpass', 'Highpass', 'Lowpass', 'Bandstop'};
        wtypes = {'rectangular', 'bartlett', 'hann', 'hamming', 'blackman', 'kaiser'};
        wtypesStr = {'Rectangular (PB dev=0.089, SB att=-21dB)', 'Bartlett (PB dev=0.056, SB att=-25dB)', 'Hann (PB dev=0.006, SB att=-44dB)', 'Hamming (PB dev=0.002, SB att=-53dB)', 'Blackman (PB dev=0.0002, SB att=-74dB)', 'Kaiser'};
        uigeom = {[1 0.75 0.75] [1 0.75 0.75] 1 [1 0.75 0.75] [1 0.75 0.75] [1 0.75 0.75] [1 1.5] 1 [1 0.75 0.75]};
        uilist = {{'Style' 'text' 'String' 'Cutoff frequency(ies) [hp lp] (-6 dB; Hz):'} ...
                  {'Style' 'edit' 'String' '' 'Tag' 'fcutoffedit'} {} ...
                  {'Style' 'text' 'String' 'Filter type:'} ...
                  {'Style' 'popupmenu' 'String' ftypesStr 'Tag' 'ftypepop'} {} ...
                  {} ...
                  {'Style' 'text' 'String' 'Window type:'} ...
                  {'Style' 'popupmenu' 'String' wtypesStr 'Tag' 'wtypepop' 'Value' 5 'Callback' 'temp = {''off'', ''on''}; set(findobj(gcbf, ''-regexp'', ''Tag'', ''^warg''), ''Enable'', temp{double(get(gcbo, ''Value'') == 6) + 1}), set(findobj(gcbf, ''Tag'', ''wargedit''), ''String'', '''')'} {} ...
                  {'Style' 'text' 'String' 'Kaiser window beta:' 'Tag' 'wargtext' 'Enable' 'off'} ...
                  {'Style' 'edit' 'String' '' 'Tag' 'wargedit' 'Enable' 'off'} ...
                  {'Style' 'pushbutton' 'String' 'Estimate' 'Tag' 'wargpush' 'Enable' 'off' 'Callback' @comwarg} ...
                  {'Style' 'text' 'String' 'Filter order (mandatory even):'} ...
                  {'Style' 'edit' 'String' '' 'Tag' 'forderedit'} ...
                  {'Style' 'pushbutton' 'String' 'Estimate' 'Callback' {@comforder, wtypes, EEG.srate}} ...
                  {} {'Style' 'checkbox', 'String', 'Use minimum-phase converted causal filter (non-linear!; beta)', 'Tag' 'minphase', 'Value', 0} ...
                  {'Style' 'edit' 'Tag' 'devedit' 'Visible' 'off'} ...
                  {} {} {'Style' 'pushbutton' 'String', 'Plot filter responses' 'Callback' {@comfresp, wtypes, ftypes, EEG.srate}}};
        result = inputgui(uigeom, uilist, 'pophelp(''pop_firws'')', 'Filter the data -- pop_firws()');
        if isempty(result), return; end

        args = {};
        if ~isempty(result{1})
            args = [args {'fcutoff'} {str2num(result{1})}];
        end
        args = [args {'ftype'} ftypes(result{2})];
        args = [args {'wtype'} wtypes(result{3})];
        if ~isempty(result{4})
            args = [args {'warg'} {str2double(result{4})}];
        end
        if ~isempty(result{5})
            args = [args {'forder'} {str2double(result{5})}];
        end
        args = [args {'minphase'} result{6}];
    else
        args = varargin;
    end

    % Convert args to structure
    args = struct(args{:});

    c = parseargs(args, EEG.srate);
    b = firws(c{:});

    % Check arguments
    if ~isfield(args, 'minphase') || isempty(args.minphase)
        args.minphase = 0;
    end

    % Filter
    disp('pop_firws() - filtering the data');
    if args.minphase
        b = minphaserceps(b);
        EEG = firfiltsplit(EEG, b, 1);
    else
        EEG = firfilt(EEG, b);
    end

    % History string
    com = sprintf('%s = pop_firws(%s', inputname(1), inputname(1));
    for c = fieldnames(args)'
        if ischar(args.(c{:}))
            com = [com sprintf(', ''%s'', ''%s''', c{:}, args.(c{:}))];
        else
            com = [com sprintf(', ''%s'', %s', c{:}, mat2str(args.(c{:})))];
        end
    end
    com = [com ');'];

% Convert structure args to cell array firws parameters
function c = parseargs(args, srate)

    % Filter order and cutoff frequencies
    if ~isfield(args, 'fcutoff') || ~isfield(args, 'forder') || isempty(args.fcutoff) || isempty(args.forder)
        error('Not enough input arguments.');
    end
    c = [{args.forder} {sort(args.fcutoff / (srate / 2))}]; % Sorting and normalization

    % Filter type
    if isfield(args, 'ftype')  && ~isempty(args.ftype)
        if (strcmpi(args.ftype, 'bandpass') || strcmpi(args.ftype, 'bandstop')) && length(args.fcutoff) ~= 2
            error('Not enough input arguments.');
        elseif (strcmpi(args.ftype, 'highpass') || strcmpi(args.ftype, 'lowpass')) && length(args.fcutoff) ~= 1
            error('Too many input arguments.');
        end
        switch args.ftype
            case 'bandstop'
                c = [c {'stop'}];
            case 'highpass'
                c = [c {'high'}];
        end
    end

    % Window type
    if isfield(args, 'wtype')  && ~isempty(args.wtype)
        if strcmpi(args.wtype, 'kaiser')
            if isfield(args, 'warg')  && ~isempty(args.warg)
                c = [c {windows(args.wtype, args.forder + 1, args.warg)'}];
            else
                error('Not enough input arguments.');
            end
        else
            c = [c {windows(args.wtype, args.forder + 1)'}];
        end
    end

% Callback estimate Kaiser beta
function comwarg(varargin)
    [warg, dev] = pop_kaiserbeta;
    set(findobj(gcbf, 'Tag', 'wargedit'), 'String', warg);
    set(findobj(gcbf, 'Tag', 'devedit'), 'String', dev);

% Callback estimate filter order
function comforder(obj, evt, wtypes, srate)
    wtype = wtypes{get(findobj(gcbf, 'Tag', 'wtypepop'), 'Value')};
    dev = get(findobj(gcbf, 'Tag', 'devedit'), 'String');
    [forder, dev] = pop_firwsord(wtype, srate, [], dev);
    set(findobj(gcbf, 'Tag', 'forderedit'), 'String', forder);
    set(findobj(gcbf, 'Tag', 'devedit'), 'String', dev);

% Callback plot filter responses
function comfresp(obj, evt, wtypes, ftypes, srate)
    args.fcutoff = str2num(get(findobj(gcbf, 'Tag', 'fcutoffedit'), 'String'));
    args.ftype = ftypes{get(findobj(gcbf, 'Tag', 'ftypepop'), 'Value')};
    args.wtype = wtypes{get(findobj(gcbf, 'Tag', 'wtypepop'), 'Value')};
    args.warg = str2num(get(findobj(gcbf, 'Tag', 'wargedit'), 'String'));
    args.forder = str2double(get(findobj(gcbf, 'Tag', 'forderedit'), 'String'));
    args.minphase = get(findobj(gcbf, 'Tag', 'minphase'), 'Value');
    causal = args.minphase;
    c = parseargs(args, srate);
    b = firws(c{:});
    if args.minphase
        b = minphaserceps(b);
    end
    H = findobj('Tag', 'filter responses', 'type', 'figure');
    if ~isempty(H)
        figure(H);
    else
        H = figure;
        set(H, 'color', [.93 .96 1], 'Tag', 'filter responses');
    end
    plotfresp(b, 1, [], srate, causal);

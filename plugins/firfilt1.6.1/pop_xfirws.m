% pop_xfirws() - Design and export xfir compatible windowed sinc FIR filter
%
% Usage:
%   >> pop_xfirws; % pop-up window mode
%   >> [b, a] = pop_xfirws; % pop-up window mode
%   >> pop_xfirws('key1', value1, 'key2', value2, 'keyn', valuen);
%   >> [b, a] = pop_xfirws('key1', value1, 'key2', value2, 'keyn', valuen);
%
% Inputs:
%   'srate'   - scalar sampling rate (Hz)
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
%   'filename'    - char array export filename
%   'pathname'    - char array export pathname {default '.'}
%
% Outputs:
%   b         - filter coefficients
%   a         - filter coefficients
%
% Note:
%   Window based filters' transition band width is defined by filter
%   order and window type/parameters. Stopband attenuation equals
%   passband ripple and is defined by the window type/parameters. Refer
%   to table below for typical parameters. (Windowed sinc) FIR filters
%   are zero phase in passband when shifted by the filters group delay
%   (what firfilt does). Pi phase jumps noticable in the phase reponse
%   reflect a negative frequency response and only occur in the
%   stopband.
%
%               Beta    Max stopband    Max passband    Max passband    Transition width    Mainlobe width
%                       attenuation     deviation       ripple (dB)     (normalized freq)   (normalized rad freq)
%                       (dB)
%   Rectangular         -21             0.0891          1.552           0.9 / m*             4 * pi / m
%   Bartlett            -25             0.0562          0.977           (2.9** / m)          8 * pi / m
%   Hann                -44             0.0063          0.109           3.1 / m              8 * pi / m
%   Hamming             -53             0.0022          0.038           3.3 / m              8 * pi / m
%   Blackman            -74             0.0002          0.003           5.5 / m             12 * pi / m
%   Kaiser      5.653   -60             0.001           0.017           3.6 / m
%   Kaiser      7.857   -80             0.0001          0.002           5.0 / m
%   * m = filter order
%   ** estimate for higher m only
%
% Example:
%   fs = 500; tbw = 2; dev = 0.001;
%   beta = pop_kaiserbeta(dev);
%   m = pop_firwsord('kaiser', fs, tbw, dev);
%   pop_xfirws('srate', fs, 'fcutoff', [1 25], 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', beta, 'forder', m, 'filename', 'foo.fir')
%
% Author: Andreas Widmann, University of Leipzig, 2011
%
% See also:
%   firfilt, firws, pop_firwsord, pop_kaiserbeta, plotfresp, windows

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2011 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function [varargout] = pop_xfirws(varargin)

% Pop-up window mode
if nargin < 1

    drawnow;
    ftypes = {'bandpass' 'highpass' 'lowpass' 'bandstop'};
    wtypes = {'rectangular' 'bartlett' 'hann' 'hamming' 'blackman' 'kaiser'};
    uigeom = {[1 0.75 0.75] 1 [1 0.75 0.75] [1 0.75 0.75] 1 [1 0.75 0.75] [1 0.75 0.75] [1 0.75 0.75] 1 [1 0.75 0.75]};
    uilist = {{'Style' 'text' 'String' 'Sampling frequency (Hz):'} ...
              {'Style' 'edit' 'String' '2' 'Tag' 'srateedit'} {} ...
              {} ...
              {'Style' 'text' 'String' 'Cutoff frequency(ies) [hp lp] (-6 dB; Hz):'} ...
              {'Style' 'edit' 'String' '' 'Tag' 'fcutoffedit'} {} ...
              {'Style' 'text' 'String' 'Filter type:'} ...
              {'Style' 'popupmenu' 'String' ftypes 'Tag' 'ftypepop'} {} ...
              {} ...
              {'Style' 'text' 'String' 'Window type:'} ...
              {'Style' 'popupmenu' 'String' wtypes 'Tag' 'wtypepop' 'Value' 5 'Callback' 'temp = {''off'', ''on''}; set(findobj(gcbf, ''-regexp'', ''Tag'', ''^warg''), ''Enable'', temp{double(get(gcbo, ''Value'') == 6) + 1}), set(findobj(gcbf, ''Tag'', ''wargedit''), ''String'', '''')'} {} ...
              {'Style' 'text' 'String' 'Kaiser window beta:' 'Tag' 'wargtext' 'Enable' 'off'} ...
              {'Style' 'edit' 'String' '' 'Tag' 'wargedit' 'Enable' 'off'} ...
              {'Style' 'pushbutton' 'String' 'Estimate' 'Tag' 'wargpush' 'Enable' 'off' 'Callback' @comwarg} ...
              {'Style' 'text' 'String' 'Filter order (mandatory even):'} ...
              {'Style' 'edit' 'String' '' 'Tag' 'forderedit'} ...
              {'Style' 'pushbutton' 'String' 'Estimate' 'Callback' {@comforder, wtypes}} ...
              {'Style' 'edit' 'Tag' 'devedit' 'Visible' 'off'} ...
              {} {} {'Style' 'pushbutton' 'String', 'Plot filter responses' 'Callback' {@comfresp, wtypes, ftypes}}};
    result = inputgui(uigeom, uilist, 'pophelp(''pop_firws'')', 'Filter the data -- pop_firws()');
    if isempty(result), return; end

    Arg = struct;
    Arg.srate = str2double(result{1});
    Arg.fcutoff = str2num(result{2});
    Arg.ftype = ftypes{result{3}};
    Arg.wtype = wtypes{result{4}};
    Arg.warg = str2num(result{5});
    Arg.forder = str2double(result{6});

% Command line mode
else
    Arg = struct(varargin{:});
end

% Sampling rate
if ~isfield(Arg, 'srate') || isempty(Arg.srate) % Use default
    Arg.srate = 2;
end

% Filter order and cutoff frequencies
if ~isfield(Arg, 'fcutoff') || ~isfield(Arg, 'forder') || isempty(Arg.fcutoff) || isempty(Arg.forder)
    error('Not enough input arguments.');
end
firwsArgArray = {Arg.forder sort(Arg.fcutoff / Arg.srate * 2)}; % Sorting and normalization

% Filter type
if ~isfield(Arg, 'ftype') || isempty(Arg.ftype) % Use default
    switch length(Arg.fcutoff)
        case 1
            Arg.ftype = 'lowpass';
        case 2
            Arg.ftype = 'bandpass';
        otherwise
            error('Wrong number of arguments.')
    end
else
    if any(strcmpi(Arg.ftype, {'bandpass' 'bandstop'})) && length(Arg.fcutoff) ~= 2
        error('Not enough input arguments.');
    elseif any(strcmpi(Arg.ftype, {'highpass' 'lowpass'})) && length(Arg.fcutoff) ~= 1
        error('Too many input arguments.');
    end
    switch Arg.ftype
        case 'bandstop'
            firwsArgArray(end + 1) = {'stop'};
        case 'highpass'
            firwsArgArray(end + 1) = {'high'};
    end
end

% Window type
if ~isfield(Arg, 'wtype') || isempty(Arg.wtype) % Use default
    Arg.wtype = 'blackman';
end

% Window parameter
if ~isfield(Arg, 'warg') || isempty(Arg.warg)
    Arg.warg = [];
    firwsArgArray(end + 1) = {windows(Arg.wtype, Arg.forder + 1)};
else
    firwsArgArray(end + 1) = {windows(Arg.wtype, Arg.forder + 1, Arg.warg)};
end

b = firws(firwsArgArray{:});
a = 1;

if nargout == 0 || isfield(Arg, 'filename')
    
    % Open file
    if ~isfield(Arg, 'filename') || isempty(Arg.filename)
        [Arg.filename Arg.pathname] = uiputfile('*.fir', 'Save filter -- pop_xfirws');
    end
    if ~isfield(Arg, 'pathname') || isempty(Arg.pathname)
        Arg.pathname = '.';
    end
    [fid message] = fopen(fullfile(Arg.pathname, Arg.filename), 'w', 'l'); 
    if fid == -1
        error(message)
    end

    % Author
    fprintf(fid, '[author]\n');
    fprintf(fid, '%s\n\n', 'pop_xfirws 1.5.1');

    % FIR design
    fprintf(fid, '[fir design]\n');
    fprintf(fid, 'method  %s\n', 'fourier');
    fprintf(fid, 'type    %s\n', Arg.ftype);
    fprintf(fid, 'fsample %f\n', Arg.srate);
    fprintf(fid, 'length  %d\n', Arg.forder + 1);
    fprintf(fid, 'fcrit%d  %f\n', [1:length(Arg.fcutoff); Arg.fcutoff]);
    fprintf(fid, 'window  %s %s\n\n', Arg.wtype, num2str(Arg.warg)); % fprintf bug
    
    % FIR
    fprintf(fid, '[fir]\n');
    fprintf(fid, '%d\n', Arg.forder + 1);
    fprintf(fid, '% 18.10e\n', b);

    % Close file
    fclose(fid);

end

if nargout > 0
    varargout = {b a};
end

% Callback estimate Kaiser beta
function comwarg(varargin)
    [warg, dev] = pop_kaiserbeta;
    set(findobj(gcbf, 'Tag', 'wargedit'), 'String', warg);
    set(findobj(gcbf, 'Tag', 'devedit'), 'String', dev);

% Callback estimate filter order
function comforder(obj, evt, wtypes)
    srate = str2double(get(findobj(gcbf, 'Tag', 'srateedit'), 'String'));
    wtype = wtypes{get(findobj(gcbf, 'Tag', 'wtypepop'), 'Value')};
    dev = str2double(get(findobj(gcbf, 'Tag', 'devedit'), 'String'));
    [forder, dev] = pop_firwsord(wtype, srate, [], dev);
    set(findobj(gcbf, 'Tag', 'forderedit'), 'String', forder);
    set(findobj(gcbf, 'Tag', 'devedit'), 'String', dev);

% Callback plot filter responses
function comfresp(obj, evt, wtypes, ftypes)
    Arg.srate = str2double(get(findobj(gcbf, 'Tag', 'srateedit'), 'String'));
    Arg.fcutoff = str2num(get(findobj(gcbf, 'Tag', 'fcutoffedit'), 'String'));
    Arg.ftype = ftypes{get(findobj(gcbf, 'Tag', 'ftypepop'), 'Value')};
    Arg.wtype = wtypes{get(findobj(gcbf, 'Tag', 'wtypepop'), 'Value')};
    Arg.warg = str2num(get(findobj(gcbf, 'Tag', 'wargedit'), 'String'));
    Arg.forder = str2double(get(findobj(gcbf, 'Tag', 'forderedit'), 'String'));
    xfirwsArgArray(1, :) = fieldnames(Arg);
    xfirwsArgArray(2, :) = struct2cell(Arg);
    [b a] = pop_xfirws(xfirwsArgArray{:});
    H = findobj('Tag', 'filter responses', 'type', 'figure');
    if ~isempty(H)
        figure(H);
    else
        H = figure;
        set(H, 'color', [.93 .96 1], 'Tag', 'filter responses');
    end
    plotfresp(b, a, [], Arg.srate);

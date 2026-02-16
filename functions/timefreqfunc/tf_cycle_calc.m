% tf_cycle_calc() - Calculate Morlet wavelet cycles from temporal or spectral (band-)width
%
% Usage:
%   >> [cycles, widths_table] = tf_cycle_calc('freqs', freqs) % pop-up window mode
%   >> [cycles, widths_table] = tf_cycle_calc('freqs', freqs, ...
%       'width', width, 'width_unit', width_unit, 'key1', value1, ...
%       'key2', value2, 'keyn', valuen);
%
% Inputs:
%   'freqs'       - vector (or scalar) wavelet center frequency/ies (Hz)
%   'width'       - vector (or scalar) temporal or spectral wavelet
%                   (band-)width in 'width_unit'. If length freqs > 2 and
%                   length width == 2, width is interpolated in linear or
%                   log-spaced steps (see 'log_spaced'). If length
%                   width == 1, constant width is used for all freqs (i.e.,
%                   STFT)
%   'width_unit'  - char array wavelet width type/unit. 'fwhm_t', 'fwhm_f',
%                   '2_sigma_t', '2_sigma_f', 'sigma_t', 'sigma_f', or
%                   'cycles'
%
% Optional inputs:
%   'log_spaced   - scalar boolean linear or log spaced width steps if
%                   length freqs > 2 and length width == 2 {default false}
%   'plot'        - scalar boolean plot spectral and temporal (band-)width
%                   in FWHM and cycles {default false}
%   'data'        - column vector data to compute TF transform
%   'fs'          - scalar sampling frequency of data vector
%   'norm'        - scalar boolean normalize tf_data output to unit energy
%                   {default false}
%
% Outputs:
%   cycles        - vector (of length freqs) width in cycles unit
%   widths_table  - array with length freqs rows and the columns 'freq',
%                   'cycles', 'fwhm_f', 'fwhm_t', '2_sigma_f', '2_sigma_t',
%                   'sigma_f', and 'sigma_t'
%   tf_data       - time x freqs array TF transformed data
%   cmorlf        - frequency x freqs array frequency domain Morlet
%                   wavelets with center frequencies freqs and width width
%
% Note:
%   For more detailed documentation, see
%   https://github.com/widmann/tf_cycle_calc. The implemented simple TF
%   transform is intended for testing and demonstration purposes only, not
%   for productive use.
%
% Author: Andreas Widmann, University of Leipzig, 2025
%
% See also:
%   pop_newtimef

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2025 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function [cycles, widths_table, tf_data, cmorlf] = tf_cycle_calc(varargin)

% Constants
sigma2fwhm = 2 * sqrt(2 * log(2));

% Parse input arguments
Args = struct(varargin{:});
if ~isfield(Args, 'freqs')
    Args.freqs = [];
end
if ~isfield(Args, 'width')
    Args.width = [];
end

% GUI
if isempty(Args.freqs) || isempty(Args.width)

    widthStr = {'Temporal width FWHM [s]', 'Spectral bandwidth FWHM [Hz]', 'Temporal width 2 * sigma_t [s]', 'Spectral bandwidth 2 * sigma_f [Hz]', 'Temporal width sigma_t [s]', 'Spectral bandwidth sigma_f [Hz]', 'Cycles'};
    widthIdx = 1;
    widthArgs = {'fwhm_t', 'fwhm_f', '2_sigma_t', '2_sigma_f', 'sigma_t', 'sigma_f', 'cycles'};
    uigeom = {[3 2 1], [3 2 1], [3 2 1], [3 2 1]};

    uilist = {{'Style', 'text', 'String', 'Specify unit:'}, ...
              {'Style', 'popupmenu', 'String', widthStr, 'Tag', 'widthpop', 'Value', widthIdx}, {}, ...
              {'Style', 'text', 'String', 'Frequency limits [min max] (Hz) or sequence:'}, ...
              {'Style', 'edit', 'String', num2str(Args.freqs(:)'), 'Tag', 'freqedit'}, {}, ...
              {'Style', 'text', 'String', '(Band-)width limits [min max] (Hz/s/cycles) or sequence:'}, ...
              {'Style', 'edit', 'String', num2str(Args.width(:)'), 'Tag', 'widthedit'}, {'Style' 'checkbox', 'String', 'Log spacing', 'Tag' 'spacingcheck', 'Value', 0}, ...
              {}, {'Style', 'pushbutton', 'String', 'Plot (band-)width', 'Callback', {@complot, widthArgs}}, {}, ...
              };
    if exist('inputgui', 'file') ~= 2
        error('GUI mode requires EEGLAB.')
    end
    result = inputgui(uigeom, uilist, 'pophelp(''tf_cycle_calc'')', 'Wavelet cycles calculator -- tf_cycle_calc()');
    if isempty(result), cycles = []; widths_table = []; return; end

    Args.width_unit = widthArgs{result{1}};
    Args.freqs = str2num(result{2}); %#ok<ST2NM>
    Args.width = str2num(result{3}); %#ok<ST2NM>
    Args.log_spaced = result{4};
    
end

% Check input arguments
if isempty(Args.freqs)
    error('Frequencies vector input argument required.')
end
Args.freqs = Args.freqs(:);
nFreq = length(Args.freqs);

if ~isfield(Args, 'width_unit') || isempty(Args.width_unit)
    error('Width unit input argument is required.')
end
if isempty(Args.width)
    error('Width input argument is required.')
end

if ~isfield(Args, 'log_spaced') || isempty(Args.log_spaced)
    Args.log_spaced = 0;
end

% Convert requested width to vector if necessary
if isscalar(Args.width)
    Args.width = repmat(Args.width, [nFreq 1]);
elseif length(Args.width) == 2 && nFreq > 2
    if Args.log_spaced
        Args.width = logspace(log10(Args.width(1)), log10(Args.width(2)), nFreq);
    else
        Args.width = linspace(Args.width(1), Args.width(2), nFreq);
    end
elseif length(Args.width) == 2 && nFreq == 2 && Args.log_spaced
    warning('Log spaced option has no effect if length frequencies is 2')
elseif length(Args.width) > 2 && length(Args.width) ~= nFreq
    error('Length width must equal length freqs if > 2.')
end
Args.width = Args.width(:);

% Convert to sigma_f
switch Args.width_unit
    case 'fwhm_f'
        sigma_f = Args.width / sigma2fwhm;
    case 'fwhm_t'
        sigma_t = Args.width / sigma2fwhm;
        sigma_f = 1 ./ (2 * pi * sigma_t);
    case '2_sigma_f'
        sigma_f = Args.width / 2;
    case '2_sigma_t'
        sigma_t = Args.width / 2;
        sigma_f = 1 ./ (2 * pi * sigma_t);
    case 'sigma_f'
        sigma_f = Args.width;
    case 'sigma_t'
        sigma_t = Args.width;
        sigma_f = 1 ./ (2 * pi * sigma_t);
    case 'cycles'
        cycles = Args.width;
        sigma_f = Args.freqs ./ cycles;
    otherwise
        error('Unkown wavelet width parameter.')
end

% Convert sigma_f to other units
sigma_t = 1 ./ (2 * pi * sigma_f);
fwhm_f = sigma_f * sigma2fwhm;
fwhm_t = sigma_t * sigma2fwhm;
cycles = Args.freqs ./ sigma_f;

widths_table = [Args.freqs cycles fwhm_f fwhm_t 2 * sigma_f 2 * sigma_t sigma_f sigma_t];

% TF transform
if isfield(Args, 'data') && ~isempty(Args.data)

    if isvector(Args.data)
        Args.data = Args.data(:);
    else
        error('Data should be a vector.')
    end
    if ~isfield(Args, 'fs') || isempty(Args.fs)
        error('Requested sampling frequency input argument required.')
    end
    if ~isfield(Args, 'norm') || isempty(Args.norm)
        Args.norm = 1;
    end

    [tf_data, cmorlf] = my_tf(Args.data, Args.fs, Args.freqs', sigma_f', Args.norm);

end

% Plot if requested
if isfield(Args, 'plot') && Args.plot

    % widths_table

    hFig = findobj('tag', 'FigCycleCalc');
    if isempty(findobj('tag', 'FigCycleCalc')) || isfield(Args, 'data')
        hFig = figure;
        set(hFig, 'Tag', 'FigCycleCalc')
    else
        figure(hFig)
    end

    % Spectral bandwidth vs. frequency
    hAx(1) = subplot(2, 2, 1, 'Tag', 'AxCycleCalcBandwidth');
    plot(Args.freqs, fwhm_f, '.', 'MarkerSize', 15)
    xlabel('Frequency [Hz]'), ylabel('FWHM spectral bandwidth [Hz]')
    title('Spectral bandwidth')

    % Temporal width vs. frequency
    hAx(2) = subplot(2, 2, 2, 'Tag', 'AxCycleCalcWidth');
    plot(Args.freqs, fwhm_t, '.', 'MarkerSize', 15)
    xlabel('Frequency [Hz]'), ylabel('FWHM width [s]')
    title('Temporal width')

    % Cycles vs. frequency
    hAx(3) = subplot(2, 2, 3, 'Tag', 'AxCycleCalcCycles');
    plot(Args.freqs, cycles, '.', 'MarkerSize', 15)
    xlabel('Frequency [Hz]'), ylabel('Cycles')
    title('Cycles')

    hold(hAx, 'on')

    if isfield(Args, 'data') && ~isempty(Args.data)
        hAx(4) = subplot(2, 2, 4);
        colormap(turbo)
        % imagesc((0:length(Args.data) - 1) / Args.fs, Args.freqs, abs(tf_data').^2)
        imagesc((0:length(Args.data) - 1) / Args.fs, Args.freqs, 2 * abs(tf_data').^2) % Correct for onesided spectrum
        colorbar
        set(hAx(4), 'YDir', 'normal')
    end

    set(hAx, 'XGrid', 'on', 'YGrid', 'on')

end

end

function [tf_data, cmorlf] = my_tf(data, fs, freqs, sigma_f, norm)

if nargin < 5
    error('Not enough input arguments.')
end
N = length(data);
if any(1 ./ (2 * pi * sigma_f) * 6 > N / fs)
    warning('Data length too short. Minimum data length of %.3f s recommended.', max(1 ./ (2 * pi * sigma_f) * 6));
end

% Prepare arguments
freqs = freqs * (N / fs);
sigma_f = sigma_f * (N / fs);
sigma_f = repmat(sigma_f, [N 1]);
f = repmat(0:N - 1, [size( freqs, 2) 1])' - repmat(freqs, [N 1]);

% Normalization
if norm
    A = sqrt(N ./ (sqrt(pi) * sigma_f));
else
    A = 1;
end

% Frequency domain complex Morlet wavelet
cmorlf = A .* exp(-f .^ 2 ./ (2 * sigma_f .^ 2));

% TF transform
tf_data = ifft(repmat(fft(data), [1 length(freqs)]) .* cmorlf);

end

function complot(obj, evt, widthArgs) %#ok<INUSD>
    widthIdx = get(findobj(gcbf, 'Tag', 'widthpop'), 'Value');
    freqs = str2num(get(findobj(gcbf, 'Tag', 'freqedit'), 'String')); %#ok<ST2NM>
    if isempty(freqs)
        error('Frequencies vector input argument is required.')
    end
    width = str2num(get(findobj(gcbf, 'Tag', 'widthedit'), 'String')); %#ok<ST2NM>
    if isempty(width)
        error('Width input argument is required.')
    end
    log_spaced = get(findobj(gcbf, 'Tag', 'spacingcheck'), 'Value');
    tf_cycle_calc('freqs', freqs, 'width_unit', widthArgs{widthIdx}, 'width', width, 'log_spaced', log_spaced, 'plot', 1);
end
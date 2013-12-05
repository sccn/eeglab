% firfilt() - Pad data with DC constant, filter data with FIR filter,
%             and shift data by the filter's group delay
%
% Usage:
%   >> EEG = firfilt(EEG, b, nFrames);
%
% Inputs:
%   EEG           - EEGLAB EEG structure
%   b             - vector of filter coefficients
%
% Optional inputs:
%   nFrames       - number of frames to filter per block {default 1000}
%
% Outputs:
%   EEG           - EEGLAB EEG structure
%
% Note:
%   Higher values for nFrames increase speed and working memory
%   requirements.
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   filter, findboundaries

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

function EEG = firfilt(EEG, b, nFrames)

if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 3 || isempty(nFrames)
    nFrames = 1000;
end

% Filter's group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;

% Find data discontinuities and reshape epoched data
if EEG.trials > 1 % Epoched data
    EEG.data = reshape(EEG.data, [EEG.nbchan EEG.pnts * EEG.trials]);
    dcArray = 1 : EEG.pnts : EEG.pnts * (EEG.trials + 1);
else % Continuous data
    dcArray = [findboundaries(EEG.event) EEG.pnts + 1];
end

% Initialize progress indicator
nSteps = 20;
step = 0;
fprintf(1, 'firfilt(): |');
strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
tic

for iDc = 1:(length(dcArray) - 1)

        % Pad beginning of data with DC constant and get initial conditions
        ziDataDur = min(groupDelay, dcArray(iDc + 1) - dcArray(iDc));
        [temp, zi] = filter(b, 1, double([EEG.data(:, ones(1, groupDelay) * dcArray(iDc)) ...
                                  EEG.data(:, dcArray(iDc):(dcArray(iDc) + ziDataDur - 1))]), [], 2);

        blockArray = [(dcArray(iDc) + groupDelay):nFrames:(dcArray(iDc + 1) - 1) dcArray(iDc + 1)];
        for iBlock = 1:(length(blockArray) - 1)

            % Filter the data
            [EEG.data(:, (blockArray(iBlock) - groupDelay):(blockArray(iBlock + 1) - groupDelay - 1)), zi] = ...
                filter(b, 1, double(EEG.data(:, blockArray(iBlock):(blockArray(iBlock + 1) - 1))), zi, 2);

            % Update progress indicator
            [step, strLength] = mywaitbar((blockArray(iBlock + 1) - groupDelay - 1), size(EEG.data, 2), step, nSteps, strLength);
        end

        % Pad end of data with DC constant
        temp = filter(b, 1, double(EEG.data(:, ones(1, groupDelay) * (dcArray(iDc + 1) - 1))), zi, 2);
        EEG.data(:, (dcArray(iDc + 1) - ziDataDur):(dcArray(iDc + 1) - 1)) = ...
            temp(:, (end - ziDataDur + 1):end);

        % Update progress indicator
        [step, strLength] = mywaitbar((dcArray(iDc + 1) - 1), size(EEG.data, 2), step, nSteps, strLength);

end

% Reshape epoched data
if EEG.trials > 1
    EEG.data = reshape(EEG.data, [EEG.nbchan EEG.pnts EEG.trials]);
end

% Deinitialize progress indicator
fprintf(1, '\n')

end

function [step, strLength] = mywaitbar(compl, total, step, nSteps, strLength)

progStrArray = '/-\|';
tmp = floor(compl / total * nSteps);
if tmp > step
    fprintf(1, [repmat('\b', 1, strLength) '%s'], repmat('=', 1, tmp - step))
    step = tmp;
    ete = ceil(toc / step * (nSteps - step));
    strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '%s %3d%%, ETE %02d:%02d'], progStrArray(mod(step - 1, 4) + 1), floor(step * 100 / nSteps), floor(ete / 60), mod(ete, 60));
end

end

function EEGfiltered = eeg_filter(EEGinput,sample_freq,lcf,hcf,order);

% eeg_filter - apply a butterworth polynomial filter
% 
% 	Usage : EEGfiltered = eeg_filter(EEGinput,sample_freq,lcf,hcf,order)
%
% 	- input arguments 
%		EEGinput    : eeg data - M samples x N channels x P epochs
%		sample_freq : sampling frequency
%		lcf         : low cutoff frequency (highpass, default 0.01)
%		hcf	        : high cutoff frequency (lowpass, default 40)
%		order       : butterworth polynomial order (default 2)
%
% 	- output argument
%		EEGfiltered : filtered EEGinput;
%
%  This function uses the filtfilt and butter functions of the matlab
%  signal processing toolbox, using the bandpass option.  The filtfilt
%  function baselines on the last point of the timeseries, so it is
%  necessary to call *eeg_baseline* after filtering.
%

% Copyright (C) 2004  Darren L. Weber
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $
% Created: 10/2004, copyright 2004 Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '[$Revision: 1.1 $]';
fprintf('\nEEG_FILTER [v%s]\n',version(12:16));

if ~exist('sample_freq','var') || isempty(sample_freq),
    error('no sample_freq defined');
end

% low cutoff frequency (default 2)
if ~exist('lcf','var') || isempty(lcf),
    lcf = 0.01;
end
% high cutoff frequency (default 20)
if ~exist('hcf','var') || isempty(hcf),
    hcf = 40;
end
% butter filter order (default 2)
if ~exist('order','var') || isempty(order),
    order = 2;
end

if hcf > (sample_freq/2),
    warning('hcf > sample_freq/2, setting hcf = sample_freq/2');
    hcf = sample_freq / 2;
end
if lcf <= 0 || lcf > (sample_freq/2) || lcf >= hcf,
    warning('lcf value is <=0 or >(sample_freq/2) or >=hcf, setting lcf = 2');
    lcf = 2;
end

% call the butter function of the signal processing toolbox
cf1 = lcf/(sample_freq/2);
cf2 = hcf/(sample_freq/2);
[B,A] = butter(order,[cf1 cf2]);

epochs = size(EEGinput,3);
for epoch = 1:epochs,
    EEGfiltered(:,:,epoch) = filtfilt(B,A,EEGinput(:,:,epoch));
end

return

function [ EEGbase ] = eeg_baseline(EEGdata,baseline_samples)

% eeg_baseline - remove mean of baseline period from all EEG epochs
%
% [ EEGbase ] = eeg_baseline(EEGdata,baseline_samples)
%
% It is assumed that EEGdata is an MxNxP matrix, with M data samples (ie,
% time), N channels and P epochs (if P=1, EEGdata is MxN).
%
% baseline_samples is an array of data sample indices (of M) that define
% the baseline period.  If baseline_samples is a scalar, the function
% assumes the baseline index array is 1:baseline_samples.
% 
% The function calculates the baseline mean for each epoch as follows:
%
% baseline_mean = mean(EEGdata(baseline_samples,:,epoch));
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

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:50 $
% Created:  10/2004, copyright 2004 Darren.Weber_at_radiology.ucsf.edu
% Modified: 05/2005, Darren.Weber_at_radiology.ucsf.edu
%                    baseline_samples can now be an array or a scalar input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


version = '[$Revision: 1.1 $]';
fprintf('\nEEG_BASELINE [v%s]\n',version(12:16)); tic

if ~exist('baseline_samples','var'),
    error('no baseline_samples specified');
end
if isempty(baseline_samples),
    error('no baseline_samples specified');
end
if length(baseline_samples) == 1,
    baseline_samples = 1:baseline_samples;
end

number_samples = size(EEGdata,1);

for epoch = 1:size(EEGdata,3),

    baseline_mean = mean(EEGdata(baseline_samples,:,epoch));

    baseline_mean = repmat(baseline_mean,number_samples,1);

    EEGbase(:,:,epoch) = EEGdata(:,:,epoch) - baseline_mean;

end

t = toc; fprintf('...done (%6.2f sec).\n\n',t);
return

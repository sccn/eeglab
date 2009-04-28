function [regions] = eeg_region_peaks(p,times,regions,missing)

% eeg_region_peaks - Find peaks in regions of electrodes
%
% Usage: [regions] = eeg_region_peaks(p,times,regions,missing)
%
% requires p.volt.data and p.volt.timeArray
%
% 'times' is 1 row x 3 columns, with values
% for start, peak and end time windows (in that order)
% in the scale of p.volt.timeArray
%
% 'regions' is optional.  When omitted or an empty matrix, 
% it is created by elec_regions.m, see that function for 
% more information.
%
% 'missing' is a boolean switch.  If missing = 0, then
% the maximal region value at the peak time (ie, times(r,2))
% is returned. This avoids missing values in the return 
% data (default).
%
% returns 'regions' - an array of structures with
% fields: name, elec, pospeaks, postimes, negpeaks &
% negtimes.  The latter are arrays of the max/min peak
% values at the time windows of 'times'.
%
% Depends on eeg_peaks & sometimes elec_regions.  Called
% in, for example, eeg_region_peaks_script.m
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:53 $

% Licence: GNU GPL, no express or implied warranties
% History: 03/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check p structure fields
if ~isfield(p.volt,'timeArray'),
    msg = sprintf('EEG_REGION_PEAKS: p.volt.timeArray does not exist.\n');
    error(msg);
end
if isempty(p.volt.timeArray),
    msg = sprintf('EEG_REGION_PEAKS: p.volt.timeArray is empty.\n');
    error(msg);
end
% check for time windows of analysis
if ~exist('times','var'),
    % Set times for beginning and end of volt.timeArray
    times = [p.volt.timeArray(1,1) 1 p.volt.timeArray(end,1)];
    missing = 1;
end
if isempty(times),
    times = [p.volt.timeArray(1,1) 1 p.volt.timeArray(end,1)];
    missing = 1;
end
% check for regions of analysis
if ~exist('regions','var'),
    regions = elec_regions;
end
if isempty(regions),
    regions = elec_regions;
end
% check for missing value instruction
if ~exist('missing','var'),
    missing = 0;
end
if isempty(missing),
    missing = 0;
end


% first get all peaks for input p.volt.data
p = eeg_peaks(p);

tic
fprintf('EEG_REGION_PEAKS: Finding regional peaks...');

for r=1:length(regions),
    
    % the regions.elec field can be empty, especially
    % when processing a control condition of a 2 condition
    % experiment where no positive or negative peak is
    % found in the experimental condition.  See
    % eeg_region_peaks_script.m
    if isnan(regions(r).elec),
        regions(r).pospeaks = NaN;
        regions(r).postimes = NaN;
        regions(r).poselecs = NaN;
        regions(r).negpeaks = NaN;
        regions(r).negtimes = NaN;
        regions(r).negelecs = NaN;
        continue;
    end
    
    % Select region voltage/time arrays
    region.data = p.volt.data(:,regions(r).elec);
    region.time = p.volt.timeArray(:,regions(r).elec);
    
    % select peak indices between start/end times
    peaks = p.volt.peaks.data(:,regions(r).elec);
    peakdata = peaks .* and(region.time >= times(1), region.time <= times(3));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select all positive region peaks between start/end times
    if isempty(find(peakdata>0)),
        if missing,
            % Create missing data
            regions(r).pospeaks = NaN;
            regions(r).postimes = NaN;
            regions(r).poselecs = NaN;
        else
            % Select region data at times(2)
            pospeaks.data = region.data .* (region.time==times(2));
            pospeaks.time = region.time .* (region.time==times(2));
            
            % find max positive peak value, timing & electrode
            [i,j,v] = find(pospeaks.data); % non-zero elements
            regions(r).pospeaks = max(v);
            index = find(v == max(v));
            if ~isequal(size(index),[1 1]),
                fprintf('\nEEG_REGION_PEAKS: Warning: Selecting first of two equal max values.\n');
            end
            regions(r).postimes = pospeaks.time(i(index(1,1)),j(index(1,1)));
            regions(r).poselecs = regions(r).elec(j(index(1,1)));
        end
    else
        pospeaks.data = region.data .* (peakdata > 0);
        pospeaks.time = region.time .* (peakdata > 0);
        
        % find max positive peak value, timing & electrode
        [i,j,v] = find(pospeaks.data); % non-zero elements
        regions(r).pospeaks = max(v);
        index = find(v == max(v));
        if ~isequal(size(index),[1 1]),
            fprintf('\nEEG_REGION_PEAKS: Warning: Selecting first of two equal max values.\n');
        end
        regions(r).postimes = pospeaks.time(i(index(1,1)),j(index(1,1)));
        regions(r).poselecs = regions(r).elec(j(index(1,1)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select all negative region peaks between start/end times
    if isempty(find(peakdata<0)),
        if missing,
            % Create missing data
            regions(r).negpeaks = NaN;
            regions(r).negtimes = NaN;
            regions(r).negelecs = NaN;
        else
            % Select region data at times(2)
            negpeaks.data = region.data .* (region.time==times(2));
            negpeaks.time = region.time .* (region.time==times(2));
            
            % find min negative peak value, timing & electrode
            [i,j,v] = find(negpeaks.data); % non-zero elements
            regions(r).negpeaks = min(v);
            index = find(v == min(v));
            if ~isequal(size(index),[1 1]),
                fprintf('\nEEG_REGION_PEAKS: Warning: Selecting first of two equal min values.\n');
            end
            regions(r).negtimes = negpeaks.time(i(index(1,1)),j(index(1,1)));
            regions(r).negelecs = regions(r).elec(j(index(1,1)));
        end
    else
        negpeaks.data = region.data .* (peakdata < 0);
        negpeaks.time = region.time .* (peakdata < 0);
        
        % find min negative peak value, timing & electrode
        [i,j,v] = find(negpeaks.data); % non-zero elements
        regions(r).negpeaks = min(v);
        index = find(v == min(v));
        if ~isequal(size(index),[1 1]),
            fprintf('\nEEG_REGION_PEAKS: Warning: Selecting first of two equal min values.\n');
        end
        regions(r).negtimes = negpeaks.time(i(index(1,1)),j(index(1,1)));
        regions(r).negelecs = regions(r).elec(j(index(1,1)));
    end
end

t = toc;
fprintf('done (%5.2f sec)\n',t);

return

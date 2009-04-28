function [p] = eeg_peaks(p)

% eeg_peaks - Find peaks in EEG data
%
% Usage: [p] = eeg_peaks(p)
%
% The p structure is explained in eeg_toolbox_defaults.
% For this function, it must contain the field called
% 'volt.data,' which is assumed to be a matrix 
% with N rows of EEG sample values from M columns 
% of electrodes.
% 
% The indices to peak values are returned to
% p.volt.peaks.data (which is size(p.volt.data)).
% The values are:
%                    0   no peak
%                    1   a +ve peak
%                   -1   a -ve peak
%
% This function finds all peaks in a waveform, where
% f(x) changes from an increasing to a decreasing 
% function and vice versa). These are local minima
% and maxima that have a first derivate of zero and
% satisfy the first derivative test.
% 
% Thus, the +ve/-ve peak refers to peaks and troughs 
% in a waveform, regardless of the +ve/-ve value of 
% the waveform at that point. The peak values are 
% commonly +ve/-ve potentials, but they may not be.
% A +ve peak can have a -ve value where the change 
% from increasing f(x) to decreasing f(x) occurs 
% entirely in the range of negative values.
% 
% All the +ve/-ve peak values are returned in:
%
%   p.volt.peaks.all
%   p.volt.peaks.pos
%   p.volt.peaks.neg
%
% If the p.volt.timeArray is defined, then the timing
% of the peaks is returned in:
%
%   p.volt.peaks.alltimes
%   p.volt.peaks.postimes
%   p.volt.peaks.negtimes
%
% eg,   p.volt.timeArray = linspace(-20,20,200)';
%       p.volt.data = sin(linspace(-20,20,200))';
%       p.volt.data(:,2) = sin(-1 .* linspace(-20,20,200))';
%      [p] = eeg_peaks(p);
%       scatter(p.volt.timeArray(p.volt.peaks.all),p.volt.data(p.volt.peaks.all))
%       hold on
%       plot(p.volt.timeArray,p.volt.data)
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  www.gnu.org GPL, no implied or express warranties
% Created:  02/2002, Darren.Weber_at_radiology.ucsf.edu
%                    - adapted a perl foreach loop to peaks_matrix
%                      to take advantage of matlab speed.
%           05/2002, Darren.Weber_at_radiology.ucsf.edu
%                    - used matlab 'diff' command
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('EEG_PEAKS\n...calculating peaks...'); tic;

p.volt.peaks.data = peaks_matrix(p.volt.data);

% Find all the peak values
p.volt.peaks.all = p.volt.data .* (p.volt.peaks.data ~= 0);
p.volt.peaks.pos = p.volt.data .* (p.volt.peaks.data  > 0);
p.volt.peaks.neg = p.volt.data .* (p.volt.peaks.data  < 0);

if ~isfield(p.volt,'timeArray'),
    fprintf('\n...p doesn''t contain ''volt.timeArray''.\n');
    return
elseif isempty(p.volt.timeArray),
    fprintf('\n...p.volt.timeArray is empty.\n');
    return
end

% make sure that volt.timeArray is same size as volt.data
if ~isequal(size(p.volt.timeArray),size(p.volt.data)),
    p.volt.timeArray = repmat(p.volt.timeArray(:,1),1,size(p.volt.data,2));
end
% Find the peak timing
p.volt.peaks.alltimes = p.volt.timeArray .* (p.volt.peaks.data ~= 0);
p.volt.peaks.postimes = p.volt.timeArray .* (p.volt.peaks.data  > 0);
p.volt.peaks.negtimes = p.volt.timeArray .* (p.volt.peaks.data  < 0);

t = toc;
fprintf('done (%5.2f sec)\n',t);
return



function [peaks] = peaks_matrix(data)
    
    % This function should be quicker than that below.
    % The code below remains as it has been proved to
    % be reliable, but this code effectively replaces it.
    % Both code sets produce the same results.
    
    [r,c] = size(data);
    % replicate data, but staggered one point
    data2 = zeros(1,c);
    data2(2:r+1,:) = data;
    % add zero values to end of data points
    data(r+1,:) = zeros(1,c);
    % generate boolean difference between two matrices
    dif = data < data2;
    
    % Now replicate this boolean matrix
    [br,bc] = size(dif);
    dif2 = zeros(1,bc);
    dif2(2:br+1,:) = dif;
    % add zeros to end of boolean dif
    dif(br+1,:) = 0;
    
    % Calculate the peaks
    peaks = dif(2:r+1,:) - dif2(2:r+1,:);
    
    % now trim peaks dif back to size of data
    peaks(1,:) = 0; % first point cannot be a peak
    peaks(r,:) = 0; % last point cannot be a peak
    
    
    
    % Could use matlab diff command, like this
    %[r,c] = size(data);
    %decdat = (diff(data)<=0) .* -1;
    %incdat =  diff(data)>=0;
    
    %zrow = zeros(1,c);
    %dirdat1 = [ zrow; decdat + incdat ];
    %dirdat2 = [ decdat + incdat; zrow ];
    
    %negpeaks = ((dirdat1 - dirdat2) < 0) .* -1;
    %pospeaks =  (dirdat1 - dirdat2) > 0;
    %peaks = negpeaks + pospeaks;
    
    %peaks(1,:)   = 0; % first point cannot be a peak
    %peaks(end,:) = 0; % last point cannot be a peak
    
    
return


% The function 'peaks_matrix' above replaced the 
% following for loop (02/2002, Darren Weber)

% foreach electrode
%for e = 1:size(p.volt.data,2),
%	v = p.volt.data(:,e);
%    n = 1;
%    volt = v(n);
%	while (n <= size(v,1)),
%        if(volt < v(n)),
%            volt = v(n); n = n + 1;
%            if(volt < v(n)), p.volt.peaks.data(n-1,e) = 0;
%            else             p.volt.peaks.data(n-1,e) = 1;
%            end
%        elseif(volt > v(n)),
%            volt = v(n); n = n + 1;
%            if(volt > v(n)), p.volt.peaks.data(n-1,e) = 0;
%            else             p.volt.peaks.data(n-1,e) = -1;
%            end
%        else
%            volt = v(n); n = n + 1;
%           p.volt.peaks.data(n-1,e) = 0;
%        end
%        if isequal(n,size(v,1)),
%            break;
%        end
%    end
%    p.volt.peaks.data(1,e) = 0; % first point cannot be a peak
%    p.volt.peaks.data(n,e) = 0; % last point cannot be a peak
%end
%p.volt.peaks.all = find(p.volt.peaks.data ~= 0);
%p.volt.peaks.pos = find(p.volt.peaks.data  > 0);
%p.volt.peaks.neg = find(p.volt.peaks.data  < 0);
%return

% timewarp()   - Given two event marker vectors, computes a matrix
%                that can be used to warp a time series so that its
%                evlatencies match newlatencies. Values of the warped
%                timeserie that falls between two frames in the original
%                timeserie will be linear interpolated.
% Usage:
%   >> warpmat = timewarp(evlatency, newlatency)
%
% Necessary inputs:
%   evlatency  - [vector] event markers in the original time series, in frames
%                Markers must be ordered by increasing frame latency. 
%                If you want to warp the entire time-series, make sure 
%                the first (1) and last frames are in the vector.
%   newlatency - [vector] desired warped event time latencies. The original
%                time series will be warped so that the frame numbers of its 
%                events (see evlatency above) match the frame numbers in 
%                newlatency. newlatency frames must be sorted by ascending 
%                latency. Both vectors must be the same length.
%   
% Optional outputs:
%      warpmat - [matrix] Multiplying this matrix with the original
%                time series (column) vectors performs the warping.
%
% Example:
%      % In 10-frame vectors, warp frames 3 and 5 to frames 4 and 8,
%      % respectively. Generate a matrix to warp data epochs in 
%      % variable 'data' of size (10,k)
%      >> warpmat = timewarp([1 3 5  10], [1 4 8 10])
%      >> warped_data = warpmat*data;
%
% Authors: Jean Hausser, SCCN/INC/UCSD, 2006
%
% See also: angtimewarp(), phasecoher(), erpimage()
%
function M=timewarp(evLatency, newLatency)
  M = [0];
  if min(sort(evLatency) == evLatency) == 0
    error('evLatency should be in ascending order');
    return;
  end
  if min(sort(newLatency) == newLatency) == 0
    error('newLatency should be in ascending order');
    return;
  end
  if length(evLatency) ~= length(newLatency)
    error('evLatency and newLatency must have the same length.');
    return;
  end
  if length(evLatency) < 2 | length(newLatency) < 2
    error(['There should be at least two events in evlatency and ' ...
          'newlatency (e.g., "begin" and "end")'] );
    return;
  end
  if evLatency(1) ~= 1
    disp(['Assuming old and new time series beginnings are synchronized.']);
    disp(['Make sure you have defined an ending event in both the old and new time series!']);
    evLatency(end+1)=1;
    newLatency(end+1)=1;
    evLatency = sort(evLatency);
    newLatency = sort(newLatency);
  end
    
  t = 1:max(evLatency);
  
  for k=1:length(evLatency)-1
    for i=evLatency(k):evLatency(k+1)-1
      tp(i) = (t(i)-evLatency(k)) * ...
              (newLatency(k+1) - newLatency(k))/...
              (evLatency(k+1) - evLatency(k)) + ...
              newLatency(k);
    end
  end
  
  
  % Check what's going on at tp(max(newLatency)), should equal t(max(evLatency))
  tp(max(evLatency)) = max(newLatency);
  ts = tp-min(newLatency)+1;
  
% $$$   M = sparse(max(newLatency)-min(newLatency)+1, max(evLatency));
  M = zeros(max(newLatency)-min(newLatency)+1, max(evLatency));
  
  k = 0;
  for i=1:size(M,1)
    while i > ts(k+1)
      k = k+1;
    end
% $$$     k = k-1;
    
    if k == 0
      % Check wether i == ts(1) and i == 1
      % In that case, M(1,1) = 1
      M(1,1) = 1;
    else
      M(i,k) = 1 - (i-ts(k))/(ts(k+1)-ts(k));
      M(i,k+1) = 1 - (ts(k+1)-i)/(ts(k+1)-ts(k));
    end
  end

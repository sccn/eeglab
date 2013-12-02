% angtimewarp() - Given two event marker vectors, computes a
%                 warping of the input angular time series so that its
%                 evlatencies match newlatencies. Values of the warped
%                 timeserie that falls between two frames in the original
%                 timeserie will be linearly interpolated under the
%                 assumption that phase change is minimal between two
%                 successive time points.
% Usage:
%   >> warpAngs = angtimewarp(evlatency, newlatency, angData)
%
% Necessary inputs:
%   evlatency  - [vector] time markers on the original time-series, in
%                frames. Markers must be ordered by increasing
%                latency. If you want to warp the entire time series, 
%                make sure frame 1 and the last frame are in the vector.
%   newlatency - [vector] desired time marker latencies. The original
%                time series will be warped so that its time markers (see
%                evlatency) match the ones in newlatency. newlatency
%                frames must be sorted by ascending latencies in frames.
%                Both vectors have to be the same length.
%   angData    - [vector] original angular time series (in radians). 
%                Angles should be between -pi and pi.
%   
% Optional outputs:
%   warpAngs   - [vector] warped angular time-course, with values between
%                -pi and pi
%
% Example:
%   >> angs = 2*pi*rand(1,10)-pi;
%   >> warpangs = angtimewarp([1 5 10], [1 6 10], angs)
%
% Authors: Jean Hausser, SCCN/INC/UCSD, 2006
%
% See also: timeWarp(), phasecoher(), erpimage(), newtimef()
%

function angdataw=angtimewarp(evLatency, newLatency, angdata)
  
  if min(sort(evLatency) == evLatency) == 0
    error('evlatency should be sorted');
    return;
  end
  if min(sort(newLatency) == newLatency) == 0
    error('newlatency should be sorted');
    return;
  end
  if length(evLatency) ~= length(newLatency)
    error('evlatency and newlatency must have the same length.');
    return;
  end
  if length(evLatency) < 2 | length(newLatency) < 2
    error(['There should be at least two events in evlatency and ' ...
          'newlatency, that is "begin" and "end"' ]);
    return;
  end
  if evLatency(1) ~= 1
    disp(['Assuming old and new time series beginnings are ' ...
          'synchronized']);
    disp(['Make sure you defined an end event for both old and new time ' ...
          'series !']);
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
  
  
  %Check what's going on at tp(max(newLatency)), should equal t(max(evLatency))
% $$$   keyboard;
  tp(max(evLatency)) = max(newLatency);
  ts = tp-min(newLatency)+1;
  
% $$$   M = sparse(max(newLatency)-min(newLatency)+1, max(evLatency));
% $$$   M = zeros(max(newLatency)-min(newLatency)+1, max(evLatency));
  angdataw = zeros(1, max(newLatency)-min(newLatency)+1);
  
  k = 0;
  for i=1:length(angdataw)
    while i > ts(k+1)
      k = k+1;
    end
% $$$     k = k-1;
    
% $$$     keyboard;
    if k == 0
      % Check wether i == ts(1) and i == 1
      % In that case, M(1,1) = 1
% $$$       keyboard;
% $$$       M(1,1) = 1;
      angdataw(1) = angdata(1);
    else
% $$$       M(i,k) = 1 - (i-ts(k))/(ts(k+1)-ts(k));
% $$$       M(i,k+1) = 1 - (ts(k+1)-i)/(ts(k+1)-ts(k));
      
      %Linear interp
      angdataw(i) = angdata(k)*(1 - (i-ts(k))/(ts(k+1)-ts(k))) + ...
          angdata(k+1)*(1 - (ts(k+1)-i)/(ts(k+1)-ts(k)));
      
      %Correction because angles have a ring structure
      theta1 = [angdata(k) angdata(k+1) angdataw(i)];
      theta2 = theta1 - min(angdata(k), angdata(k+1));
      theta2max = max(theta2(1), theta2(2));
% $$$       if i == 13
% $$$         keyboard; %so we get a chance to check wether 0<theta2(3)<2pi
% $$$       end
      if ~ ( (theta2max <= pi & theta2(3) <= theta2max) | ...
             (theta2max >= pi & theta2(3) >= theta2max) | ...
             theta2(3) == theta2(1) | theta2(3) == theta2(2) )
        angdataw(i) = angdataw(i) + pi;
      end
      if angdataw(i) > pi %Make sure we're still on [-pi, pi]
        angdataw(i) = angdataw(i) - 2*pi;
      end
    end
  end

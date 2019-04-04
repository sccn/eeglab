function outvec = gauss(frames,sds)
% gauss() - return a smooth Gaussian window
%
% Usage:
%   >> outvector = gauss(frames,sds);
%
% Inputs:
%   frames = window length
%   sds    = number of +/-std. deviations = steepness 
%            (~0+ -> flat; >>10 -> spike)

outvec = [];
if nargin < 2
  help gauss
  return
end
if sds <=0 || frames < 1
  help gauss
  return
end

incr = 2*sds/(frames-1);
outvec = exp(-(-sds:incr:sds).^2);

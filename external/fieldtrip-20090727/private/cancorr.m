function [px, py, wx, wy, powx, powy] = cancorr(C,x,y,flag,trunc)

% CANCORR computes the canonical correlation between multiple variables
%
% Canonical correlation analysis (CCA) is a way of measuring the linear
% relationship between two multidimensional variables. It finds two bases,
% one for each variable, that are optimal with respect to correlations and,
% at the same time, it finds the corresponding correlations.
%
% Use as
%   [px, py, wx, wy] = cancorr(C,x,y)

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.2  2007/07/17 09:55:41  jansch
% added possibility to truncate the values on the diagonal, and to output the
% rotated powers
%
% Revision 1.1  2005/05/09 14:17:36  roboos
% new implementation, used in freqdescriptives
%

if nargin<3,
  error('to compute canonical coherence, you need to specify the indices to the cross-spectral density making up the dependent and independent variable');
elseif nargin==3,
  flag  = 0;
  trunc = 0;
elseif nargin==4,
  trunc = 0;
end

Cxx = C(x,x);
Cyy = C(y,y);
Cxy = C(x,y);
Cyx = C(y,x);

[wx, px] = eig(pinv(Cxx)*Cxy*pinv(Cyy)*Cyx);
[wy, py] = eig(pinv(Cyy)*Cyx*pinv(Cxx)*Cxy);

[srtx,indx] = sort(diag(px), 'descend');
[srty,indy] = sort(diag(py), 'descend');

px = px(indx,indx);
wx = wx(:,   indx);
py = py(indy,indy);
wy = wy(:,   indy);

if flag,
  powx = wx'*Cxx*wx;
  powy = wy'*Cyy*wy;
end

if trunc>0 && trunc<1,
  px = px.*double(px>trunc);
  py = py.*double(py>trunc);
elseif trunc>=1,
  px(trunc+1:end,trunc+1:end) = 0;
  py(trunc+1:end,trunc+1:end) = 0;
end

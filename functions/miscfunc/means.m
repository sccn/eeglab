function [M,stderr,V,grpids] = means(X,grps)
% MEANS:  Means, standard errors and variances.  For column vectors, means(x) 
%         returns the mean value.  For matrices or row vectors, means(x) is a 
%         row vector containing the mean value of each column.  The basic 
%         difference from the Matlab functions mean() and var() is for a row vector, 
%         where means() returns the row vector instead of the mean value of the 
%         elements of the row.  Also allows for missing data, passed as NaNs.
%
%         If an optional grouping vector is supplied, returns a vector of means
%         for each group in collating sequence.
%
%     Usage: [M,stderr,V,grpids] = means(X,{grps})
%
%         X =       [n x p] data matrix.
%         grps =    optional [n x 1] grouping vector for k groups.
%         -----------------------------------------------------------------------
%         M =       [k x p] matrix of group means.
%         stderr =  corresponding [k x 1] vector of standard errors of the means.
%         V =       corresponding vector of variances.
%         grpids =  corresponding vector of group identifiers, for multiple 
%                     groups.
%

% RE Strauss, 2/2/99
%   12/26/99 - added standard errors.
%   10/19/00 - sort group-identifiers into collating sequence.
%   12/21/00 - corrected problem with uninitialized 'fg' for single group.
%    2/24/02 - added estimation of variances.
%   10/15/02 - corrected documentation.

  if nargin == 0, help means; return; end

  if (nargin < 2) grps = []; end

  [n,p] = size(X);

  if (isempty(grps))
    grps = ones(n,1);
    ug = 1;
    fg = n;
    k = 1;
  else
    [ug,fg] = uniquef(grps,1);  
    k = length(ug);
  end

  grpids = ug;
  M = zeros(k,p);
  stderr = zeros(k,p);

  for ik = 1:k
    ir = find(grps==ug(ik));
    for c = 1:p
      x = X(ir,c);
      ic = find(isfinite(x));
      if (isempty(ic))
        M(ik,c) = NaN;
        stderr(ik,c) = NaN;
      else
        M(ik,c) = mean(x(ic));
        s = std(x(ic));
        stderr(ik,c) = s./sqrt(fg(ik));
        V(ik,c) = s.^2;
      end
    end
  end

  return;

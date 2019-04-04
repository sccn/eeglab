function [value,freq,index] = uniquef(grp,sortflag)
% UNIQUEF: Given a matrix containing group labels, returns a vector containing 
%         a list of unique group labels, in the sequence found, and a 
%         vector of corresponding frequencies of each group.  
%         Optionally sorts the indices into ascending sequence.
%
%         Note: it might be necessary to truncate ('de-fuzz') real numbers to 
%         some arbitrary number of decimal positions (see TRUNCATE) before 
%         finding unique values.
%
%     Syntax: [value,freq,index] = uniquef(grp,sortflag)
%
%         grp -      matrix of a set of labels.
%         sortflag - boolean flag indicating that list of labels, and 
%                    corresponding frequencies, are to be so sorted 
%                    [default=0, =FALSE].
%         ------------------------------------------------------------
%         value -    column vector of unique labels.
%         freq -     corresponding absolute frequencies.
%         index -    indices of the first observation having each value.
%

% RE Strauss, 6/5/95
%   6/29/98 - modified to return indices.
%   1/25/00 - changed name from unique to uniquef to avoid conflict with 
%             Matlab v5 function.

  if (nargin < 2) sortflag = []; end

  get_index = 0;
  if (nargout > 2)
    get_index = 1;
  end

  if (isempty(sortflag))
    sortflag = 0;
  end

  tol = eps * 10.^4;
  grp = grp(:);                           % Convert input matrix to vector

  if (get_index)                          % Create vector of indices
    ind = [1:length(grp)]';
  end

  if (any([~isfinite(grp)]))                % Remove NaN's and infinite values
    i = find(~isfinite(grp));               %   from input vector and index vector
    grp(i) = [];
    if (get_index)                     
      ind(i) = [];
    end
  end

  value = [];
  freq = [];

  for i = 1:length(grp)                   % For each element of grp,
    b = (abs(value-grp(i)) < tol);        %   check if already in value list
    if (sum(b) > 0)                       % If so,
      freq(b) = freq(b) + 1;              %   increment frequency counter
    else                                  % If not,
      value = [value; grp(i)];            %   add to value list
      freq =  [freq; 1];                  %   and initialize frequency counter
    end
  end

  if (sortflag)
    [value,i] = sort(value);
    freq = freq(i);
  end

  if (get_index)
    nval = length(value);                 % Number of unique values
    index = zeros(nval,1);                % Allocate vector of indices
    for v = 1:nval                        % For each unique value,
      i = find(grp == value(v));          %   Find observations having value
      index(v) = ind(i(1));               %   Save first
    end
  end

  return;


function dists = eucl(crds1,crds2)
% EUCL: Calculates the euclidean distances among a set of points, or between a
%       reference point and a set of points, or among all possible pairs of two 
%       sets of points, in P dimensions.  Returns a single distance for two points.
%
%     Syntax: dists = eucl(crds1,crds2)
%
%        crds1 = [N1 x P] matrix of point coordinates.  If N=1, it is taken to
%                   be the reference point.
%        crds2 = [N2 x P] matrix of point coordinates.  If N=1, it is taken to
%                   be the reference point.
%        -----------------------------------------------------------------------
%        dists = [N1 x N1] symmetric matrix of pairwise distances (if only crds1
%                   is specified);
%                [N1 x 1]  col vector of euclidean distances (if crds1 & ref
%                   are specified);
%                [1 x N2]  row vector of euclidean distances (if ref & crds2
%                   are specified);
%                [N1 x N2] rectangular matrix of pairwise distances (if crds1
%                   & crds2 are specified);
%                [1 x 1]   scalar (if crds1 is a [2 x P] matrix or ref1 & ref2
%                   are specified);
%

% RE Strauss, 5/4/94
%   10/28/95 - output row (rather than column) vector for the (reference
%               point)-(set of points) case; still outputs column vector for the
%               (set of points)-(reference point) case.
%   10/30/95 - for double for-loops, put one matrix-col access in outer loop
%               to increase speed.
%   10/12/96 - vectorize inner loop to increase speed.
%    6/12/98 - allow for P=1.
%   11/11/03 - initialize dists to NaN for error return.

  if nargin == 0, help eucl; return; end
  dists = NaN;

  if (nargin < 2)                     % If only crds1 provided,
    [N,P] = size(crds1);
    if (N<2)
      error('  EUCL: need at least two points');
    end

    crds1 = crds1';                     % Transpose crds
    dists = zeros(N,N);                 % Calculate pairwise distances

    for i = 1:N-1
      c1 = crds1(:,i) * ones(1,N-i);
      if (P>1)
        d = sqrt(sum((c1-crds1(:,(i+1:N))).^2));
      else
        d = abs(c1-crds1(:,(i+1:N)));
      end
      dists(i,(i+1:N)) = d;
      dists((i+1:N),i) = d';
    end
    if (N==2)                            % Single distance for two points
      dists = dists(1,2);
    end

  else                                % If crds1 & crds2 provided,
    [N1,P1] = size(crds1);
    [N2,P2] = size(crds2);
    if (P1~=P2)
      error('  EUCL: sets of coordinates must be of same dimension');
    else
      P = P1;
    end

    crds1 = crds1';                     % Transpose crds
    crds2 = crds2';

    if (N1>1 && N2>1)                    % If two matrices provided,
      dists = zeros(N1,N2);               % Calc all pairwise distances between them
      for i = 1:N1
        c1 = crds1(:,i) * ones(1,N2);
        if (P>1)
          d = sqrt(sum((c1-crds2).^2));
        else
          d = abs(c1-crds2);
        end
        dists(i,:) = d;
      end
    end

    if (N1==1 && N2==1)                  % If two vectors provided,
      dists = sqrt(sum((crds1-crds2).^2));  % Calc scalar
    end

    if (N1>1 && N2==1)                   % If matrix && reference point provided,
      crds1 = crds1 - (ones(N1,1)*crds2')'; % Center points on reference point
      if (P>1)                              % Calc euclidean distances in P-space
         dists = sqrt(sum(crds1.^2))';
      else
         dists = abs(crds1)';
      end
    end;                                    % Return column vector

    if (N1==1 && N2>1)                   % If reference point && matrix provided,
      crds2 = crds2 - (ones(N2,1)*crds1')'; % Center points on reference point
      if (P>1)                              % Calc euclidean distances in P-space
        dists = sqrt(sum(crds2.^2));
      else
        dists = abs(crds2);
      end
    end;                                    % Return row vector
  end

  return;


% UPGMA: Unweighted pair-group hierarchical cluster analysis of a distance
%        matrix.  Produces plot of dendrogram.  To bootstrap cluster support,
%        see cluster().
%
%     Usage: [topology,support] = upgma(dist,{labels},{doplot},{fontsize})
%
%         dist =     [n x n] symmetric distance matrix.
%         labels =   optional [n x q] matrix of group labels for dendrogram.
%         doplot =   optional boolean flag indicating, if present and true, that
%                      graphical dendrogram output is to be produced
%                      [default = true].
%         fontsize = optional font size for labels [default = 10].
%         -----------------------------------------------------------------------
%         topology = [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         support =  [(n-2) x (n-1)] matrix, with one row for all but the base 
%                       node, specifying group membership (support) at each node.
%

% To boostrap cluster support, see cluster().

% RE Strauss, 5/27/96
%   9/7/99  -  miscellaneous changes for Matlab v5.
%   9/24/01 - check diagonal elements against eps rather than zero.
%   3/14/04 - changed 'suppress' flag to 'doplot'.

function [topology,support] = upgma(dist,labels,doplot,fontsize)
  if (nargin < 2) labels = []; end;
  if (nargin < 3) doplot = []; end;
  if (nargin < 4) fontsize = []; end;

  suprt = 0;
  if (nargout > 1) suprt = 1; end;
  if (isempty(doplot)) doplot = 1; end;

  [n,p] = size(dist);
  if (n~=p | any(diag(dist)>eps))
    dist
    error('  UPGMA: input matrix is not a distance matrix.');
  end;

  if (~isempty(labels))
    if (size(labels,1)~=n)
      error('  UPGMA: numbers of taxa and taxon labels do not match.');
    end;
  end;

  clstsize = ones(1,n);               % Number of elements in clusters/otus
  id = 1:n;                           % Cluster IDs
  topology = zeros(n-1,4);            % Output dendrogram-topology matrix

  plug = 10e6;
  dist = dist + eye(n)*plug;          % Replace diagonal with plugs

  for step = 1:(n-1)                  % Clustering steps
    min_dist = min(dist(:));            % Find minimum pairwise distance
    [ii,jj] = find(dist==min_dist);     % Find location of minimum
    k = 1;                              % Use first identified minimum
    while  (ii(k)>jj(k))                 %   for which i<j
      k = k+1;
      if k > length(ii) ,keyboard;end;
    end;
    i = ii(k);
    j = jj(k);
    if (id(i)<id(j))
      topology(step,:) = [id(i) id(j) n+step min_dist];
    else
      topology(step,:) = [id(j) id(i) n+step min_dist];
    end;
    id(i) = n+step;
    dist(i,j) = plug;
    dist(j,i) = plug;

    new_clstsize = clstsize(i) + clstsize(j);
    alpha_i = clstsize(i) / new_clstsize;
    alpha_j = clstsize(j) / new_clstsize;
    clstsize(i) = new_clstsize;

    for k = 1:n                         % For all other clusters/OTUs,
      if (k~=i & k~=j)                  %   adjust distances to new cluster
        dist(k,i) = alpha_i * dist(k,i) + alpha_j * dist(k,j);
        dist(i,k) = alpha_i * dist(i,k) + alpha_j * dist(j,k);
        dist(k,j) = plug;
        dist(j,k) = plug;
      end;
    end;
  end; % for step

  if (doplot)                         % Plot dendrogram
    dendplot(topology,labels,fontsize);
  end;

  if (suprt)                          % Specify group membership at nodes
    support = clstsupt(topology);
  end;

  return;

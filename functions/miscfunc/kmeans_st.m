% KMEANS: K-means clustering of n points into k clusters so that the
%         within-cluster sum of squares is minimized.  Based on Algorithm
%         AS136, which seeks a local optimum such that no movement of a
%         point from one cluster to another will reduced the within-cluster
%         sum of squares.  Tends to find spherical clusters.
%         Not globally optimal against all partitions, which is an NP-complete 
%         problem; thus allows for multiple restarts.
%
%     Usage: [centr,clst,sse] = kmeans(X,k,restarts)
%
%         X =        [n x p] data matrix.
%         k =        number of desired clusters.
%         restarts = optional number of restarts, after finding a new minimum 
%                     sse, needed to end search [default=0].
%         --------------------------------------------------------------------
%         centr =    [k x p] matrix of cluster centroids.
%         clst =     [n x 1] vector of cluster memberships.
%         sse =      total within-cluster sum of squared deviations.
%

% Hartigan,JA & MA Wong. 1979.  Algorithm AS 136: A k-means clustering
%   algorithm.  Appl. Stat. 200:100-108.
% Milligan,G.W. & L. Sokol. 1980. A two-stage clustering algorithm with
%   robustness recovery characteristics.  Educational and Psychological
%   Measurement 40:755-759.

% RE Strauss, 8/26/98
%   8/21/99 -  changed misc statements for Matlab v5.
%   10/14/00 - use means() rather than mean() to update cluster centers;
%                remove references to TRUE and FALSE.

function [centr,clst,sse] = kmeans_st(X,k,restarts)
  if (nargin<3) restarts = []; end;
  
  if (isempty(restarts))
    restarts = 0;
  end;

  [n,p] = size(X);
  if (k<1 | k>n)
    error('KMEANS: k out of range 1-N');
  end;

  max_iter = 50;
  highval = 10e8;
  upgma_flag = 1;                     % Do UPGMA first time thru

  % Repeat optimization until have 'restart' attempts with no better solution
  restart_iter = restarts + 1;
  best_sse = highval;

  while (restart_iter > 0)
    restart_iter = restart_iter - 1;

    c1 = zeros(n,1);
    c2 = c1;
    clst = c1;
    nclst = zeros(k,1);

    % Estimate initial cluster centers via UPGMA the first time, and randomly
    % thereafter.

    if (upgma_flag)                     % If UPGMA,
      upgma_flag = 0;
      dist = eucl(X);                     % All interpoint distances
      topo = upgma(dist,[],0);            % UPGMA dendrogram topology

      grp = [1:n]';                       % Identify obs in each of k grps
      for step = 1:(n-k)
        t1 = topo(step,1);
        t2 = topo(step,2);
        t3 = topo(step,3);
        indx = find(grp==t1);
        grp(indx) = t3 * ones(length(indx),1);
        indx = find(grp==t2);
        grp(indx) = t3 * ones(length(indx),1);
      end;

      grpid = uniquef(grp);                % Unique group identifiers
      centr = zeros(k,p);
      for kk = 1:length(grpid);           % Calculate group centroids
        g = grpid(kk);
        indx = find(grp == g);
        Xk = X(indx,:);
        if (size(Xk,1)>1)
          centr(kk,:) = mean(Xk);
        else
          centr(kk,:) = Xk;
        end;
      end;

    else                                % Else if randomized,
      rndprm = randperm(n);               % Pull out random set of points
      centr = X(rndprm(1:k),:);
    end;

    % For each point i, find its two closest centers, c1(i) and c2(i), and
    % assign it to c1(i)

    dist = zeros(n,k);
    for kk = 1:k                        % Dists to centers of all clusters
      dist(:,kk) = eucl(X,centr(kk,:));
    end;
    [d,indx] = sort(dist');             % Sort points separately

    c1 = indx(1,:)';
    c2 = indx(2,:)';

    % Update cluster centers to be the centroids of points contained within them

    for kk = 1:k
      indx = find(c1==kk);
      len_indx = length(indx);
%      if (len_indx>1)
%        centr(kk,:) = mean(X(indx,:));
%      else
%        centr(kk,:) = X(indx,:);
%      end;
      centr(kk,:) = means(X(indx,:));
      nclst(kk) = len_indx;
    end;

    % Initialize working matrices

    live = ones(k,1);
    R1 = zeros(n,1);
    R2 = zeros(k,1);
    adj1 = nclst ./ (nclst+1);
    adj2 = nclst ./ (nclst-1+eps);
    update = n * ones(k,1);

    for i = 1:n               % Adjusted dist from each pt to own cluster center
      R1(i) = adj2(c1(i)) * eucl(X(i,:),centr(c1(i),:))^2;
    end;

    % Iterate

    iter = 0;
    while ((iter < max_iter) & (any(live)))
      iter = iter+1;

      % Optimal-transfer stage

      i = 0;
      while ((i<n) & any(live))         % Cycle thru points
        i = i+1;
        update = update-1;                % Decrement cluster update counters

        L1 = c1(i);                       % Pt i is in cluster L1
        if (nclst(L1)>1)                  % If point i is not sole member of L1,
          R2 = adj1(L1) .* eucl(centr,X(i,:)).^2; % Dists from pt to cluster centers

          if (live(L1))                     % If L1 is in live set,
            R2(L1) = highval;               %   find min R2 over all clusters
            [r2min,L2] = min(R2);
          else                              % If L1 is not in live set,
            indx = find(~live);             %   find min R2 over live clusters only
            R2([L1;indx]) = highval * ones(length(indx)+1,1);
            [r2min,L2] = min(R2);
          end;

          if (r2min >= R1(i))             % No reallocation necessary
            c2(i) = L2;
          else                            % Else reallocate pt to new cluster
            c1(i) = L2;
            c2(i) = L1;

            for kk = [L1,L2]                % Recalc cluster centers and sizes
              indx = find(c1==kk);
              len_indx = length(indx);
              if (len_indx>1)
                centr(kk,:) = mean(X(indx,:));
              else
                centr(kk,:) = X(indx,:);
              end;
              nclst(kk) = len_indx;
            end;

            adj1([L1,L2]) = nclst([L1,L2]) ./ (nclst([L1,L2])+1);
            adj2([L1,L2]) = nclst([L1,L2]) ./ (nclst([L1,L2])-1+eps);
            R1(i) = adj2(c1(i)) * eucl(X(i,:),centr(c1(i),:))^2;

            live([L1,L2]) = [1,1];          % Put into live set
            update([L1,L2]) = [n,n];        % Reinitialize update counters
          end;
        end;  % if nclst(L1)>1

        indx = find(update<1);          % Remove clusters from live set
        if (~isempty(indx))             %   that haven't been recently updated
          live(indx) = zeros(length(indx),1);
        end;
      end;  % while (i<n & any(live))

      if (any(live))

        % Quick transfer stage

        for i = 1:n                         % Cycle thru points
          L1 = c1(i);
          L2 = c2(i);
          if (any(live([L1,L2])) & nclst(L1)>1)
            R2 = adj1(L2) .* eucl(centr(L2,:),X(i,:)).^2; % Dists from pt to L2 center
            if (R1>=R2)
              temp = c1(i);                   % Switch L1 & L2
              c1(i) = c2(i);
              c2(i) = temp;

              for kk = [L1,L2]                % Recalc cluster centers and sizes
                indx = find(c1==kk);
                centr(kk,:) = mean(X(indx,:));
                nclst(kk) = length(indx);
              end;

              adj1([L1,L2]) = nclst([L1,L2]) ./ (nclst([L1,L2])+1);
              adj2([L1,L2]) = nclst([L1,L2]) ./ (nclst([L1,L2])-1+eps);
              R1(i) = adj2(c1(i)) * eucl(X(i,:),centr(c1(i),:))^2;

              live([L1,L2]) = [1,1];          % Put into live set
              update([L1,L2]) = [n,n];        % Reinitialize update counters
            end;
          end;
        end;
      end;  % if any(live)
    end;  % while iter<max_iter & any(live)

    sse = 0;                          % Calc final total sum-of-squares
    for i=1:n
      sse = sse + eucl(X(i,:),centr(c1(i),:))^2;
    end;

    if (sse < best_sse)               % Update best solution so far
      best_sse = sse;
      best_c1 = c1;
      best_centr = centr;
      restart_iter = restarts;
    end;
  end;  % while (restart_iter > 0)

  clst =  best_c1;
  centr = best_centr;
  sse =   best_sse;

  return;

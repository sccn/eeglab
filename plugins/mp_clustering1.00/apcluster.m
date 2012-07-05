%
% [idx,netsim,dpsim,expref]=apcluster(s,p)
%
% APCLUSTER uses affinity propagation (Frey and Dueck, Science,
% 2007) to identify data clusters, using a set of real-valued
% pair-wise data point similarities as input. Each cluster is
% represented by a data point called a cluster center, and the
% method searches for clusters so as to maximize a fitness
% function called net similarity. The method is iterative and
% stops after maxits iterations (default of 500 - see below for
% how to change this value) or when the cluster centers stay
% constant for convits iterations (default of 50). The command
% apcluster(s,p,'plot') can be used to plot the net similarity
% during operation of the algorithm.
%
% For N data points, there may be as many as N^2-N pair-wise
% similarities (note that the similarity of data point i to k
% need not be equal to the similarity of data point k to i).
% These may be passed to APCLUSTER in an NxN matrix s, where
% s(i,k) is the similarity of point i to point k. In fact, only
% a smaller number of relevant similarities are needed for
% APCLUSTER to work. If only M similarity values are known,
% where M < N^2-N, they can be passed to APCLUSTER in an Mx3
% matrix s, where each row of s contains a pair of data point
% indices and a corresponding similarity value: s(j,3) is the
% similarity of data point s(j,1) to data point s(j,2).
%
% APCLUSTER automatically determines the number of clusters,
% based on the input p, which is an Nx1 matrix of real numbers
% called preferences. p(i) indicates the preference that data
% point i be chosen as a cluster center. A good choice is to 
% set all preference values to the median of the similarity
% values. The number of identified clusters can be increased or
% decreased  by changing this value accordingly. If p is a
% scalar, APCLUSTER assumes all preferences are equal to p.
%
% The fitness function (net similarity) used to search for
% solutions equals the sum of the preferences of the the data
% centers plus the sum of the similarities of the other data
% points to their data centers.
%
% The identified cluster centers and the assignments of other
% data points to these centers are returned in idx. idx(j) is
% the index of the data point that is the cluster center for
% data point j. If idx(j) equals j, then point j is itself a
% cluster center. The sum of the similarities of the data
% points to their cluster centers is returned in dpsim, the
% sum of the preferences of the identified cluster centers is
% returned in expref and the net similarity (sum of the data
% point similarities and preferences) is returned in netsim.
%
% EXAMPLE
%
% N=100; x=rand(N,2); % Create N, 2-D data points
% M=N*N-N; s=zeros(M,3); % Make ALL N^2-N similarities
% j=1;
% for i=1:N
%   for k=[1:i-1,i+1:N]
%     s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);
%     j=j+1;
%   end;
% end;
% p=median(s(:,3)); % Set preference to median similarity
% [idx,netsim,dpsim,expref]=apcluster(s,p,'plot');
% fprintf('Number of clusters: %d\n',length(unique(idx)));
% fprintf('Fitness (net similarity): %f\n',netsim);
% figure; % Make a figures showing the data and the clusters
% for i=unique(idx)'
%   ii=find(idx==i); h=plot(x(ii,1),x(ii,2),'o'); hold on;
%   col=rand(1,3); set(h,'Color',col,'MarkerFaceColor',col);
%   xi1=x(i,1)*ones(size(ii)); xi2=x(i,2)*ones(size(ii)); 
%   line([x(ii,1),xi1]',[x(ii,2),xi2]','Color',col);
% end;
% axis equal tight;
%
% PARAMETERS
% 
% [idx,netsim,dpsim,expref]=apcluster(s,p,'NAME',VALUE,...)
% 
% The following parameters can be set by providing name-value
% pairs, eg, apcluster(s,p,'maxits',1000):
%
%   Parameter    Value
%   'sparse'     No value needed. Use when the number of data
%                points is large (eg, >3000). Normally,
%                APCLUSTER passes messages between every pair
%                of data points. This flag causes APCLUSTER
%                to pass messages between pairs of points only
%                if their input similarity is provided and
%                is not equal to -Inf.
%   'maxits'     Any positive integer. This specifies the
%                maximum number of iterations performed by
%                affinity propagation. Default: 500.
%   'convits'    Any positive integer. APCLUSTER decides that
%                the algorithm has converged if the estimated
%                cluster centers stay fixed for convits
%                iterations. Increase this value to apply a
%                more stringent convergence test. Default: 50.
%   'dampfact'   A real number that is less than 1 and
%                greater than or equal to 0.5. This sets the
%                damping level of the message-passing method,
%                where values close to 1 correspond to heavy
%                damping which may be needed if oscillations
%                occur.
%   'plot'       No value needed. This creates a figure that
%                plots the net similarity after each iteration
%                of the method. If the net similarity fails to
%                converge, consider increasing the values of
%                dampfact and maxits.
%   'details'    No value needed. This causes idx, netsim,
%                dpsim and expref to be stored after each
%                iteration.
%   'nonoise'    No value needed. Degenerate input similarities
%                (eg, where the similarity of i to k equals the
%                similarity of k to i) can prevent convergence.
%                To avoid this, APCLUSTER adds a small amount
%                of noise to the input similarities. This flag
%                turns off the addition of noise.
%
% Copyright (c) Brendan J. Frey and Delbert Dueck (2006). This
% software may be freely used and distributed for
% non-commercial purposes.

function [idx,netsim,dpsim,expref]=apcluster(s,p,varargin);

% Handle arguments to function
if nargin<2 error('Too few input arguments');
else
    maxits=500; convits=50; lam=0.5; plt=0; details=0; nonoise=0;
    i=1;
    while i<=length(varargin)
        if strcmp(varargin{i},'plot')
            plt=1; i=i+1;
        elseif strcmp(varargin{i},'details')
            details=1; i=i+1;
		elseif strcmp(varargin{i},'sparse')
			[idx,netsim,dpsim,expref]=apcluster_sparse(s,p,varargin{:});
			return;
        elseif strcmp(varargin{i},'nonoise')
            nonoise=1; i=i+1;
        elseif strcmp(varargin{i},'maxits')
            maxits=varargin{i+1};
            i=i+2;
            if maxits<=0 error('maxits must be a positive integer'); end;
        elseif strcmp(varargin{i},'convits')
            convits=varargin{i+1};
            i=i+2;
            if convits<=0 error('convits must be a positive integer'); end;
        elseif strcmp(varargin{i},'dampfact')
            lam=varargin{i+1};
            i=i+2;
            if (lam<0.5)||(lam>=1)
                error('dampfact must be >= 0.5 and < 1');
            end;
        else i=i+1;
        end;
    end;
end;
if lam>0.9
    fprintf('\n*** Warning: Large damping factor in use. Turn on plotting\n');
    fprintf('    to monitor the net similarity. The algorithm will\n');
    fprintf('    change decisions slowly, so consider using a larger value\n');
    fprintf('    of convits.\n\n');
end;

% Check that standard arguments are consistent in size
if length(size(s))~=2 error('s should be a 2D matrix');
elseif length(size(p))>2 error('p should be a vector or a scalar');
elseif size(s,2)==3
    tmp=max(max(s(:,1)),max(s(:,2)));
    if length(p)==1 N=tmp; else N=length(p); end;
    if tmp>N
        error('data point index exceeds number of data points');
    elseif min(min(s(:,1)),min(s(:,2)))<=0
        error('data point indices must be >= 1');
    end;
elseif size(s,1)==size(s,2)
    N=size(s,1);
    if (length(p)~=N)&&(length(p)~=1)
        error('p should be scalar or a vector of size N');
    end;
else error('s must have 3 columns or be square'); end;

% Construct similarity matrix
if N>3000
    fprintf('\n*** Warning: Large memory request. Consider activating\n');
    fprintf('    the sparse version of APCLUSTER.\n\n');
end;
if size(s,2)==3
    S=-Inf*ones(N,N); 
    for j=1:size(s,1) S(s(j,1),s(j,2))=s(j,3); end;
else S=s;
end;

% In case user did not remove degeneracies from the input similarities,
% avoid degenerate solutions by adding a small amount of noise to the
% input similarities
if ~nonoise
    rns=randn('state'); randn('state',0);
    S=S+(eps*S+realmin*100).*rand(N,N);
    randn('state',rns);
end;

% Place preferences on the diagonal of S
if length(p)==1 for i=1:N S(i,i)=p; end;
else for i=1:N S(i,i)=p(i); end;
end;

% Allocate space for messages, etc
dS=diag(S); A=zeros(N,N); R=zeros(N,N); t=1;
if plt netsim=zeros(1,maxits+1); end;
if details
    idx=zeros(N,maxits+1);
    netsim=zeros(1,maxits+1); 
    dpsim=zeros(1,maxits+1); 
    expref=zeros(1,maxits+1); 
end;

% Execute parallel affinity propagation updates
e=zeros(N,convits); dn=0; i=0;
while ~dn
    i=i+1; 

    % Compute responsibilities
    Rold=R;
    AS=A+S; [Y,I]=max(AS,[],2); for k=1:N AS(k,I(k))=-realmax; end;
    [Y2,I2]=max(AS,[],2);
    R=S-repmat(Y,[1,N]);
    for k=1:N R(k,I(k))=S(k,I(k))-Y2(k); end;
    R=(1-lam)*R+lam*Rold; % Damping

    % Compute availabilities
    Aold=A;
    Rp=max(R,0);
    for k=1:N Rp(k,k)=R(k,k); end;
    A=repmat(sum(Rp,1),[N,1])-Rp;
    dA=diag(A); A=min(A,0); for k=1:N A(k,k)=dA(k); end;
    A=(1-lam)*A+lam*Aold; % Damping

    % Check for convergence
    E=((diag(A)+diag(R))>0); e(:,mod(i-1,convits)+1)=E; K=sum(E);
    if i>=convits || i>=maxits
        se=sum(e,2);
        unconverged=(sum((se==convits)+(se==0))~=N);
        if (~unconverged&&(K>0))||(i==maxits) dn=1; end;
    end;

    % Handle plotting and storage of details, if requested
    if plt||details
        if K==0
            tmpnetsim=nan; tmpdpsim=nan; tmpexpref=nan; tmpidx=nan;
        else
            I=find(E); [tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c);
            tmpnetsim=sum(S((tmpidx-1)*N+[1:N]'));
            tmpexpref=sum(dS(I)); tmpdpsim=tmpnetsim-tmpexpref;
        end;
    end;
    if details
        netsim(i)=tmpnetsim; dpsim(i)=tmpdpsim; expref(i)=tmpexpref;
        idx(:,i)=tmpidx;
    end;
    if plt
        netsim(i)=tmpnetsim;
        figure(234); 
        tmp=1:i; tmpi=find(~isnan(netsim(1:i)));
        plot(tmp(tmpi),netsim(tmpi),'r-');
        xlabel('# Iterations');
        ylabel('Fitness (net similarity) of quantized intermediate solution');
        drawnow; 
    end;
end;
I=find(diag(A+R)>0); K=length(I); % Identify exemplars
if K>0
    [tmp c]=max(S(:,I),[],2); c(I)=1:K; % Identify clusters
    % Refine the final set of exemplars and clusters and return results
    for k=1:K ii=find(c==k); [y j]=max(sum(S(ii,ii),1)); I(k)=ii(j(1)); end;
    [tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c);
    tmpnetsim=sum(S((tmpidx-1)*N+[1:N]')); tmpexpref=sum(dS(I));
else
    tmpidx=nan*ones(N,1); tmpnetsim=nan; tmpexpref=nan;
end;
if details
    netsim(i+1)=tmpnetsim; netsim=netsim(1:i+1);
    dpsim(i+1)=tmpnetsim-tmpexpref; dpsim=dpsim(1:i+1);
    expref(i+1)=tmpexpref; expref=expref(1:i+1);
    idx(:,i+1)=tmpidx; idx=idx(:,1:i+1);
else
    netsim=tmpnetsim; dpsim=tmpnetsim-tmpexpref;
    expref=tmpexpref; idx=tmpidx;
end;
if plt||details
    fprintf('\nNumber of identified clusters: %d\n',K);
    fprintf('Fitness (net similarity): %f\n',tmpnetsim);
    fprintf('  Similarities of data points to exemplars: %f\n',dpsim(end));
    fprintf('  Preferences of selected exemplars: %f\n',tmpexpref);
    fprintf('Number of iterations: %d\n\n',i);
end;
if unconverged
    fprintf('\n*** Warning: Algorithm did not converge. The similarities\n');
    fprintf('    may contain degeneracies - add noise to the similarities\n');
    fprintf('    to remove degeneracies. To monitor the net similarity,\n');
    fprintf('    activate plotting. Also, consider increasing maxits and\n');
    fprintf('    if necessary dampfact.\n\n');
end;

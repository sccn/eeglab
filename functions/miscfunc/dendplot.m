function dendplot(topology,labels,fontsize)
% DENDPLOT: Plots a dendrogram given a topology matrix.
%
%     Usage: dendplot(topology,{labels},{fontsize})
%
%         topology = [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         labels =   optional [n x q] matrix of group labels for dendrogram.
%         fontsize = optional font size for labels [default = 10].
%

% RE Strauss, 5/27/98
%   8/20/99 - miscellaneous changes for Matlab v5

  if (nargin<2)
    labels = [];
  end
  if (nargin < 3)
    fontsize = [];
  end

  if (isempty(fontsize))              % Default font size for labels
    fontsize = 10;
  end

  r = size(topology,1);
  n = r+1;                            % Number of taxa

  links = dendhier([],topology,n-1);  % Find dendrogram links (branches)
  otu_indx = find(links(:,1)<=n);     % Get sequence of OTUs
  otus = links(otu_indx,1);
  y = zeros(2*n-1,1);                 % Y-coords for plot
  y(otus) = 0.5:(n-0.5);
  for i = 1:(n-1)
    y(topology(i,3)) = mean([y(topology(i,1)),y(topology(i,2))]);
  end

  clf;                                % Begin plot
  hold on;

  for i = 1:(2*n-2)                   % Horizontal lines
    desc = links(i,1);
    anc =  links(i,2);
    X = [links(i,3) links(i,4)];
    Y = [y(desc) y(desc)];
    plot(X,Y,'k');
  end

  for i = (n+1):(2*n-1)               % Vertical lines
    indx = find(links(:,2)==i);
    X = [links(indx,4)];
    Y = [y(links(indx(1),1)) y(links(indx(2),1))];
    plot(X,Y,'k');
  end

  maxdist = max(links(:,4));
  for i = 1:n                         % OTU labels
    if (~isempty(labels))
      h = text(-.02*maxdist,y(i),labels(i,:));   % For OTUs on right
      set(h,'fontsize',fontsize);
    else
      h = text(-.02*maxdist,y(i),num2str(i));  % For OTUs on right
      set(h,'fontsize',fontsize);
%     text(-.06*maxdist,y(i),num2str(i)); % For UTOs on left
    end
  end

  axis([0 maxdist+0.03*maxdist 0 n]); % Axes
  axis('square');
  set(gca,'Ytick',[]);                % Suppress y-axis labels and tick marks
  set(gca,'Ycolor','w');              % Make y-axes invisible
  set(gca,'Xdir','reverse');          % For OTUs on right

  hold off;
  return;

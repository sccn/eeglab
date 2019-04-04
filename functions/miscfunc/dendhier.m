function [links,topology,node] = dendhier(links,topology,node)
% DENDHIER: Recursive algorithm to find links and distance coordinates on a 
%             dendrogram, given the topology matrix.
%
%     Usage: [links,topology,node] = dendhier(links,topology,node)
%
%         links =     4-col matrix of descendants, ancestors, descendant 
%                       distances, and ancestor distances; pass to 
%                       function as null vector []
%         topology =  dendrogram topology matrix
%         node =      current node; pass as N-1
%

% RE Strauss, 7/13/95 

  n = size(topology,1)+1;             % Number of OTUs

  c1 =   topology(node,1);
  c2 =   topology(node,2);
  clst = topology(node,3);
  dist = topology(node,4);

  if (c1 <= n)
    links = [links; c1 clst 0 dist];
  else                                
    prevnode = find(topology(:,3)==c1);
    prevdist = topology(prevnode,4);
    links = [links; c1 clst prevdist dist];
    [links,topology,node] = dendhier(links,topology,prevnode);
  end

  if (c2 <= n)
    links = [links; c2 clst 0 dist];
  else
    prevnode = find(topology(:,3)==c2);
    prevdist = topology(prevnode,4);
    links = [links; c2 clst prevdist dist];
    [links,topology,node] = dendhier(links,topology,prevnode);
  end

  return;

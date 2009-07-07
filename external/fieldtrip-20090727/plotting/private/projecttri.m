function [tri] = projecttri(pnt);

% PROJECTTRI makes a closed triangulation of a list of vertices by
% projecting them onto a unit sphere and subsequently by constructing
% a convex hull triangulation.
%
% Use as
%   [tri] = projecttri(pnt)

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.2  2007/05/08 07:36:42  roboos
% updated help
%
% Revision 1.1  2006/12/12 11:27:45  roboos
% created subfunction into a seperate function
%

ori = (min(pnt) + max(pnt))./2;
pnt(:,1) = pnt(:,1) - ori(1);
pnt(:,2) = pnt(:,2) - ori(2);
pnt(:,3) = pnt(:,3) - ori(3);
nrm = sqrt(sum(pnt.^2, 2));
pnt(:,1) = pnt(:,1)./nrm;
pnt(:,2) = pnt(:,2)./nrm;
pnt(:,3) = pnt(:,3)./nrm;
tri = convhulln(pnt);


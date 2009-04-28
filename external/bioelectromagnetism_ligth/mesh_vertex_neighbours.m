function [nI,nXYZ,nD] = mesh_vertex_neighbours(FV,index),

% mesh_vertex_neighbours - find all immediate neighbours of a vertex
%
% [nI,nXYZ,nD] = mesh_vertex_neighbours(FV,index)
% 
% FV is a struct with fields:
%
% FV.vertices   - Nx3 matrix of Cartesian vertex coordindates (X,Y,Z)
% FV.faces      - Mx3 matrix of triangulation of FV.vertices
% FV.edge       - NxN sparse matrix of edge lengths (see mesh_edges)
%
% [nI,nXYZ,nD]  - vertex neighbour Indices, XYZ coordinates and Distance
%                 from vertex
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% Licence:  GNU GPL, no implied or express warranties
% History:  04/2004, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('index','var'),
    fprintf('...no index specified, using vertex 1\n');
    index = 1;
end
if isempty(index),
    fprintf('...no index specified, using vertex 1\n');
    index = 1;
end

if ~isfield(FV,'edge'),
    FV.edge = mesh_edges(FV);
end

Nvert = size(FV.vertices,1);

fprintf('...finding vertex neighbours...'); tic;

% the indices of the neighbours
nI = find(FV.edge(index,:))';

% the coordinates of the neighbours
nXYZ = FV.vertices(nI,:);

% only return unique vertex coordinates
%[nXYZ, i, j] = unique(FV.vertices(nI,:),'rows');

% extract neighbour edge lengths
nD = full(FV.edge(index,nI))';


t = toc; fprintf('done (%5.2f sec).\n',t);

return

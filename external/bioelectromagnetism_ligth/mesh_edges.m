function edge = mesh_edges(FV)

% mesh_edges - Calculate edge lengths of triangulation
% 
% edge = mesh_edges(FV)
% 
% FV.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces      - triangulation of vertices
% 
% edge          - edge lengths, indexed by vertex 
%                 number (sparse NxN matrix)
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  07/2002, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...searching for mesh edges...');

nvertex = size(FV.vertices,1);
nface   = size(FV.faces,1);

% the 'edge' matrix is the connectivity of all vertices
edge = sparse(nvertex,nvertex);

for f = 1:nface,
    
    % compute the length of all triangle edges (Diff is [3x3])
    Diff = [FV.vertices(FV.faces(f,[1 2 3]),:) - FV.vertices(FV.faces(f,[2 3 1]),:)];
    Norm = sqrt( sum(Diff.^2, 2) );
    
    edge(FV.faces(f,1),FV.faces(f,2)) = Norm(1);
    edge(FV.faces(f,2),FV.faces(f,3)) = Norm(2);
    edge(FV.faces(f,3),FV.faces(f,1)) = Norm(3);
    
    % make sure that all edges are symmetric
    edge(FV.faces(f,2),FV.faces(f,1)) = Norm(1);
    edge(FV.faces(f,3),FV.faces(f,2)) = Norm(2);
    edge(FV.faces(f,1),FV.faces(f,3)) = Norm(3);
end

t=toc;
fprintf('done (%5.2f sec).\n',t);

return

function [faceNormals,unit_faceNormals,centroids] = mesh_face_normals(FV)

% mesh_face_normals - Calculate face surface normals
% 
% [normals,unit_normals,origins] = mesh_face_normals(FV)
% 
% FV.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces      - triangulation of vertices (Mx3 matrix)
% 
% normals       - face surface normals (Mx3 matrix)
% unit_normals  - normalized normals!
%
% origins       - the point at the center of each face
%
% When calculating the surface normal, this function assumes the convention
% that the vertices are given in FV.faces in clockwise order:
%
%        1
%       /\
%      /  \
%   e2/    \e1
%    /      \
%   /        \
%  /__________\
% 3	    	   2
%
% We then define edge1 (e1) from vertex1 to vertex2 and edge2 (e2) from
% vertex1 to vertex3.  The cross product of these two edge vectors gives
% the surface normal.  The direction of the normal is either into the
% page or out of the page, depending on the order of the cross product,
%
% edge1 x edge2 = -1 * ( edge2 x edge1 )
%
% So, the order of the vertex points, the direction of the edge vectors
% and their cross products are very important if you want a particular
% direction.
%
% In this function, we assume that the vertices of each face are given in
% clockwise order, when viewed from the "outside".  The resulting surface
% normals are oriented "outward".  Here is an example of how to view them:
%
% figure
% Hp = patch('faces',FV.faces,'vertices',FV.vertices,...
%     'facecolor',[.7 .7 .7],'facealpha',0.8,'edgecolor',[.8 .8 .8]); 
% camlight('headlight','infinite'); daspect([1 1 1]); axis vis3d; axis off
% material dull; hold on
% [faceNormals,faceNormalsUnit,centroids] = mesh_face_normals(FV);
% quiver3(centroids(f,1),centroids(f,2),centroids(f,3),...
%     faceNormals(f,1),faceNormals(f,2),faceNormals(f,3));
%


% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  04/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...calculating face surface normals...');

nface   = size(FV.faces,1);

faceNormals = zeros(nface,3);
unit_faceNormals = faceNormals;
centroids = faceNormals;

for f = 1:nface,
    
    vertex_index1 = FV.faces(f,1);
    vertex_index2 = FV.faces(f,2);
    vertex_index3 = FV.faces(f,3);
    
    vertex1 = FV.vertices(vertex_index1,:);
    vertex2 = FV.vertices(vertex_index2,:);
    vertex3 = FV.vertices(vertex_index3,:);
    
    % If the vertices are given in clockwise order, when viewed from the
    % outside, then following calculates the "outward" surface normals.
    
    edge_vector1 = vertex2 - vertex1;
    edge_vector2 = vertex3 - vertex1;
    
    faceNormals(f,:) = cross( edge_vector2, edge_vector1 );
    
    magnitude = vector_magnitude(faceNormals(f,:));
    
    unit_faceNormals(f,:) = faceNormals(f,:) / magnitude;
    
    
    % Now find the midpoint between all vertices
    a = (vertex1 + vertex2) ./ 2;
    b = (vertex2 + vertex3) ./ 2;
    c = (vertex3 + vertex1) ./ 2;
    
    % Now find the centroid length of the medians
    centroids(f,:) = b + ( (vertex1 - b) ./3 );
    
end

t=toc;
fprintf('done (%5.2f sec).\n',t);

return

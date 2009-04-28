function [vertNormals,vertNormalsUnit] = mesh_vertex_normals(FV,index)

% mesh_vertex_normals - Calculate vertex surface normals
% 
% [normals,unit_normals] = mesh_vertex_normals(FV,index)
% 
% FV.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces      - triangulation of vertices (Mx3 matrix)
%
% index is an array of vertex indices (default is all vertices)
% 
% normals       - vertex surface normals (Nx3 matrix)
% unit_normals  - normalized normals!
%
% This routine first calculates the surface normals of each face.  It then
% finds all the faces that contain a specific vertex and sums across
% the face normals (weighted by the face area).  
%
% If the faces are defined
% according to the right-hand rule, all their normals will be "outward"
% normals and the average should be a sensible value for the vertex normal.
%
% See also the VertexNormals property of the patch and surface commands.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  04/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;
fprintf('...calculating vertex surface normals...');

Nvert = size(FV.vertices,1);

if ~exist('index','var'), index = 1:Nvert; end
if isempty(index), index = 1:Nvert; end

vertNormals = zeros(length(index),3);
vertNormalsUnit = vertNormals;

for i = 1:length(index),
    
    vi = index(i);
    
    % To calculate a vertex normal you need to perform an iteration
    % beginning with a triangle that holds the vertex and perform a cross
    % product on the two vectors of that triangle that extend from that
    % vertex. This will result in a normal for that triangle. Then go
    % around to each triangle containing that vertex and perform a cross
    % product for that triangle. You will end up with a set of normals for
    % all the triangles. To compute the vertex normal, sum all the face
    % normal vectors.
    
    % get all the faces that contain this vertex
    [faceIndexI,faceIndexJ] = find(FV.faces == vi);
    
    Nface = length(faceIndexI);
    
    faceNormals = zeros(Nface,3);
    
    for face = 1:Nface,
        
        f = faceIndexI(face);
        
        vi1 = FV.faces(f,1);
        vi2 = FV.faces(f,2);
        vi3 = FV.faces(f,3);
        
        v1 = FV.vertices(vi1,:);
        v2 = FV.vertices(vi2,:);
        v3 = FV.vertices(vi3,:);
        
        % If the vertices are given in clockwise order, when viewed from the
        % outside, then following calculates the "outward" surface normals.
        
        edge1 = v2 - v1;
        edge2 = v3 - v1;
        
        faceNormals(face,:) = cross( edge2, edge1 );
        
    end
    
    %faceNormalsMagnitude = vector_magnitude(faceNormals);
    %faceNormalsUnit = faceNormals ./ repmat(faceNormalsMagnitude,1,3);
    
    % Area of Triangle = || AB x AC|| / 2 ; ie, the absolute value of
    % the length of the cross product of AB and AC divided by 2
    %faceArea = abs(faceNormalsMagnitude) / 2;
    
    % Weight the faceNormals by the faceArea
    %vertNormals(v,:) = sum( faceNormals .* repmat(faceArea,1,3) ) / sum( faceArea );
    
    vertNormals(i,:) = sum(faceNormals) / Nface;
    
    vertNormalsMagnitude = norm(vertNormals(i,:));
    vertNormalsUnit(i,:) = vertNormals(i,:) / vertNormalsMagnitude;
    
end

t=toc;
fprintf('done (%5.2f sec).\n',t);

return

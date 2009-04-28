function [faceArea] = mesh_face_area(FV)

% mesh_face_area - Calculate face area
% 
% [area] = mesh_face_area(FV)
% 
% FV.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces      - triangulation of vertices (Mx3 matrix)
% 
% area          - face area (Mx1 array)
%
% If the lengths of the three sides are known then Heron's formula
% can be used: sqrt(s*(s-a)*(s-b)*(s-c)) (where a, b, c are the sides
% of the triangle, and s = (a + b + c)/2 is half of its perimeter).
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  05/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...calculating face area...');

nface   = size(FV.faces,1);

faceArea = zeros(nface,1);

for f = 1:nface,
    
    % get the vertex index numbers into FV.vertices
    vertex_index1 = FV.faces(f,1);
    vertex_index2 = FV.faces(f,2);
    vertex_index3 = FV.faces(f,3);
    
    % get the vertex coordinates in Cartesian XYZ
    vertex1 = FV.vertices(vertex_index1,:);
    vertex2 = FV.vertices(vertex_index2,:);
    vertex3 = FV.vertices(vertex_index3,:);
    
    % face surface area
    % If the lengths of the three sides are known then Heron's formula
    % can be used: sqrt(s*(s-a)*(s-b)*(s-c)) (where a, b, c are the sides
    % of the triangle, and s = (a + b + c)/2 is half of its perimeter).
    
    % a is the distance from vertex1 to vertex2
    a = sqrt ( sum( ( vertex2 - vertex1 ).^2 ) );
    % b is the distance from vertex2 to vertex3
    b = sqrt ( sum( ( vertex3 - vertex2 ).^2 ) );
    % c is the distance from vertex3 to vertex1
    c = sqrt ( sum( ( vertex1 - vertex3 ).^2 ) );
    
    s = (a + b + c)/2;
    
    faceArea(f) = sqrt(s*(s-a)*(s-b)*(s-c));
    
end

t=toc;
fprintf('done (%5.2f sec).\n',t);

return

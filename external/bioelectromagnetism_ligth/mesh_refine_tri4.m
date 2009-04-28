function [ FV ] = mesh_refine_tri4(FV)

% mesh_refine_tri4 - creates 4 triangle from each triangle of a mesh
%
% [ FV ] = mesh_refine_tri4( FV )
%
% FV.vertices   - mesh vertices (Nx3 matrix)
% FV.faces      - faces with indices into 3 rows
%                 of FV.vertices (Mx3 matrix)
% 
% For each face, 3 new vertices are created at the 
% triangle edge midpoints.  Each face is divided into 4
% faces and returned in FV.
%
%        B
%       /\
%      /  \
%    a/____\b       Construct new triangles
%    /\    /\       [A,a,c]
%   /  \  /  \      [a,B,b]
%  /____\/____\     [c,b,C]
% A	     c	   C    [a,b,c]
% 
% It is assumed that the vertices are listed in clockwise order in
% FV.faces (A,B,C above), as viewed from the outside in a RHS coordinate
% system.
% 
% See also: mesh_refine, sphere_tri, sphere_project
% 


% ---this method is not implemented, but the idea here remains...
% This can be done until some minimal distance (D) of the mean 
% distance between vertices of all triangles is achieved.  If
% no D argument is given, the function refines the mesh once.
% Alternatively, it could be done until some minimum mean 
% area of faces is achieved.  As is, it just refines once.


% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  05/2002, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...refining mesh (tri4)...')

% NOTE
% The centroid is located one third of the way from each vertex to 
% the midpoint of the opposite side. Each median divides the triangle 
% into two equal areas; all the medians together divide it into six 
% equal parts, and the lines from the median point to the vertices 
% divide the whole into three equivalent triangles.

% Each input triangle with vertices labelled [A,B,C] as shown
% below will be turned into four new triangles:
%
% Make new midpoints
% a = (A+B)/2
% b = (B+C)/2
% c = (C+A)/2
%
%        B
%       /\
%      /  \
%    a/____\b       Construct new triangles
%    /\    /\       [A,a,c]
%   /  \  /  \      [a,B,b]
%  /____\/____\     [c,b,C]
% A	     c	   C    [a,b,c]
%

% Initialise a new vertices and faces matrix
Nvert = size(FV.vertices,1);
Nface = size(FV.faces,1);
V2 = zeros(Nface*3,3);
F2 = zeros(Nface*4,3);

for f = 1:Nface,
    
    % Get the triangle vertex indices
    NA = FV.faces(f,1);
    NB = FV.faces(f,2);
    NC = FV.faces(f,3);
    
    % Get the triangle vertex coordinates
    A = FV.vertices(NA,:);
    B = FV.vertices(NB,:);
    C = FV.vertices(NC,:);
    
    % Now find the midpoints between vertices
    a = (A + B) ./ 2;
    b = (B + C) ./ 2;
    c = (C + A) ./ 2;
    
    % Find the length of each median
    %A2blen = sqrt ( sum( (A - b).^2, 2 ) );
    %B2clen = sqrt ( sum( (B - c).^2, 2 ) );
    %C2alen = sqrt ( sum( (C - a).^2, 2 ) );
    
    % Store the midpoint vertices, while
    % checking if midpoint vertex already exists
    [FV, Na] = mesh_find_vertex(FV,a);
    [FV, Nb] = mesh_find_vertex(FV,b);
    [FV, Nc] = mesh_find_vertex(FV,c);
    
    % Create new faces with orig vertices plus midpoints
    F2(f*4-3,:) = [ NA, Na, Nc ];
    F2(f*4-2,:) = [ Na, NB, Nb ];
    F2(f*4-1,:) = [ Nc, Nb, NC ];
    F2(f*4-0,:) = [ Na, Nb, Nc ];
    
end

% Replace the faces matrix
FV.faces = F2;

t=toc; fprintf('done (%5.2f sec)\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV, N] = mesh_find_vertex(FV,vertex)

    Vn = size(FV.vertices,1);
    Va = repmat(vertex,Vn,1);
    Vexist = find( FV.vertices(:,1) == Va(:,1) & ...
                   FV.vertices(:,2) == Va(:,2) & ...
                   FV.vertices(:,3) == Va(:,3) );
    if Vexist,
        if size(Vexist) == [1,1],
            N = Vexist;
        else,
            msg = sprintf('replicated vertices');
            error(msg);
        end
    else
        FV.vertices(end+1,:) = vertex;
        N = size(FV.vertices,1);
    end

return

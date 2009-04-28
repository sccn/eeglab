function [ FV ] = mesh_refine_tri6(FV)

% mesh_refine_tri6 - creates 6 smaller triangles from a triangle mesh
%
% [ FV ] = mesh_refine_tri6( FV )
%
% FV.vertices is vertices (Nx3 matrix)
% FV.faces is faces with indices into rows of V (Mx3 matrix)
%
% For each face of F, 4 new vertices are created at the 
% triangle edge midpoints and the triangle centroid.  Each
% face is divided into 6 faces and returned in FV.faces.
% 
%        B
%       /+\
%      / + \
%    a/__X__\b       Construct new triangles
%    /  +|+  \       [A,a,X], [A,X,c]
%   / +  |  + \      [B,X,a], [B,b,X]
%  /+___ | ___+\     [C,X,b], [C,c,X]
% A	     c	   C
%
% It is assumed that the vertices are listed in clockwise order in
% FV.faces (A,B,C above), as viewed from the outside in a RHS coordinate
% system.
% 
% See also: mesh_refine, mesh_refine_tri4, 
%           sphere_tri, sphere_project
%


% This can be done until some minimal distance (D) of the mean 
% distance between vertices of all triangles is achieved.  If
% no D argument is given, the function refines the mesh once.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  05/2002, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic; fprintf('...refining mesh (tri6)...');

% NOTE
% The centroid is located one third of the way from each vertex to 
% the midpoint of the opposite side. Each median divides the triangle 
% into two equal areas; all the medians together divide it into six 
% equal parts, and the lines from the median point to the vertices 
% divide the whole into three equivalent triangles.

% Each input triangle with vertices labelled [A,B,C] as shown
% below will be turned into six new triangles:
%
% Make new midpoints
% a = (A+B)/2
% b = (B+C)/2
% c = (C+A)/2
%
% Make triangle centroid
% X = b + ( (A - b) ./3 );
% 
%        B
%       /+\
%      / + \
%    a/__X__\b       Construct new triangles
%    /  +|+  \       [A,a,X], [A,X,c]
%   / +  |  + \      [B,X,a], [B,b,X]
%  /+___ | ___+\     [C,X,b], [C,c,X]
% A	     c	   C
%




% Initialise a new vertices and faces matrix
Nvert = size(FV.vertices,1);
Nface = size(FV.faces,1);
V2 = zeros(Nface*4,3);
F2 = zeros(Nface*6,3);

for f=1:Nface,
    
    % Get the triangle vertex indices
    NA = FV.faces(f,1);
    NB = FV.faces(f,2);
    NC = FV.faces(f,3);
    
    % Get the triangle vertex coordinates
    A = FV.vertices(NA,:);
    B = FV.vertices(NB,:);
    C = FV.vertices(NC,:);
    
    % Now find the midpoint between all vertices
    a = (A + B) ./ 2;
    b = (B + C) ./ 2;
    c = (C + A) ./ 2;
    
    % Now find the centroid length of the medians
    X = b + ( (A - b) ./3 );
    %Bc = c + ( (B - c) ./3 );  % Bc = X
    %Cc = a + ( (C - a) ./3 );  % Cc = X
    
    % Store the midpoint and the centroid vertices,
    % checking if vertex already exists
    [FV, Na] = mesh_find_vertex(FV,a);
    [FV, Nb] = mesh_find_vertex(FV,b);
    [FV, Nc] = mesh_find_vertex(FV,c);
    [FV, NX] = mesh_find_vertex(FV,X);
    
    % Create new faces with centroid
    F2(f*6-5,:) = [NA, Na, NX ];
    F2(f*6-4,:) = [NA, NX, Nc ];
    F2(f*6-3,:) = [NB, NX, Na ];
    F2(f*6-2,:) = [NB, Nb, NX ];
    F2(f*6-1,:) = [NC, NX, Nb ];
    F2(f*6-0,:) = [NC, Nc, NX ];
    
    
    %figure; patch('vertices',[A;B;C],'faces',[1 2 3],'facecolor',[.7 .7 .7]); hold on;
    %plot3(A(1),A(2),A(3),'ro');
    %plot3(b(1),b(2),b(3),'ro');
    %plot3(Ac(1),Ac(2),Ac(3),'bo')
    %if isequal(r,2), return; end
    
end

% Replace the faces matrix
FV.faces = F2;


t=toc; fprintf('done (%5.2f sec)\n',t);

return


% Find the length of each median
%A2bLength = sqrt ( sum( (A - b).^2, 2 ) );
%B2cLength = sqrt ( sum( (B - c).^2, 2 ) );
%C2aLength = sqrt ( sum( (C - a).^2, 2 ) );


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

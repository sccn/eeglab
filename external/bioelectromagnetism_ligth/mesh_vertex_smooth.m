function [FV] = mesh_vertex_smooth(FV, index),

% [FV] = mesh_vertex_smooth(FV, index)
%
% FV.vertices - Nx3
% FV.faces - Nx3
% FV.edges - NxN, see mesh_edges
%
% index is the indices of FV.vertices to smooth (all by default)
%
% This function shifts vertex location toward the center of mass of its
% neighbouring vertices.  It is useful to maintain surface smoothness.
%
% This code is developed on the basis of Smith (2002, fig 5).
%
% Smith, S. (2002). Fast robust automated brain extraction.
%    Human Brain Mapping, 17(3): 143-155.


% This function adapts Smith, S. (2002), Fast robust automated brain
% extraction.  Human Brain Mapping, 17: 143-155.  This function
% corresponds to update component 2: surface smoothness control.


if isfield(FV,'edges'),
    if isempty(FV.edges),
        FV.edges = mesh_edges(FV);
    end
else
    FV.edges = mesh_edges(FV);
end

% Define limits for radius of curvature, empirically optimized per
% Smith (2002), see figure 6.
radius_min =  3.33; % mm
radius_max = 10.00; % mm
% Sigmoid function parameters,
% "where E and F control the scale and offset of the sigmoid"
E = mean([(1 / radius_min),  (1 / radius_max)]);
F = 6 * ( (1 / radius_min) - (1 / radius_max) );

Nvert = size(FV.vertices,1);

if ~exist('index','var'), index = 1:Nvert; end
if isempty(index), index = 1:Nvert; end

fprintf('...smoothing %d vertices of surface triangulation\n',length(index));

for i = 1:length(index),

    vi = index(i);
    
    v = FV.vertices(vi,:);

    % Find neighbouring vertex coordinates
    neighbour_index = find(FV.edges(vi,:));  % the indices
    
    % Find distances between vertex and neighbours, using edge lengths.
    % The mean value is lower(L) in Smith (2002, eq. 4)
    edge_distance = full(FV.edges(vi,neighbour_index));
    L = mean(edge_distance);

    % Find neighbour mean location;
    % this is 'mean position of A and B'
    % in figure 4 of Smith (2002)
    neighbour_vertices = FV.vertices(neighbour_index,:);
    Nmean = mean(neighbour_vertices);

    % Find difference in distance between the central vertex and its
    % neighbours; this value is 's' in Smith (2002, fig 4); we take the
    % direction of this vector to be inward toward mean location of
    % neighbours
    s = Nmean - v;

    % see function below
    vertNormal = surf_normal(FV, vi);
    
    unit_normal = vertNormal / norm(vertNormal);

    % Find the vector sn
    % the projection of s in the direction of the surface normal
    sn = dot( s, unit_normal ) * unit_normal;

    % Find the vector st, the component of s orthogonal to the
    % surface normal vector.
    st = s - sn;
    
    % Calculate radius of local curvature, solve Smith (2002, eq. 4)
    r = L^2 / (2 * norm(sn));
    
    % calculate the fraction of Sn to move the vertex
    f2 = (1 + tanh( F * (1 / r - E))) / 2;

    % If we move the vertices at this point, all subsequent calculations
    % for the neighbours will account for the current move.  So, we
    % calculate the surface normals during this loop, and we apply the
    % movements incrementally.  If we accumulate the movement vectors and
    % apply them after the loop, the risk is that vertices will move too
    % far toward one another.
    
    % sum([ v; st / 2; f2 * sn ])
    
    FV.vertices(vi,:) = v + (st / 2) + ( f2 * sn );
    
    plotface = 0;
    if plotface,
        figure; hold on
        %patch('faces',FV.faces,'vertices',FV.vertices,'facecolor',[.2 .2 .2],'facealpha',0.6)
        [faceIndexI,faceIndexJ] = find(FV.faces == vi);
        patch('faces',FV.faces(faceIndexI,:),'vertices',FV.vertices,'facecolor',[.6 .6 .6])        
        scatter3(v(1), v(2), v(3), 40, 'r', 'filled')
        tmp = neighbour_vertices;
        scatter3(tmp(:,1), tmp(:,2), tmp(:,3), 40, 'b', 'filled')
        vn = vertNormal;
        quiver3(v(1), v(2), v(3), vn(1), vn(2), vn(3), 0)
        tmpS = 10 * s;
        quiver3(v(1), v(2), v(3), tmpS(1), tmpS(2), tmpS(3), 0)
        tmpSn = 10 * sn;
        quiver3(v(1), v(2), v(3), tmpSn(1), tmpSn(2), tmpSn(3), 0)
        tmpSt = 10 * st;
        quiver3(v(1)+tmpSn(1), v(2)+tmpSn(2), v(3)+tmpSn(3), tmpSt(1), tmpSt(2), tmpSt(3), 0)
        
        
        motion = -10;
        tmpVn = v + (st / 2) + ( f2 * motion * sn);
        
        
        scatter3(tmpVn(1), tmpVn(2), tmpVn(3), 40, 'g', 'filled')
        rotate3d on
    end

end

% We could use St and Sn to move the vertex toward the center of mass
% defined by the neighbour vertices, but this is best done inside the loop
% above
%F2 = repmat(f2,1,3);
%FV.vertices = FV.vertices + (St / 2) + ( F2 .* Sn );

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vertNormal = surf_normal(FV, vi),

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

nface = length(faceIndexI);

faceNormals = zeros(nface,3);

for i = 1:nface,

    f = faceIndexI(i);

    i1 = FV.faces(f,1);
    i2 = FV.faces(f,2);
    i3 = FV.faces(f,3);

    v1 = FV.vertices(i1,:);
    v2 = FV.vertices(i2,:);
    v3 = FV.vertices(i3,:);

    edge1 = v2 - v1;
    edge2 = v3 - v1;
    
    % If the vertices are given in clockwise order, when viewed from the
    % outside, the following calculates the "outward" surface normal. This
    % should be consistent with the matlab patch command.
    faceNormals(i,:) = cross( edge2, edge1 );

    % If the vertices are given in counterclockwise order, when viewed from
    % the outside, the following calculates the "outward" surface normal.
    % This should be consistent with the right hand rule.
    %faceNormals(i,:) = cross( edge1, edge2 );

end

%faceNormalsMagnitude = vector_magnitude(faceNormals);
%faceNormalsUnit = faceNormals ./ repmat(faceNormalsMagnitude,1,3);

% Area of Triangle = || AB x AC|| / 2 ; ie, the absolute value of
% the length of the cross product of AB and AC divided by 2
%faceArea = abs(faceNormalsMagnitude) / 2;

% Weight the faceNormals by the faceArea
%vertNormals(v,:) = sum( faceNormals .* repmat(faceArea,1,3) ) / sum( faceArea );

vertNormal = sum(faceNormals) / nface;

return

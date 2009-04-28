function [FV] = mesh_vertex_spacing(FV),

% [FV] = mesh_vertex_spacing(FV)
%
% FV.vertices - Nx3
% FV.faces - Nx3
% FV.edges - NxN, see mesh_edges
%
% This function shifts vertex location toward the center of its
% neighbouring vertices.  It is useful to maintain vertex spacing.
%
% This code is developed on the basis of Smith (2002, fig 4).
%
% Smith, S. (2002). Fast robust automated brain extraction. 
%    Human Brain Mapping, 17(3): 143-155.


if isfield(FV,'edges'),
    if isempty(FV.edges),
        FV.edges = mesh_edges(FV);
    end
else
    FV.edges = mesh_edges(FV);
end

%[normals,unit_normals] = mesh_vertex_normals(FV);

% get surface normals
hf = figure('Visible','off');
hp = patch('faces',FV.faces,'vertices',FV.vertices);
normals = get(hp,'VertexNormals');
close(hf); clear hf hp
%Convert to unit normals
[normals,unit_normals] = colnorm(normals');
unit_normals = unit_normals';
clear normals


Nvertices = size(FV.vertices,1);

for index = 1:Nvertices,

    v = FV.vertices(index,:);
    x = FV.vertices(index,1);
    y = FV.vertices(index,2);
    z = FV.vertices(index,3);

    unit_normal = unit_normals(index,:);

    % Find neighbouring vertex coordinates
    vi = find(FV.edges(index,:));  % the indices
    neighbour_vertices = FV.vertices(vi,:);
    X = neighbour_vertices(:,1);
    Y = neighbour_vertices(:,2);
    Z = neighbour_vertices(:,3);

    % Find neighbour mean location; this is 'mean position of A and B' in
    % figure 4 of Smith (2002)
    Xmean = mean(X);
    Ymean = mean(Y);
    Zmean = mean(Z);

    % Find difference in distance between the vertex of interest and its
    % neighbours; this value is 's' and 'sn' in figure 4 of
    % Smith (2002, eq. 1 to 4)
    s = [ Xmean - x, Ymean - y, Zmean - z]; % inward toward mean

    % Find the vector sn
    % the projection of s in the direction of the surface normal
    sn = dot( s, unit_normal ) * unit_normal;

    % Find the vector st, the component of s orthogonal to the
    % surface normal vector.
    st = s - sn; % absolute value

    
    S(index,:) = s;
    Sn(index,:) = sn;
    St(index,:) = st;

end

% We can now use St to move the vertex toward the mean location of the
% neighbour vertices

FV.vertices = FV.vertices + St;

return

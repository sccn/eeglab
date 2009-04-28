function [FV] = mesh_vertex_smooth(FV,index,origin),

% This function adapts Smith, S. (2002), Fast robust automated brain
% extraction.  Human Brain Mapping, 17: 143-155.  This function
% corresponds to update component 2: surface smoothness control.


fprintf('this is in development\n');
return


xo = origin(1); yo = origin(2); zo = origin(3);

v = FV.vertices(index,:);
x = FV.vertices(index,1);
y = FV.vertices(index,2);
z = FV.vertices(index,3);

% Find radial distance of vertex from origin
r = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );

% Calculate unit vector
v_unit_vector = ( v - origin ) / r;

% Find direction cosines for line from center to vertex
l = (x-xo)/r; % cos alpha
m = (y-yo)/r; % cos beta
n = (z-zo)/r; % cos gamma

% Find neighbouring vertex coordinates
vi = find(FV.edge(index,:));  % the indices
neighbour_vertices = FV.vertices(vi,:);
X = neighbour_vertices(:,1);
Y = neighbour_vertices(:,2);
Z = neighbour_vertices(:,3);

% Find neighbour radial distances
r_neighbours = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
r_neighbours_mean = mean(r_neighbours);

% Find difference in radial distance between the vertex of interest and its
% neighbours; this value approximates the magnitude of sn in 
% Smith (2002, eq. 1 to 4)
r_diff = r - r_neighbours_mean;

% Find the vector sn, in the direction of the vertex of interest, given the
% difference in radial distance between vertex and mean of neighbours
sn = r_diff * v_unit_vector;

% Find distances between vertex and neighbours, using edge lengths.
% The mean value is l in Smith (2002, eq. 4)
edge_distance = FV.edge(index,vi);
edge_distance_mean = mean(edge_distance);

% Calculate radius of local curvature, solve Smith (2002, eq. 4)
if r_diff,
  radius_of_curvature = (edge_distance_mean ^ 2) / (2 * r_diff);
else
  radius_of_curvature = 10000;
end

% Define limits for radius of curvature
radius_min =  3.33; % mm
radius_max = 10.00; % mm

% Sigmoid function parameters,
% "where E and F control the scale and offset of the sigmoid"
E = mean([(1 / radius_min),  (1 / radius_max)]);
F = 6 * ( (1 / radius_min) - (1 / radius_max) );

Fsigmoid = (1 + tanh( F * (1 / radius_of_curvature - E))) / 2;

% multiply sigmoid function by sn
move_vector = Fsigmoid * sn;

FV.vertices(index,:) = v + move_vector;

return

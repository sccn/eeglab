function [distances] = mesh_vertex_distance(vertices)

% mesh_vertex_distance - Calculates spherical intervertex distances (arc length).
%
% Usage: [distances] = mesh_vertex_distance(vertices)
%
% Notes:    vertices is an Nx3 matrix of rectangular Cartesian coordinates
%           with a centroid at (0,0,0).  If the centroid is elsewhere, all
%           vertex points are internally shifted to adjust the centroid 
%           because the arc length method employed assumes a 
%           centroid at (0,0,0).
%           
%           Sphere radius is estimated by the average radius from
%           (xo,yo,zo) to any 2 pairs of vertices.  It will vary
%           from one pair to another, but this method works well for
%           small theta (eg, nearest neighbours).
%           
%           Returns a matrix for paired electrode arc length estimates.
%           The matrix is a full matrix, but the upper triange is zeros.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% Author:   Darren.Weber_at_radiology.ucsf.edu
% Created:  18/05/00 - linear distance
% Modified: 24/06/01 - spherical arc length (should be OK for small theta
%                      but for large theta, elliptical arc length may be
%                      preferable).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s\n', 'Calculating inter-vertex spherical arc length.');

% check centroid, unless input parameters defined
mean = mean(vertices);
xo = mean(1);
yo = mean(2);
zo = mean(3);

% define vertices A,B, given centroid (xo,yo,zo)
A = [ (vertices(:,1)-xo) (vertices(:,2)-yo) (vertices(:,3)-zo) ];

A_len = sqrt( sum(A.^2,2) );   % length of vertex vector A

% Initialise distances matrix
distances = zeros(length(vertices));

progress_bar('init');
for a=1:length(vertices),
    
    progress = a / length(vertices);
    progress_bar('set',progress);
    
    Aa = A(a,:);
    Al = A_len(a);  % length of vertex vector Aa
    
    r = (Al + A_len(a:end))/2;  % estimate sphere radius from Al and B_len
    dotAB = zeros(length(A)-a,1);
    for i=a:length(A), dotAB(i-(a-1),1) = dot(Aa,A(i,:)); end
    theta = acos( dotAB ./ (Al * A_len(a:end)) );
    arc_len = r .* theta;  % arc length = radius * theta
    distances(a:length(A),a) = arc_len;
    
    %for b=a:length(vertices),
    %    Bb = A(b,:);
    %    Bl = A_len(b); % length of vertex vector Bb
    %    if( Aa == Bb ),
    %        arc_len = 0;  % no distance from vertex to itself
    %    else
    %        r = (Al + Bl)/2;  % estimate sphere radius from A_len and B_len
    %        theta = acos( dot(Aa,Bb) / (Al * Bl) );  % Angle between A & B, in radians
    %        arc_len = r * theta;  % arc length = radius * theta
    %    end
    %    distances(a,b) = arc_len;
    %end
end
distances = real(distances);

progress_bar('clear');

return
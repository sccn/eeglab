% elec_distance_testscript - test electrode distance calculations


% Create 2 points a & b in Cartesian 3space
elec_labels = cellstr(['a' 'b']');
x = [ 2 -2]';
y = [ 3  2]';
z = [ 2  0]';

xd = x(1) - x(2); yd = y(1) - y(2); zd = z(1) - z(2);
dist = sqrt ( xd^2 + yd^2 + zd^2 );
fprintf('Cartesian distance = %6.4f\n', dist);

[labels, dist] = elec_distance(elec_labels,x,y,z);
fprintf('Spherical arc length = %6.4f\n', dist(2,1));
function [T] = fiducial_coregister(P,Q)

% fiducial_coregister - Determine fiducial points coregistration
% 
% [T] = fiducial_coregister(P,Q)
% 
% P is [3x3], three points in Cartesian XYZ
% Q is [3x3], three points in Cartesian XYZ
% 
% This function 3D rotates and translates P into
% the space of Q.  The method is described in
% this m file, developed from Mortenson, M. (1985),
% Geometric Modelling, New York: John Wiley & Sons,
% pp. 353-355.
% 
% The return 4x4 matrix (T) contains a 3x3 rotation 
% matrix in the upper left and a 1x3 translation
% row vector in the bottom row. The last column contains
% a projection vector in T(1:3,4)(all zeros here) and 
% the scalar at T(4,4) is a homogeneous coordinate scale 
% factor, usually 1.
% 
% Using this T matrix:
% to get P values in the space of Q,
% P2Q = [P,ones(size(P,1),1)] * T;
% to get Q values in the space of P,
% Q2P = [Q,ones(size(Q,1),1)] * inv(T);
% 
% In this manner, the rotation and translation is
% applied to any Nx3 points in relation to P space.
% 
% For use in coregistration of electrodes to a scalp
% mesh, P & Q are the fiducial points:
% P(1,:) is the electrode nasion XYZ
% P(2,:) is the electrode right preauricular XYZ
% P(3,:) is the electrode left preauricular XYZ
% Similarly:
% Q(1,:) is the scalp mesh nasion XYZ
% Q(2,:) is the scalp mesh right preauricular XYZ
% Q(3,:) is the scalp mesh left preauricular XYZ
% The order of these points is essential!
% 
% For example, using emse_read_reg to obtain reg, we have:
% T = fiducial_coregister(reg.fiducials.head,reg.fiducials.mri);
% where T is almost equivalent to reg.elec2mri.  This
% function is a hard fiducial transformation, whereas 
% the EMSE coregistration is a bit more sophisticated,
% providing weights to the fiducial points and some MRI
% surface points.
%
% See also: MESH_FIT_ELEC
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no express or implied warranties
% History:  09/2002, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DB = 0; % debug, 0 off, 1 on


if ~exist('P','var') | ~exist('Q','var'),
    msg = sprintf('ELEC_COREGISTER: no input P or Q\n');
    error(msg);
elseif isempty(P) | isempty(Q),
    msg = sprintf('ELEC_COREGISTER: empty input P or Q\n');
    error(msg);
elseif ~(isequal(size(P),size(ones(3,3)))),
    msg = sprintf('ELEC_COREGISTER: P must be 3x3 matrix\n');
    error(msg);
elseif ~(isequal(size(Q),size(ones(3,3)))),
    msg = sprintf('ELEC_COREGISTER: Q must be 3x3 matrix\n');
    error(msg);
end


% Given a geometric model that contains 3D points p1,p2,p3
% and three other points q1,q2,q3, find the total rigid
% body transformation that:
% (1) transforms p1 into q1;
% (2) transforms the vector (p2 - p1) into the 
%     vector (q2 - q1), direction only;
% (3) transforms the plane containing the three points
%     p1,p2,p3 into the plane containing q1,q2,q3

% Note that when p1,p2,p3 are reference points among
% many points of a geometric model, the transformations
% here can be applied to all points to bring the model
% into coincidence with the points in q1,q2,q3

if DB,
	% These 2 point sets can be coregistered with just
	% a rotation of 90 degrees around Z
	P = [0 0 0; 1 0 0;  0 1 0];
	Q = [0 0 0; 0 1 0; -1 0 0];
    P = P + 10;
end

% Step 2 of Mortenson (1985, p. 354)
V1 = P(2,:) - P(1,:);
W1 = Q(2,:) - Q(1,:);

% Step 3 of Mortenson (1985, p. 354)
% V3 & W3 are orthogonal to the plane containing
% the points p1,p3 and q1,q3 respectively
V3 = cross( V1, [P(3,:) - P(1,:)] );
W3 = cross( W1, [Q(3,:) - Q(1,:)] );

% Step 4 of Mortenson (1985, p. 354)
% V2 & W2 are orthogonal to the vectors
% V1,V3 and W1,W3 respectively
V2 = cross( V3, V1 );
W2 = cross( W3, W1 );

% Thus, V1,V2,V3 form a right hand orthogonal system,
% as do W1,W2,W3

% Step 5 of Mortenson (1985, p. 354)
% Calculate unit vectors
v1 = V1 / sqrt( V1(1)^2 + V1(2)^2 + V1(3)^2 );
v2 = V2 / sqrt( V2(1)^2 + V2(2)^2 + V2(3)^2 );
v3 = V3 / sqrt( V3(1)^2 + V3(2)^2 + V3(3)^2 );
w1 = W1 / sqrt( W1(1)^2 + W1(2)^2 + W1(3)^2 );
w2 = W2 / sqrt( W2(1)^2 + W2(2)^2 + W2(3)^2 );
w3 = W3 / sqrt( W3(1)^2 + W3(2)^2 + W3(3)^2 );

% Step 6 of Mortenson (1985, p. 354)
% To transform any point p in the v system into the
% w system, use the transformation relationship
% p* = pR + T
% Obviously we need to find R and T first

% Step 7 of Mortenson (1985, p. 354)
% [w1 w2 w3] = [v1 v2 v3]R, since [w] and [v] are
% the unit vector matrices.
w = [w1;w2;w3];
v = [v1;v2;v3];
% Then the required rotation matrix with respect 
% to the w system is simply:
R = v\w;
% Note in matlab A\B is the matrix division 
% of A into B, which is roughly the same as INV(A)*B,
% except it is computed by Gaussian elimination.

% Step 8 of Mortenson (1985, p. 354)
% Obtain the translation matrix by substituting
% R back into step 6 and solving for T, using
% p1 and q1
T = Q(1,:) - [ P(1,:) * R ];

% Now construct the homogeneous coordinate
% transformation matrix that encorporates both
% of the above R & T transforms
temp1 = [ [ R; zeros(1,3) ], zeros(4,1) ];
temp1(4,4) = 1;
temp1(4,1:3) = T;

T = temp1;

if DB,
    P2Q = [P,ones(3,1)] * T;
	plot3(P(:,1),P(:,2),P(:,3),'bo'); hold on; view(2)
	plot3(Q(:,1),Q(:,2),Q(:,3),'ro');
	plot3(P2Q(:,1),P2Q(:,2),P2Q(:,3),'b.');
    Q2P = [Q,ones(3,1)] * inv(T);
    plot3(Q2P(:,1),Q2P(:,2),Q2P(:,3),'r.');
end


return

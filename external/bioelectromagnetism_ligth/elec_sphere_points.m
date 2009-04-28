function [x,y,z] = elec_sphere_points(Epoints,Rpoints,R)

% elec_sphere_points - Generate points on a hemisphere, radius R
%
% Useage: [X,Y,Z] = elec_sphere_points(Epoints, Rpoints, R)
%
%           Epoints     number of elevations (phi)
%           Rpoints     number of points in rotation (theta)
%                       for each level of elevation. Works best
%                       when a multiple of 4 or 2, asymmetrical 
%                       when odd.
%           R           sphere radius
%
%           [X,Y,Z]     cartesian points
%                       size = (Epoints * Rpoints) - (Rpoints - 1)
%   `                   because the vertex point has no rotation.
%
% Example:  [x,y,z] = elec_sphere_points(10,16,1); plot3(x,y,z,'o');
%
% Note:     The points are more numerous closer to the vertex.  It
%           is tricky to create points with constant arc length
%           separation.
%
%           This function cannot create standard 10-20 electrode
%           positions.  These are given in spherical coordinates
%           in the text files 'elec_10-20*.txt'.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/01 Darren.Weber_at_radiology.ucsf.edu
%                 It would be nice to decrease number of Rotation
%                 points for every increase in elevation, but it
%                 is tricky to maintain symmetry.  If possible,
%                 this could maintain near constant arc length
%                 between points.
%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


eegversion = '$Revision: 1.1 $';
fprintf('ELEC_SPHERE_POINTS [v %s]\n',eegversion(11:15));
fprintf('...generating spherical points\n');
tic;

phi    = linspace(0,pi/2,Epoints);     % elevation, vertex = 0

theta1 = linspace(0,2*pi,Rpoints+1);   % 360 degree rotation
theta1 = theta1(1:end-1);              % remove 2*pi (== 0)

thetad = (theta1(2) - theta1(1))/2;    % provide staggered rotation

n = 0; odd = 1;
for i = 1:length(phi)
    
    % every 2nd level of elevation, stagger rotation by thetad
    if(odd == 1), odd = 0;  theta = theta1;
    else,         odd = 1;  theta = theta1 + thetad;
    end
    
    if isequal(phi(i),0),  % vertex point
        n = n + 1;  x(n) = 0;  y(n) = 0;  z(n) = R;
    elseif(phi(i) < (pi/8)),
        for j = 1:length(theta)
            
            A(1) = x(n);
            A(2) = y(n);
            A(3) = z(n);
            B(1) = R * sin(phi(i)) * cos(theta(j));
            B(2) = R * sin(phi(i)) * sin(theta(j));
            B(3) = R * cos(phi(i));
            
            A_len = sqrt ( sum(A.^2) );
            B_len = sqrt ( sum(B.^2) );
            
            angle = acos( dot(A,B) / (A_len * B_len) );
            
            if and((phi(i) < (pi/16)),(angle>(pi/90))),
                n = n + 1;
                x(n) = B(1); y(n) = B(2); z(n) = B(3);
            elseif(angle>(pi/45)),
                n = n + 1;
                x(n) = B(1); y(n) = B(2); z(n) = B(3);
            end
        end
    else
        for j = 1:length(theta)
            n = n + 1;
            x(n) = R * sin(phi(i)) * cos(theta(j));
            y(n) = R * sin(phi(i)) * sin(theta(j));
            z(n) = R * cos(phi(i));
        end
    end
end
x = x';
y = y';
z = z';

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return

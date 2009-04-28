function [labels, distances] = elec_distance(elec_labels,X,Y,Z,xo,yo,zo)

% elec_distance - Calculates spherical interelectrode distances (arc length).
%
% Usage: [labels, distances] = elec_distance(elec_labels,X,Y,Z,xo,yo,zo)
%
% Notes:    Arc length method assumes a centroid at (xo,yo,zo) = (0,0,0)
%           and input values are in rectangular Cartesian coordinates.
%           
%           Sphere radius is estimated by the average radius from
%           (xo,yo,zo) to any 2 pairs of electrodes.  It will vary
%           from one pair to another, but this method works well for
%           small theta (eg, nearest neighbours).
%
%           Returns 2 matrices, one for paired electrode labels and 
%           another for spherical arc length estimates.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Author:   Darren.Weber_at_radiology.ucsf.edu
% Created:  18/05/00 - linear distance
% Modified: 24/06/01 - spherical arc length (should be OK for small theta
%                      but for large theta, elliptical arc length may be
%                      preferable).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialise centroid, unless input parameters defined
    if ~exist('xo','var')  xo = 0;   end
    if ~exist('yo','var')  yo = 0;   end
    if ~exist('zo','var')  zo = 0;   end

    fprintf('%s\n', 'Calculating inter-electrode spherical arc length.');

    labels = cell(1);
    distances = [];
 
    % convert numerical electrode labels to cell strings
    if isnumeric(elec_labels) > 0
        if size(elec_labels,1) > size(elec_labels,2)
            elec_labels = cellstr(strcat(num2str(elec_labels)));
        else
            elec_labels = cellstr(strcat(num2str(elec_labels')));
        end
    end
    
    A = [ (X-xo) (Y-yo) (Z-zo) ];  B = A;              % define electrode vectors A,B, given centroid (xo,yo,zo)
    rows = 50;                                         % progress indicator
    for a =   1:length(X)
        fprintf('.'); if( a == rows ) fprintf('\n'); rows = rows + 50; end  % progress indicator
        
        Aa = A(a,:);
        A_len = sqrt( sum( Aa.^2 ) );                        % length of electrode vector A
        
        for b = 1:length(X)
            
            Bb = B(b,:);
            B_len = sqrt ( sum(Bb.^2) );                     % length of electrode vector B
            
            if( Aa == Bb )
                arc_len = 0;                                    % no distance from electrode to itself
            else
                r = (A_len + B_len)/2;                          % estimate sphere radius from A_len and B_len
                theta = acos( dot(Aa,Bb) / (A_len * B_len) );   % Angle between A & B, in radians
                arc_len = r * theta;                            % arc length = radius * theta
            end
            
            distances(a, b) = arc_len;
            labels(a, b) = strcat(elec_labels(a), ', ', elec_labels(b));
        end
    end
    fprintf('\n');

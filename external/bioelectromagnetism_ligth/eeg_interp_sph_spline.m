function [FV,C] = eeg_interp_sph_spline(Zi,Ei)

% eeg_interp_sph_spline - Spherical Spline Interpolation of Potential
%
% Useage: [FV,C] = eeg_interp_sph_spline(Zi,Ei)
%
% where:    Zi is Nelec x 1, an EEG/ERP measurement at time t
%           Ei is Nelec x 3, [X Y Z] electrode positions.
%           The origin of Ei is assumed (0,0,0).
%
% FV => interpolated spherical surface (see sphere_tri)
%
% FV.faces    => triangulation of FV.vertices 
% FV.vertices => cartesian coordinates (Nx3)
% FV.Cdata    => spherical spline potential at FV.vertices
% 
% C => interpolation coefficients of Ei (includes co = C(1))
% 
% Notes:    This function calculates the spherical spline of 
%           Perrin et al (1989).  Electroenceph. & Clin. 
%             Neurophysiology, 72: 184-187. Corrigenda (1990),
%             Electroenceph. & Clin. Neurophysiology, 76: 565.
%             (see comments in the .m file for details).

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2001 Darren.Weber_at_radiology.ucsf.edu, with
%                   mathematical advice from
%                   Dr. Murk Bottema (Flinders University of SA)
%           10/2003 Darren.Weber_at_radiology.ucsf.edu, with
%                   mathematical advice and LegendreP function from 
%                   Dr. Tom.Ferree_at_radiology.ucsf.edu
%                   revised, tested & initial verification complete
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for correct size & orientation of Ei & Zi
[e1,e2] = size(Ei);
[v1,v2] = size(Zi);
if e1 < e2, Ei = Ei'; [e1,e2] = size(Ei); end
if v1 < v2, Zi = Zi'; [v1,v2] = size(Zi); end
if ~and(isequal(e1,v1),and(isequal(e2,3),isequal(v2,1)))
  error('...Ei must be Nx3 & Zi must be Nx1');
end
nElectrodes = e1; % The number of electrodes
clear e1 e2 v1 v2;


% -------------------------------------------------------------------------
% estimate spherical radius of the electrodes and
% obtain spherical projections of Ei
[r,x,y,z] = elec_sphere_project(Ei(:,1),Ei(:,2),Ei(:,3));
%Ei = [ x y z ]; clear x y z;

% create spherical interpolation surface
FV = sphere_tri('ico',4,r);

% -------------------------------------------------------------------------
% Calculate the cosines, if Ei is Nx3, COS is NxN matrix
% We use (Ei,Ei) here because it gives the cosines between
% each electrode and every other electrode.  This is required
% here because we solve a linear system of equations below
% that will find the interpolated value at a given electrode
% location, which must be equal to the measured 
% potential at that location.
EiCOS = elec_cosines(Ei,Ei);

% create zeros on the diagonal elements
% [not sure why this works, but it does.]
for i = 1:length(EiCOS), EiCOS(i,i) = 0; end

% -------------------------------------------------------------------------
% Calculate g(x), nElectrodes x nElectrodes
Gx = eeg_interp_sph_spline_g(EiCOS);

% -------------------------------------------------------------------------
% calculate the spherical interpolation coefficients
C = eeg_interp_sph_spline_c(Zi,Gx); clear Gx

% -------------------------------------------------------------------------
% cosines between electrodes and interpolation points
FvCOS = elec_cosines(Ei,FV.vertices);

% Calculate g(x), nElectrodes x NinterpolationPoints
Gx = eeg_interp_sph_spline_g(FvCOS);

eegversion = '$Revision: 1.1 $';
fprintf('EEG_INTERP_SPH_SPLINE [v %s]\n',eegversion(11:15)); tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that Ci is solved, we can obtain interpolated potentials at Ej (Eq. 1)
% U(Ej) = c(0) + ( for i=1:n, sum = (sum + (c(i) * g(cos(Ei,Ej)))) )
% U(Ej) = c(0) + sum( Ci * g(x) )

% Solve Eq 1. (where FV.Cdata = U)
Co = C(1);
Ci = C(2:end);
FV.Cdata = Co + ( Ci' * Gx );

t=toc; fprintf('...done (%6.2f sec)\n',t);

return

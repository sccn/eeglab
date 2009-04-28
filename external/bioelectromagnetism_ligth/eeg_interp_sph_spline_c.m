function [C] = eeg_interp_sph_spline_c(Zi,G)

% eeg_interp_sph_spline_c - Spherical Spline Coefficients
%
% Useage: [C] = eeg_interp_sph_spline_c(Zi,G)
%
% Zi - Nx1 is an EEG/ERP measurement at time t from N sensors
% G  - NxN Legendre function of the cosine matrix for N sensor positions
%      eg, COS = elec_cosines(Ei,Ei); where Ei is Nx3 [X Y Z]
%          G = eeg_interp_sph_spline_g(COS)
%
% C - spherical spline coefficients for Ei (includes co = C(1))
% 
% Notes:    This function calculates the spherical spline coefficients of 
%           Perrin et al (1989).  Electroenceph. & Clin. 
%             Neurophysiology, 72: 184-187. Corrigenda (1990),
%             Electroenceph. & Clin. Neurophysiology, 76: 565.
%             (see comments in the .m file for details).
%
% see elec_cosines, eeg_interp_sph_spline, eeg_interp_sph_spline_g

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eegversion = '$Revision: 1.1 $';
fprintf('EEG_INTERP_SPH_SPLINE_C [v %s]\n',eegversion(11:15)); tic

%--------------------------------------------------------------------------
% Check for correct size & orientation of Zi & G
[g1,g2] = size(G);
[z1,z2] = size(Zi);
if g1 ~= g2, error('G matrix must be square'); end
if z1 < z2, Zi = Zi'; [z1,z2] = size(Zi); end
if ~and(isequal(g1,z1),isequal(z2,1)),
  error('...G must be NxN & Zi must be Nx1');
end
nElectrodes = z1; % The number of electrodes
clear g1 g2 z1 z2;


%--------------------------------------------------------------------------
% Computations described by Perrin et al. (1989)
%
% Let z(i) be the potential value measured at the ith electrode whose
% spherical projection is denoted Ei.  The value at any point E of the
% sphere, of the spherical spline U that interpolates from zi at all Ei,
% is given by (Eq. 1):
% 
% U(E) = Co + ( for i=1:nElectrodes, sum = (sum + (C(i) * g(cos(E,E(i)))) )
% 
% Where the column array C(i) contains the spherical spline coefficients.
% These are solutions of the simultaneous equations (Eq. 2):
% 
% a) GC + TCo = Z
% b) T'C = 0
% 
% with:
% T' = ( 1, 1,..., 1),
% C' = (c1,c2,...,cn),  - unknown
% Z' = (z1,z2,...,zn),  - potentials at Ei electrodes
% G  = (Gij) = (g(cos(Ei,Ej))) - cosine matrix betweeen E and Ei
% 
% That is, cos(Ei,Ej) denotes the cosine of the angle between the
% electrode projections Ei and the sphere interpolation points Ej.
% 
% The function g(x) is the sum of the series (Eq. 3),
%
% g(x) = (1/(4*pi) * ...
%        for n=1:inf, sum = sum + ( ( (2*n+1)/(n^m * (n+1)^m) ) * Pn(x) )
%
% where m is a constant > 1 and Pn(x) is the nth degree Legendre
% polynomial.  Perrin et al. (1989) evaluated m=1:6 and recommend m=4,
% for which the first 7 terms of Pn(x) are sufficient to obtain
% a precision of 10^-6 for g(x) (ie, the above for loop is n=1:7).
% Perrin et al. (1989) recommend tabulation of g(x) for
% x = linspace(-1,1,2000) to be used as a lookup given actual
% values for cos(E(i),E(j)).

% Initialise T and C
T =  ones(nElectrodes,1);
C = zeros(nElectrodes,1); % so Eq. 2b, T' * C = 0

% G is NxN, where N = nElectrodes
% C is Nx1 coefficients that will solve Eq. 2a
% Z is Nx1 potential values

% GC is Nx1, a linear combination of the columns of G
% TCo is Nx1, a constant column vector
%
% We can rearrange Eq. 2a so we have
% GC = Z, effectively incorporating TCo into G,C and Z.
% First, we consider Co at the top of C, so we get C' = (co,c1,c2,...cn) 
% and C becomes an (nElectrodes+1)x1 column array (ie [0;C]).  (Note that
% we solve for C below, so this concatenation is not actually done.)
% Then we have to concatenate [0 T'] onto the top of G (ie, [ [0 T']; G]),
% so that [0 T'] * C = 0 (this is the first line of Gx multiplied by C).
% Lastly, we have to concatenate a column of ones into the first column
% of G and place a zero into G(1,1), so that effectively we add co to the
% sum of multiplying G by C.  The resulting linear system is:
%
% [ 0 1 .. 1] [co]   [ 0]
% [ 1       ] [c1]   [z1]
% [ .       ] [c2] = [z2]
% [ .   G   ] [..]   [..]
% [ 1       ] [cn]   [zn]
%
% so...

% add ones to first row of G, call it Gx
Gx = [T'; G];
% now add ones to first column of Gx, but
% with a zero at Gx(1,1)
Gx = [[0;T] Gx];

% Append 0 to the top of Zi (Z becomes (nElectrodes+1)x1 )
Z = [0;Zi];

% Find the inverse of Gx
%invGx = inv(Gx);

% or
invGx = Gx\eye(size(Gx));

% or
%lambda = 0;
%tolerance = 10^-7 * 10^(lambda*-1);
%invGx = pinv(Gx,tolerance);

% Calculate inverse of Gx
% Ginv = inv(GT);
% EMSE NOTE: The inversion of G is sensitive to errors because 
% electrode locations are not really on a sphere. Therefore, 
% the solution should be regularized.  In EMSE v4.2, the 
% regularization involves Singluar Value Decomposition:
% G = UWV'
% inv(G) approx = V W- U'
%
% W- = [w-(i,j)], w-(i,j) { 1/w(ij) for w(i,j) > cutoff
%                           0       otherwise
%
% cutoff is determined by a smoothing parameter:
%
% cutoff = k * 10 ^ -lambda
%
% where k = { 10^-5 for laplacian
%             10^-7 for voltage
% and lambda is a the smoothing parameter.
% Smaller values of lambda (<0), fewer terms are
% truncated and the solution is less smooth.
%
% [U,W,S] = svd(Gx);  %Gx = U*W*S';
% invW = inv(W);
% lambda = 0;
% cutoff = 10^-7 * 10^(lambda*-1);
% Winv(find(invW<cutoff)) = 0;
% invGx = S * invW * U.';



% now we solve for C = inv(Gx)*Z
C = invGx * Z;

Co = C(1);
Ci = C(2:end);

t=toc; fprintf('...done (%6.2f sec)\n\n',t);

return

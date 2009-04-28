function [Cos] = elec_cosines(A,B)

% elec_cosines - Compute cosine of angle between electrode position vectors
%
% Usage: [Cos] = elec_cosines(A,B)
%
%         A,B  both Nx3 (X,Y,Z) electrode coordinates
%              assuming a common origin at (0,0,0)
%         COS  cosine angle between each row of A & B
%              If A is Nx3 and B is Mx3, COS is NxM.
%
% Note:       For a complete cosine matrix, let A = B.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2001, Darren.Weber_at_radiology.ucsf.edu
%           01/2004, Darren.Weber_at_radiology.ucsf.edu
%                    cos(0) = 1, not cos(0) = 0 at line 43 ???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eegversion = '$Revision: 1.1 $';
fprintf('ELEC_COSINES [v %s]\n',eegversion(11:15));

if isempty(A), error('Input matrix A is empty\n'); end
if isempty(B), error('Input matrix B is empty\n'); end

As = size(A);
if ~isequal(As(2),3),
  if isequal(As(1),3), A = A';
  else error('Input matrix A should be Nx3 [x y z] coordinates\n');
  end
end
Bs = size(B);
if ~isequal(Bs(2),3),
  if isequal(Bs(1),3), B = B';
  else error('Input matrix B should be Nx3 [x y z] coordinates\n');
  end
end

fprintf('...calculating cosines\n'); tic

for   a = 1:size(A,1),  Aa = A(a,:); A_len = norm(Aa);
  for b = 1:size(B,1),  Bb = B(b,:); B_len = norm(Bb);
    
    Cos(a,b) = dot(Aa,Bb) / (A_len * B_len);
  end
end

t=toc; fprintf('...done (%6.2f sec)\n\n',t);

return

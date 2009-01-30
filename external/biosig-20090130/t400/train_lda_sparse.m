function [CC] = train_lda_sparse(X,G,par,tol)
% Linear Discriminant Analysis for the Small Sample Size Problem as described in
% Algorithm 1 of J. Duintjer Tebbens, P. Schlesinger: 'Improving
% Implementation of Linear Discriminant Analysis for the High Dimension/Small Sample Size
% Problem', to appear in Computational Statistics and Data Analysis in 2007. 
% Input:
%               X                 ......       (sparse) training data matrix
%               G                 ......       group coding matrix of the training data
%               test              ......       (sparse) test data matrix
%               Gtest             ......       group coding matrix of the test data
%               par               ......       if par = 0 then classification exploits sparsity too
%               tol               ......       tolerance to distinguish zero eigenvalues
% Output:
%               err               ......       Wrong classification rate (in %)
%               trafo             ......       LDA transformation vectors
%
% Reference(s): 
% J. Duintjer Tebbens, P. Schlesinger: 'Improving
% Implementation of Linear Discriminant Analysis for the High Dimension/Small Sample Size
% Problem', to appear in Computational Statistics and Data Analysis in
% 2007. 
%
% Copyright (C) by J. Duintjer Tebbens, Institute of Computer Science of the Academy of Sciences of the Czech Republic,
% Pod Vodarenskou vezi 2, 182 07 Praha 8 Liben, 18.July.2006. 
% This work was supported by the Program Information Society under project
% 1ET400300415.
%
%
% Modified for the use with Matlab6.5 by A. Schlögl, 22.Aug.2006
%
%	$Id: train_lda_sparse.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step (1)
%p = length(X(1,:));n = length(X(:,1));g = length(G(1,:));
G = sparse(G);
[n,p]=size(X); 
g = size(G,2);

for j=1:g
        nj(j) = norm(G(:,j))^2;
end
Dtild = spdiags(nj'.^(-1),[0],g,g);
Xtild = X*X';
Xtild1 = Xtild*ones(n,1);
help = ones(n,1)*Xtild1'/n - (ones(1,n)*Xtild'*ones(n,1))/(n^2); 
matrix = Xtild - Xtild1*ones(1,n)/n - help;
% eliminate non-symmetry of matrix due to rounding error:
matrix = (matrix+matrix')/2;
[V0,S] = eig(matrix);
% [s,I] = sort(diag(S),'descend');
[s,I] = sort(-diag(S)); s = -s; 

cc = sum(s<tol);

count = n-cc;
V1 = V0(:,I(1:count));
D1inv = diag(s(1:count).^(-1.0));
Dhalfinv = diag(s(1:count).^(-0.5));
Dhalf = diag(s(1:count).^(0.5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step (2)
help2 = V1*D1inv;
M1 = Dtild*G'*Xtild;
B1 = (G*(M1*(speye(n)-1/n))-help)*help2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step (3)
opts.issym = 1;opts.isreal = 1;opts.disp = 0;
%if 0, 
try,
        [V0,S,flag] = eigs(B1'*B1,g-1,'lm',opts);
        EV = Dhalfinv*V0;
        [s,I] = sort(-diag(S)); s = -s; 
        %else
catch
        % needed as long as eigs is not supported by Octave
        [V0,S] = eig(B1'*B1);
        flag   = 0;
        [s,I]  = sort(-diag(S)); s = -s(I(1:g-1));
        EV = Dhalfinv * V0(:,I(1:g-1));
        I = 1:g-1;
end;
%EV = Dhalfinv*V0;
%[s,I] = sort((diag(S)),'descend');
%[s,I] = sort(-diag(S)); s = -s; 
if flag ~= 0
        'eigs did not converge'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step (4)
for j=1:g-1,
        C(:,j) = EV(:,I(j))/norm(EV(:,I(j)));
end
cc = 0;
for j=1:g-1,
        if (1-s(j))<tol
                cc = cc+1;
                V2(:,j) = EV(:,I(j));
        else
                break
        end
end
if cc > 0
        [Q,R] = qr(V2,0);
        matrix = B1*Dhalf*Q;
        [V0,S] = eig(matrix'*matrix);
        %[s,I] = sort(diag(S),'descend');
        [s,I] = sort(-diag(S)); s = -s; 
        for j=1:cc
                C(:,j) = Q*V0(:,I(j));
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step (5)
C1 = help2*Dhalf*C;
trafo(:,1:g-1) = X'*C1 - (X'*ones(n,1))*(ones(1,n)*C1/n);
for j=1:g-1
        trafo(:,j) = trafo(:,j)/norm(trafo(:,j));
end
CC.trafo = trafo; 

if par == 0
%    X2 = full(test*X');
%    [pred] = classifs(C1,M1,X2);
        CC.C1 = C1;
        CC.M1 = M1;
        CC.X  = X;
else
%    M = Dtild*G'*X;
%    [pred] = classifs(trafo,M,test);
        CC.C1 = trafo; 
        CC.M1 = Dtild*G'*X;
end


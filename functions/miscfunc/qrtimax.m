% qrtimax() -  perform Quartimax rotation of rows of a data matrix.
%
% Usage: >> [Q,B] = qrtimax(data);   
%        >> [Q,B] = qrtimax(data,tol,'[no]reorder');
%
% Inputs:
%        data      - input matrix
%        tol       - the termination tolerance {default: 1e-4}
%        noreorder - rotate without negation/reordering
%
% Outputs:
%        B         - B=Q*A the Quartimax rotation of A
%        Q         - the orthogonal rotation matrix
%
% Author: Sigurd Enghoff, CNL / Salk Institute, 6/18/98

% Copyright (C) Sigurd Enghoff - CNL / Salk Institute, La Jolla 6/18/98
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% Reference: Jack O. Nehaus and Charles Wrigley (1954) 
% The Quartimax Method: an analytic approach to orthogonal 
% simple structure, Br J Stat Psychol, 7:81-91.

% 01-25-02 reformated help & license -ad 

function [Q,B] = qrtimax(A,tol,reorder)

if nargin < 1
	help qrtimax
	return
end

DEFAULT_TOL = 1e-4;
MAX_ITERATIONS = 50;

if nargin < 3
	reorder = 1;
elseif isempty(reorder) || reorder == 0
	reorder = 1; % set default
else
	reorder = strcmp('reorder',reorder);
end

if nargin < 2
	eps1 = DEFAULT_TOL;
	eps2 = DEFAULT_TOL;
else
	eps1 = tol;
	eps2 = tol;
end

% Do unto 'Q' what is done to A

Q = eye(size(A,1));

% Compute the cross-products of the rows of the squared loadings,
% i.e. the cost function.
%
%  ---  ---
%  \    \     2    2
%  /    /    f    f
%  ---  ---   ij   ik
%   i   j<k
%
%  See reference, p. 85, line 3

B = tril((A.^2)*(A.^2)');
crit = [sum(sum(B)) - trace(B) , 0];

% Initialize variables

inoim = 0;
iflip = 1;
ict = 0;

% Main iterative loop: keep looping while no two consecutive trials 
% satisfy tolerance constraint AND less than MAX_ITERATIONS trials 
% AND one or more rotations were performed during last trial.

while inoim < 2 & ict < MAX_ITERATIONS & iflip,
	iflip = 0;

% Run through all combinations of j and k

	for j = 1:size(A,1)-1,
		for k = j+1:size(A,1),
%
%          ---                                   ---
%          \                2     2              \       2     2  2      2   2
% fnum = 4 /    [f   f    (f   - f  )]  , fden = /    [(f   - f  )  - 4 f   f  ]
%          ---    ki  kj    ki    kj             ---     ki    kj        ki  kj
%           k                                     k
%
%  See equation (5)

			u = A(j,:) .^ 2 - A(k,:) .^2;
			v = 2 * A(j,:) .* A(k,:);
			c = sum(u .^ 2 - v .^ 2);
			d = sum(u .* v);

			fden = c;
			fnum = 2 * d;

% Skip rotation if angle is too small

			if abs(fnum) > eps1 * abs(fden)
				iflip = 1;
				angl = atan2(fnum, fden);

% Set angle of rotation according to Table I

				if fnum > 0
					if fden > 0
						angl = .25 * angl;
					else
						angl = .25 * (pi - angl);
					end
				else
					if fden > 0
						angl = .25 * (2 * pi - angl);
					else
						angl = .25 * (pi + angl);
					end
				end

% Perform rotation

				tmp    =  cos(angl) * Q(j,:) + sin(angl) * Q(k,:);
				Q(k,:) = -sin(angl) * Q(j,:) + cos(angl) * Q(k,:);
				Q(j,:) = tmp;

				tmp    =  cos(angl) * A(j,:) + sin(angl) * A(k,:);
				A(k,:) = -sin(angl) * A(j,:) + cos(angl) * A(k,:);
				A(j,:) = tmp;
			end
		end
	end

% Compute cost function.

	B = tril((A.^2)*(A.^2)');
	crit = [sum(sum(B)) - trace(B) , crit(1)];
	
	inoim = inoim + 1;
	ict = ict + 1;

	fprintf('#%d - crit = %g\n',ict,(crit(1)-crit(2))/crit(1));

% Check relative change of cost function (termination criterion).

	if (crit(1) - crit(2)) / crit(1) > eps2
		inoim = 0;
	end
end

% Reorder and negate if required. Determine new row order based on
% row norms and reorder accordingly. Negate those rows in which the
% accumulated sum is negative.

if reorder
	fprintf('Reordering rows...');
	[fnorm index] = sort(sum(A'.^2));
	Q = Q .* ((2 * (sum(A') > 0) - 1)' * ones(1, size(Q,2)));
	A = A .* ((2 * (sum(A') > 0) - 1)' * ones(1, size(A,2)));
	Q = Q(fliplr(index),:);
	A = A(fliplr(index),:);
	fprintf('\n');
else
	fprintf('Not reordering rows.\n');
end

B=A;

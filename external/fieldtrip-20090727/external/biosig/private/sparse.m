%% Copyright (C) 2000  Pascal Fleury
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%% Fake sparse function: represents sparse matrices using full matrices
%% sparse (A) [A matrix] 
%%              return A
%% sparse (m, n)
%%              return an mxn matrix of zeros
%% sparse (i, j, s [, m, n [, maxnz]])
%%              return an mxn matrix with values of vector s at
%%              locations given by the corresponding indices in
%%              vector i and j.  Any of i, j or s may be given as
%%              scalars, in which case the same value is used for
%%              all sparse entries, otherwise i, j and s must
%%              be the same length.  The maximum number of non-zero
%%              entries allowed in the sparse matrix, maxnz, is
%%              accepted for compatibility but otherwise ignored.
%%              m and n default to max(i) and max(j) respectively.
%% Example:
%%              A = sparse( [1,2,4], [2,1,3], [0.5,0.1,0.2], 4, 3)
%%              A = [  0   0.5   0
%%                    0.1   0    0
%%                     0    0    0
%%                     0    0   0.2 ]

%% Author: Pascal Fleury
%% Modified-by: Paul Kienzle (for speed and further compatibility)

function A = sparse(i,j,s,m,n,maxnz)
  if (nargin < 1 || nargin > 6 || nargin == 4)
    usage ('sparse(A) or sparse(m, n) or sparse(i, j, s [, m, n [, maxnz]])');
  end

  if ( nargin == 1 )
    %% We get a full matrix to sparsify
    A = i;

  elseif ( nargin == 2 )
    %% We get only the size of the matrix
    A = zeros(i,j);

  else
    %% ignore original shapes for indices and values
    i = i(:); j=j(:); s=s(:);

    %% assign defaults for m, n, maxnz
    if nargin < 5
      m = max(i);
      n = max(j);
    end
    sizes = [length(i) length(j) length(s)];
    nnz = max(sizes);
    if nargin < 6
      maxnz = nnz;
    end

    %% Verify that the index and value vectors are the same shape
    if ( any (sizes ~= nnz & sizes ~= 1) )
      error('sparse: index and value vectors i,j, and s must be of same size'); 
    end

    %% Verify that the indices lie within A
    if ( any ( i < 1 | i > m | j < 1 | j > n ) )
      error('sparse: index [i,j] must lie within the matrix [m,n]');
    end
      
    %% Force indices to be integers
    if ( any( i ~= floor(i) | j ~= floor(j) ) )
      warning('sparse: index [i,j] should be integers')
      i = floor(i); j = floor(j);
    end
      
    %% Ok, set the values!
    A = zeros(m,n);
    try dfi = do_fortran_indexing;
    catch dfi = 0;
    end
    try wfi = warn_fortran_indexing;
    catch wfi = 0;
    end
    unwind_protect
      do_fortran_indexing = 1;
      warn_fortran_indexing = 0;
      A((j-1)*m + i) = s;
    unwind_protect_cleanup
      do_fortran_indexing = dfi;
      warn_fortran_indexing = wfi;
    end_unwind_protect
  end

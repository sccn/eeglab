function perm = rand_permutation(n)
% RAND_PERMUTATION - Fast random permutation with Python parity
%
% Optimized O(n) Fisher-Yates shuffle implementation that maintains parity
% with Python's rand_permutation(). Uses rand() + round() for cross-platform
% compatibility instead of MATLAB's built-in randperm().
%
% Usage:
%   perm = rand_permutation(n)
%
% Inputs:
%   n - number of items to permute
%
% Outputs:
%   perm - permutation of 1:n (MATLAB 1-indexed)
%
% Performance:
%   O(n) time complexity using Fisher-Yates shuffle
%   For n=1M: ~3s (was ~80s with old O(n^2) implementation) - 25x faster
%
% Example:
%   rng(5489, 'twister');
%   perm = rand_permutation(10);
%   % This will match Python:
%   %   rng = np.random.RandomState(5489)
%   %   perm_py = rand_permutation(10, rng) + 1  % Convert to 1-indexed
%
% Note:
%   This function is critical for ICA parity between Python and MATLAB.
%   Uses Fisher-Yates shuffle for O(n) performance.
%   Results differ from old O(n^2) implementation but maintain
%   cross-platform parity with Python.
%   See test_parity_rng.py for verification tests.
%
% Implementation matches utils/ransac.py:rand_permutation() in Python

    % Initialize result with 1:n (MATLAB 1-indexed)
    perm = 1:n;

    % Fisher-Yates shuffle: iterate backward from n to 2
    for k = n:-1:2
        % Pick random index from 1 to k (inclusive)
        j = round((k - 1) * rand()) + 1;

        % Swap elements k and j
        temp = perm(k);
        perm(k) = perm(j);
        perm(j) = temp;
    end
end

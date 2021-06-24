function res = ismatlab
% true if called from Matlab; false if called from Octave
res = exist('OCTAVE_VERSION','builtin') == 0;

function [dat] = preproc_rectify(dat)

% PREPROC_RECTIFY rectifies the data, i.e. converts all samples with a
% negative value into the similar magnitude positive value
%
% Use as
%   [dat] = preproc_rectify(dat)
% where
%   dat        data matrix (Nchans X Ntime)
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%

dat = abs(dat);

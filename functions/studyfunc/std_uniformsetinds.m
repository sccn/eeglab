% std_uniformsetinds() - Check uniform channel distribution accross datasets
%
% Usage:    
%                >> boolval = std_uniformsetinds(STUDY);   
% Inputs:
%   STUDY      - EEGLAB STUDY
%
% Outputs:
%   boolval    - [0|1] 1 if uniform
%
% Authors: Arnaud Delorme, SCCN/UCSD, CERCO/CNRS, 2010-

% Copyright (C) Arnaud Delorme
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function uniformchannels = std_uniformsetinds( STUDY );

uniformchannels = 1;
for c = 1:length(STUDY.changrp(1).setinds(:))
    tmpind = cellfun(@(x)(length(x{c})), { STUDY.changrp(:).setinds });   
    if length(unique(tmpind)) ~= 1
        uniformchannels = 0;
    end;
end;

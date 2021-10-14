% scanfold() - scan folder content
%
% Usage:    
%    >> [cellres textres] = scanfold(foldname);
%    >> [cellres textres] = scanfold(foldname, ignorelist, maxdepth);
%
% Inputs:
%   foldname   - [string] name of the folder
%   ignorelist - [cell] list of folders to ignore
%   maxdepth   - [integer] maximum folder depth
%
% Outputs:
%   cellres   - cell array containing all the files
%   textres   - string array containing all the names preceded by "-a"
% 
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2009

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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

function [ cellres, textres ] = scanfold(foldname, ignorelist, maxdepth)

if nargin < 2, ignorelist = {}; end
if nargin < 3, maxdepth = 100; end
foldcontent = dir(foldname);
textres = '';
cellres = {};
if maxdepth == 0, return; end
for i = 1:length(foldcontent)
    if (exist(foldcontent(i).name) == 7 || strcmpi(foldcontent(i).name, 'functions')) && ~ismember(foldcontent(i).name, ignorelist)
        if ~strcmpi(foldcontent(i).name, '..') && ~strcmpi(foldcontent(i).name, '.')
            disp(fullfile(foldname, foldcontent(i).name));
            [tmpcellres tmpres] = scanfold(fullfile(foldname, foldcontent(i).name), ignorelist, maxdepth-1);
            textres = [ textres tmpres ];
            cellres = { cellres{:} tmpcellres{:} };
        end
    elseif length(foldcontent(i).name) > 2
        if strcmpi(foldcontent(i).name(end-1:end), '.m')
            textres = [ textres ' -a ' foldcontent(i).name ];
            cellres = { cellres{:} foldcontent(i).name };
        end
    else 
        disp( [ 'Skipping ' fullfile(foldname, foldcontent(i).name) ]);
    end
end
return;

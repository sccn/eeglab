% readeetraklocs() - read 3-D location files saved using the EETrak
%                    digitizing software.
% Usage:
%   >> CHANLOCS = readeetraklocs( filename );
%
% Inputs:
%   filename       - [string] file name
%
% Outputs:
%   CHANLOCS       - EEGLAB channel location data structure. 
%                    See help readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, Nov 2003
%
% See also: readlocs()

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function chanlocs = readeetraklocs( filename )
    
    if nargin < 1
        help readeetraklocs;
        return;
    end
    
    % read location file
    % ------------------
    locs  = loadtxt( filename );
        
    % get label names
    % ---------------
    indlabels = [];
    indpos    = [];
    for ind = 1:size(locs,1)
        if ischar(locs{ind,1}) 
            if strcmpi(locs{ind,1}, 'Labels')
                indlabels = ind;
            end
            if strcmpi(locs{ind,1}, 'Positions')
                indpos = ind;
            end
        end
    end
    if isempty(indpos) || isempty(indlabels)
        error('Could not find ''Labels'' or ''Position'' tag in electrode file');
    end
    
    % get positions
    % -------------
    if strcmp(locs(indpos+1,2),':')
        positions = locs(indpos+1:indlabels-1,3:5);
    else
        positions = locs(indpos+1:indlabels-1,1:3);
    end
    labels    = locs(indlabels+1:end,:);
        
    % create structure
    % ----------------
    for index = 1:length(labels)
        chanlocs(index).labels = labels{index};
        chanlocs(index).X      = positions{index,1};
        chanlocs(index).Y      = positions{index,2};
        chanlocs(index).Z      = positions{index,3};
    end
        
    chanlocs = convertlocs(chanlocs, 'cart2all');

% readneurolocs() - read neuroscan polar location file (.asc)
%
% Usage:
%   >> CHANLOCS = readneurolocs( filename );
%   >> CHANLOCS = readneurolocs( filename, 'key1', val1, ...);
%
% Inputs:
%   filename       - file name or matlab cell array { names x_coord y_coord }
%
% Optional inputs:
%   same as caliblocs()
%   note that if no optional input are provided, re-centering will be
%   performed automatically and re-scaling of coordinates will be
%   performed for '.asc' files (not '.dat' files).
%
% Outputs:
%   CHANLOCS       - EEGLAB channel location data structure. See
%                    help readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 4 March 2003
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

function chanlocs = readneurolocs( filename, varargin)

if nargin < 1
    help readneurolocs;
    return;
end
if nargin < 2
    plottag = 0;
end

% read location file
% ------------------
if ischar(filename)
    locs  = loadtxt( filename, 'delim', 9 );
end

if ~ischar(filename) || locs{1,1}(1) == ';' || size(locs,2) < 5
    if ~ischar(filename)
        names = filename{1};
        x     = filename{2};
        y     = filename{3};
    else
        if locs{1,1}(1) == ';'
            % remove trailing control channels
            % --------------------------------
            while isnumeric( locs{end,1} ) & locs{end,1} ~= 0
                locs  = locs(1:end-1,:);
            end

            % find first numerical index
            % --------------------------
            index = 1;
            while ischar( locs{index,1} )
                index = index + 1;
            end

            % extract location array
            % ----------------------
            nchans = size( locs, 1 ) - index +1;
            chans  = [locs{end-nchans+1:end, 1:5}];
            chans  = reshape(chans,nchans,5);               %% Added this line in order to get x = chans(:,3)
            names  = locs(end-nchans*2+1: end-nchans, 2);
            for index = 1:length(names)
                if ~ischar(names{index})
                    names{index} = int2str(names{index});
                end
            end
            x = chans(:,3);
            y = -chans(:,4);
        else
            [tmp2 tmpind] = sort( [ locs{:,1} ]);
            locs = locs(tmpind,:);
            y      = [ locs{:,end} ];
            x      = [ locs{:,end-1} ];
            x      = x/513.1617*44;
            y      = y/513.1617*44;
            names = locs(:,end-2);
        end
    end

    % second solution using angle
    % ---------------------------
    [phi,theta] = cart2pol(x, y);
    phi = phi/pi*180;

    % convert to other types of coordinates
    % -------------------------------------
    labels = names';
    chanlocs = struct('labels', labels, 'sph_theta_besa', mattocell(theta)', 'sph_phi_besa', mattocell(phi)');      %% labels instead of labels(:) 
    chanlocs = convertlocs( chanlocs, 'sphbesa2all');

    for index = 1:length(chanlocs)
        chanlocs(index).labels = num2str(chanlocs(index).labels);
    end

    % re-calibration
    % --------------
    chanlocs = adjustlocs(chanlocs, 'autoscale', 'on', 'autorotate', 'off', varargin{:});

else % 5 rows, xyz positions
    try
        for index = 1:size(locs,1)
            locs{index,3} = - locs{index,3};
        end
        chanlocs = struct('labels', locs(:,1), 'type', locs(:,2), 'X', locs(:,4), 'Y', locs(:,3), 'Z', locs(:,5));
        chanlocs = convertlocs( chanlocs, 'cart2all');
    catch
        chanlocs = readlocs(filename, 'filetype', 'custom', 'format', { 'labels' 'ignore' '-Y' 'X' 'Z' });
    end
end
    

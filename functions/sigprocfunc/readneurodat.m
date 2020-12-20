% readneurodat() - read neuroscan location file (.dat)
%
% Usage:
%   >> [ CHANLOCS labels theta theta ] = readneurodat( filename );
%
% Inputs:
%   filename       - file name or matlab cell array { names x_coord y_coord }
%
% Outputs:
%   CHANLOCS       - [structure] EEGLAB channel location data structure. See
%                    help readlocs()
%   labels         - [cell arrya] channel labels
%   theta          - [float array]array containing 3-D theta angle electrode
%                    position (in degree)
%   phi            - [float array]array containing 3-D phi angle electrode
%                    position (in degree)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 28 Nov 2003
%
% See also: readlocs(), readneurolocs()

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

function [chanlocs, labels, theta, phi] = readneurodat(filename);
    
    if nargin < 1
        help readneurodat;
        return;
    end
    
    % enter file name here
    % --------------------
    %tmp = loadtxt('/home/ftp/pub/locfiles/neuroscan/cap128.dat');
    %tmp = loadtxt('/home/arno/temp/quik128.DAT');
    tmp = loadtxt(filename);

    % resort electrodes
    % -----------------
    [tmp2 tmpind] = sort(celltomat(tmp(:,1))');
    tmp = tmp(tmpind,:);

    % convert to polar coordinates
    % ----------------------------
    %figure; plot(celltomat(tmp(:,2)), celltomat(tmp(:,3)), '.');
    [phi,theta] = cart2pol(celltomat(tmp(:,end-1)), celltomat(tmp(:,end)));
    theta = theta/513.1617*44;
    phi   = phi/pi*180;

    % convert to other types of coordinates
    % -------------------------------------
    labels = tmp(:,end-2)';
    chanlocs = struct('labels', labels, 'sph_theta_besa', mattocell(theta)', 'sph_phi_besa', mattocell(phi)');
    chanlocs = convertlocs( chanlocs, 'sphbesa2all');
    for index = 1:length(chanlocs)
        chanlocs(index).labels = num2str(chanlocs(index).labels);
    end
    theta = theta/pi*180;
    fprintf('Note: .dat file contain polar 2-D coordinates. EEGLAB will use these coordinates\n');
    fprintf('      to recreated 3-D coordinates.\n');
    fprintf('\n');
    fprintf('Warning: if the channels are clustered in the middle or oddly spaced removed\n');
    fprintf('         the peripheral channels and then used the optimize head function\n');
    fprintf('         in the Channel Locations GUI. It seems that sometimes the peripherals\n');
    fprintf('         are throwing off the spacing.\n');
    

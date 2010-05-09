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

function [chanlocs, labels, theta, phi] = readneurodat(filename);
    
    if nargin < 1
        help readneurodat;
        return;
    end;
    
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
    end;
    theta = theta/pi*180;

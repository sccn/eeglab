% readneurolocs() - read neuroscan polar location file (.asc)
%
% Usage:
%   >> CHANLOCS = readneurolocs( filename );
%   >> CHANLOCS = readneurolocs( filename, 'key1', val1, 'key2', val2, ...);
%
% Inputs:
%   filename       - file name or matlab cell array { names x_coord y_coord }
%
% Optional inputs:
%   'hcalib'       - [cell array] name and phi angle of two electrodes with 
%                    left/right simetry. Ex: { 'c3' 44 'c4' 44 }.
%   'vcalib'       - [cell array] name and phi angle of two electrodes with 
%                    bottom/up simetry. Ex: { 'fz' 44 'pz' 44 }.
%   'autocalib'    - ['on'|'off'] attempt to automatically detect 'c3, 'c4'
%                    'fz', 'pz' to calibrate the coordinates. Default is 'on'
%                    for '.asc' files and 'off' for '.dat' files.
%   'plot'         - ['on'|'off'] if 'on', plot the electrode locations. 
%                    Default 'off'.
%
% Outputs:
%   CHANLOCS       - EEGLAB channel location data structure. See
%                    help readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 4 March 2003
%
% See also: readlocs()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.7  2003/12/01 02:31:01  arno
% correct slight innacuracy when reading file
%
% Revision 1.6  2003/08/05 18:24:04  arno
% header
%
% Revision 1.5  2003/07/29 21:25:20  arno
% debug C3-C4
%
% Revision 1.4  2003/07/29 21:05:53  arno
% handle cell array for channel locations
%
% Revision 1.3  2003/06/13 00:53:47  arno
% debuging
%
% Revision 1.2  2003/03/04 20:04:17  arno
% debuging
%
% Revision 1.1  2003/03/04 19:18:35  arno
% Initial revision
%

function chanlocs = readneurolocs( filename, varargin)
    
    if nargin < 1
        help readneurolocs;
        return;
    end;
    if nargin < 2
        plottag = 0;
    end;
    
    % check input parameters
    % ----------------------
    g = finputcheck( varargin, { 'hcalib'    'cell'  []   {};
                                 'vcalib'    'cell'  []   {};
                                 'autocalib' 'string'  { 'on' 'off' 'auto' }   'auto';
                                 'plot'      'string'  { 'on' 'off' }   'off' });
    if isstr(g), error(g); end;
    
    % read location file
    % ------------------
    if isstr( filename )
        locs  = loadtxt( filename );
    
        if locs{1,1}(1) == ';'
            % remove trailing control channels
            % --------------------------------
            while isnumeric( locs{end,1} ) & locs{end,1} ~= 0
                locs  = locs(1:end-1,:);
            end;    
            
            % find first numerical index
            % --------------------------
            index = 1;
            while isstr( locs{index,1} )
                index = index + 1;
            end;
            
            % extract location array
            % ----------------------   
            nchans = size( locs, 1 ) - index +1;
            chans  = cell2mat(locs(end-nchans+1:end, 1:5));
            names  = locs(end-nchans*2+1: end-nchans, 2);
            for index = 1:length(names)
                if ~isstr(names{index})
                    names{index} = int2str(names{index});
                end;
            end;
            x = chans(:,3);
            y = chans(:,4);
            if strcmpi(g.autocalib, 'auto')
                g.autocalib = 'on';
            end;
        else
            [tmp2 tmpind] = sort(cell2mat(locs(:,1))');
            locs = locs(tmpind,:);
            y      = cell2mat(locs(:,end));
            x      = cell2mat(locs(:,end-1));
            x      = x/513.1617*44;
            y      = y/513.1617*44;
            names = locs(:,end-2);
            if strcmpi(g.autocalib, 'auto')
                g.autocalib = 'off';
            end;
        end;
    else
        names = filename{1};
        x     = filename{2};
        y     = filename{3};
        if strcmpi(g.autocalib, 'auto')
            g.autocalib = 'on';
        end;
    end;
    
    % auto calibration
    % ----------------
    if strcmpi(g.autocalib, 'on')
        indexc3 = strmatch( 'c3', lower(names'), 'exact' );
        indexc4 = strmatch( 'c4', lower(names'), 'exact' );
        if ~isempty(indexc3) & ~isempty(indexc4)
            g.hcalib = { 'c3' 45 'c4' 45 };
        end;
        indexfz = strmatch( 'fz', lower(names'), 'exact' );
        indexpz = strmatch( 'pz', lower(names'), 'exact' );
        if ~isempty(indexfz) & ~isempty(indexpz)
            g.vcalib = { 'fz' 45 'pz' 45 };
        end;
        if isempty(g.vcalib) & isempty(g.hcalib)
            disp('No electrodes found for position re-calibration')
        else
            disp('Landmark electrodes found for position re-calibration')
        end;
    end;
    
    % calibrate
    % ---------
    if ~isempty(g.hcalib)
        index1 = strmatch(lower( g.hcalib{1}), lower(names'), 'exact' );
        index2 = strmatch( lower(g.hcalib{3}), lower(names'), 'exact' );
        if isempty(index1) | isempty(index2)
            error('Electrode not found for horizontal calibration');
        end;
        disp([ 'Readneurolocs: electrode location re-calibrated along x axis using ''' ...
               g.hcalib{1} ''' and '''  g.hcalib{3} '''']);
        centerx = (x(index1)+x(index2))/2;
        centery = (y(index1)+y(index2))/2;
        x = - x + centerx;
        y = - y + centery;
        x = x/x(index1)*g.hcalib{2};
    end;
    if ~isempty(g.vcalib)
        index1 = strmatch( lower(g.vcalib{1}), lower(names'), 'exact' );
        index2 = strmatch( lower(g.vcalib{3}), lower(names'), 'exact' );
        if isempty(index1) | isempty(index2)
            error('Electrode not found for horizontal calibration');
        end;
        disp([ 'Readneurolocs: electrode location re-calibrated along y axis using ''' ...
               g.vcalib{1} ''' and '''  g.vcalib{3} '''']);
        if isempty(g.hcalib)
            centerx = (x(index1)+x(index2))/2;
            centery = (y(index1)+y(index2))/2;
            x = - x + centerx;
            y = - y + centery;
        end;
        y = y/y(index1)*g.vcalib{2};
    end;
    
    % plot all channels
    % -----------------
    if strcmpi(g.plot, 'on')
        figure;
        for index = 1:length(x)
            plot( x(index), y(index), '+');
            hold on; % do not erase
            text( x(index)+0.01, y(index), int2str(index));
        end;
    end;
       
    % second solution using angle
    % ---------------------------
    [phi,theta] = cart2pol(-x, y);
    phi = phi/pi*180;
    
    % convert to other types of coordinates
    % -------------------------------------
    labels = names';
    chanlocs = struct('labels', labels, 'sph_theta_besa', mat2cell(theta)', 'sph_phi_besa', mat2cell(phi)');
    chanlocs = convertlocs( chanlocs, 'sphbesa2all');
    for index = 1:length(chanlocs)
        chanlocs(index).labels = num2str(chanlocs(index).labels);
    end;

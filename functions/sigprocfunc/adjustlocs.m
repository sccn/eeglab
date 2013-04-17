% adjustlocs() - read neuroscan polar location file (.asc)
%
% Usage:
%   >> chanlocs = adjustlocs( chanlocs );
%   >> chanlocs = adjustlocs( chanlocs, 'key1', val1, 'key2', val2, ...);
%
% Inputs:
%   chanlocs       - EEGLAB channel location data structure. See
%                    help readlocs()
%
% Optional inputs:
%   'center'       - [cell array] one or several electrode names to compute
%                    center location (if several electrodes are provided, the
%                    function use their iso-barycenter). Ex: { 'cz' } or 
%                    { 'fz' 'pz' }.
%   'autocenter'   - ['on'|'off']  attempt to automatically detect all symetrical
%                    electrode in the 10-20 system to compute the center. 
%                    Default: 'on'.
%   'rotate'       - [cell array] name and planar angle of landmark electrodes
%                    to use as template planar rotation. Ex: { 'c3', 90 }.
%   'autorotate'   - ['on'|'off'] attempt to automatically detect 
%                    electrode in the 10-20 system to compute the average
%                    planar rotation. Default 'on'.
%   'scale'        - [cell array] name and phi angle of electrodes along
%                    horizontal central line. Ex: { 'c3' 44 }. Scale uniformly
%                    all direction. Use 'hscale' and 'vscale' for scaling x and
%                    y axis independently.
%   'hscale'       - [cell array] name and phi angle of electrodes along
%                    honrizontal central line. Ex: { 'c3' 44 }.
%   'vscale'       - [cell array] name and phi angle of one electrodes along
%                    central vertical axis. Ex: { 'fz' 44 }.
%   'autoscale'    - ['on'|'off'] automatic scaling with 10-20 system used as
%                    reference. Default is 'on'.
%   'uniform'      - ['on'|'off'] force the scaling to be uniform along the X
%                    and the Y axis. Default is 'on'.
%   'coordinates'  - ['pol'|'sph'|'cart'] use polar coordinates ('pol'), sperical
%                    coordinates ('sph') or cartesian ('cart'). Default is 'sph'.
%                    (Note that using polar and spherical coordinates have the 
%                    same effect, except that units are different). 
%                    Default is 'sph'.
%
% Outputs:
%   chanlocs       - EEGLAB channel location data structure. See
%                    help readlocs()
%
% Note: operations are performed in the following order, first re-centering
%       then planar rotation and finally sperical re-scaling.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 1 Dec 2003
%
% See also: readlocs()

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

function chanlocs = adjustlocs( chanlocs, varargin)
    
    if nargin < 1
        help adjustlocs;
        return;
    end;
    
    % check input parameters
    % ----------------------
    g = finputcheck( varargin, { 'hscale'     'cell'  []   {};
                                 'vscale'     'cell'  []   {};
                                 'scale'      'cell'  []   {};
                                 'center'     'cell'  []   {};
                                 'rotate'     'cell'  []   {};
                                 'autoscale'  'string'  { 'on','off' }    'off';
                                 'autocenter' 'string'  { 'on','off' }    'on';
                                 'autorotate' 'string'  { 'on','off' }    'on';
                                 'uniform'    'string'  { 'on','off' }    'on';
                                 'coordinates' 'string' { 'pol','sph','cart' } 'sph' });
    if ischar(g), error(g); end;
    
    names = { chanlocs.labels };
    
    % auto center
    % -----------
    if strcmpi(g.autocenter, 'on') & isempty(g.center)
        disp('Reading template 10-20 file');
        locs1020 = readlocs('eeglab1020.ced');
        
        % scan electrodes for horiz pos
        % -----------------------------
        tmpnames      = lower(names);
        tmpnames1020  = lower({ locs1020.labels });
        [tmp indelec] = intersect_bc(tmpnames1020, tmpnames);

        % remove non-symetrical electrodes
        % --------------------------------
        if ~isempty(indelec)
            if find(indelec == 79), indelec(end+1) = 80; end; % for Cz
            ind2remove = [];
            for index = 1:length(indelec)
                if mod(indelec(index),2)
                    if ~ismember(indelec(index)+1, indelec)
                        ind2remove = [ ind2remove index ];
                    end;
                else
                    if ~ismember(indelec(index)-1, indelec)
                        ind2remove = [ ind2remove index ];
                    end;
                end;
            end;
            indelec(ind2remove) = [];
            if find(indelec == 80), indelec(end) = []; end; % for Cz
        
            g.center = tmpnames1020(indelec);
        end;
        if isempty(g.center)
            disp('No electrodes found for auto-centering')
        else
            disp([ num2str(length(g.center)) ' landmark electrodes found for position auto-centering' ])
        end;
    end;

    % auto rotate
    % -----------
    if strcmpi(g.autorotate, 'on') & isempty(g.rotate)
        if exist('locs1020') ~= 1
            disp('Reading template 10-20 file');
            locs1020 = readlocs('eeglab1020.ced');
        end;
                
        % scan electrodes for horiz pos
        % -----------------------------
        tmpnames      = lower(names);
        tmpnames1020  = lower({ locs1020.labels });
        tmptheta1020 = { locs1020.theta };
        [tmp indelec] = intersect_bc(tmpnames1020(1:end-1), tmpnames); % do not use cz

        g.rotate(1:2:2*length(indelec))   = tmpnames1020 (indelec);
        g.rotate(2:2:2*length(indelec)+1) = tmptheta1020(indelec);
        
        if isempty(g.rotate)
            disp('No electrodes found for auto planar rotation')
        else
            disp([ num2str(length(g.rotate)/2) ' landmark electrodes found for auto planar rotation' ])
        end;
    end;

    % auto scale
    % ----------
    if strcmpi(g.autoscale, 'on') & isempty(g.hscale) & isempty(g.vscale)
        if exist('locs1020') ~= 1
            disp('Reading template 10-20 file');
            locs1020 = readlocs('eeglab1020.ced');
        end;
        
        if strcmpi(g.uniform, 'off')
            % remove all vertical electrodes for horizontal scaling
            % -----------------------------------------------------
            theta  = [ locs1020.theta ];
            indxh  = find(abs(theta) == 90);
            indxv  = union_bc(find(theta == 0) , find(theta == 180));
            locs1020horiz = locs1020(indxh);
            locs1020vert  = locs1020(indxv);
            
            % scan electrodes for horiz pos
            % -----------------------------
            tmpnames      = lower(names);
            tmpnames1020  = lower({ locs1020horiz.labels });
            tmpradius1020 = { locs1020horiz.radius };
            [tmp indelec] = intersect_bc(tmpnames1020, tmpnames);
            
            if isempty(indelec)
                disp('No electrodes found for horiz. position spherical re-scaling')
            else
                disp([ num2str(length(indelec)) ' landmark electrodes found for horiz. position spherical re-scaling' ]);
                if strcmpi(g.coordinates, 'cart')
                    g.hscale(2:2:2*length(indelec)+1) = [ locs1020horiz.Y ];
                else
                    g.hscale(2:2:2*length(indelec)+1) = tmpradius1020(indelec);
                end;
            end;
            
            % scan electrodes for vert. pos
            % -----------------------------
            tmpnames1020  = lower({ locs1020vert.labels });
            tmpradius1020 = { locs1020vert.radius };
            [tmp indelec] = intersect_bc(tmpnames1020, tmpnames);
            
            if isempty(indelec)
                disp('No electrodes found for vertical position spherical re-scaling')
            else
                disp([ num2str(length(indelec)) ' landmark electrodes found for vertical spherical re-scaling' ]);
                g.vscale(1:2:2*length(indelec))   = tmpnames1020 (indelec);
                if strcmpi(g.coordinates, 'cart')
                    g.vscale(2:2:2*length(indelec)+1) = [ locs1020vert.X ];
                else
                    g.vscale(2:2:2*length(indelec)+1) = tmpradius1020(indelec);
                end;
            end;    
        else
            % uniform scaling
            % ---------------
            tmpnames      = lower(names);
            tmpnames1020  = lower({ locs1020.labels });
            tmpradius1020 = { locs1020.radius };
            [tmp indelec] = intersect_bc(tmpnames1020, tmpnames);
            
            if isempty(indelec)
                disp('No electrodes found for uniform spherical re-scaling')
            else
                disp([ num2str(length(indelec)) ' landmark electrodes found for uniform spherical re-scaling' ]);
                g.scale(1:2:2*length(indelec))   = tmpnames1020 (indelec);
                if strcmpi(g.coordinates, 'cart')
                    tmpabsxyz = mattocell(abs([ locs1020.X ]+j*[ locs1020.Y ]));
                    g.scale(2:2:2*length(indelec)+1) = tmpabsxyz(indelec);
                else
                    g.scale(2:2:2*length(indelec)+1) = tmpradius1020(indelec);
                end;
            end;
        end;
        if strcmpi(g.coordinates, 'sph')
            g.coordinates = 'pol'; % use polar coordinates for scaling
        end;
    end;
    
    % get X and Y coordinates
    % -----------------------
    if strcmpi(g.coordinates, 'sph') | strcmpi(g.coordinates, 'pol')
        [X Y] = pol2cart( [ chanlocs.theta ]/180*pi, [ chanlocs.radius ]); Z = 1;
        if strcmpi(g.coordinates, 'sph')
            X = X/0.25*46;
            Y = Y/0.25*46;
        end;
    else
        X = [ chanlocs.X ];
        Y = [ chanlocs.Y ];
        Z = [ chanlocs.Z ];
    end;
    
    % recenter
    % --------
    if ~isempty(g.center)
        for index = 1:length(g.center)
            tmpindex = strmatch( lower(g.center{index}), lower(names), 'exact' );
            if isempty(tmpindex)
                error(['Electrode ''' g.center{index} ''' not found for re-centering']);
            end;
            indexelec(index) = tmpindex;
        end;
        showmsg('Using electrode', 'for re-centering', g.center);
        centerx = mean(X(indexelec));
        centery = mean(Y(indexelec));
        X = X - centerx;
        Y = Y - centery;
    end;

    % planar rotation
    % ---------------    
    if ~isempty(g.rotate)
        % find electrodes
        % ---------------
        clear elec;
        for index = 1:2:length(g.rotate)
            tmpindex = strmatch( lower(g.rotate{index}), lower(names), 'exact' );
            if isempty(tmpindex)
                error(['Electrode ''' g.rotate{index} ''' not found for left-right scaling']);
            end;
            elec((index+1)/2) = tmpindex;
        end;
        vals = [ g.rotate{2:2:end} ];
        
        % compute average scaling factor
        % ------------------------------
        [ allangles tmp ] = cart2pol(X(elec), Y(elec));
        allangles = allangles/pi*180;
        diffangle = allangles - vals;
        %diffangle2 = allangles + vals;
        %if abs(diffangle1) > abs(diffangle2), diffangle = diffangle2;
        %else                                  diffangle = diffangle1;
        %end;
        tmpind    = find(diffangle >  180); diffangle(tmpind) = diffangle(tmpind)-360;
        tmpind    = find(diffangle < -180); diffangle(tmpind) = diffangle(tmpind)+360;
        anglerot  = mean(diffangle);
        tmpcplx = (X+j*Y)*exp(-j*anglerot/180*pi);
        X = real(tmpcplx);
        Y = imag(tmpcplx);
        showmsg('Using electrode', ['for planar rotation (' num2str(anglerot,2) ' degrees)'], g.rotate(1:2:end));
    end;
    
    % computing scaling factors
    % -------------------------
    if ~isempty(g.scale)
        % find electrodes
        % ---------------
        clear elec;
        for index = 1:2:length(g.scale)
            tmpindex = strmatch( lower(g.scale{index}), lower(names), 'exact' );
            if isempty(tmpindex)
                error(['Electrode ''' g.scale{index} ''' not found for left-right scaling']);
            end;
            elec((index+1)/2) = tmpindex;
        end;
        vals = [ g.scale{2:2:end} ];
        
        % compute average scaling factor
        % ------------------------------
        nonzero = find(vals > 0);
        hscalefact = mean(abs(Y(elec(nonzero))+j*X(elec(nonzero)))./vals(nonzero)); % *46/0.25; %/44/0.25;
        vscalefact = hscalefact;
        showmsg('Using electrode', ['for uniform spherical re-scaling (x' num2str(1/hscalefact,4) ')'], g.scale(1:2:end));
    else
        if ~isempty(g.hscale)
            % find electrodes
            % ---------------
            clear elec;
            for index = 1:2:length(g.hscale)
                tmpindex = strmatch( lower(g.hscale{index}), lower(names), 'exact' );
                if isempty(tmpindex)
                    error(['Electrode ''' g.hscale{index} ''' not found for left-right scaling']);
                end;
                elec((index+1)/2) = tmpindex;
            end;
            vals = [ g.hscale{2:2:end} ];
            showmsg('Using electrode', [ 'for left-right spherical re-scaling (x' num2str(1/hscalefact,4) ')'], g.hscale(1:2:end));
            
            % compute average scaling factor
            % ------------------------------
            hscalefact = mean(abs(Y(elec))./vals); % *46/0.25; %/44/0.25;
            if isempty(g.vscale)
                vscalefact =  hscalefact;
            end;
        end;
        if ~isempty(g.vscale)
            % find electrodes
            % ---------------
            clear elec;
            for index = 1:2:length(g.vscale)
                tmpindex = strmatch( lower(g.vscale{index}), lower(names), 'exact' );
                if isempty(tmpindex)
                    error(['Electrode ''' g.vscale{index} ''' not found for rear-front scaling']);
                end;
                elec((index+1)/2) = tmpindex;
            end;
            vals = [ g.vscale{2:2:end} ];
            showmsg('Using electrode', ['for rear-front spherical re-scaling (x' num2str(1/vscalefact,4) ')'], g.vscale(1:2:end));
            
            % compute average scaling factor
            % ------------------------------
            vscalefact = mean(abs(X(elec))./vals); % *46/0.25; %/44/0.25;
            if isempty(g.vscale)
                hscalefact =  vscalefact;
            end;
        end;
    end;
    
    % uniform?
    % --------
    if strcmpi(g.uniform, 'on') & ( ~isempty(g.vscale) | ~isempty(g.hscale))
        disp('uniform scaling: averaging left-right and rear-front scaling factor');
        hscalefact = mean([hscalefact vscalefact]);
        vscalefact = hscalefact;
    end;
    
    % scaling data
    % ------------
    if ~isempty(g.vscale) | ~isempty(g.hscale) | ~isempty(g.scale)
        Y = Y/hscalefact;
        X = X/vscalefact;
        Z = Z/((hscalefact+vscalefact)/2);
    end;
    
    % updating structure
    % ------------------
    if strcmpi(g.coordinates, 'sph') |  strcmpi(g.coordinates, 'pol') 
        [phi,theta] = cart2pol(Y, X);
        phi = phi/pi*180;
        if strcmpi(g.coordinates, 'pol')
            theta = theta/0.25*46;
        end;
        
        % convert to other types of coordinates
        % -------------------------------------
        labels = names';
        chanlocs = struct('labels', names, 'sph_theta_besa', mattocell(theta), ...
                'sph_phi_besa', mattocell(phi), 'sph_radius', { chanlocs.sph_radius });
        chanlocs = convertlocs( chanlocs, 'sphbesa2all');
    else
        for index = 1:length(chanlocs)
            chanlocs(index).X = X(index);
            chanlocs(index).Y = Y(index);
            chanlocs(index).Z = Z(index);
        end;
        chanlocs = convertlocs(chanlocs, 'cart2all');
    end;
    
function showmsg(begmsg, endmsg, struct);
    if length(struct) <= 1
        disp([ begmsg ' ''' struct{1} ''' ' endmsg]);
    elseif length(struct) <= 2
        disp([ begmsg ' ''' struct{1} ''' and ''' struct{2}  ''' ' endmsg]);
    elseif length(struct) <= 3
        disp([ begmsg ' ''' struct{1} ''', ''' struct{2} ''' and ''' struct{3}  ''' ' endmsg]);
    else
        disp([ begmsg ' ''' struct{1} ''', ''' struct{2} ''', ''' struct{3}  ''' ... ' endmsg]);
    end;
    

% readneurolocs() - read neuroscan polar location file (.asc)
%
% Usage:
%   >> CHANLOCS = readneurolocs( filename );
%   >> CHANLOCS = readneurolocs( filename, plot, czindex, fzindex, c3index );
%
% Inputs:
%   filename       - file name or matlab cell array { names x_coord y_coord }
%
% Optional inputs:
%   plot           - [0|1] if 1, plot the electrode locations. Default 0.
%   czindex        - index of electrode Cz if the label is not defined 
%                    in the file. Default is [] (find electrode label
%                    automatically).
%   fzindex        - index of electrode Fz if the label is not defined 
%                    in the file. (for 3D accurate conversion of electrode 
%                    positions). Set it to [] so that the data will not
%                    be rescaled along the y axis.
%   c3index        - index of electrode c3 if the label is not defined 
%                    in the file. (for 3D accurate conversion of electrode 
%                    positions). Set it to [] so that the data will not
%                    be rescaled along the x axis.
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

function chanlocs = readneurolocs( filename, plottag, indexcz, indexfz)
    
    if nargin < 1
        help readneurolocs;
        return;
    end;
    if nargin < 2
        plottag = 0;
    end;
    
    % read location file
    % ------------------
    if isstr( filename )
        locs  = loadtxt( filename );
    
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
    else
        names = filename{1};
        x     = filename{2};
        y     = filename{3};
    end;
    
    % find Cz and Fz
    % --------------
    if exist('indexcz') ~= 1
        indexcz = strmatch( 'cz', lower(names'), 'exact' );
        if ~isempty(indexcz), 
            centerx = x(indexcz);
            centery = y(indexcz);
        end;
    end;
    if exist('indexc3') ~= 1
        indexc3 = strmatch( 'c3', lower(names'), 'exact' );
        if exist('indexc4') ~= 1
            indexc4 = strmatch( 'c4', lower(names'), 'exact' );
            if exist('centerx') ~= 1 & ~isempty(indexc4) & ~isempty(indexc3)
                centerx = (x(indexc3)+x(indexc4))/2;
                centery = (y(indexc3)+y(indexc4))/2;
            end;
        end;
    end;
    if exist('indexfz') ~= 1
        indexfz = strmatch( 'fz', lower(names'), 'exact' );
        if exist('indexpz') ~= 1
            indexpz = strmatch( 'pz', lower(names'), 'exact' );
            if exist('centerx') ~= 1 & ~isempty(indexpz) & ~isempty(indexfz)
                centerx = (x(indexfz)+x(indexpz))/2;
                centery = (y(indexfz)+y(indexpz))/2;
            end;
        end;
    end;
    if exist('centerx') ~= 1 
        error('Unable to find electrode names for rescaling, define the electrode index manually in readneurolocs()');
    end;
    
    % plot all channels
    % -----------------
    if plottag
        figure;
        for index = 1:length(x)
            plot( x(index), y(index), '+');
            hold on; % do not erase
            text( x(index)+0.01, y(index), int2str(index));
        end;
    end;
    
    x = - x + centerx;
    y = - y + centery;

    % shrink coordinates
    % ------------------
    if ~isempty(indexfz)
        yy = y/y(indexfz)*45;
        y  = y/y(indexfz)*0.25556;
        disp('Readneurolocs: electrode location scaled along the y axis for standard Fz location');
        if isempty(indexc3), indexc3 = indexfz; end;
    end;
    if ~isempty(indexc3)
        xx = x/x(indexc3)*45;
        x  = x/x(indexc3)*0.25556;
        disp('Readneurolocs: electrode location scaled along the x axis for standard C3 location');
    end;

    % convert to polar
    % ----------------
    [theta r] = cart2pol(x, y);    
    
    % create structure
    % ----------------
    for index = 1:length(names)
        chanlocs(index).labels = names{index};
        chanlocs(index).theta  = theta(index)/pi*180-90;
        chanlocs(index).radius = r(index);
    end;
    % return;
   
    % second solution using angle
    % ---------------------------
    [phi,theta] = cart2pol(-xx, yy);
    phi = phi/pi*180;
    
    % convert to other types of coordinates
    % -------------------------------------
    labels = names';
    chanlocs = struct('labels', labels, 'sph_theta_besa', mat2cell(theta)', 'sph_phi_besa', mat2cell(phi)');
    chanlocs = convertlocs( chanlocs, 'sphbesa2all');
    for index = 1:length(chanlocs)
        chanlocs(index).labels = num2str(chanlocs(index).labels);
    end;

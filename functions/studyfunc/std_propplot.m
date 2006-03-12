% std_propplot() - Commandline function to plot cluster properties. 
%                  Displays either mean cluster ERSPs, or all cluster component ERSPs 
%                  plut the mean cluster ERSP in one figure per condition.
%                  The ERSPs can be plotted only when component ERSPs have been computed 
%                  and saved in the EEG datasets of the STUDY.
%                  These may be computed during pre-clustering using the gui-based function
%                  pop_preclust() or via the equivalent commandline functions eeg_createdata() 
%                  and eeg_preclust(). This function is called by pop_clustedit().
% Usage:    
%              >> [STUDY] = std_propplot(STUDY, ALLEEG, clusters);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets in the STUDY. 
%                ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Optional inputs:
%   clusters   - [numeric vector]  -> cluster numbers to plot.
%                            'all' -> plot all clusters in STUDY {default: 'all'}.
% Outputs:
%   STUDY      - the input STUDY set structure modified with the plotted cluster 
%                mean mproperties, to allow quick replotting (unless cluster means 
%                already existed in the STUDY).  
% Example:
%              >> [STUDY] = std_propplot(STUDY,ALLEEG, 5);
%                 % Plot mean properties for cluster 5 in one figure. 
%
% See also:  pop_clustedit() 
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, July, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, July 12, 2005, hilit@sccn.ucsd.edu
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
% Revision 1.8  2006/03/09 18:58:25  scott
% help msg -sm
%
% Revision 1.7  2006/03/08 21:03:37  arno
% rename func
%
% Revision 1.6  2006/03/08 20:56:31  arno
% rename func
%
% Revision 1.5  2006/03/07 19:41:41  arno
% warning text and message
%
% Revision 1.4  2006/03/07 19:35:31  arno
% ploting
%

function STUDY = std_propplot(STUDY, ALLEEG,  varargin)
icadefs;

% Set default values
cls = 1:length(STUDY.cluster); % plot all clusters in STUDY

if length(varargin) > 0
    if length(varargin) == 1, varargin{2} = varargin{1}; end; % backward compatibility
    if isnumeric(varargin{2})
        cls = varargin{2};
    elseif isstr(varargin{2}) & strcmpi(varargin{2}, 'all')
        cls = 1:length(STUDY.cluster);
    else
        error('std_propplot: clusters input is either specific clusters (numeric vector) or keyword ''all''.');
    end
end

len = length(cls);
% Plot clusters mean properties
for k = 1: len
    if k == 1
        try 
            % optional 'CreateCancelBtn', 'delete(gcbf); error(''USER ABORT'');', 
            h_wait = waitbar(0,['Computing cluster properties ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Computing cluster properties ...'],'position', [300, 200, 300, 48]);
        end
    end  
    warningon = 0;
    figure
    orient tall
    set(gcf,'Color', BACKCOLOR);
    subplot(2,3,1),
    try,
        STUDY = std_topoplot(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off');
    catch
        axis off; text(0.5, 0.5, 'No scalp information', 'horizontalalignment', 'center');
        warningon = 1;
    end
    waitbar(k/(len*6),h_wait)
    subplot(2,3,2),
    try,
        STUDY = std_erpplot(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off');
    catch
        axis off; text(0.5, 0.5, 'No ERP information', 'horizontalalignment', 'center');
        warningon = 1;
    end
    waitbar((k*2)/(len*6),h_wait)
    subplot(2,3,3),
    try,
        [STUDY] = std_erspplot(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off' );
    catch, 
        axis off; text(0.5, 0.5, 'No ERSP information', 'horizontalalignment', 'center');
        warningon = 1;
    end
    waitbar((k*3)/(len*6),h_wait)
    axes('unit', 'normalized', 'position', [0.1 0.16 0.2 0.28]); %subplot(2,3,4),
    try,
        STUDY = std_dipplot(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'apart', 'figure', 'off'); set(gcf,'Color', BACKCOLOR);
    catch
        axis off; text(0.5, 0.5, 'No dipole information', 'horizontalalignment', 'center');
        warningon = 1;
    end
    waitbar((k*4)/(len*6),h_wait)
    subplot(2,3,5),
    try,
        STUDY = std_specplot(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off');
    catch
        axis off; text(0.5, 0.5, 'No spectrum information', 'horizontalalignment', 'center');
        warningon = 1;
    end
    waitbar((k*5)/(len*6),h_wait)
    subplot(2,3,6),
    try,
        [STUDY] = std_itcplot(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off' );
    catch, 
        axis off; text(0.5, 0.5, 'No ITC information', 'horizontalalignment', 'center');
        warningon = 1;
    end
    waitbar((k*6)/(len*6),h_wait);
    %subplot('position', [0.77 0.16 0.15 0.28]),
    maintitle = ['Cluster '''  STUDY.cluster(cls(k)).name ''' average properties (' num2str(length(STUDY.cluster(cls(k)).comps)) ' comps).' ];
    a = textsc(maintitle, 'title'); 
    set(a, 'fontweight', 'bold');     

    if warningon
        disp('Some properties could not be plotted. To plot these properties, first');
        disp('include them in pre-clustering. There, you may spcify 0 dimensions if you');
        disp('want the property (scalp_map, ERSP, etc...) to be computed but not included');
        disp('in the clustering procedure - see the clustering tutorial).');
    end;

end  % Finished all conditions
delete(h_wait)
    

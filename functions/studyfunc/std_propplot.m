% cls_plotclust() - Commandline function, to visualizing cluster/s measures. 
%                   Displays either mean cluster/s ERSP/s, or all cluster/s component 
%                   ERSPs with the mean cluster/s ERSP in one figure (per condition).
%                   The ERSPs can be visualized only if component ERSPs     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = cls_plotclust(STUDY, ALLEEG, clusters);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Optional inputs:
%   clusters   - [numeric vector]  -> specific cluster numbers to plot.
%                     'all'                         -> plot all clusters in STUDY.
%                     {default: 'all'}.
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster 
%                     mean measures, to allow quick replotting (unless cluster means 
%                     already exists in the STUDY).  
%
%   Example:
%                         >> [STUDY] = cls_plotclust(STUDY,ALLEEG, 5);
%                    Plots the mean measures of cluster 5 on the same figure. 
%
%  See also  pop_clustedit, cls_plotcomp        
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

function STUDY = cls_plotclust(STUDY, ALLEEG,  varargin)
icadefs;

% Set default values
cls = 1:length(STUDY.cluster); % plot all clusters in STUDY

if ~isempty(varargin)
    if isnumeric(varargin{1})
        cls = varargin{1};
    elseif isstr(varargin{1}) & strcmpi(varargin{1}, 'all')
        cls = 1:length(STUDY.cluster);
    else
        error('cls_plotclust: clusters input is either specific clusters (numeric vector) or keyword ''all''.');
    end
end

len = length(cls);
% Plot clusters mean measures
for k = 1: len
    if k == 1
        try 
            % optional 'CreateCancelBtn', 'delete(gcbf); error(''USER ABORT'');', 
            h_wait = waitbar(0,['Computing cluster measures ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Computing cluster measures ...'],'position', [300, 200, 300, 48]);
        end
    end    
    figure
    orient tall
    set(gcf,'Color', BACKCOLOR);
    subplot(2,3,1),
    try,
        STUDY = cls_plotclustmap(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off');
    catch
    end
    waitbar(k/(len*6),h_wait)
    subplot(2,3,2),
     try,
         [STUDY] = cls_plotclustersp(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off' );
     catch, end
    waitbar((k*2)/(len*6),h_wait)
    subplot(2,3,3),
    try,
        STUDY = cls_plotclusterp(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off');
    catch
    end
    waitbar((k*3)/(len*6),h_wait)
    subplot(2,3,5),
    try,
        [STUDY] = cls_plotclustitc(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off' );
    catch, end
    waitbar((k*4)/(len*6),h_wait);
    subplot(2,3,4),
    try,
        STUDY = cls_plotclustspec(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'centroid', 'figure', 'off');
    catch
    end
    waitbar((k*5)/(len*6),h_wait)
    subplot('position', [0.77 0.16 0.15 0.28]),
    %subplot(2,3,6),
    try,
        STUDY = cls_plotclustdip(STUDY,ALLEEG, 'clusters', cls(k), 'mode', 'apart', 'figure', 'off');
        set(gcf,'Color', BACKCOLOR);
        %title([ STUDY.cluster(cls(k)).name ' dipoles']);
    catch
    end
    waitbar((k*6)/(len*6),h_wait)
    maintitle = ['Cluster '''  STUDY.cluster(cls(k)).name ''' average measures (' num2str(length(STUDY.cluster(cls(k)).comps)) ' comps).' ];
    a = textsc(maintitle, 'title'); 
    set(a, 'fontweight', 'bold');     
end  % Finished all conditions
delete(h_wait)
    

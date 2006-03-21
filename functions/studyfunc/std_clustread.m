% std_clustread() - load data for one or more requested component measures 
%                      ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map']
%                   for all components of a specified STUDY cluster.  
%                   Can be used in scripts to load cluster data.
%                   Called by cluster plotting functions: std_envtopo(), 
%                   std_erpplot(), std_erspplot(), etc.
% Usage:
%         >> clustinfo = std_clustread(STUDY, ALLEEG, ...
%                                         cluster, infotype, condition);
% Inputs:
%         STUDY - studyset structure containing some or all files in ALLEEG
%        ALLEEG - vector of loaded EEG datasets
%       cluster - cluster number in STUDY
%      infotype - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map'] type of stored
%                 cluster information to read. May also be a cell array of
%                 these types, for example: { 'erp' 'map' 'dipole' }
% Optional input:
%     condition - ['string'] STUDY condition name to read {default: [] if only
%                 one condition. If more than one condition, the name of the 
%                 first condition in the STUDY} UNIMPLEMENTED !!
%
% Output:
%      clustinfo - structure of specified cluster information with fields:
%         .clustname     % STUDY cluster name
%         .clustnum      % STUDY cluster index
%         .condition     % the STUDY condition asked for
%
%         .comp          % [integer array] of component indices 
%                        %  in their respective STUDY datasets
%         .subject       % {cell array} of component subject codes UNIMPLMENTED!
%         .group         % {cell array} of component group codes   UNIMPLMENTED!
%
%         .erp           % [(ncomps, ntimes) array] of component ERPs
%           .erp_times   % [(1,ntimes) array] of ERP epoch latencies
%
%         .spec          % [(ncomps, nfreqs) array] of component spectra
%           .spec_freqs  % [(1,nfreqs) array] of spectral frequencies 
%
%         .ersp          % [(ncomps,ntimes,nfreqs) array] of comp. ERSPs
%           .ersp_times  % [(1,ntimes)] of ERSP latencies
%           .ersp_freqs  % [(1,ntimes)] of ERSP frequencies
%
%         .itc           % [(ncomps,ntimes,nfreqs) array] of comp. ITCs
%           .itc_times   % [(1,ntimes) array] of ITC latencies
%           .itc_freqs   % [(1,nfreqs) array] of ITC frequencies
%
%         .scalp         % [(ncomps, ngrid, ngrid)] comp. scalp map grids
%           .xi          % [(1, ngrid) array] of abscissa values of grid rows 
%           .yi          % [(1, ngrid) array] of ordinate values of grid cols 
%                        % {default ngrid: 65}
%
%         .dipole        % [struct array] of component dipole structures
%                        % Same format as EEG.dipfit.model
% Example:
%         % To plot the ERPs for all Cluster-3 components in a loaded STUDY
%         %
%         clsinfo = std_clustread(STUDY, ALLEEG, 3, 'erp');
%         figure; plot(clsinfo.erp_times, clsinfo.erp);
% 
% Author: Hilit Serby, Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, 2005-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, July 22, 2005, smakeig@ucsd.edu
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

function clustinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, condition);

if nargin < 4
    help std_clustread;
    return;
end
if nargin < 5
  condition = []; % default
end

if ~iscell(infotype), 
    infotype = { infotype }; 
end;

clustinfo = [];
clustinfo.clustname  = STUDY.cluster(cluster).name;
clustinfo.clustnum   = cluster;
clustinfo.comp       = STUDY.cluster(cluster).comps;
clustinfo.condition  = condition;

ncomps = length(STUDY.cluster(cluster).comps);

for k = 1:ncomps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for each cluster component %%%%%%%%%%%%%%
    
    abset = [STUDY.datasetinfo(STUDY.setind(conition,STUDY.cluster(cluster).sets(1,k))).index];
    comp  = STUDY.cluster(cluster).comps(k);
    % clustinfo.subject{k} = ??? UNIMPLEMENTED 
    % clustinfo.group{k} = ??? UNIMPLEMENTED BECAUSE OF CLUSTER.SETS PROBLEM
    
    for index = 1:length(infotype) %%%%%%%%%%%%%%%% for each information type %%%%%%%%%%%%%%%%
        switch infotype{index}        
            case 'erp'
                [erp, t] = std_readerp(ALLEEG, abset, comp, STUDY.preclust.erpclusttimes);
                if  k == 1
                    clustinfo.erp       = zeros(ncomps,length(erp));
                    clustinfo.erp_times = t;
                end
                clustinfo.erp(k,:) = erp;

            case 'spec'
                [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
                if  k == 1
                    clustinfo.spec       = zeros(ncomps,length(spec));
                    clustinfo.spec_freqs = f;
                end
                clustinfo.spec(k,:) = spec;

            case 'ersp'
                [ersp, logfreqs, timevals] = std_readersp(ALLEEG, abset, comp, ...
                          STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
                if  k == 1
                 clustinfo.ersp = zeros(ncomps,length(STUDY.preclust.erspclusttimes), ...
                                       length(STUDY.preclust.erspclustfreqs));
                 clustinfo.ersp_freqs = STUDY.preclust.erspclustfreqs;
                 clustinfo.ersp_times = STUDY.preclust.erspclusttimes;
                end
                clustinfo.ersp(k,:,:)  = ersp;

            case 'itc'
                [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp, ...
			STUDY.preclust.erspclusttimes, STUDY.preclust.erspclustfreqs );
                if  k == 1
                 clustinfo.itc = zeros(ncomps,length(STUDY.preclust.itcclusttimes), ...
                                       length(STUDY.preclust.itcclustfreqs));
                 clustinfo.itc_freqs = STUDY.preclust.itcclustfreqs;
                 clustinfo.itc_times = STUDY.preclust.itcclusttimes;
                end
                clustinfo.itc[k,:,:]   = itc;

            case 'dipole'
                clustinfo.dipole(k) = ALLEEG(abset).dipfit.model(comp);            

            case { 'map' 'scalp' 'topo' }
                [grid, yi, xi] = std_readtopo(ALLEEG, abset, comp); 
                if  k == 1
                 clustinfo.scalp = zeros(ncomps,length(xi),length(yi));
                 clustinfo.xi = xi;
                 clustinfo.yi = yi;
                end
                clustinfo.scalp[k,:,:] = grid;

            otherwise, error('Unrecognized ''infotype'' entry');
        end; % switch
    end; % infotype
end % comp

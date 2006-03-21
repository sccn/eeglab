% std_clustread() - load one or more requested component data measures 
%                     ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map']
%                   for all components of a specified cluster.  
%                   Useful for accessing cluster compoent data in scripts.
%                   Called by cluster plotting functions: std_envtopo(), 
%                   std_erpplot(), std_erspplot(), etc.
% Usage:
%         >> clustinfo = std_clustread(STUDY,ALLEEG, cluster, ...
%                                             infotype(s), condition(s));
% Inputs:
%         STUDY - studyset structure containing some or all files in ALLEEG
%        ALLEEG - vector of loaded EEG datasets including STUDY datasets
%       cluster - cluster number in STUDY
%      infotype(s) - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map'] type(s) of 
%                    cluster information to read. May also be a cell array of
%                    these types. For example: { 'erp' 'map' 'dipole' }
% Optional input:
%     condition(s) - ['string'] or {cell array of 'strings'} STUDY condition 
%                    name(s) to return info for {default|[]: all conditions}
% Output:
%      clustinfo   - structure of specified cluster information 
%
% With fields:
%         .clustname     % cluster name
%         .clustnum      % cluster index 
%         .condition     % {dell array of 'strings'} condition(s) returned
%         .subject       % {cell array} of component subject codes 
%         .group         % {cell array} of component group codes  UNIMPLMENTED!
%         .comp          % [integer array] of comp. dataset indices 
%
% And optional fields (defined for the requested measures):
%         .erp           % [(ncomps, ntimes) array] component ERPs
%           .erp_times   % [(1,ntimes) array] ERP epoch latencies (ms)
%
%         .spec          % [(ncomps, nfreqs) array] component spectra
%           .spec_freqs  % [(1,nfreqs) array] spectral frequencies (Hz)
%
%         .ersp          % [(ncomps, ntimes, nfreqs) array] component ERSPs
%           .ersp_times  % [(1,ntimes) array] ERSP latencies (ms)
%           .ersp_freqs  % [(1,nfreqs) array] ERSP frequencies (Hz)
%
%         .itc           % [(ncomps, ntimes, nfreqs) array] component ITCs
%           .itc_times   % [(1,ntimes) array] ERSP latencies (ms)
%           .itc_freqs   % [(1,nfreqs) array] ERSP frequencies (Hz)
%
%         .scalp         % [(ncomps, ngrid, ngrid) array] comp. map grids
%           .xi          % [(1, ngrid) array] abscissa values 
%           .yi          % [(1, ngrid) array] oridnate values 
%                        % {default ngrid: 65}
%
%         .dipole        % [(1,ncomps) struct array] comp. dipole information 
%                        % Same format as EEG.dipfit.model ( see >> help dipfit)
% Example:
%         % To overplot the ERPs for all Cluster 3 components in a STUDY
%         %
%         info = std_clustread(STUDY, ALLEEG, 3, 'erp');
%         figure; plot(info.erp_times, info.erp);
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

clustinfo = [];      % initialize

if isempty(condition)
  clustinfo.condition   = STUDY.condition; % all conditions
else % condition(s) specified
  if ~iscell(condition), 
    condition = { condition }; 
  end;
  if ~ismember(condition,STUDY.condition) % ???? WRONG ???
    error('Named condition not found in STUDY')
  else
    clustinfo.condition  = condition;
  end
end

if ~iscell(infotype), 
    infotype = { infotype }; 
end;

clustinfo.clustname  = STUDY.cluster(cluster).name;
clustinfo.clustnum   = cluster;
clustinfo.subject    = [];  % UNIMPLEMENTED! ???
clustinfo.group      = [];  % UNIMPLEMENTED! ???
clustinfo.comp       = STUDY.cluster(cluster).comps;

ncomps = length(STUDY.cluster(cluster).comps);
nconditions = length(condition);
condnums = zeros(1,nconditions);

% for c=1:nconditions % find condition indices
%    condnums(c) = ismember(condition(c),STUDY.condition); % ??? WRONG ??? should return the condition numbers
% end

for k = 1:ncomps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for each cluster component %%%%%%%%%%%%%%
    
 for c = 1:nconditions
   
    % NB: BELOW, CONDITION NEEDS TO BE THE INDICES OF THE REQUESTED CONDITIONS, NO???

    abset = [STUDY.datasetinfo(STUDY.cluster(cluster).sets(condnums(c),k)).index]; % ??? condition ???
    comp  = STUDY.cluster(cluster).comps(k);
    % clustinfo(c).clustsubj{k} = ??? UNIMPLEMENTED 
    % clustinfo(c).clustgrp{k}  = ??? UNIMPLEMENTED - NOTE CLUSTER.SETS PROBLEM
    
    for index = 1:length(infotype) %%%%%%%%%%%%%%%% for each information type %%%%%%%%%%%%%%%%
        switch infotype{index}        
            case 'erp'
                [erp, t] = std_readerp(ALLEEG, abset, comp, STUDY.preclust.erpclusttimes);
                if  k == 1
                    clustinfo.erp       = zeros(ncomps,length(erp));
                    clustinfo.erp_times = t;
                end
                clustinfo.erp(k,:) = erp;
                if k==ncomps & c==conditions, fprintf('Cluster component ERPs read.\n'); end;
            case 'spec'
                [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
                if  k == 1
                    clustinfo.spec       = zeros(ncomps,length(spec));
                    clustinfo.spec_freqs = f;
                end
                clustinfo.spec(k,:) = spec;
                if k==ncomps & c==conditions, fprintf('Cluster component spectra read.\n'); end;

            case 'ersp'
                [ersp, logfreqs, timevals] = std_readersp(ALLEEG, abset, comp, ...
                          STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
                clustinfo.ersp_freqs{k} = logfreqs;
                clustinfo.ersp_times{k} = timevals;
                clustinfo.ersp{k}       = ersp;
                if k==ncomps & c==conditions, fprintf('Cluster component ERSPs read.\n'); end;

            case 'itc'
                [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp, ...
                         STUDY.preclust.erspclusttimes, STUDY.preclust.erspclustfreqs );
                clustinfo.itc_freqs{k} = logfreqs;
                clustinfo.itc_times{k} = timevals;
                clustinfo.itc{k}       = itc;
                if k==ncomps & c==conditions, fprintf('Cluster component ITCs read.\n'); end;

            case 'dipole'
                clustinfo.dipole(k) = ALLEEG(abset).dipfit.model(comp);            
                if k==ncomps & c==conditions, fprintf('Cluster component dipole models read.\n'); end;

            case { 'map' 'scalp' 'topo' }
                [grid, yi, xi] = std_readtopo(ALLEEG, abset, comp); 
                if  k == 1
                 clustinfo.scalp = zeros(ncomps,length(xi),length(yi));
                 clustinfo.xi = xi;
                 clustinfo.yi = yi;
                end
                clustinfo.scalp(k,:,:) = grid;
                if k==ncomps & c==conditions, fprintf('Cluster component scalp maps grids read.\n'); end;

            otherwise, error('Unrecognized ''infotype'' entry');
        end; % switch
    end; % infotype
end % comp

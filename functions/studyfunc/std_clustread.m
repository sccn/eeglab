% std_clustread() - load one or more requested measures 
%                   ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map']
%                   for all components of a specified cluster.  
%                   Called by cluster plotting functions: std_envtopo(), 
%                   std_erpplot(), std_erspplot(), ...
% Usage:
%         >> clustinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, condition);
% Inputs:
%         STUDY - studyset structure containing some or all files in ALLEEG
%        ALLEEG - vector of loaded EEG datasets
%       cluster - cluster number in STUDY
%      infotype - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map'] type of stored
%                 cluster information to read. May also be a cell array of
%                 these types, for example: { 'erp' 'map' 'dipole' }
%     condition - STUDY condition number to read {default: 1}
%
% Output:
%      clustinfo - structure of specified cluster information:
%
%         clustinfo.name          % cluster name
%         clustinfo.clusternum    % cluster index
%         clustinfo.condition     % index of the condition asked for
%
%         clustinfo.comp[]        % array of component indices 
%         clustinfo.subject{}     % cell array of component subject codes UNIMPLMENETED
%         clustinfo.group{}       % cell array of component group codes  UNIMPLMENETED!
%
%         clustinfo.erp[]         % (ncomps, ntimes) array of component ERPs
%           clustinfo.erp_times[] % vector of ERP epoch latencies
%
%         clustinfo.spec[]        % (ncomps, nfreqs) array of component spectra
%           clustinfo.spec_freqs[]% vector of spectral frequencies 
%
%         clustinfo.ersp[]        % (ncomps,ntimes,nfreqs) array of component ERSPs
%           clustinfo.ersp_times[]% vector of ERSP latencies
%           clustinfo.ersp_freqs[]% vector of ERSP frequencies
%
%         clustinfo.itc[]         % (ncomps,ntimes,nfreqs) array of component ITCs
%           clustinfo.itc_times[] % vector of ITC latencies
%           clustinfo.itc_freqs[] % vector of ITC frequencies
%
%         clustinfo.scalp[]       % (ncomps,65,65) array of component scalp map grids
%           clustinfo.xi[]        % abscissa values for columns of the scalp maps
%           clustinfo.yi[]        % ordinate values for rows of the scalp maps
%
%         clustinfo.dipole        % array of component dipole information structs
%                                 % with same format as EEG.dipfit.model
% Example:
%         % To plot the ERPs for all Cluster-3 components from a STUDY
%         %
%         clustinfo = std_clustread(STUDY, ALLEEG, 3, 'erp');
%         figure; plot(clustinfo.erp_times, clustinfo.erp);
% 
% Author: Hilit Serby, Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, 2005-

function clustinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, condition);

if nargin < 4
    help std_clustread;
    return;
end
if nargin < 5
  condition = 1; % default
end

if ~iscell(infotype), 
    infotype = { infotype }; 
end;

clustinfo = [];
clustinfo.name       = STUDY.cluster(cluster).name;
clustinfo.clusternum = cluster;
clustinfo.comps      = STUDY.cluster(cluster).comps;
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

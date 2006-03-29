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
  condition = [1:length(STUDY.condition)]; % default
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
for k = 1:ncomps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for each channel component %%%%%%%%%%%%%%
    
    for n = 1:length(condition)

        abset = [STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).index];
        comp  = STUDY.cluster(cluster).comps(k);
        clustinfo.subject{k} = STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).subject;
        clustinfo.group{k} = STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).group;

        for index = 1:length(infotype) %%%%%%%%%%%%%%%% for each information type %%%%%%%%%%%%%%%%
            switch infotype{index}        
                case 'erp'
                    [erp, t] = std_readerp(ALLEEG, abset, comp, STUDY.preclust.erpclusttimes);
                    clustinfo.erp_times = t;
                    clustinfo.erp{n,k} = erp;

                case 'spec'
                    [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
                    clustinfo.spec_freqs = f;
                    clustinfo.spec{n,k} = spec;

                case 'ersp'
                    if n == 1
                        abset = [STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition,k)).index];
                        [ersp, logfreqs, timevals] = std_readersp(ALLEEG, abset, comp, ...
                                  STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
                        for cc = 1:length(condition)
                            clustinfo.ersp{condition(cc), k} = ersp(:,:,cc);
                        end;
                        clustinfo.ersp_freqs = logfreqs;
                        clustinfo.ersp_times = timevals;
                    end;
                case 'itc'
                    [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp, ...
                             STUDY.preclust.erspclusttimes, STUDY.preclust.erspclustfreqs );
                    clustinfo.itc_freqs = logfreqs;
                    clustinfo.itc_times = timevals;
                    clustinfo.itc{n,k}       = itc;

                case 'dipole'
                    if n == 1, clustinfo.dipole(k) = ALLEEG(abset).dipfit.model(comp); end;

                case { 'map' 'scalp' 'topo' }
                    if n == 1
                        [grid, yi, xi] = std_readtopo(ALLEEG, abset, comp); 
                        if  k == 1
                         clustinfo.xi = xi;
                         clustinfo.yi = yi;
                        end
                        clustinfo.scalp{k} = grid;
                    end;
                otherwise, error('Unrecognized ''infotype'' entry');
            end; % switch
        end;
    end; % infotype
end % comp

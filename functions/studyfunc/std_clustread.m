% std_clustread() - this function has been replaced by std_readdata() for
%                   consistency. Please use std_readdata() instead.
%                   load one or more requested measures 
%                   ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map']
%                   for all components of a specified cluster.  
%                   Calls std_readerp(), std_readersp(), etc.
% Usage:
%         >> clustinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, condition);
% Inputs:
%         STUDY - studyset structure containing some or all files in ALLEEG
%        ALLEEG - vector of loaded EEG datasets
%       cluster - cluster number in STUDY
%      infotype - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map'] type of stored
%                 cluster information to read. May also be a cell array of
%                 these types, for example: { 'erp' 'map' 'dipole' }
%     condition - STUDY condition number to read {default: all}
%
% Output:
%      clustinfo - structure of specified cluster information:
%         clustinfo.name          % cluster name
%         clustinfo.clusternum    % cluster index
%         clustinfo.condition     % index of the condition asked for
%
%         clustinfo.comp[]        % array of component indices 
%         clustinfo.subject{}     % cell array of component subject codes
%         clustinfo.group{}       % cell array of component group codes
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
%         % To plot the ERPs for all components in cluster 3 of a loaded STUDY
%         >> clustinfo = std_clustread(STUDY, ALLEEG, 3, 'erp');
%            figure; plot(clustinfo.erp_times, clustinfo.erp);
% 
% See also: std_readerp(), std_readspec(), std_readersp(), std_readitc(), std_readtopo()
%
% Authors: Hilit Serby, Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, 2005-
%
% RCS-recorded version number, date, editor and comments
function clustinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, condition);

    help std_clustread;
    return;
    
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
for k = 1:ncomps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for each channel | component %%%%%%%%%%%%%%
    
    for n = 1:length(condition) %%%%%%%%%%%%%%%%%%%%%%%% for each STUDY condition %%%%%%%%%%%%%%

        abset = [STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).index];
        comp  = STUDY.cluster(cluster).comps(k);
        clustinfo.subject{k} = STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).subject;
        clustinfo.group{k} = STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).group;

        for index = 1:length(infotype) %%%%%%%%%%%%%% for each information type %%%%%%%%%%%%%%%%
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

% FUTURE HEADER
% std_clustread() - return detailed information and (any) requested component 
%                   measures for all components of a specified cluster. Restrict 
%                   component info to components from specified subjects, groups,
%                   sessions, and/or conditions. Use in scripts handling results 
%                   of component clustering. Called by cluster plotting 
%                   functions: std_envtopo(), std_erpplot(), std_erspplot(), ...
% Usage:
%         >> clsinfo = std_clustread(STUDY, ALLEEG, cluster);  % use defaults
%         >> clsinfo = std_clustread(STUDY, ALLEEG, cluster, ...
%                                             'keyword1', keyval1, ...
%                                                   'keyword2', keyval2, ...);
% Inputs:
%         STUDY     - studyset structure containing some or all files in ALLEEG
%         ALLEEG    - vector of loaded EEG datasets including those in STUDY
%         cluster   - cluster number in STUDY to return information for 
%
% Optional keywords - and values:                                               IMPLEMENT !!
%      'measure'    - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'topo'] stored
%                     cluster measure(s) to read. May also be a cell array of
%                     these, for example: { 'erp' 'map' 'dipole' }. 
%                     Else 'all', meaning all measures available.
%                     Else 'cls', meaning all measures clustered on.
%                     Else [], for none {default: 'cls'}                                            
%     'condition'   - STUDY condition 'name' or {'names'} to read,
%                     Else 'all' {default: 'all'}
%     'condnum''    - STUDY condition [number(s)] to read,
%                     Else 0 -> all {default: 0}
%     'subject'     - STUDY subjects 'name' or {'names'} to read, 
%                     Else 'all' {default: 'all'}
%     'subjnum'     - STUDY subjects [number(s)] to read, 
%                     Else 0 -> all {default: 0}
%     'group'       - STUDY subject group 'name' or {'names'} to read, 
%                     Else 'all' {default: 'all'}
%     'groupnum'    - STUDY subject group [number(s)] to read, 
%                     Else 0 -> all {default: 0}
%     'session'     - STUDY session [number(s)] to read, 
%                     Else 0 -> all {default: 0}
%
% Output:
%      clsinfo - structure containing information about cluster components in fields:
%
%         .clustername   % cluster name
%           .clusternum  % cluster index
%
%         .dataset       % {(conditions,components) cell array} component dataset indices 
%                                                           in the input ALLEEG array
%         .component     % [(1,components) int array]  component decomposition indices 
%         .subject       % {(1,components) cell array} component subject codes 
%           .subjectnum  % [(1,components) int array] component subject indices
%         .group         % {(1,components) cell array} component group codes  
%           .groupnum    % [(1,components) int array] component group indices
%
%         .condition     % {(1,components) cell array} component condition codes
%           .conditionnum % [(1,components) int array] component condition indices
%         .session       % [(1,components) int array] component session indices
%
%         .erp           % {(conditions, components) cell array} 
%                                       (1, latencies) component ERPs             CHECK DIM
%           .erp_times   % [num vector] ERP epoch latencies (s)
%
%         .spec          % {(conditions, components) cell array} 
%                                       (1,frequencies) component spectra         CHECK DIM
%           .spec_freqs  % [num vector] spectral frequencies (Hz)
%
%         .ersp          % {(conditions, components) cell array} 
%                                       (freqs,latencies) component ERSPs         CHECK DIM
%           .ersp_times  % [num vector] ERSP epoch latencies (s)
%           .ersp_freqs  % [num vector] ERSP frequencies (Hz)
%
%         .itc           % {(conditions, components) cell array} 
%                                       (freqs,latencies) component ITCs           CHECK DIM
%           .itc_times   % [num vector] ITC epoch latencies (s)
%           .itc_freqs   % [num vector] ITC frequencies (Hz)
%
%         .topo          % {(1,components) cell array} 
%                                              (65,65) component topo map grids    CHECK DIM
%           .xi          % [vector] topo grid abscissa values 
%           .yi          % [vector] topo grid ordinate values
%
%         .dipole        % [struct array] component dipole information 
%                        % structures with same format as "EEG.dipfit.model"
%                        % See >> help dipfit                                      CHECK HELP
% Example:
%         % To plot the ERPs for all Cluster-3 components in a one-condition STUDY
%         %
%         clsinfo3 = std_clustread(STUDY, ALLEEG, 3, 'measure', 'erp'); % assume 1 condition
%         times = clsinfo3.erp_times; figure; plot(times, clsinfo3.erp');
%         %
%         % To plot ERPs for only those Cluster-3 components from subjects in group 'female'
%         %
%         feminfo3 = std_clustread(STUDY, ALLEEG, 3, 'measure', 'erp', 'group', 'female');
%         figure; plot(times, feminfo3.erp');
%         %
%         % Alternatively, to extract 'female' subject components from clsinfo3 above
%         %
%         femidx = find(strcmp({clsinfo3.group},'female'));
%         figure; plot(times, clustinfo.erp(femidx,:)');                            CHECK EXAMPLE
% 
% Authors: Hilit Serby, Scott Makeig, Toby Fernsler & Arnaud Delorme, SCCN/INC/UCSD, 2005-

% std_clustread() - return one or more requested measures 
%                    ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map'|'all']
%                   for all components of a specified cluster in specified 
%                   conditions. Useful for scripts handling results of 
%                   component clustering. Called by cluster plotting 
%                   functions std_envtopo(), std_erpplot(), std_erspplot(), ...
% Usage:
%         >> clustinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, condition);
% Inputs:
%         STUDY - studyset structure containing some or all files in ALLEEG
%        ALLEEG - vector of loaded EEG datasets
%       cluster - cluster number in STUDY to return information for 
%
% Optional inputs:
%      infotype - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'topo'] type of stored
%                 cluster information to read. May also be a cell array of
%                 these types, for example: { 'erp' 'map' 'dipole' } or 'all'
%                 for all available {default: 'all'}                             IMPLEMENT 'all'
%     condition - STUDY condition number(s) to read {default: 1}
%
% Output:
%      clustinfo - structure containing information about the cluster in fields:
%
%         .clustername   % cluster name
%         .clusternum    % cluster index
%         .conditionname % name of the specified condition(s)                    IMPLEMENT .conditionname
%         .conditionnum  % index of the specified condition(s) 
%
%         .comp          % [(1,components) int array]  component indices 
%         .set           % {(conditions,components) cell array} component dataset 
%                                      indices in ALLEEG 
%         .subject       % {(1,components) cell array} component subject codes 
%         .group         % {(1,components) cell array} component group codes  
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
%         % To plot the (condition 1) ERPs for all (cluster 3) components 
%         % from a study on the same axis
%         %
%         clustinfo = std_clustread(STUDY, ALLEEG, 3, 'erp');
%         times = clustinfo.erp_times;
%         figure; plot(times, clustinfo.erp');
%         %
%         % To restrict the plot to subjects from group 'female'
%         %
%         femidx = find(strcmp({clustinfo.group},'female'));
%         figure; plot(times, clustinfo.erp(femidx,:)');                            CHECK EXAMPLE
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
clustinfo.clustername   = STUDY.cluster(cluster).name;
clustinfo.clusternum    = cluster;
clustinfo.comps         = STUDY.cluster(cluster).comps;
clustinfo.conditionname = [];                                                     % IMPLEMENT
clustinfo.conditionnum  = condition;

ncomps = length(STUDY.cluster(cluster).comps);

for k = 1:ncomps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for each cluster component %%%%%%%%%%%%%%
    for n = 1:length(condition) %%%%%%%%%%%%%%%%%%%%%% for each cluster condition %%%%%%%%%%%%%

        abset = [STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).index];
        if n == 1
           comp  = STUDY.cluster(cluster).comps(k);
           clustinfo.subject{k} = STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).subject;
           clustinfo.group{k} = STUDY.datasetinfo(STUDY.cluster(cluster).sets(condition(n),k)).group;
        end
        clustinfo.set{n,k} = abset;

        for index = 1:length(infotype) %%%%%%%%%%%%%%%% for each information type %%%%%%%%%%%%%%%%
            switch infotype{index}        
                case 'erp'
                    [erp, t] = std_readerp(ALLEEG, abset, comp, STUDY.preclust.erpclusttimes);
                    if n == 1 & k == 1, 
                        clustinfo.erp_times = t; 
                    end;
                    clustinfo.erp{n,k} = erp;

                case 'spec'
                    [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
                    if (n == 1 & k == 1) 
                        clustinfo.spec_freqs = f; 
                    end
                    clustinfo.spec{n,k} = spec;

                case 'ersp'
                    [ersp, logfreqs, timevals] = std_readersp(ALLEEG, abset, comp, ...
                                  STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
                    if k==1 & n==1, 
                        clustinfo.ersp_freqs = logfreqs;
                        clustinfo.ersp_times = timevals;
                    end
                    clustinfo.ersp{n, k} = ersp;
                case 'itc'
                    [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp, ...
                                  STUDY.preclust.erspclusttimes, STUDY.preclust.erspclustfreqs );
                    if k==1 & n==1 
                        clustinfo.itc_freqs = logfreqs;
                        clustinfo.itc_times = timevals;
                    end
                    clustinfo.itc{n, k} = itc;
                case 'dipole'
                    if n==1, 
                        clustinfo.dipole(k) = ALLEEG(abset).dipfit.model(comp); 
                    end;

                case { 'topo' 'map' 'scalp' }
                    if n==1
                        [grid, yi, xi] = std_readtopo(ALLEEG, abset, comp); 
                        if  k==1
                            clustinfo.xi = xi;
                            clustinfo.yi = yi;
                        end
                        clustinfo.topo{k} = grid;
                    end;
                case 'all'
                   error('Type ''all'' UNIMPLEMENTED.');                       % IMPLEMENT
                otherwise, error('Unrecognized ''infotype'' entry');
            end; % switch
        end;
    end; % infotype
end % comp

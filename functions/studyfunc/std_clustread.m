% std_clustread() - load a requested measure (i.e. ['erp'|'spec'|'ersp'|'itc'|...
%                   'dipole'|'scalp']), for all components of a specified cluster.  
%                   Called by cluster plotting functions: std_envtopo(), 
%                   std_erpplot(), std_erspplot(), ...
% Usage:
%         >> clustinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, cond);
%
% Inputs:
%         STUDY - studyset structure containing some or all files in ALLEEG
%        ALLEEG - vector of loaded EEG datasets
%       cluster - cluster number in STUDY
%      infotype - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'scalp'] type of stored
%                 cluster information to read. May also be a cell array of
%                 these type, for instance { 'erp' 'scalp' 'dipole' }
%          cond - STUDY condition number to read {default: 1}
%
% Output:
%      clustinfo - structure containing relevant information fields:
%                 clustinfo.name         % cluster name
%                 clustinfo.clusternum   % cluster index
%                 clustinfo.comps        % component indices of the component
%                 clustinfo.sets         % dataset indices (one row per condition)
%                 clustinfo.condition    % condition read
%                 clustinfo.erp          % vector of erp value for cluster
%                 ->clustinfo.erp_times  % vector of erp latencies
%                 clustinfo.spec         % vector of spectral values for cluster
%                 ->clustinfo.spec_freqs % vector of spec frequencies
%                 clustinfo.ersp         % 2-D ERSP array for cluster
%                 ->clustinfo.ersp_times % vector of ersp latencies
%                 ->clustinfo.ersp_freqs % vector of ersp frequencies
%                 clustinfo.itc          % 2-D ITC array for cluster
%                 ->clustinfo.itc_times  % vector of itc latencies
%                 ->clustinfo.itc_freqs  % vector of itc frequencies
%                 clustinfo.scalp        % scalp topomap array for cluster centroid
%                 ->clustinfo.xi         % abscicia of columns in array above
%                 ->clustinfo.yi         % ordinate of rows in array above
%                 clustinfo.dipole       % when dipole information is returned, 
%                                         clustinfo.dipole is has the same format 
%                                         as EEG.dipfit.model
%                                         for the cluster components.
% Example:
%   % assuming some cluster have been computed
%   % to plot the ERP of all the member of the cluster
%   clustinfo = std_clustread(STUDY, ALLEEG, 3, 'erp');
%   figure; plot(clustinfo.erp_times, clustinfo.erp);
% 
% Author: Hilit Serby & Arnaud Delorme, SCCN/INC/UCSD, 2005

function clustinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, cond);

if nargin < 5
  cond = 1; % default
end

if nargin < 4
    help std_clustread;
    return;
end
if ~iscell(infotype), infotype = { infotype }; end;
clustinfo = [];
clustinfo.name       = STUDY.cluster(cluster).name;
clustinfo.clusternum = cluster;
clustinfo.comps      = STUDY.cluster(cluster).comps;
clustinfo.condition  = cond;

len = length(STUDY.cluster(cluster).comps);
for k = 1:len
    if nargin < 5
        abset = STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cluster).sets(1,k))).index;
    else
        abset = [STUDY.datasetinfo(STUDY.setind(cond,STUDY.cluster(cluster).sets(1,k))).index];
    end
    comp = STUDY.cluster(cluster).comps(k);
    
    for index = 1:length(infotype)
        switch infotype{index}        
            case 'erp'
                [erp, t] = std_readerp(ALLEEG, abset, comp, STUDY.preclust.erpclusttimes);
                if  k == 1
                    clustinfo.erp       = zeros(len,length(erp));
                    clustinfo.erp_times = t;
                end
                clustinfo.erp(k,:) = erp;

            case 'spec'
                [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
                if  k == 1
                    clustinfo.spec       = zeros(len,length(spec));
                    clustinfo.spec_freqs = f;
                end
                clustinfo.spec(k,:) = spec;

            case 'ersp'
                [ersp, logfreqs, timevals] = std_readersp(ALLEEG, abset, comp, ...
                          STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
                clustinfo.ersp{k}  = ersp;
                clustinfo.ersp_freqs{k}  = logfreqs;
                clustinfo.ersp_times{k} = timevals;

            case 'itc'
                [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp, STUDY.preclust.erspclusttimes, ...
                                                                        STUDY.preclust.erspclustfreqs );
                clustinfo.itc{k}   = itc;
                clustinfo.itc_freqs{k}  = logfreqs;
                clustinfo.itc_times{k} = timevals;

            case 'dipole'
                clustinfo.dipole(k) = ALLEEG(abset).dipfit.model(comp);            

            case { 'scalp' 'topo' }
                [grid, yi, xi] = std_readtopo(ALLEEG, abset, comp); 
                clustinfo.scalp{k} = grid;
                clustinfo.yi{k} = yi;
                clustinfo.xi{k} = xi;
            otherwise, error('Unrecognized ''infotype'' entry');
        end
    end;     
end

        

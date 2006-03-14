% std_clustread() - load a requested measure (i.e. ['erp'|'spec'|'ersp'|'itc'|...
%                   'dipole'|'scalp']), for all components of a specified cluster.  
%                   Called by cluster plotting functions: std_envtopo(), std_erpplot(),
%                   std_erspplot(), ...

function clusinfo = std_clustread(STUDY,ALLEEG, cluster, infotype, cond);

% infotype - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'scalp']

if nargin < 4
    help std_clustread;
    return;
end
clusinfo = [];

len = length(STUDY.cluster(cluster).comps);
for k = 1:len
    if nargin < 5
        abset = STUDY.datasetinfo(STUDY.setind(1,STUDY.cluster(cluster).sets(1,k))).index;
    else
        abset = [STUDY.datasetinfo(STUDY.setind(cond,STUDY.cluster(cluster).sets(1,k))).index];
    end
    comp = STUDY.cluster(cluster).comps(k);
    
	switch infotype        
        case 'erp'
            if nargin < 5
                error('Reading cluster ERP requires condition number');
            end
            [erp, t] = std_readerp(ALLEEG, abset, comp, STUDY.preclust.erpclusttimes);
            if  k == 1
                clusinfo.erp = zeros(len,length(erp));
                clusinfo.t = t;
            end
            clusinfo.erp(k,:) = erp;
            
        case 'spec'
            if nargin < 5
                error('Reading cluster spectrum requires condition number');
            end
            [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
            if  k == 1
                clusinfo.spec = zeros(len,length(spec));
                clusinfo.f = f;
            end
            clusinfo.spec(k,:) = spec;
            
        case 'ersp'
            if nargin < 5
                error('Reading cluster ERSP requires condition number');
            end
            [ersp, logfreqs, timevals] = std_readersp(ALLEEG, abset, comp, ...
                      STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
            clusinfo.ersp{k}  = ersp;
            clusinfo.logf{k}  = logfreqs;
            clusinfo.times{k} = timevals;
            
        case 'itc'
            if nargin < 5
                error('Reading cluster ITC requires condition number');
            end
            [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp, STUDY.preclust.erspclusttimes, ...
                                                                    STUDY.preclust.erspclustfreqs );
            clusinfo.itc{k}   = itc;
            clusinfo.logf{k}  = logfreqs;
            clusinfo.times{k} = timevals;
            
        case 'dipole'
            clusinfo(k) = ALLEEG(abset).dipfit.model(comp);            
            
        case { 'scalp' 'topo' }
            [grid, yi, xi] = std_readtopo(ALLEEG, abset, comp); 
            clusinfo.grid{k} = grid;
            clusinfo.yi{k} = yi;
            clusinfo.xi{k} = xi;
    end
            
end

        

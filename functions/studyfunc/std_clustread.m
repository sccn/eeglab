% std_clustread() - load a requested measure (i.e. ['erp'|'spec'|'ersp'|'itc'|...
%        'dipole'|'scalp']), for all the components of a specified cluster.  
%        Used mostly by cluster plotting functions. 
%        This helper function is called by std_envtopo(), std_erpplot(),
%        std_erspplot(), ...

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
            if ~isfield(ALLEEG(abset).etc, 'icaerpparams')
                error([ 'No ERP information available in dataset ' num2str(abset) ]);   
            end
            [erp, t] = std_readerp(ALLEEG, abset, comp);
            if  k == 1
                clusinfo.erp = zeros(len,length(erp));
                clusinfo.t = t;
            end
            clusinfo.erp(k,:) = erp;
            
        case 'spec'
            if nargin < 5
                error('Reading cluster spectrum requires condition number');
            end
            if ~isfield(ALLEEG(abset).etc, 'icaspecparams')
                error([ 'No spectrum information available in dataset ' num2str(abset) ]);   
            end
            [spec, f] = std_readspec(ALLEEG, abset, comp, 'm');
            if  k == 1
                clusinfo.spec = zeros(len,length(spec));
                clusinfo.f = f;
            end
            clusinfo.spec(k,:) = spec;
            
        case 'ersp'
            if nargin < 5
                error('Reading cluster ERSP requires condition number');
            end
            if ~isfield(ALLEEG(abset(1)).etc, 'icaerspparams')
                error([ 'No ERSP information available in dataset ' num2str(abset(1)) ]);   
            end
            [ersp, logfreqs, timevals] = std_readersp(ALLEEG, abset, comp, ...
                      STUDY.preclust.erspclusttimes,  STUDY.preclust.erspclustfreqs);
            clusinfo.ersp{k}  = ersp;
            clusinfo.logf{k}  = log(logfreqs);
            clusinfo.times{k} = timevals;
            
        case 'itc'
            if nargin < 5
                error('Reading cluster ITC requires condition number');
            end
            if ~isfield(ALLEEG(abset).etc, 'icaitcparams')
                error([ 'No ITC information available in dataset ' num2str(abset) ]);   
            end
            [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp);
            clusinfo.itc{k}   = itc;
            clusinfo.logf{k}  = log(logfreqs);
            clusinfo.times{k} = logfreqs;
            
        case 'dipole'
            if ~isfield(ALLEEG(abset).etc, 'dipfit')
                error([ 'No dipole information available in dataset ' num2str(abset) ]);   
            end
            clusinfo(k) = ALLEEG(abset).dipfit.model(comp);            
            
        case { 'scalp' 'topo' }
            if ~isfield(ALLEEG(abset).etc,'icascalpparams')
                error([ 'Dataset ' num2str(abset) ' has no topoplot image information']);   
            end
            [grid, yi, xi] = std_readtopo(ALLEEG, abset, comp); 
            clusinfo.grid{k} = grid;
            clusinfo.yi{k} = yi;
            clusinfo.xi{k} = xi;
    end
            
end

        

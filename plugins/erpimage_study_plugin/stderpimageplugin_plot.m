function stderpimageplugin_plot(com, hdl);

userdat  = get(hdl, 'userdat');
ALLEEG   = userdat{1}{1};
STUDY    = userdat{1}{2};
cls      = userdat{1}{3};
allchans = userdat{1}{4};
plotting_option = 'erpimageplot';

if strcmpi(com, 'params')
    [STUDY com] = pop_erpimparams(STUDY);
    if ~isempty(com)
        STUDY.history =  sprintf('%s\n%s',  STUDY.history, com);
    end;
    userdat{1}{2} = STUDY;
    set(hdl, 'userdat',userdat); %update information (STUDY)
end;

if strcmpi(com, 'changrp') || strcmpi(com, 'onechan')
    changrp = get(findobj('parent', hdl, 'tag', 'chan_list')   , 'value');
    onechan = get(findobj('parent', hdl, 'tag', 'chan_onechan'), 'value');
    
    if strcmpi(com, 'changrp')
        changrpstr = allchans(changrp);
        a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ');' ];
        % update Study history
        eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
        userdat{1}{2} = STUDY;
        set(hdl, 'userdat',userdat);
    else
        changrpstr    = allchans(changrp);
        if onechan(1) ~= 1  % check that not all onechan in channel are requested
            subject = STUDY.design(STUDY.currentdesign).cases.value{onechan-1};
            a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''subject'', ''' subject ''' );' ];
            eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
        else
            a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''plotsubjects'', ''on'' );' ];
            eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
        end;
        userdat{1}{2} = STUDY;
        set(hdl, 'userdat',userdat);
    end;
    
else
    clus     = get(findobj('parent', hdl, 'tag', 'clus_list'), 'value');
    comp_ind = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'Value');
    
    if strcmpi(com, 'onecomp')
        if (clus ~= 1 ) % specific cluster option
            if ~isempty(STUDY.cluster(cls(clus-1)).comps)
                a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(clus-1)) ');' ];
                eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
            end
        else % all clusters
            % All clusters does not include 'Notclust' 'ParentCluster' and 'Outliers' clusters.
            tmpcls = [];
            for k = 1:length(cls)
                if ~strncmpi(STUDY.cluster(cls(k)).name,'Notclust',8) & ~strncmpi(STUDY.cluster(cls(k)).name,'Outliers',8) & ...
                        (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)) & ~isempty(STUDY.cluster(cls(k)).comps)
                    tmpcls = [ tmpcls cls(k)];
                end
            end
            a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'',['  num2str(tmpcls) ']);' ];
            %if strcmpi(plotting_option, 'dipplot'), a = [a(1:end-2) ',''mode'', ''together'');' ]; end;
            eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
        end
        userdat{1}{2} = STUDY;
        set(hdl, 'userdat',userdat);
        
    else
        if (clus ~= 1 ) %specific cluster
            if comp_ind(1) ~= 1  % check that not all comps in cluster are requested
                subject = STUDY.datasetinfo( STUDY.cluster(cls(clus-1)).sets(1,comp_ind-1)).subject;
                a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(clus-1)) ', ''comps'', ' num2str(comp_ind-1) ' );' ];
                eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
            else
                a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(clus-1)) ', ''plotsubjects'', ''on'' );' ];
                eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
            end
        else
            comp_list = get(findobj('parent', hdl, 'tag', 'clust_comp'), 'String');
            comp_name = comp_list(comp_ind);
            for ci = 1:length(comp_name)
                num_comps = 0;
                tmp = strfind(comp_name{ci},'''');
                clust_name = comp_name{ci}(tmp(1)+1:tmp(end)-1);
                for k = 1:length(cls)
                    if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) & ~strncmpi('Outliers',STUDY.cluster(cls(k)).name,8) & ...
                            (~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13))
                        if strcmpi(STUDY.cluster(cls(k)).name, clust_name)
                            cind = comp_ind(ci) - num_comps; % component index in the cluster
                            subject = STUDY.datasetinfo( STUDY.cluster(cls(k)).sets(1,cind)).subject;
                            a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''clusters'','  num2str(cls(k)) ', ''comps'',' num2str(cind) ' );' ];
                            eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                            break;
                        else
                            num_comps = num_comps + length(STUDY.cluster(cls(k)).comps);
                        end
                    end
                end
            end
        end
        userdat{1}{2} = STUDY;
        set(hdl, 'userdat',userdat);
    end;
end;

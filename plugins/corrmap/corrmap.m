% corrmap() -  compute correlations between IC template from a single dataset and all ICs
%              from current STUDY, cluster ICs above correlation threshold, and display
%              summary plot containing topographical maps for clustered ICs, average map,
%              correlation plot and similarity index plot.
%
% Usage:
%          >> [CORRMAP,STUDY,ALLEEG] = corrmap(STUDY, ALLEEG, n_tmp,index);
%
% Inputs:
%  STUDY - input STUDY structure
%  ALLEEG - input ALLEEG structure
%  n_tmp - number of the template dataset
%  index - index for component of template dataset that is going to be
%          correlated with all ICs from all datasets stored in ALLEEG
%
% Optional inputs: ('key', 'value')
%   'chanlocs'  - [string] path for a channel locations file. It is used to interpolate
%               missing channels ONLY if datasets in ALLEEG have different
%               number of channels. {default: ''}
%   'th'        - [string] absolute correlation threshold (0<th<1).
%               Only ICs with a correlation value above "th" are plotted.
%               In "auto" mode a suggested threshold  is calculated by corrmap and only the second summary
%               plot is displayed. {default: 'auto'}
%   'ics'       - [integer] number of ICs with correlation above "th" that will be selected from each dataset.
%                Maximum value allowed: ics=3. {default: 'ics=2'}
%   'pl'        - [string] diplay options.  {default:'2nd',when th='auto'; default:'both',when 0<th<1}
%               'none' - plots are not displayed - only acessible from command line;
%               '2nd'  - only second summary plot is displayed (template: average topo map);
%               'both' - both plots are displayed;
%   'title'     - [string] title for the summary plots created when pl~='none'. {default:''}
%   'clname'    - [string] name for the cluster that will be stored in the STUDY structure.
%                 ATTENTION: STUDY needs to be preclustered before running
%                 corrmap(). {default:''}
%   'badcomps'  - [string] storing clustered ICs in each dataset in a field called EEG.badcomps. {default:'no'}
%                 'yes' - EEG.badcomps are created
%                 'no'  - EEG.badcomps are not stored
%   'plot'      - ['on'|'off'] plot results
%   'resetcluster' - ['on'|'off'] reset the parent cluster to include all
%                 components. If this is not performed only component which
%                 where preselected by residual variace are included in the
%                 final clusters (although ALL components are always used
%                 during the calculation phase).
%
% Outputs:
% CORRMAP - structure that contains information about template,input
%           parameters,correlation values and clustered ICs.
% STUDY - STUDY structure updated
% ALLEEG - ALLEEG structure
%
% Examples:
%
% Running in "auto" mode:
%
% [CORRMAP,STUDY,ALLEEG]=corrmap(STUDY,ALLEEG,1,1,'chanlocs','Q:\myfolder\mydata\68_channel_besa_sphere.elp',...
%                                 ...'th','auto','ics',2,'title','plot','clname','name');
%
% Running in "manual" mode:
%
% [CORRMAP,STUDY,ALLEEG]=corrmap(STUDY,ALLEEG,1,1,'th','0.90','ics',3,'pl','none','title','plot','clname','name');
%
%
% See also:  pop_corrmap(), corrmap_plot_v1(),corrmap_plot_v2()
%
% Author: F. Campos-Viola, MRC-IHR, 16/07/2007 (f.viola@soton.ac.uk)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) F. Campos-Viola, MRC-IHR
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

% revised by F Campos-Viola - corrmap1.01 (30/01/2009)

% revised by F Campos-Viola - adding extra feature to save EEG.badcomps. (28/01/2009)

% revised by F Campos-Viola - preclustering warning. (27/01/2009);

% revised by F Campos-Viola - added a second set of iterations for auto mode
% option. (29/02/2008)

% revised by F Campos-Viola - corrmap1.00 (06/02/2008)


function [CORRMAP,STUDY,ALLEEG]=corrmap(STUDY,ALLEEG,n_tmp,index,varargin)
tic
CORRMAP=[];


%%%%% BLOCK 1 - checking if mandatory input parameters exist %%%%%%%%%%%%%
if nargin < 3
    fprintf('ERROR: Not enough input parameters. Please, check: help corrmap. \n');
    return
end

% dataset indices (only unique ICAs)
% ----------------------------------
sameicas = std_findsameica(ALLEEG);
datindsori = [];
flagtmp = 0;
for tmpi = 1:length(sameicas)
    datindsori = [ datindsori sameicas{tmpi}(1) ];
    if any(sameicas{tmpi} == n_tmp) && flagtmp == 0 n_tmp = length(datindsori); flagtmp = 1; end; % convert dataset index
end;
SELEEG = ALLEEG(datindsori);
n=length(SELEEG); % number of total datasets
CORRMAP.datasetindices = datindsori;

%checking if mandatory input parameters have "impossible values"
if n_tmp>length(SELEEG);
    fprintf('Error "n_tmp": Dataset %11.4g does not exist.Only %11.4g were loaded. \n',n_tmp,n);
    fprintf('see >> help corrmap. \n');
    return
end

tot_ics=length(SELEEG(n_tmp).icawinv); %total number of ICs in template dataset

if index>tot_ics
    fprintf('Error "index": IC %11.4g  does not exist. There are only %11.4g ICs. \n',index,tot_ics);
    fprintf('see >> help corrmap. \n');
    return
end

%%%%%%%%%%%%%% BLOCK 2 - decoding input parameters %%%%%%%%%%%%
g = finputcheck(varargin, { 'chanlocs' { 'string','struct' } {[] []} '';...
                            'th'     'string'    []     'auto' ;...
                            'ics'    'integer'  [1 2 3]    2 ;...
                            'pl'     'string'   {'none','2nd', 'both'}     '2nd';...
                            'resetclusters'     'string'   {'on','off'}    'off';...
                            'plot'              'string'   {'on','off'}    'on';...
                            'title'  'string'   []    '';...
                            'clname' 'string'   []    '';...
                            'badcomps' 'string' {'yes', 'no'} 'no'});

if isstr(g), error(g); end;
opt = g;

if ~isempty(g.chanlocs)
    chanlocs = g.chanlocs;
else
    chanlocs='';
end

if ~isempty(g.th)
    th = g.th;
else
    th='auto';
end

if ~isempty(g.ics)
    
    if g.ics>3
        fprintf('Error "number of ICs": Maximum number allowed is 3.\n');
        fprintf('see >> help corrmap. \n');
        return
    else
        ics = g.ics;
    end
    
else
    ics=2;
end

if ~isempty(g.pl)
    pl=g.pl;
else
    pl='2nd';
end

if ~isempty(g.title)
    title=g.title;
else
    title='';
end

if ~isempty(g.clname)
    clname=g.clname;
else
    clname='';
end

if ~isempty(g.badcomps)
    badcomps=g.badcomps;
else
    badcomps='no';
end

%%%%%%%%% BLOCK 3 - creating variables %%%%%%%%%%%
chan=SELEEG(n_tmp).chanlocs; % channel information

% aux=1:n;

for ii=1:n; %auxiliary variable
    
    comp{ii}=SELEEG(ii).icawinv; %creating cell with all datasets
    ics_set(ii)=size(SELEEG(ii).icawinv,2); %number of ICs in each dataset
    names{ii}=SELEEG(ii).setname; %creating cell with all dataset names
    
    %checking if all datasets have the same number of channels
    ch(ii)=SELEEG(ii).nbchan;
    
end

tot_ics=sum(ics_set); %total number of ICs

%checking if all datasets have the same number of ICs
ch_aux=find(ch(2:end)==ch(1));

if length(ch_aux)==0 || length(ch_aux)<n-1
    
    if length(chanlocs)==0
        chanlocs=eeg_mergelocs(ALLEEG.chanlocs);
    end;
    
    % code to interpolate missing channels - calls interp_chan()
    fprintf('Warning: Datasets do not have the same number of channels.\nMissing channels will be interpolated using channels location file.\n');
    fprintf('>> \n');
    clear comp
    %         comp1={SELEEG(aux).icawinv}; %creating cell with all datasets
    CHANS = readlocs(chanlocs);

    %%%%%%%%%%% Added by Romain 20 Aug, 2010 %%%%%%%%%%

    % Get labels for channels used in all datasets
    nbDatasets = length(SELEEG);
    chanLabels = {};
    nbChans = 0;
    for i=1:nbDatasets
        for j=1:length({SELEEG(i).chanlocs.labels})
            nbChans = nbChans+1;
            chanLabels(nbChans) = {SELEEG(i).chanlocs(j).labels};
        end
    end
    % Get unique labels
    uniChanLabels =  unique(chanLabels);

    % Find matching channels indexes in loc file
    chanIndexes = [];
    for i=1:length(uniChanLabels)
        chanIndex = find(strcmp({CHANS.labels}, uniChanLabels(i)));
        if length(chanIndex)==1
            chanIndexes(i) = chanIndex;
        elseif length(chanIndex)==0
            fprintf('>');
            fprintf('ERROR: No channel in location file match channel %s.\n', uniChanLabels{i});
            fprintf('> \n');
            return
        else
            fprintf('>');
            fprintf('ERROR: Several channels in location file match channel %s.\n', uniChanLabels{i});
            fprintf('> \n');
            return
        end
    end

    % Use only channels from loc file which are present in datasets
    CHANS = CHANS(chanIndexes);

    % Use only channels with location
    % CHANS = CHANS(~cellfun(@isempty, { chanlocs.theta })); % Modified by Romain 20 Aug. 2010
    CHANS = CHANS(~cellfun(@isempty, { CHANS.theta }));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    chan  = CHANS;
    for ind = 1:length(SELEEG)
        TMPEEG = eeg_emptyset;
        TMPEEG.data       = SELEEG(ind).icawinv;
        TMPEEG.chanlocs   = SELEEG(ind).chanlocs(SELEEG(ind).icachansind);
        emptychans        = find(cellfun(@isempty, { TMPEEG.chanlocs.theta }));
        TMPEEG.nbchan     = length(SELEEG(ind).icachansind);
        TMPEEG.trials     = 1;
        TMPEEG.pnts       = size(TMPEEG.data,2);
        TMPEEG.xmin       = 0;
        TMPEEG.srate      = SELEEG(ind).srate;
        TMPEEG.xmax       = (TMPEEG.pnts-1)/TMPEEG.srate;
        TMPEEG.icawinv    = [];
        if ~isempty(emptychans), TMPEEG = pop_select(TMPEEG, 'nochannel', emptychans); end;
        TMPEEG            = eeg_interp(TMPEEG, CHANS);
        comp{ind} = TMPEEG.data;
        tmpchans{ind} = TMPEEG.chanlocs;
    end;
    %[comp,CHANS]=interp_chan(SELEEG,chanlocs); %function to create cell with new datasets
    %chan=CHANS; % new channel locations
    %tot_ics=size(comp{1},2); %number of ICs can be different from number of channels
end

% check electrode order
% ---------------------
%tmplab = cellfun(@(x)({x.labels}), tmpchans, 'uniformoutput', false);
%tmplab2 = cellfun(@(x)[x{:}], tmplab, 'uniformoutput', false);
%strvcat(tmplab2)

%%%%%%%%% BLOCK 4 - storing OUTPUT info in CORRMAP structure %%%%%%%
CORRMAP.template.setname=ALLEEG(n_tmp).setname;
CORRMAP.template.index=n_tmp;
CORRMAP.template.ic=index;

CORRMAP.datasets.setnames=names;
CORRMAP.datasets.index=1:n;
CORRMAP.datasets.ics=ics_set;

CORRMAP.input.chanlocs=chanlocs;
CORRMAP.input.corr_th=th;
CORRMAP.input.ics_sel=ics;
CORRMAP.input.plots=pl;
CORRMAP.input.title=title;
CORRMAP.input.clname=clname;
CORRMAP.input.badcomps=badcomps;

%%%%%%%%%%%%%%%%%%%%% BLOCK 5 - loop to calculate correlations %%%%%%%%%%%%%%%
%it works in two steps: 1st - using template IC: 2nd
%using average map (result from previous calculation)

rep=cell(1,2); %used as auxiliary variable to store info to check ICs contributing from the same dataset and adding same colour to label

%checking if will operate in manual mode or auto mode
%when working in auto mode an extra loop is needed (loop "a")

if  strcmp(th,'auto')
    %parameters for auto threshold
    tres=0.95:-0.01:0.80; % Default
    leng=length(tres);
    if  strcmp(pl,'none')
         pl='none';
    else pl='2nd';
    end
    dum=1; % used as reference for plotting automatic threshold
else
    %parameters for normal threshold
    leng=1;
    tres=str2num(th);
    if  strcmp(pl,'none')
         pl='none';
    else pl='both';
    end
    dum=0;  % used as reference for plotting automatic threshold
end

for a=1:leng % first loop: only runs more than once in auto mode
    
    for m=1:2; %loop for 2 main steps
        
        h=1;
        if m==1
            tmp=comp{n_tmp}(:,index); % first step: template = selected by user
        elseif m==2
            tmp=reshape(mediaplot(a,1,:),length(chan),1); %second step: template = average map previous step
        end %closes if
        
        %%%%%%%% calculating correlation
        
        for i=1:n %loop for all datasets
            
            for j=1:size(comp{i},2)%loop for all ICs
                if i==n_tmp && j==index && m==1
                    res(j)=0; %eliminating template
                elseif exist('corr') == 2 && license('checkout', 'statistics_toolbox')
                    res(j)=corr(tmp,comp{i}(:,j));
                else
                    tmpmat = corrcoef(tmp,comp{i}(:,j));
                    if length(tmpmat) == 1
                        res(j)=tmpmat; % Octave
                    else
                        res(j)=tmpmat(2);
                    end;
                end
            end %closes "j" loop
            
            %saving correlation for selected ICs: ics==1 or ics==2 or ics==3
            for f=h:h+ics-1
                
                aux(f)=max(abs(res)); %abs max corr
                all_info(f,2)=find(abs(res)==aux(f)); %corr_index
                all_info(f,1)=res(all_info(f,2)); % max_corr
                
                warning off MATLAB:divideByZero %Warning: Divide by zero added by FViola on 04/06/08
                all_info(f,3)=0.5*log((1+all_info(f,1))/(1.00000015-all_info(f,1))); %fisher_z
                all_info(f,4)=CORRMAP.datasets.index(i); %dataset number
                res(all_info(f,2))=0;
                warning on MATLAB:divideByZero %Warning: Divide by zero added by FViola on 04/06/08
            end
            
            h=h+ics;
            
        end %closes "i" loop
        clear res
        
        %sorting correlation in descending order
        corr_dec=sortrows((all_info),-1);
        corr_dec_abs=sortrows(abs(all_info),-1);
        
        %%%%% storing OUTPUT in temporary CORRMAP structure - done this way because
        %%%%% when using auto mode
        CORRMAP_temp.corr{a}.abs_values{m}=corr_dec_abs(:,1);
        CORRMAP_temp.corr{a}.sets{m}=corr_dec_abs(:,4);
        CORRMAP_temp.corr{a}.ics{m}=corr_dec_abs(:,2);
        CORRMAP_temp.temp{a}.fz{m}=corr_dec_abs(:,3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % finding ICs above threshold
        varaux=find(corr_dec_abs(:,1)>tres(a)); %this variable is cleared in each iteration
        
        loop=isempty(varaux);
        ind=1;
        
        %error message for user
        if leng==1 && loop==1;
            fprintf('>');
            fprintf('ERROR: There are no ICs above the selected correlation threshold (th=%11.4g). \n',tres(a));
            fprintf('> \n');
            return
        else
            while loop==1;
                if ind >length(tres)-1
                    fprintf('Warning: There are no ICs above the lowest correlation threshold (th=%11.4g). \n',tres(a+ind-1));
                    if tres(a+ind-1)==0.8
                        tres=0.79:-0.01:0.55; %changed by FViola 08/10/2008
                        ind=1;
                        fprintf('CORRMAP is going to be run with a starting correlation threshold of th=%11.4g. \n',tres(a));
                        fprintf('> \n');
                    else
                        fprintf('Warning: There are no ICs above the lowest correlation threshold (th=%11.4g). \n',tres(a+ind-1));
                        fprintf('> \n');
                        return
                    end
                    
                else
                    varaux=find(corr_dec_abs(:,1)>tres(a+ind));
                    ind=ind+1;
                    loop=isempty(varaux);
                end
            end
        end
        
        % keeping ICs to use to calculate average map
        sel_ic=zeros(length(chan),length(varaux)+1);
        plot_ics_temp{a}(m)=max(varaux);
        dataset=zeros(1,n);
        
        for g=1:plot_ics_temp{a}(m)
            
            ind_c(g)=find(CORRMAP.datasets.index==corr_dec_abs(g,4)); % finding datasets above threshold
            dataset(ind_c(g))=dataset(ind_c(g))+1; %auxiliary to count number of total ICs and total sets contributing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%root mean square normalization %%%%%%%%%%%%%%%%%%%%%%
            rms=sqrt( mean( ( comp{ind_c(g)}(:,corr_dec_abs(g,2)) ).^2 ) );
            sel_ic(:,g)=(comp{ind_c(g)}(:,corr_dec_abs(g,2)))./rms;
            ttt=find(-corr_dec_abs(g,1)==all_info(:,1));
            
            % inverting polarities for ICs with negative correlations - important
            % for average map
            
            if isempty(ttt)
                sel_ic(:,g)=comp{ind_c(g)}(:,corr_dec_abs(g,2));
            else
                sel_ic(:,g)=-comp{ind_c(g)}(:,corr_dec_abs(g,2));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        rep{a}.temp{m}=dataset; %used to check ICs contributing from the same dataset and adding same colour to label
        
        %mean correlation - using fisher-z value
        corr_dec_abs(1:max(varaux),3); %using (max(varaux)) means that only ICs in "cluster" contribute to mean
        med_z=mean(corr_dec_abs(1:max(varaux),3)); %mean fisher z values
        med(a,m)=(exp(2*(med_z)) - 1)/((exp(2*(med_z)) + 1)); %mean correlation
        
        %mean plot
        mediaplot(a,m,:)=mean(sel_ic,2);
        
        clear sel_ic
        
        %finding total sets contributing
        if ics==1
            N=length(varaux); %this variable is cleared in each iteration
        else
            N=length(find(dataset>0)); %this variable is cleared in each iteration
        end
        
        %finding total sets that contribute with two ics
        number=length(find(dataset>1)); %this variable is cleared in each iteration
        
        %%%%%OUTPUT info - storing in CORRMAP structure
        CORRMAP_temp.clust.ics{a}(m)=length(varaux);
        CORRMAP_temp.clust.sets{a}.number(m)=N;
        CORRMAP_temp.clust.sets{a}.more_oneIC(m)=number;
        
        % finding datasets not contributing
        not_data=find(dataset==0);  %this variable is cleared in each iteration
        if isempty(not_data)
            CORRMAP_temp.clust.sets{a}.absent{m}=0;
        else
            CORRMAP_temp.clust.sets{a}.absent{m}=not_data;
        end
        
        CORRMAP_temp.clust.mean_corr{a}(m)=med(a,m);
        clear varaux
        clear N
        clear number
        clear not_data
        clear ind_c
        clear dataset
        clear all_info
        clear corr_dec
        clear corr_dec_abs
        
    end %closes "m" loop
    
    %correlation between average maps
    cofcor=corrcoef(squeeze(mediaplot(a,1,:)),squeeze(mediaplot(a,2,:)));
    if length(cofcor) == 1
        dif_corr=abs(cofcor);
    else
        dif_corr=abs(cofcor(2,1));
    end;
    %%%%%OUTPUT info - storing in CORRMAP structure
    CORRMAP_temp.clust.similarity(a)=dif_corr;
    
end %closes "a" loop

%%%%%%%%%%%%%% BLOCK 6 - when auto, finding the best threshold; updating CORRMAP %%%%%%%%%%%%%%%


if leng==1
    
    %updating CORRMAP structure
    CORRMAP.corr=CORRMAP_temp.corr{1};
    
    %CORRMAP.temp=CORRMAP_temp.temp;
    
    CORRMAP.clust.ics=cell2mat(CORRMAP_temp.clust.ics);
    CORRMAP.clust.sets=cell2mat(CORRMAP_temp.clust.sets);
    CORRMAP.clust.mean_corr=cell2mat(CORRMAP_temp.clust.mean_corr);
    CORRMAP.clust.similarity=CORRMAP_temp.clust.similarity;
    plot_ics=cell2mat(plot_ics_temp);
    rep=rep{1}.temp;
    tt=1; % when operating in manual mode position for threshold is always the first ("loop a" runs only once)
    flg=0;%used as flag parameter inside corrmap_plot() to identify if is operating in manual or auto mode
    
else
    
    format('short', 'e'); %using this to see if function was picking correct max
    tt=find(CORRMAP_temp.clust.similarity==max(CORRMAP_temp.clust.similarity)); %finding position for suggested threshold
    th=tres(tt(1)); %saving suggested threshold
    
    %updating CORRMAP structure
    CORRMAP.clust.best_th=th;
    CORRMAP.corr=CORRMAP_temp.corr{tt(1)};
    %CORRMAP.temp=CORRMAP_temp.temp{tt(1)};
    CORRMAP.clust.ics=CORRMAP_temp.clust.ics{tt(1)};
    CORRMAP.clust.sets=CORRMAP_temp.clust.sets{tt(1)};
    CORRMAP.clust.mean_corr= CORRMAP_temp.clust.mean_corr{tt(1)};
    CORRMAP.clust.similarity=CORRMAP_temp.clust.similarity(tt(1));
    plot_ics=plot_ics_temp{tt(1)};
    rep=rep{tt(1)}.temp;
    
    flg=1; %used as flag parameter inside corrmap_plot() to identify if is operating in manual or auto mode
end

%%%% Warning: before plotting CORRMAP structure must have all information stored


%%%% BLOCK 7 - PLOTING OUTPUT - calls corrmap_plot()
% options for plots
switch pl
    
    case 'none'
        disp('> Calculation is finished. Plots were not displayed.')
        return
        
    case '2nd'
        in=2;
        fin=in;
        disp('> Calculation is finished. Second plot is going to be displayed.')
        
    case 'both'
        in=1;
        fin=2;
        disp('> Calculation is finished. Both plots are going to be displayed.')
end

CORRMAP.output.chanlocs=chan;
CORRMAP.output.average_plot{1}=reshape(mediaplot(tt(1),1,:),length(chan),1);
CORRMAP.output.average_plot{2}=reshape(mediaplot(tt(1),2,:),length(chan),1);
CORRMAP.output.sets{1}= CORRMAP.corr.sets{1}(1:CORRMAP.clust.ics(1));
CORRMAP.output.sets{2}= CORRMAP.corr.sets{2}(1:CORRMAP.clust.ics(2));
CORRMAP.output.ics{1}= CORRMAP.corr.ics{1}(1:CORRMAP.clust.ics(1));
CORRMAP.output.ics{2}= CORRMAP.corr.ics{2}(1:CORRMAP.clust.ics(2));


if strcmpi(opt.plot, 'on')
    %%%% checking Matlab version to use correct plotting options
    vers=version; %7.3.0.267 (R2006b)or 7.5.0.342 (R2007b)
    kvers1=strfind(vers, '7.3.');
    if ~isempty(kvers1)
        corrmap_plot_v1(CORRMAP,CORRMAP_temp,comp,chan,rep,m,in,fin,plot_ics, mediaplot,tt(1),flg);
    else
        corrmap_plot_v2(CORRMAP,CORRMAP_temp,comp,chan,rep,m,in,fin,plot_ics, mediaplot,tt(1),flg);
    end
end;

%%% BLOCK 8 - saving EEG:badcomps %%%%%

if strcmp(badcomps,'yes')
    display('>EEG.badcomps are going to be stored')
    
    %storing EEG.badcomps
    for i=1:length(CORRMAP.output.sets{2})
        a=CORRMAP.output.sets{2}(i);
        if isempty(SELEEG(a).badcomps) % EEG.badcomponents is empty
            SELEEG(a).badcomps(i)=CORRMAP.output.ics{2}(i);
        else
            nbc=length(SELEEG(a).badcomps); %number of bad components already stored
            SELEEG(a).badcomps(nbc+1)=CORRMAP.output.ics{2}(i);
        end
        
    end
    
    for jj=1:length(CORRMAP.datasets.index)
        aux_ind=find(SELEEG(jj).badcomps==0);
        SELEEG(jj).badcomps(aux_ind)=[];
        CURRENTSTUDY = 0;EEG=SELEEG(jj); [SELEEG EEG CURRENTSET] = pop_newset(SELEEG, EEG, jj , 'overwrite','on', 'study',1);
        
        pop_saveset( EEG,  'filename', EEG.filename, 'filepath', EEG.filepath);
    end
    
end

%%% BLOCK 9 - creating cluster in STUDY structure %%%%%
val_cl=strcmp(CORRMAP.input.clname,''); %evaluating if cluster is going to be saved

% creating and saving cluster in STUDY structure
if val_cl==0
    
    if strcmpi(opt.resetclusters, 'on')
        STUDY.cluster = [];
        for index = 1:length(ALLEEG)
            STUDY.datasetinfo(index).comps = 1:size(ALLEEG(index).icaweights,1);
        end;
        STUDY = std_checkset(STUDY, ALLEEG);
    end;
    
    % create cluster indices
    datinds  = datindsori(CORRMAP.output.sets{2});
    compinds = CORRMAP.output.ics{2};
    ncomps = length(STUDY.cluster(1).comps);
    clust = zeros(1,ncomps);
    for index = 1:length(datinds)
        allsets = STUDY.cluster(1).sets';
        tmpdatind  = mod(find(allsets(:) == datinds(index))-1, ncomps)+1;
        tmpcompind = find( STUDY.cluster(1).comps(tmpdatind) == compinds(index) );
        oriind = tmpdatind(tmpcompind);
        clust(oriind) = 1;
    end;
    
    STUDY = std_createclust(STUDY,ALLEEG, 'name', CORRMAP.input.clname, 'clusterind', clust, 'ignore0', 'on');
    %sname= CORRMAP.input.clname;
    %namestud=strcat(sname,'.study');
    %[STUDY ALLEEG] = pop_savestudy(STUDY, ALLEEG,'filename',namestud,'filepath',STUDY.filepath);
    
end

fprintf('\n')
fprintf('done \n')
toc

%Useful info

% Warning: Divide by zero.
%   (Type "warning off MATLAB:divideByZero" to suppress this warning.)

% replacement function for the corr function
function corrmat = replacement_corr(a, b);

corrmat = zeros(size(a,2), size(b,2));
for inda = 1:size(a,2)
    for indb = 1:size(b,2)
        tmpmat = corrcoef(a(:,inda), b(:,indb));
        corrmat(inda, indb) = tmpmat(2);
    end;
end;

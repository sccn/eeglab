% eeg_region_peaks_script - analyse regional ERP peaks

eeg_load_data_script;

% Define ERP time windows of analysis,
% corresponding to ERP component peak(s)
times = [300 700];

clear exp con

for g=1:length(groups),
	for s=1:length(subs),
        
        for t=1:size(times,1),
            
            % Find regional peaks for experimental condition
            
            sub = sprintf('%s%02d%s%s',groups{g},subs(s),cond{1},ext);
            fprintf('Selecting %s\n',sub);
            p.volt.data = eval(sub);
		    p.volt.timeArray = timeArray;
            
            regions = eeg_region_peaks(p,times(t,:));
            
            exp(t,g).sub(s,:)      = {sub};
            exp(t,g).group(s,:)    = groups(g);
            exp(t,g).cond(s,:)     = cond(1);
            exp(t,g).times(s,:)    = times(t,:);
            exp(t,g).names(s,:)    = {regions.name};
            exp(t,g).pospeaks(s,:) = [regions.pospeaks];
            exp(t,g).postimes(s,:) = [regions.postimes];
            exp(t,g).poselecs(s,:) = [regions.poselecs];
            exp(t,g).negpeaks(s,:) = [regions.negpeaks];
            exp(t,g).negtimes(s,:) = [regions.negtimes];
            exp(t,g).negelecs(s,:) = [regions.negelecs];
            
            
            % Run control condition, being careful to select the
            % same electrodes in each region where peaks are
            % detected in the above experimental condition
            
            sub = sprintf('%s%02d%s%s',groups{g},subs(s),cond{2},ext);
            fprintf('Selecting %s\n',sub);
            p.volt.data = eval(sub);
		    p.volt.timeArray = timeArray;
	
            % obtain all +ve/-ve peak electrodes from the experimental
            % condition to redefine the region structure array and
            % thereby control exp/con electrode selection
            for r=1:length(regions),
                posregion(r).name = regions(r).name;
                posregion(r).elec = regions(r).poselecs;
                negregion(r).name = regions(r).name;
                negregion(r).elec = regions(r).negelecs;
            end
            
            clear regions
            
            regions = eeg_region_peaks(p,times(t,:),posregion);
            
            con(t,g).sub(s,:)      = {sub};
            con(t,g).group(s,:)    = groups(g);
            con(t,g).cond(s,:)     = cond(2);
            con(t,g).times(s,:)    = times(t,:);
            con(t,g).names(s,:)    = {regions.name};
            con(t,g).pospeaks(s,:) = [regions.pospeaks];
            con(t,g).postimes(s,:) = [regions.postimes];
            con(t,g).poselecs(s,:) = [regions.poselecs];
            
            clear regions
            
            regions = eeg_region_peaks(p,times(t,:),negregion);
            
            con(t,g).negpeaks(s,:) = [regions.negpeaks];
            con(t,g).negtimes(s,:) = [regions.negtimes];
            con(t,g).negelecs(s,:) = [regions.negelecs];
            
            clear regions posregion negregion
            
        end
	end
end

clear g s r sub p

regions = elec_regions;

% For export to SPSS, use this format
Sub = char([exp(1,1).sub; exp(1,2).sub]);
Gp  = [ones(length(exp(1,1).sub),1); ones(length(exp(1,1).sub),1) * 2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to repeat for loop below for postimes/negpeaks/negtimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t = 1:size(times,1),
    
	[path,name,ext] = fileparts(strcat(eegpath,'regions_',num2str(times(t)),'.txt'));
	file = fullfile(path,[name ext]);
    
    fid = fopen(file,'w');
	
	fprintf(fid,'%-15s\t%-8s\t','Subjects','Group');
	ExpCond = strcat(exp(1,1).names(1,:),'exp');
	fprintf(fid,'%-10s\t',ExpCond{:});
	ConCond = strcat(con(1,1).names(1,:),'con');
	fprintf(fid,'%-10s\t',ConCond{:});
	fprintf(fid,'\n');
	
	Exp = [exp(t,1).pospeaks; exp(t,2).pospeaks];
	Con = [con(t,1).pospeaks; con(t,2).pospeaks];
	
	Data = [Exp,Con];
	
	format = '%-15s\t%-8d\t';
	for i=1:length(regions)*2,format = strcat(format,'%10.4f\t'); end
	format = strcat(format,'\n');
	
	for i=1:length(Sub),
        fprintf(fid,format,Sub(i,:),Gp(i),Data(i,:));
	end
    fclose(fid);
    
    fprintf('\n\nRegion peaks results in file:\n... %s\n',file);
    fprintf('\nThis file can be imported into SPSS as text.\n');
end




% Cannot run full MANOVA with matlab (04/2002)
return




% Run stats in matlab
MODEL = 'full';
SSTYPE = 3;
DISPLAYOPT = 'off';
for t = 1:size(times,1),
	for r = 1:length(regions),
        
        cont = [exp(t,1).pospeaks(:,r); con(t,1).pospeaks(:,r)];
        ptsd = [exp(t,2).pospeaks(:,r); con(t,2).pospeaks(:,r)];
        
        Y = [cont; ptsd];
        
        GROUP1 = {'CONT'};
        GROUP2 = {'PTSD'};
        GROUPS = [ repmat(GROUP1,size(exp(t,1).pospeaks(:,r))); ...
                   repmat(GROUP1,size(con(t,1).pospeaks(:,r))); ...
                   repmat(GROUP2,size(exp(t,2).pospeaks(:,r))); ...
                   repmat(GROUP2,size(con(t,2).pospeaks(:,r))) ];
        
        EXP = {'tac'};
        CON = {'oac'};
        CONDITIONS = [ repmat(EXP,size(exp(t,1).pospeaks(:,r))); ...
                       repmat(CON,size(con(t,1).pospeaks(:,r))); ...
                       repmat(EXP,size(exp(t,2).pospeaks(:,r))); ...
                       repmat(CON,size(con(t,2).pospeaks(:,r))) ];
        
        GROUP = {GROUPS CONDITIONS};
        GNAME = strvcat('Group', 'Task');
        
        [anova(t,r).p,anova(t,r).table] = ANOVAN(Y,GROUP,MODEL,SSTYPE,GNAME,DISPLAYOPT) ;
        
	end
end
clear Y GROUP MODEL SSTYPE GNAME DISPLAYOPT;

% eeg_region_dif_peaks_script - analyse regional ERP peaks

eeg_load_difdata_script;

% Define ERP time windows of analysis,
% corresponding to ERP component peak(s)
% Each row is start time, peak time, end time (msec)
times = [150 300 450; 300 550 700];

% Replace missing data with peak time potential
missing = 0;

clear dif

for g=1:length(groups),
	for s=1:length(subs),
        
        for t=1:size(times,1),
            
            % Find regional peaks for difference waves
            
            sub = sprintf('%s%02d%s%s',groups{g},subs(s),cond{1},ext);
            fprintf('Selecting %s\n',sub);
            p.volt.data = eval(sub);
		    p.volt.timeArray = timeArray;
            
            regions = eeg_region_peaks(p,times(t,:),[],missing);
            
            dif(t,g).sub(s,:)      = {sub};
            dif(t,g).group(s,:)    = groups(g);
            dif(t,g).cond(s,:)     = cond(1);
            dif(t,g).times(s,:)    = times(t,:);
            dif(t,g).names(s,:)    = {regions.name};
            dif(t,g).pospeaks(s,:) = [regions.pospeaks];
            dif(t,g).postimes(s,:) = [regions.postimes];
            dif(t,g).poselecs(s,:) = [regions.poselecs];
            dif(t,g).negpeaks(s,:) = [regions.negpeaks];
            dif(t,g).negtimes(s,:) = [regions.negtimes];
            dif(t,g).negelecs(s,:) = [regions.negelecs];
            
            clear regions posregion negregion
            
        end
	end
end

clear g s r sub p


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TO SPSS TEXT FILES

regions = elec_regions;

% For export to SPSS, use this format
Sub = char([dif(1,1).sub; dif(1,2).sub]);
Gp  = [ones(length(dif(1,1).sub),1); ones(length(dif(1,1).sub),1) * 2];

status = mkdir(eegpath,'regions');
if status > 0,
    [path,name,ext] = fileparts(strcat(eegpath,filesep,'regions',filesep,''));
    outpath = fullfile(path);
end

for t = 1:size(times,1),
    
    % Define output file path/names
    [path,name,ext] = fileparts(strcat(outpath,filesep,'eeg_regions_dif_p',sprintf('%04d',times(t,2)),'_amp.txt'));
    pampfile = fullfile(path,[name ext]);
    [path,name,ext] = fileparts(strcat(outpath,filesep,'eeg_regions_dif_p',sprintf('%04d',times(t,2)),'_lat.txt'));
    platfile = fullfile(path,[name ext]);
    [path,name,ext] = fileparts(strcat(outpath,filesep,'eeg_regions_dif_n',sprintf('%04d',times(t,2)),'_amp.txt'));
    nampfile = fullfile(path,[name ext]);
    [path,name,ext] = fileparts(strcat(outpath,filesep,'eeg_regions_dif_n',sprintf('%04d',times(t,2)),'_lat.txt'));
    nlatfile = fullfile(path,[name ext]);
    
    % Open the output files for writing
    pampfid = fopen(pampfile,'w');
    platfid = fopen(platfile,'w');
    nampfid = fopen(nampfile,'w');
    nlatfid = fopen(nlatfile,'w');
	
    % Create headings
    fprintf(pampfid,'%-15s\t%-8s\t','Subjects','Group');
    fprintf(platfid,'%-15s\t%-8s\t','Subjects','Group');
    fprintf(nampfid,'%-15s\t%-8s\t','Subjects','Group');
    fprintf(nlatfid,'%-15s\t%-8s\t','Subjects','Group');
    
    DifCond = strcat(dif(1,1).names(1,:),'dif');
    fprintf(pampfid,'%-10s\t',DifCond{:}); fprintf(pampfid,'\n');
    fprintf(platfid,'%-10s\t',DifCond{:}); fprintf(platfid,'\n');
    fprintf(nampfid,'%-10s\t',DifCond{:}); fprintf(nampfid,'\n');
    fprintf(nlatfid,'%-10s\t',DifCond{:}); fprintf(nlatfid,'\n');
    
    % Create data output format string
    format = '%-15s\t%-8d\t';
	for i=1:length(regions),format = strcat(format,'%10.4f\t'); end
	format = strcat(format,'\n');
	
    % Pos peaks output
    Dif = [dif(t,1).pospeaks; dif(t,2).pospeaks];
    Data = [Dif];
	for i=1:length(Sub),
        fprintf(pampfid,format,Sub(i,:),Gp(i),Data(i,:));
	end
    % Pos latency output
    Dif = [dif(t,1).postimes; dif(t,2).postimes];
    Data = [Dif];
	for i=1:length(Sub),
        fprintf(platfid,format,Sub(i,:),Gp(i),Data(i,:));
	end
    % Neg peaks output
    Dif = [dif(t,1).negpeaks; dif(t,2).negpeaks];
    Data = [Dif];
	for i=1:length(Sub),
        fprintf(nampfid,format,Sub(i,:),Gp(i),Data(i,:));
	end
    % Neg latency output
    Dif = [dif(t,1).negtimes; dif(t,2).negtimes];
    Data = [Dif];
	for i=1:length(Sub),
        fprintf(nlatfid,format,Sub(i,:),Gp(i),Data(i,:));
	end
    
    fclose('all');
    
    fprintf('\n\nRegion peaks results in files:\n');
    fprintf('... %s\n',pampfile);
    fprintf('... %s\n',platfile);
    fprintf('... %s\n',nampfile);
    fprintf('... %s\n',nlatfile);
end
fprintf('\nThese files can be imported into SPSS as tab delimited text.\n');




return

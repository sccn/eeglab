
% eeg_anova_script - example task ERP analysis

%  Analysis script for ERP components in
%  fixed and variable task waveforms for PTSD-PET project

% Define sample rate & epoch parameters
fprintf('Defining ERP parameters\n');
sample_rate = 2.5;
epoch_start = -200;
epoch_end = 1500;
points = 680;
times = linspace(epoch_start,(epoch_end - sample_rate),points);

% Define component timing & rows
comp = [145 195 250 385 550 750 950 1150 1350];
comp_rows = (comp - epoch_start) / sample_rate;

% Load data (680 rows x 124 cols)
path = 'd:\data_emse\ptsdpet\link14hz\';
for i=1:10,
    file = sprintf('%sc%02do_link14hz.txt',path,i);
    var  = sprintf('c%02do_link14hz',i);
    if ~exist(var,'var')
        fprintf('loading %s\n',file);
        load(file);
    end
    file = sprintf('%sc%02dt_link14hz.txt',path,i);
    var  = sprintf('c%02dt_link14hz',i);
    if ~exist(var,'var')
        fprintf('loading %s\n',file);
        load(file);
    end
    file = sprintf('%sp%02do_link14hz.txt',path,i);
    var  = sprintf('p%02do_link14hz',i);
    if ~exist(var,'var')
        fprintf('loading %s\n',file);
        load(file);
    end
    file = sprintf('%sp%02dt_link14hz.txt',path,i);
    var  = sprintf('p%02dt_link14hz',i);
    if ~exist(var,'var')
        fprintf('loading %s\n',file);
        load(file);
    end
end

% Select data rows for each component

cont = zeros(1,2,124);
ptsd = zeros(1,2,124);

for i = 1:max(size(comp))
    c = comp(i);
    t = comp_rows(i);
    for j=1:10,
        sub = sprintf('c%02do_link14hz',j);
        fprintf('Component %d, selecting subject %s\n',c,sub);
        data = eval(sub);
        cont(j,1,:) = data(t,:);

        sub = sprintf('c%02dt_link14hz',j);
        fprintf('Component %d, selecting subject %s\n',c,sub);
        data = eval(sub);
        cont(j,2,:) = data(t,:);
        
        sub = sprintf('p%02do_link14hz',j);
        fprintf('Component %d, selecting subject %s\n',c,sub);
        data = eval(sub);
        ptsd(j,1,:) = data(t,:);
        
        sub = sprintf('p%02dt_link14hz',j);
        fprintf('Component %d, selecting subject %s\n',c,sub);
        data = eval(sub);
        ptsd(j,2,:) = data(t,:);
    end
    
    % Combine component data into single matrix
    data = [cont;ptsd]; % 20x2x124
    
    % Output ANOVA analysis.
    file = sprintf('%scomp%d_anova.txt',path,c);
    ANOVA = fopen(file,'w');
    file = sprintf('%scomp%d_p_task.txt',path,c);
    TASK = fopen(file,'w');
    file = sprintf('%scomp%d_p_group.txt',path,c);
    GROUP = fopen(file,'w');
    file = sprintf('%scomp%d_p_taskxgroup.txt',path,c);
    TASKxGROUP = fopen(file,'w');
    
    formatstr = '';
    
    fprintf('Doing ANOVA on electrode:\n');
    
    x = 0;
    for k=1:124,
        
        % Extract data for analysis
        anova_data = data(:,:,k);
        
        if (x==24),
            x=0;
            fprintf('%3d\n',k);
        else,
            x = x + 1;
            fprintf('%3d ',k);
        end
        fprintf(ANOVA,'\n\nELECTRODE %3d\n',k);
        
        % Calculate & output summary stats
        avg.cont = mean(anova_data( 1:10,:));
        avg.ptsd = mean(anova_data(10:20,:));
        sd.cont = std(anova_data( 1:10,:));
        sd.ptsd = std(anova_data(10:20,:));
        se.cont = sd.cont / sqrt(10);
        se.ptsd = sd.ptsd / sqrt(10);
        fprintf(ANOVA,'%-12s\t%-12s\t%12s\t%12s\t%12s\n','Group','Task','Mean','StDev','SE');
        fprintf(ANOVA,'%-12s\t%-12s\t%12s\t%12s\t%12s\n','-----','----','----','-----','--');
        fprintf(ANOVA,'%-12s\t%-12s\t%12.6f\t%12.6f\t%12.6f\n','Cont','Fixed',   avg.cont(1),sd.cont(1),se.cont(1));
        fprintf(ANOVA,'%-12s\t%-12s\t%12.6f\t%12.6f\t%12.6f\n','Cont','Variable',avg.cont(2),sd.cont(2),se.cont(2));
        fprintf(ANOVA,'%-12s\t%-12s\t%12.6f\t%12.6f\t%12.6f\n','PTSD','Fixed',   avg.ptsd(1),sd.ptsd(1),se.ptsd(1));
        fprintf(ANOVA,'%-12s\t%-12s\t%12.6f\t%12.6f\t%12.6f\n\n','PTSD','Variable',avg.ptsd(2),sd.ptsd(2),se.ptsd(2));
        
        % Run 2way ANOVA
        [p, tab] = anova2(anova_data,10,'off');
        
        % Collate the p values
        task(k) = p(1);
        group(k) = p(2);
        taskxgroup(k) = p(3);
        
        % Output ANOVA table
        tab(2,1) = {'Task'};
        tab(3,1) = {'Group'};
        tab(4,1) = {'Task x Group'};
        for m=1:6,
            for n=1:6,
                if n==1, fprintf(ANOVA,'%-12s\t',tab{m,n});
                else
                    if m==1, fprintf(ANOVA,'%12s\t',tab{m,n});
                    else     fprintf(ANOVA,'%12.6f\t',tab{m,n});
                    end
                end
            end
            fprintf(ANOVA,'\n');
        end
        
        % This string is used to output all the p values below
        formatstr = strcat(formatstr,'%12.6f\t');
    
    end
    fclose(ANOVA);
    
    % Output the p values
    fprintf(TASK,formatstr,task); fclose(TASK);
    fprintf(GROUP,formatstr,group); fclose(GROUP);
    fprintf(TASKxGROUP,formatstr,taskxgroup); fclose(TASKxGROUP);
    
end

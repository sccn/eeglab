% eeg_grand_mean_load_asc_script - load grand mean ERPs

clear; close all;

inpath = 'd:\MyDocuments\emse_data\ptsd-pet\grand_mean\';

file_ext = '_scd14hz.txt';
%bin_ext  = '_scd14hz.mat';

groups = {'c' 'p'};
conditions = {'o' 'oac' 'oat' 'ouc' 'out' 't' 'tac' 'tad' 'tat' 'tuc' 'tud' 'tut'};

for c = 1:length(conditions),
    
    cond = char(conditions(c));
    
    for g = 1:length(groups),
        
        group = char(groups(g));

        infile = strcat(inpath,group,cond,file_ext);
        
        data = eeg_load_ascii(infile);
        
        %figure('Name',infile); plot(data)
    end
end

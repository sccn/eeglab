% eeg_grand_mean_ascii_script - generate grand mean ERP

clear; close all;

inpath  = 'd:\MyDocuments\emse_data\ptsd-pet\scd14hz\';
outpath = 'd:\MyDocuments\emse_data\ptsd-pet\grand_mean\';

file_ext = '_scd14hz.txt';
bin_ext  = '_scd14hz.mat';

format = '%12.4f ';

groups = {'c' 'p'};
subjects = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10'};
conditions = {'o' 'oac' 'oat' 'ouc' 'out' 't' 'tac' 'tad' 'tat' 'tuc' 'tud' 'tut'};

for c = 1:length(conditions)        cond = char(conditions(c));
    for g = 1:length(groups)        group = char(groups(g));
        for s = 1:length(subjects)  subj = strcat(group,char(subjects(s)));

            infile = strcat(inpath,subj,cond,file_ext);
            input = eeg_load_ascii(infile);
            
            if(s == 1)   add = zeros(size(input));
            else         add = add + input;
            end
        end
        
        avg = add / s;
        
        % output grand mean
        outfile = strcat(outpath,group,cond,file_ext);
        eeg_write_ascii(outfile,avg,format);
        
    end
end

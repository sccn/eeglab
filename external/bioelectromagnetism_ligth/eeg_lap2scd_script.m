% eeg_lap2scd_script - load ERP Laplacian data and convert to SCD
% Author: Darren.Weber_at_radiology.ucsf.edu
% calls eeg_load_ascii, eeg_write_ascii, & eeg_lap2scd

clear; close all;

inpath  = 'd:\MyDocuments\emse_data\ptsd-pet\lap14hz\';
outpath = 'd:\MyDocuments\emse_data\ptsd-pet\scd14hz\';

infile_ext  = '_lap14hz.txt';
outfile_ext = '_scd14hz.txt';

format = '%12.4f ';


groups = {'c' 'p'};
subjects = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10'};
conditions = {'o' 'oac' 'oat' 'ouc' 'out' 't' 'tac' 'tad' 'tat' 'tuc' 'tud' 'tut'};

for g = 1:length(groups)
    for s = 1:length(subjects)        subj = strcat(char(groups(g)),char(subjects(s)));
        for c = 1:length(conditions)  cond = char(conditions(c));
            
            infile = strcat(inpath,subj,cond,infile_ext);    % load input Laplacian ascii data
            lap = eeg_load_ascii(infile);                    % for this data, lap in uV/m^2
            
            scd = eeg_lap2scd(lap);                          % calculate scalp current density
                                                             % given lap uV/m^2, scd in uA/m^3
            
            outfile = strcat(outpath,subj,cond,outfile_ext); % output scalp current density data
            eeg_write_ascii(outfile,scd,format);
        end
    end
end

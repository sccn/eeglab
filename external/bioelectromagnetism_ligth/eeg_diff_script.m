% eeg_diff_script - calculate ERP difference waves
% Author: Darren.Weber_at_radiology.ucsf.edu

path = 'd:\MyDocuments\emse_data\ptsd-pet\link14hz\';
% path = '/floyd3/psdlw/ptsd-pet/original-scan/';

infile_ext = '_link14hz.txt';
outfile_ext = 'sa_dif_link14hz.txt';

format = '%12.4f ';


groups = {'c' 'p'};
subjects = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10'};
conditions = {'oac' 'ouc'}; % 'out' 't' 'tac' 'tad' 'tat' 'tuc' 'tud' 'tut'};

for g = 1:length(groups)
    for s = 1:length(subjects)
        subj = strcat(char(groups(g)),char(subjects(s)));
        
        % load input data
        infile1 = strcat(path,subj,conditions(1),infile_ext);
        infile2 = strcat(path,subj,conditions(2),infile_ext);
        
        data1 = eeg_ascii_load(infile1);
        data2 = eeg_ascii_load(infile2);
            
        % calculate difference
        dif = data1 - data2;
            
        % output difference data
        outfile = strcat(path,subj,outfile_ext);
        eeg_ascii_write(outfile,dif,format);
        
        %return

    end
end

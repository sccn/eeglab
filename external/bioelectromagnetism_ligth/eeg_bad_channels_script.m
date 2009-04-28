
% eeg_bad_channels_script - process files with eeg_bad_channels


clear all;
close all;

path = '/floyd3/psdlw/ptsd-pet/original-scan/';
volt_ext = '.avg';
%var_ext = 'var.asc'; % Not used, assumed to be voltfile.var

elec_path = '/floyd3/psdlw/ptsd-pet/3dd_mean/corrected_elp/';
elec_ext = 'txt';


groups = {'c' 'p'};
subjects = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10'};
conditions = {'o' 'oac' 'ouc' 'out' 't' 'tac' 'tad' 'tat' 'tuc' 'tud' 'tut'};

for g=1:length(groups),
    
    for s=1:length(subjects),
        
        for c = 1:length(conditions),
            
            subj = strcat(char(groups(g)),char(subjects(s)));
            cond = char(conditions(c));
            
           [p] = eeg_toolbox_defaults('create');
            
            p.volt.path = path;
            p.volt.file = strcat(subj,cond,volt_ext);
            p.volt.type = 'Scan3x';
            
           [p] = eeg_open(p);
            
            p.elec.path = elec_path;
            p.elec.file = strcat(subj,'_124fit',elec_ext);
            p.elec.type = 'Cartesian';
            
           [p] = elec_open(p);
            
           [p] = eeg_bad_channels(p);
        end
    end
end

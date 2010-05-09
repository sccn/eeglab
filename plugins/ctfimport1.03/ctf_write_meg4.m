function writeTrials(datfile,savefile,types,bottop_perc)

% writeTrials(datfile,savefile,types,bottop_perc)
% datfile = *.dat file used to create the appropriate types
% savefile = file to write trials to
% types = trial type codes as discussed in groupDatHiN


DATDIR = '/data/dnl3/HighN_cuebase/dat_session_files/';
SAVEDIR = '/data/dnl3/HighN_cuebase/trialLists/';

[alltypes,resps] = groupDatHiN([DATDIR datfile]);
vals = selectDatHiN(alltypes,resps,types);

if (nargin == 4)
    [tmp,sortorder] = sort(vals.Latencies);
    trials = vals.ETrials(sortorder);
    
    bottom_trials = trials(end:-1:end-round(bottop_perc*length(trials)/100)+1);
    top_trials = trials(1:round(bottop_perc*length(trials)/100));
    
    fid = fopen([SAVEDIR savefile 'top' num2str(bottop_perc)],'w');
    for i = 1:length(top_trials)
        fprintf(fid,[num2str(top_trials(i)) ' ']);
    end
    fclose(fid);
    
    fid = fopen([SAVEDIR savefile 'bottom' num2str(bottop_perc)],'w');
    for i = 1:length(bottom_trials)
        fprintf(fid,[num2str(bottom_trials(i)) ' ']);
    end
    fclose(fid);
else
    fid = fopen([SAVEDIR savefile],'w');
    for i = 1:length(vals.ETrials)
        fprintf(fid,[num2str(vals.ETrials(i)) ' ']);
    end
    fclose(fid);
end

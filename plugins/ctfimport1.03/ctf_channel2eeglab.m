indir = '/data/dnl4/ctf_data/hin/cuebase/';
outdir = '/home/awu9/eeglab_tutorial/';
outputFile = 'data.mat';
outputCHANFile = 'CHAN.xyz';
CHAN1020 = {'Fp1','Fpz','Fp2','AF7','AFz','AF8','F9','F7','F5','F3','Fz','F2','F4','F6','F8','F10','FT9','FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8','FT10','T9','T7','C5','C3','C1','Cz','C2','C4','C6','T8','T10','TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP8','P9','P7','P5','P3','P1','Pz','P2','P4','P6','P8','P8','P10','PO9','PO7','POz','PO8','PO10','O1','Oz','O2'};

cd(indir)



ctf = ctf_read;

CHAN = ctf_channel_select1020(ctf,CHAN1020);
CHANXYZ = ctf.sensor.location(:,CHAN);

fid = fopen([outdir outputCHANFile],'w');
for i = 1:length(CHAN1020);
    fprintf(fid,'%d\t%g\t%g\t%g\t%s\n',i,CHANXYZ(1,i),CHANXYZ(2,i),CHANXYZ(3,i),CHAN1020{i});
end

eegdata = ctf.data(:,CHAN,1);
eegdata = eegdata';

save([outdir outputFile],'eegdata');


function out = eeglab2cnt(eegin, outputfile)

% eeglab2cnt(eegin,outputfile)
% 
% writes eeglab ANT matrix to cnt 
% 
% eegin = input eeg matrix, which is the result of reading an asa .cnt file
% with eeglab
% outputfile = name of output file including path, including .cnt extension
%  
% 
% Define variables to put in msr file
    
charlabels = char(eegin.chanlocs.labels);
eeglabels = cellstr(charlabels);
eegrate = eegin.srate;
eegnpnt = eegin.pnts;
eegnchan = eegin.nbchan;
eegnsample = eegnpnt;
eegtime = [eegin.xmin:1/eegrate:eegin.xmax];
eegdata = double(eegin.data);

eeg.label = eeglabels;
eeg.rate = eegrate;
eeg.npnt = eegnpnt;
eeg.nchan = eegnchan;
eeg.nsample = eegnsample;
eeg.time = eegtime;
eeg.data = eegdata;

trgfile_dest = strrep(outputfile,'.cnt','.trg');

trgfile_source = strrep([eegin.filepath eegin.filename],'.cnt','.trg');
copyfile(trgfile_source,trgfile_dest);


write_eep_cnt(outputfile,eeg);
 
out = true;
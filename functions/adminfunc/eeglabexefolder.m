% eeglabexefolder - return the exe folder for EEGLAB

function str = eeglabexefolder;

str = ctfroot;
inds = find(str == filesep);
str = str(1:inds(end)-1);

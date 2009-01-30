function [data] = readmarkerfile(folder);
% function reads a CTF data file (.ds), Read 'markerfile.mrk' 
% marker is the read-only 'markerfile.mrk'
% creates structure 'data' which contains number_markers, number_samples,
%       marker_names, trial_times, marker.

%
%
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                      > %  
%      <                    DISCLAIMER:                       > %
%      <                                                      > %
%      <  THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. > %
%      <  THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   > %
%      <                    OFFICIAL USE.                     > %
%      <                                                      > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%
%
%


a = exist(fullfile(folder,'MarkerFile.mrk'));
if a == 0
  disp(' ');
  disp('Error: Your data could not be found. MarkerFile.mrk does not exist in this directory.'); 
  disp(' ');
  data = [];
  return
end   % error message if incorrect directory is displayed

markr =  textread(fullfile(folder,'MarkerFile.mrk'),'%s','delimiter','\n');
% creates markerfile.mrk in FOLDER directory. Type 'marker' to view.

% define markers
number_markers_id = strmatch('NUMBER OF MARKERS:',markr,'exact');
markers = str2num(markr{number_markers_id+1});

% define samples
number_samples_id = strmatch('NUMBER OF SAMPLES:',markr,'exact');
samples = str2num(char(markr(number_samples_id+1)));

for i = 1:length(samples)
  if samples(i) == 0
    disp('Error: one of the markers is bad.');
    data = [];
    return
  end
end

name_id = strmatch('NAME:',markr,'exact');
names = markr(name_id+1);
% defines names

trial = strmatch('LIST OF SAMPLES',markr)+2;

for i = 1:markers
  trials{i} = str2num(char(markr(trial(i):trial(i)+samples(i)))); 
  trials{i}(:,1) = trials{i}(:,1) + 1;
end

data = struct('number_markers',{markers},'number_samples',{samples},'marker_names',{names},'trial_times',{trials});   
% defines data

return

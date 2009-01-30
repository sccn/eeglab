function [data,ctf] = readepochs(folder,varargin);
%
%
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                       >
%      <                      DISCLAIMER:                      >
%      <                                                       >
%      <  THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      <  THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                     OFFICIAL USE.                     >
%      <                                                       >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
%
% function to read set time windows (epochs) around event markers.
%
[marker_info] = readmarkerfile(folder);
ctf = ctf_read_res4(folder,1);
CHAN = ctf.sensor.index.megsens;
trials = [1:ctf.setup.number_trials];
win = 'alltimes';
for i = 1:size(varargin,2)
  
  if ~iscell(varargin{i})
    
    switch num2str(varargin{i})
      
      case 'megsens',
        CHAN = ctf.sensor.index.megsens(varargin{i+1});
        
      case 'refsens',
        CHAN = ctf.sensor.index.refsens(varargin{i+1});
        
      case 'eegsens',
        CHAN = ctf.sensor.index.eegsens(varargin{i+1});
        
      case 'othersens',
        CHAN = ctf.sensor.index.othersens(varargin{i+1});
        
      case 'vc',
        CHAN = ctf.sensor.index.vc(varargin{i+1});
        
      case 'allchans',
        CHAN = [1:ctf.setup.number_channels];
        
      case 'trials',
        trials = varargin{i+1};
        
      case 'markers',
        marks = varargin{i+1};
        
      case 'window',
        win = varargin{i+1};
        
    end
  end
end

if ~exist('marks','var'),
  
  ctf = ctf_read_meg4(folder,ctf,CHAN,win,trials);
  
  epochs = cell(1,1);
  
  epochs{1} = zeros(size(ctf.data{1},1),size(ctf.data{1},2),size(ctf.data,1));
  
  for i = 1:size(ctf.data,1),
    
    epochs{1}(:,:,i) = ctf.data{i};
    
  end
  
else
  
  nm = size(marks,2);
  
  epochs = cell(1,nm);
  
  for mkr = 1:nm,
    
    mk=ismember(marker_info.marker_names,marks(mkr));
    
    marker_info.marker_names(mk);
    
    nsamp=marker_info.number_samples(mk)
    
    nss=0;
    
    for ns=1:nsamp
      
      tr = marker_info.trial_times{mk}(ns,1);
      
      if ismember(tr,trials)
        
        nss=nss+1;
        
        tim=marker_info.trial_times{mk}(ns,2);
        
        times=win+tim;
        
        ctf = ctf_read_ds(folder,ctf,CHAN,times,tr);
        
        temp=ctf.data{1};
        
        if nss==1,
          epochs{mkr}=zeros(size(temp,1),size(temp,2),1);
        end
        
        epochs{mkr}(:,:,nss)=temp;
        
      end
    end
  end
end
data = struct('epochs',{epochs},'setup',{ctf.setup},'sensor',{ctf.sensor});
return

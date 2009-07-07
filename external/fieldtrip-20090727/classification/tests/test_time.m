function [acc,sig,cv] = test_time(myproc,cvfolds)
% test classification procedure in the time domain and return accuracy plus significance
%
% good classification performance on subject 1 only:
% [a,b,c] = test_time({cspprocessor('numchan',38,'numpatterns',1,'outputdatatype','powcsp') rfda});
% witout SG filtering; 1.5 - 2.5 s; MLO + MRO
%
% a = 0.8431 (subject 1)
% a = 0.5    (subject 2)
%
%

  if nargin == 1
    cvfolds = 10;
  end
  
  %%
  % load data
  load ~/data/fieldtrip/covattpre1
  
  % subset of channels
  cvec = ismember(left.label,channelselection({'MLO' 'MRO'},left.label));

  % subset of time segment
  tvec = (left.time{1} >= 2.0 & left.time{1} <= 2.5);
  
  data = zeros([length(left.trial) + length(right.trial) sum(cvec) sum(tvec)]);
  
  for c=1:length(left.trial)
      
      % apply savitzky-golay filter; helps when not using CSP 
      x = transpose(sgolayfilt(left.trial{c}',3,41));     
      %x = left.trial{c};
      
      data(c,:,:) = x(cvec,tvec);
      
  end
  
  for c=1:length(right.trial)
     
      % apply savitzky-golay filter; helps when not using CSP  
      x = transpose(sgolayfilt(right.trial{c}',3,41));           
      %x = right.trial{c};

      data(c+length(left.trial),:,:) =  x(cvec,tvec);

  end
  
  design  = [ones(length(left.trial),1); 2*ones(length(right.trial),1)];
  
  cv = crossvalidator('procedure',myproc,'cvfolds',cvfolds,'randomize',true,'verbose',true);

  cv = cv.validate(data,design);
  
  acc = cv.evaluate;
  sig = cv.significance;
  
end
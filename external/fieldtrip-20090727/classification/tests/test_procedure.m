function [acc,sig,cv] = test_procedure(myproc,cvfolds,data,design)
% test classification procedure and return accuracy plus significance

  if nargin == 1
    cvfolds = 10;
  end

  if nargin < 4
    
    % load data
    load covattfrq1
    
    cvec = ismember(left.label,channelselection({'MLO' 'MRO'},left.label)); % subset of channels
    fvec = (left.freq >= 8 & left.freq <= 14); % subset of frequencies
    tvec = (left.time >= 1.5 & left.time <= 2.5); % subset of time segment
    
    data    = [squeeze(mean(left.powspctrm(:,cvec,fvec,tvec),4)); squeeze(mean(right.powspctrm(:,cvec,fvec,tvec),4))];
    design  = [ones(size(left.powspctrm,1),1); 2*ones(size(right.powspctrm,1),1)];
  
  end
  
  cv = crossvalidator('procedure',myproc,'cvfolds',cvfolds,'randomize',true,'verbose',true,'compact',true,'balanced',false);

  if isa(myproc{end},'transfer_learner') && ~iscell(data)
    cv = cv.validate({data data},{design design});
  else
    cv = cv.validate(data,design);
  end
  
  acc = cv.evaluate;
  sig = cv.significance;
  
end
function stat = test_fieldtrip(myproc)
% test classification procedure via fieldtrip and return stat

  % load data
  load covattfrq1

  %cfg.compact = true;
  cfg.method  = 'crossvalidate';
  cfg.cproc   = myproc;
  cfg.design  = [ones(size(left.powspctrm,1),1); 2*ones(size(right.powspctrm,1),1)];
  
  cfg.channel = {'MLO' 'MRO'};
  cfg.latency = [2.0 2.5];
  cfg.frequency = [8 14];
  cfg.avgovertime = 'yes';
  
  stat = freqstatistics(cfg,left,right);
      
end
ds_path='/export/data/test_oscill/disc04_DiscAud_Passive_correct-fEOG.ds'
info = ds2brainstorm(ds_path,0,2);
no_channels = info.gSetUp.no_channels;            
no_trials=info.gSetUp.no_trials;
no_samples = info.gSetUp.no_samples;
sample_rate=info.gSetUp.sample_rate;
preTrigPts = info.gSetUp.preTrigPts;
imegsens = info.imegsens;
ieegsens = info.ieegsens;
Time = [-preTrigPts no_samples-preTrigPts]/sample_rate;
%F=ds2brainstorm(ds_path,0,0,irefsens,[min(Time) max(Time)],1:10,2)
 [Ftmp,Channel,imegsens,ieegsens,iothersens,irefsens,grad_order_no,no_trials,filter,Time, Comment] = ...
                            ds2brainstorm(ds_path,0,0,[],[],1,2);

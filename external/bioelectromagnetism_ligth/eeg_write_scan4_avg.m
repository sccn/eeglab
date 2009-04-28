function eeg_write_scan4_avg(f,filename)

% eeg_write_scan4_avg - Write a Neuroscan 'scan4.x' average file
%
% Useage: eeg_write_scan4_avg(f,filename)
%
% where:    filename is a complete 'path\fileprefix.ext' string;
%           the file will be a neuroscan 4.x avg file.
%
%           f is a structure containing:
%           
%               f.header
%               f.electloc
%               f.data
%               f.variance
%               f.tag
%
%           See companion function: eeg_load_scan4_avg.m
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no implied or express warranties
% Created:  08/2000, Sean.Fitzgibbon_at_flinders.edu.au
% Modified: 08/2001, Darren.Weber_at_radiology.ucsf.edu
%                    distribute under GPL, with Sean's permission
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('EEG_WRITE_SCAN4_AVG [v %s]\n',ver(11:15)); tic;

fid = fopen(filename,'w');

fwrite(fid,f.header.rev,'char');
fwrite(fid,f.header.nextfile,'long');
fwrite(fid,f.header.prevfile,'long');
fwrite(fid,f.header.type,'char');
fwrite(fid,f.header.id,'char');
fwrite(fid,f.header.oper,'char');
fwrite(fid,f.header.doctor,'char');
fwrite(fid,f.header.referral,'char');
fwrite(fid,f.header.hospital,'char');
fwrite(fid,f.header.patient,'char');
fwrite(fid,f.header.age,'short');
fwrite(fid,f.header.sex,'char');
fwrite(fid,f.header.hand,'char');
fwrite(fid,f.header.med,'char');
fwrite(fid,f.header.category,'char');
fwrite(fid,f.header.state,'char');
fwrite(fid,f.header.label,'char');
fwrite(fid,f.header.date,'char');
fwrite(fid,f.header.time,'char');
fwrite(fid,f.header.mean_age,'float');
fwrite(fid,f.header.stdev,'float');
fwrite(fid,f.header.n,'short');
fwrite(fid,f.header.compfile,'char');
fwrite(fid,f.header.spectwincomp,'float');
fwrite(fid,f.header.meanaccuracy,'float');
fwrite(fid,f.header.meanlatency,'float');
fwrite(fid,f.header.sortfile,'char');
fwrite(fid,f.header.numevents,'int');
fwrite(fid,f.header.compoper,'char');
fwrite(fid,f.header.avgmode,'char');
fwrite(fid,f.header.review,'char');
fwrite(fid,f.header.nsweeps,'ushort');
fwrite(fid,f.header.compsweeps,'ushort');
fwrite(fid,f.header.acceptcnt,'ushort');
fwrite(fid,f.header.rejectcnt,'ushort');
fwrite(fid,f.header.pnts,'ushort');
fwrite(fid,f.header.nchannels,'ushort');
fwrite(fid,f.header.avgupdate,'ushort');
fwrite(fid,f.header.domain,'char');
fwrite(fid,f.header.variance,'char');
fwrite(fid,f.header.rate,'ushort');
fwrite(fid,f.header.scale,'double');
fwrite(fid,f.header.veogcorrect,'char');
fwrite(fid,f.header.heogcorrect,'char');
fwrite(fid,f.header.aux1correct,'char');
fwrite(fid,f.header.aux2correct,'char');
fwrite(fid,f.header.veogtrig,'float');
fwrite(fid,f.header.heogtrig,'float');
fwrite(fid,f.header.aux1trig,'float');
fwrite(fid,f.header.aux2trig,'float');
fwrite(fid,f.header.heogchnl,'short');
fwrite(fid,f.header.veogchnl,'short');
fwrite(fid,f.header.aux1chnl,'short');
fwrite(fid,f.header.aux2chnl,'short');
fwrite(fid,f.header.veogdir,'char');
fwrite(fid,f.header.heogdir,'char');
fwrite(fid,f.header.aux1dir,'char');
fwrite(fid,f.header.aux2dir,'char');
fwrite(fid,f.header.veog_n,'short');
fwrite(fid,f.header.heog_n,'short');
fwrite(fid,f.header.aux1_n,'short');
fwrite(fid,f.header.aux2_n,'short');
fwrite(fid,f.header.veogmaxcnt,'short');
fwrite(fid,f.header.heogmaxcnt,'short');
fwrite(fid,f.header.aux1maxcnt,'short');
fwrite(fid,f.header.aux2maxcnt,'short');
fwrite(fid,f.header.veogmethod,'char');
fwrite(fid,f.header.heogmethod,'char');
fwrite(fid,f.header.aux1method,'char');
fwrite(fid,f.header.aux2method,'char');
fwrite(fid,f.header.ampsensitivity,'float');
fwrite(fid,f.header.lowpass,'char');
fwrite(fid,f.header.highpass,'char');
fwrite(fid,f.header.notch,'char');
fwrite(fid,f.header.autoclipadd,'char');
fwrite(fid,f.header.baseline,'char');
fwrite(fid,f.header.offstart,'float');
fwrite(fid,f.header.offstop,'float');
fwrite(fid,f.header.reject,'char');
fwrite(fid,f.header.rejstart,'float');
fwrite(fid,f.header.rejstop,'float');
fwrite(fid,f.header.rejmin,'float');
fwrite(fid,f.header.rejmax,'float');
fwrite(fid,f.header.trigtype,'char');
fwrite(fid,f.header.trigval,'float');
fwrite(fid,f.header.trigchnl,'char');
fwrite(fid,f.header.trigmask,'short');
fwrite(fid,f.header.trigisi,'float');
fwrite(fid,f.header.trigmin,'float');
fwrite(fid,f.header.trigmax,'float');
fwrite(fid,f.header.trigdir,'char');
fwrite(fid,f.header.autoscale,'char');
fwrite(fid,f.header.n2,'short');
fwrite(fid,f.header.dir,'char');
fwrite(fid,f.header.dispmin,'float');
fwrite(fid,f.header.dispmax,'float');
fwrite(fid,f.header.xmin,'float');
fwrite(fid,f.header.xmax,'float');
fwrite(fid,f.header.automin,'float');
fwrite(fid,f.header.automax,'float');
fwrite(fid,f.header.zmin,'float');
fwrite(fid,f.header.zmax,'float');
fwrite(fid,f.header.lowcut,'float');
fwrite(fid,f.header.highcut,'float');
fwrite(fid,f.header.common,'char');
fwrite(fid,f.header.savemode,'char');
fwrite(fid,f.header.manmode,'char');
fwrite(fid,f.header.ref,'char');
fwrite(fid,f.header.rectify,'char');
fwrite(fid,f.header.displayxmin,'float');
fwrite(fid,f.header.displayxmax,'float');
fwrite(fid,f.header.phase,'char');
fwrite(fid,f.header.screen,'char');
fwrite(fid,f.header.calmode,'short');
fwrite(fid,f.header.calmethod,'short');
fwrite(fid,f.header.calupdate,'short');
fwrite(fid,f.header.calbaseline,'short');
fwrite(fid,f.header.calsweeps,'short');
fwrite(fid,f.header.calattenuator,'float');
fwrite(fid,f.header.calpulsevolt,'float');
fwrite(fid,f.header.calpulsestart,'float');
fwrite(fid,f.header.calpulsestop,'float');
fwrite(fid,f.header.calfreq,'float');
fwrite(fid,f.header.taskfile,'char');
fwrite(fid,f.header.seqfile,'char');
fwrite(fid,f.header.spectmethod,'char');
fwrite(fid,f.header.spectscaling,'char');
fwrite(fid,f.header.spectwindow,'char');
fwrite(fid,f.header.spectwinlength,'float');
fwrite(fid,f.header.spectorder,'char');
fwrite(fid,f.header.notchfilter,'char');
fwrite(fid,f.header.headgain,'short');
fwrite(fid,f.header.additionalfiles,'int');
fwrite(fid,f.header.unused,'char');
fwrite(fid,f.header.fspstopmethod,'short');
fwrite(fid,f.header.fspstopmode,'short');
fwrite(fid,f.header.fspfvalue,'float');
fwrite(fid,f.header.fsppoint,'short');
fwrite(fid,f.header.fspblocksize,'short');
fwrite(fid,f.header.fspp1,'ushort');
fwrite(fid,f.header.fspp2,'ushort');
fwrite(fid,f.header.fspalpha,'float');
fwrite(fid,f.header.fspnoise,'float');
fwrite(fid,f.header.fspv1,'short');
fwrite(fid,f.header.montage,'char');
fwrite(fid,f.header.eventfile,'char');
fwrite(fid,f.header.fratio,'float');
fwrite(fid,f.header.minor_rev,'char');
fwrite(fid,f.header.eegupdate,'short');
fwrite(fid,f.header.compressed,'char');
fwrite(fid,f.header.xscale,'float');
fwrite(fid,f.header.yscale,'float');
fwrite(fid,f.header.xsize,'float');
fwrite(fid,f.header.ysize,'float');
fwrite(fid,f.header.acmode,'char');
fwrite(fid,f.header.commonchnl,'uchar');
fwrite(fid,f.header.xtics,'char');
fwrite(fid,f.header.xrange,'char');
fwrite(fid,f.header.ytics,'char');
fwrite(fid,f.header.yrange,'char');
fwrite(fid,f.header.xscalevalue,'float');
fwrite(fid,f.header.xscaleinterval,'float');
fwrite(fid,f.header.yscalevalue,'float');
fwrite(fid,f.header.yscaleinterval,'float');
fwrite(fid,f.header.scaletoolx1,'float');
fwrite(fid,f.header.scaletooly1,'float');
fwrite(fid,f.header.scaletoolx2,'float');
fwrite(fid,f.header.scaletooly2,'float');
fwrite(fid,f.header.port,'short');
fwrite(fid,f.header.numsamples,'long');
fwrite(fid,f.header.filterflag,'char');
fwrite(fid,f.header.lowcutoff,'float');
fwrite(fid,f.header.lowpoles,'short');
fwrite(fid,f.header.highcutoff,'float');
fwrite(fid,f.header.highpoles,'short');
fwrite(fid,f.header.filtertype,'char');
fwrite(fid,f.header.filterdomain,'char');
fwrite(fid,f.header.snrflag,'char');
fwrite(fid,f.header.coherenceflag,'char');
fwrite(fid,f.header.continuoustype,'char');
fwrite(fid,f.header.eventtablepos,'long');
fwrite(fid,f.header.continuousseconds,'float');
fwrite(fid,f.header.channeloffset,'long');
fwrite(fid,f.header.autocorrectflag,'char');
fwrite(fid,f.header.dcthreshold,'uchar');

for n = 1:f.header.nchannels
fwrite(fid,f.electloc(n).lab,'char');
fwrite(fid,f.electloc(n).reference,'char');
fwrite(fid,f.electloc(n).skip,'char');
fwrite(fid,f.electloc(n).reject,'char');
fwrite(fid,f.electloc(n).display,'char');
fwrite(fid,f.electloc(n).bad,'char');
fwrite(fid,f.electloc(n).n,'ushort');
fwrite(fid,f.electloc(n).avg_reference,'char');
fwrite(fid,f.electloc(n).clipadd,'char');
fwrite(fid,f.electloc(n).x_coord,'float');
fwrite(fid,f.electloc(n).y_coord,'float');
fwrite(fid,f.electloc(n).veog_wt,'float');
fwrite(fid,f.electloc(n).veog_std,'float');
fwrite(fid,f.electloc(n).snr,'float');
fwrite(fid,f.electloc(n).heog_wt,'float');
fwrite(fid,f.electloc(n).heog_std,'float');
fwrite(fid,f.electloc(n).baseline,'short');
fwrite(fid,f.electloc(n).filtered,'char');
fwrite(fid,f.electloc(n).fsp,'char');
fwrite(fid,f.electloc(n).aux1_wt,'float');
fwrite(fid,f.electloc(n).aux1_std,'float');
fwrite(fid,f.electloc(n).senstivity,'float');
fwrite(fid,f.electloc(n).gain,'char');
fwrite(fid,f.electloc(n).hipass,'char');
fwrite(fid,f.electloc(n).lopass,'char');
fwrite(fid,f.electloc(n).page,'uchar');
fwrite(fid,f.electloc(n).size,'uchar');
fwrite(fid,f.electloc(n).impedance,'uchar');
fwrite(fid,f.electloc(n).physicalchnl,'uchar');
fwrite(fid,f.electloc(n).rectify,'char');
fwrite(fid,f.electloc(n).calib,'float');
end

for i = 1:f.header.nchannels
   fwrite(fid,f.data(i).header,'char');
   fwrite(fid,f.data(i).samples,'float');
end

for j = 1:f.header.nchannels
   fwrite(fid,f.variance(j).samples,'float');
end

fwrite(fid,f.tag,'char');

frewind(fid);
fclose(fid);

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return

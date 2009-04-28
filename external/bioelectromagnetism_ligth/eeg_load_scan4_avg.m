function [f,fid] = eeg_load_scan4_avg(filename)

% eeg_load_scan4_avg - Load a Neuroscan 'scan4.x' average file
%
% Useage: [f,fid] = eeg_load_scan4_avg(filename)
%
% where:    filename is a complete 'path\fileprefix.ext' string;
%           the file must be a neuroscan 4.x avg file.
%
%           f is a structure containing:
%           
%               f.header        - general header parameters
%               f.electloc      - channel specific parameters
%               f.data.header   - small channel data header
%               f.data.samples  - channel data (not uV)
%               f.variance      - channel variance
%               f.tag           - scan4.1 file tags
%
%           fid is a file handle to 'filename'
%
%           See the .m file for more details on structure fields
%           or use the fieldnames command.
%
% Note:     The f.data.samples are not in uV, they are in
%           raw amplifier units.
%           To scale data to microvolts, multiply by the channel-specific 
%           calibration factor ([f.electloc.calib]) and divide by the 
%           number of sweeps in the average ([f.electloc.n]).
%           Also, subtract baseline value ([f.electloc.baseline])
%           to be sure to obtain baselined averaged values.  eg:
%           ampData = [f.data.samples]; % elect in col, samples in row
%           baseline = repmat([f.electloc.baseline],[f.header.pnts],1);
%           calibration = repmat([f.electloc.calib],[f.header.pnts],1);
%           n = repmat([f.electloc.n],[f.header.pnts],1);
%           uV = ( ampData - baseline ) .* calibration ./ n;
%
% See also  eeg_load_scan_avg
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no implied or express warranties
% Created:  08/2000, Sean Fitzgibbon <psspf@id.psy.flinders.edu.au>
% Modified: 07/2001, Darren.Weber_at_radiology.ucsf.edu
%                    tab formatted this file and added help, esp.
%                    info about uV conversion
%           04/2002, Darren.Weber_at_radiology.ucsf.edu
%                    added file exist check
%                    added 'ieee-le' to fopen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('EEG_LOAD_SCAN4_AVG [v %s]\n',ver(11:15)); tic;

if ~isequal(exist(filename),2),
    lookfile = which(filename);
    if isempty(lookfile),
        msg = sprintf('...cannot locate %s\n', filename);
        error(msg);
    else
        filename = lookfile;
    end
end
fprintf('...reading file:\n...%s\n',filename);

fid = fopen(filename,'r','ieee-le');

h.rev               = fread(fid,12,'char');
h.nextfile          = fread(fid,1,'long');
h.prevfile          = fread(fid,1,'long');
h.type              = fread(fid,1,'char');
h.id                = fread(fid,20,'char');
h.oper              = fread(fid,20,'char');
h.doctor            = fread(fid,20,'char');
h.referral          = fread(fid,20,'char');
h.hospital          = fread(fid,20,'char');
h.patient           = fread(fid,20,'char');
h.age               = fread(fid,1,'short');
h.sex               = fread(fid,1,'char');
h.hand              = fread(fid,1,'char');
h.med               = fread(fid,20, 'char');
h.category          = fread(fid,20, 'char');
h.state             = fread(fid,20, 'char');
h.label             = fread(fid,20, 'char');
h.date              = fread(fid,10, 'char');
h.time              = fread(fid,12, 'char');
h.mean_age          = fread(fid,1,'float');
h.stdev             = fread(fid,1,'float');
h.n                 = fread(fid,1,'short');
h.compfile          = fread(fid,38,'char');
h.spectwincomp      = fread(fid,1,'float');
h.meanaccuracy      = fread(fid,1,'float');
h.meanlatency       = fread(fid,1,'float');
h.sortfile          = fread(fid,46,'char');
h.numevents         = fread(fid,1,'int');
h.compoper          = fread(fid,1,'char');
h.avgmode           = fread(fid,1,'char');
h.review            = fread(fid,1,'char');
h.nsweeps           = fread(fid,1,'ushort');
h.compsweeps        = fread(fid,1,'ushort');
h.acceptcnt         = fread(fid,1,'ushort');
h.rejectcnt         = fread(fid,1,'ushort');
h.pnts              = fread(fid,1,'ushort');
h.nchannels         = fread(fid,1,'ushort');
h.avgupdate         = fread(fid,1,'ushort');
h.domain            = fread(fid,1,'char');
h.variance          = fread(fid,1,'char');
h.rate              = fread(fid,1,'ushort');
h.scale             = fread(fid,1,'double');
h.veogcorrect       = fread(fid,1,'char');
h.heogcorrect       = fread(fid,1,'char');
h.aux1correct       = fread(fid,1,'char');
h.aux2correct       = fread(fid,1,'char');
h.veogtrig          = fread(fid,1,'float');
h.heogtrig          = fread(fid,1,'float');
h.aux1trig          = fread(fid,1,'float');
h.aux2trig          = fread(fid,1,'float');
h.heogchnl          = fread(fid,1,'short');
h.veogchnl          = fread(fid,1,'short');
h.aux1chnl          = fread(fid,1,'short');
h.aux2chnl          = fread(fid,1,'short');
h.veogdir           = fread(fid,1,'char');
h.heogdir           = fread(fid,1,'char');
h.aux1dir           = fread(fid,1,'char');
h.aux2dir           = fread(fid,1,'char');
h.veog_n            = fread(fid,1,'short');
h.heog_n            = fread(fid,1,'short');
h.aux1_n            = fread(fid,1,'short');
h.aux2_n            = fread(fid,1,'short');
h.veogmaxcnt        = fread(fid,1,'short');
h.heogmaxcnt        = fread(fid,1,'short');
h.aux1maxcnt        = fread(fid,1,'short');
h.aux2maxcnt        = fread(fid,1,'short');
h.veogmethod        = fread(fid,1,'char');
h.heogmethod        = fread(fid,1,'char');
h.aux1method        = fread(fid,1,'char');
h.aux2method        = fread(fid,1,'char');
h.ampsensitivity    = fread(fid,1,'float');
h.lowpass           = fread(fid,1,'char');
h.highpass          = fread(fid,1,'char');
h.notch             = fread(fid,1,'char');
h.autoclipadd       = fread(fid,1,'char');
h.baseline          = fread(fid,1,'char');
h.offstart          = fread(fid,1,'float');
h.offstop           = fread(fid,1,'float');
h.reject            = fread(fid,1,'char');
h.rejstart          = fread(fid,1,'float');
h.rejstop           = fread(fid,1,'float');
h.rejmin            = fread(fid,1,'float');
h.rejmax            = fread(fid,1,'float');
h.trigtype          = fread(fid,1,'char');
h.trigval           = fread(fid,1,'float');
h.trigchnl          = fread(fid,1,'char');
h.trigmask          = fread(fid,1,'short');
h.trigisi           = fread(fid,1,'float');
h.trigmin           = fread(fid,1,'float');
h.trigmax           = fread(fid,1,'float');
h.trigdir           = fread(fid,1,'char');
h.autoscale         = fread(fid,1,'char');
h.n2                = fread(fid,1,'short');
h.dir               = fread(fid,1,'char');
h.dispmin           = fread(fid,1,'float');
h.dispmax           = fread(fid,1,'float');
h.xmin              = fread(fid,1,'float');
h.xmax              = fread(fid,1,'float');
h.automin           = fread(fid,1,'float');
h.automax           = fread(fid,1,'float');
h.zmin              = fread(fid,1,'float');
h.zmax              = fread(fid,1,'float');
h.lowcut            = fread(fid,1,'float');
h.highcut           = fread(fid,1,'float');
h.common            = fread(fid,1,'char');
h.savemode          = fread(fid,1,'char');
h.manmode           = fread(fid,1,'char');
h.ref               = fread(fid,10,'char');
h.rectify           = fread(fid,1,'char');
h.displayxmin       = fread(fid,1,'float');
h.displayxmax       = fread(fid,1,'float');
h.phase             = fread(fid,1,'char');
h.screen            = fread(fid,16,'char');
h.calmode           = fread(fid,1,'short');
h.calmethod         = fread(fid,1,'short');
h.calupdate         = fread(fid,1,'short');
h.calbaseline       = fread(fid,1,'short');
h.calsweeps         = fread(fid,1,'short');
h.calattenuator     = fread(fid,1,'float');
h.calpulsevolt      = fread(fid,1,'float');
h.calpulsestart     = fread(fid,1,'float');
h.calpulsestop      = fread(fid,1,'float');
h.calfreq           = fread(fid,1,'float');
h.taskfile          = fread(fid,34,'char');
h.seqfile           = fread(fid,34,'char');
h.spectmethod       = fread(fid,1,'char');
h.spectscaling      = fread(fid,1,'char');
h.spectwindow       = fread(fid,1,'char');
h.spectwinlength    = fread(fid,1,'float');
h.spectorder        = fread(fid,1,'char');
h.notchfilter       = fread(fid,1,'char');
h.headgain          = fread(fid,1,'short');
h.additionalfiles   = fread(fid,1,'int');
h.unused            = fread(fid,5,'char');
h.fspstopmethod     = fread(fid,1,'short');
h.fspstopmode       = fread(fid,1,'short');
h.fspfvalue         = fread(fid,1,'float');
h.fsppoint          = fread(fid,1,'short');
h.fspblocksize      = fread(fid,1,'short');
h.fspp1             = fread(fid,1,'ushort');
h.fspp2             = fread(fid,1,'ushort');
h.fspalpha          = fread(fid,1,'float');
h.fspnoise          = fread(fid,1,'float');
h.fspv1             = fread(fid,1,'short');
h.montage           = fread(fid,40,'char');
h.eventfile         = fread(fid,40,'char');
h.fratio            = fread(fid,1,'float');
h.minor_rev         = fread(fid,1,'char');
h.eegupdate         = fread(fid,1,'short');
h.compressed        = fread(fid,1,'char');
h.xscale            = fread(fid,1,'float');
h.yscale            = fread(fid,1,'float');
h.xsize             = fread(fid,1,'float');
h.ysize             = fread(fid,1,'float');
h.acmode            = fread(fid,1,'char');
h.commonchnl        = fread(fid,1,'uchar');
h.xtics             = fread(fid,1,'char');
h.xrange            = fread(fid,1,'char');
h.ytics             = fread(fid,1,'char');
h.yrange            = fread(fid,1,'char');
h.xscalevalue       = fread(fid,1,'float');
h.xscaleinterval    = fread(fid,1,'float');
h.yscalevalue       = fread(fid,1,'float');
h.yscaleinterval    = fread(fid,1,'float');
h.scaletoolx1       = fread(fid,1,'float');
h.scaletooly1       = fread(fid,1,'float');
h.scaletoolx2       = fread(fid,1,'float');
h.scaletooly2       = fread(fid,1,'float');
h.port              = fread(fid,1,'short');
h.numsamples        = fread(fid,1,'long');
h.filterflag        = fread(fid,1,'char');
h.lowcutoff         = fread(fid,1,'float');
h.lowpoles          = fread(fid,1,'short');
h.highcutoff        = fread(fid,1,'float');
h.highpoles         = fread(fid,1,'short');
h.filtertype        = fread(fid,1,'char');
h.filterdomain      = fread(fid,1,'char');
h.snrflag           = fread(fid,1,'char');
h.coherenceflag     = fread(fid,1,'char');
h.continuoustype    = fread(fid,1,'char');
h.eventtablepos     = fread(fid,1,'long');
h.continuousseconds = fread(fid,1,'float');
h.channeloffset     = fread(fid,1,'long');
h.autocorrectflag   = fread(fid,1,'char');
h.dcthreshold       = fread(fid,1,'uchar');

for n = 1:h.nchannels
e(n).lab            = fread(fid,10,'char');
e(n).reference      = fread(fid,1,'char');
e(n).skip           = fread(fid,1,'char');
e(n).reject         = fread(fid,1,'char');
e(n).display        = fread(fid,1,'char');
e(n).bad            = fread(fid,1,'char');
e(n).n              = fread(fid,1,'ushort');
e(n).avg_reference  = fread(fid,1,'char');
e(n).clipadd        = fread(fid,1,'char');
e(n).x_coord        = fread(fid,1,'float');
e(n).y_coord        = fread(fid,1,'float');
e(n).veog_wt        = fread(fid,1,'float');
e(n).veog_std       = fread(fid,1,'float');
e(n).snr            = fread(fid,1,'float');
e(n).heog_wt        = fread(fid,1,'float');
e(n).heog_std       = fread(fid,1,'float');
e(n).baseline       = fread(fid,1,'short');
e(n).filtered       = fread(fid,1,'char');
e(n).fsp            = fread(fid,1,'char');
e(n).aux1_wt        = fread(fid,1,'float');
e(n).aux1_std       = fread(fid,1,'float');
e(n).sensitivity    = fread(fid,1,'float');
e(n).gain           = fread(fid,1,'char');
e(n).hipass         = fread(fid,1,'char');
e(n).lopass         = fread(fid,1,'char');
e(n).page           = fread(fid,1,'uchar');
e(n).size           = fread(fid,1,'uchar');
e(n).impedance      = fread(fid,1,'uchar');
e(n).physicalchnl   = fread(fid,1,'uchar');
e(n).rectify        = fread(fid,1,'char');
e(n).calib          = fread(fid,1,'float');
end

for i = 1:h.nchannels
   d(i).header      = fread(fid,5,'char');
   d(i).samples     = fread(fid,h.pnts,'float');
end

for j = 1:h.nchannels
   v(j).samples     = fread(fid,h.pnts,'float');
end

t = fread(fid,'char');

f.header = h;
f.electloc = e;
f.data = d;
f.variance = v;
f.tag = t;

frewind(fid);
fclose(fid);

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return

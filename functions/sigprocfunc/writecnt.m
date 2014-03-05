% writecnt() - Write a Neuroscan continuous signal file.
%
% Usage:
%   >> writecnt(filename, CNT-dataset, varargin) 
%
% Inputs:
%   filename - name of the file with extension
%   dataset -  name of the CNT-dataset, a structure with the following fields
%                   cntdataset.header
%                   cntdataset.electloc
%                   cntdataset.data
%                   cntdataset.Teeg
%                   cntdataset.event
%                   cntdataset.tag
%                   cntdataset.endtag
%                   
%                   optional fields
%                   cntdataset.dataformat:  'int32' or 'int16'
%
% Optional inputs:
%   'header':       bool, write header. (default=true)
%   'electrodes':   bool, write electrode information. (default=true)
%   'data':         bool, write data (default=true)
%   'eventtable':   bool, write event table (default=true)
%   'endtag':       bool, write end tag (default=true)
%   'append':   bool, append requested information to end of file.
%               (default=false) Requires that file exists. 
%  't1'         - start at time t1, default 0; ERROR WITH NONZERO VALUES
%  'sample1'    - start at sample1, default 0, overrides t1 ERROR WITH NONZERO VALUES
%  'lddur'      - duration of segment to load, default = whole file
%  'ldnsamples' - number of samples to load, default = whole file, 
%                 overrides lddur
%  'scale'      - ['on'|'off'] scale data to microvolt (default:'on')
%  'dataformat' - ['int16'|'int32'] default is 'int16' for 16-bit data.
%                 Use 'int32' for 32-bit data.
%   
% Outputs:
%  file         - file containing NeuroScan CNT file with the continuous
%                 data and other information
%
% Authors:   Sean Fitzgibbon, Arnaud Delorme, Michiel Vestjens and
%            and Chris Bishop 2000-2014. Maintained by Chris Bishop
%            (cwbishop_at_ucdavis.edu).
%
% Known limitations: 
%  For more see http://www.cnl.salk.edu/~arno/cntload/index.html    

% Copyright (C) 2000 Sean Fitzgibbon, <psspf@id.psy.flinders.edu.au>
% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function f = writecnt(filename,cntdataset,varargin)

if ~isempty(varargin)
	 WriteOptions=struct(varargin{:});
else WriteOptions = []; 
end;

try, WriteOptions.t1;         catch, WriteOptions.t1=0; end
try, WriteOptions.sample1;    catch, WriteOptions.sample1=[]; end
try, WriteOptions.lddur;      catch, WriteOptions.lddur=[]; end
try, WriteOptions.ldnsamples; catch, WriteOptions.ldnsamples=[]; end
try, WriteOptions.scale;      catch, WriteOptions.scale='on'; end
try, WriteOptions.dataformat;  catch, WriteOptions.dataformat='int16'; end

%% CATCH FOR KNOWN ERRORS WITH NON ZERO T1 and SAMPLE1
if WriteOptions.t1~=0 || (~isempty(WriteOptions.sample1) && WriteOptions.sample1~=0)
    error('writecnt:NonzeroStartPosition', 'writecnt cannot deal with nonzero write positions.');
end % if WriteOptions

%% ADDITIONAL CHECKS BY CWB
%   Set defaults for optional inputs
if ~isfield(WriteOptions, 'header') || isempty(WriteOptions.header), WriteOptions.header=true; end
if ~isfield(WriteOptions, 'electrodes') || isempty(WriteOptions.electrodes), WriteOptions.electrodes=true; end
if ~isfield(WriteOptions, 'data') || isempty(WriteOptions.data), WriteOptions.data=true; end
if ~isfield(WriteOptions, 'eventtable') || isempty(WriteOptions.eventtable), WriteOptions.eventtable=true; end
if ~isfield(WriteOptions, 'append') || isempty(WriteOptions.append), WriteOptions.append=false; end
if ~isfield(WriteOptions, 'endtag') || isempty(WriteOptions.endtag), WriteOptions.endtag=true; end

%% SIZE OF EVENT DATA (differs based on version of SCAN)
sizeEvent1 = 8  ; %%% 8  bytes for Event1  
sizeEvent2 = 19 ; %%% 19 bytes for Event2 

type='cnt';
if nargin ==1 
    scan=0;
end     

h   = cntdataset.header;
e   = cntdataset.electloc;
dat = cntdataset.data;
eT  = cntdataset.Teeg;
ev2 = cntdataset.event;
t   = cntdataset.tag;
if ~isfield(cntdataset,'endtag')
    endtag=[];
else
    endtag=cntdataset.endtag;
end

%% APPEND DATA?
%   If not, then open data for general writing
if WriteOptions.append
    fid=fopen(filename,'a+');
%     display('Appending requested information ...');
else
    fid = fopen(filename,'w');
end % if WriteOptions.append


% disp(['Writing file ' filename ' ...'])

% HEADER : 900 bytes => Starts at 0h, finishes with 383h.
if WriteOptions.header
    fwrite(fid,h.rev,'char');
    fwrite(fid,h.nextfile,'long');
    fwrite(fid,h.prevfile,'long');
    fwrite(fid,h.type,'char');
    fwrite(fid,h.id,'char');
    fwrite(fid,h.oper,'char');
    fwrite(fid,h.doctor,'char');
    fwrite(fid,h.referral,'char');
    fwrite(fid,h.hospital,'char');
    fwrite(fid,h.patient,'char');
    fwrite(fid,h.age,'short');
    fwrite(fid,h.sex,'char');
    fwrite(fid,h.hand,'char');
    fwrite(fid,h.med,'char');
    fwrite(fid,h.category,'char');
    fwrite(fid,h.state,'char');
    fwrite(fid,h.label,'char');
    fwrite(fid,h.date,'char');
    fwrite(fid,h.time,'char');
    fwrite(fid,h.mean_age,'float');
    fwrite(fid,h.stdev,'float');
    fwrite(fid,h.n,'short');
    fwrite(fid,h.compfile,'char');
    fwrite(fid,h.spectwincomp,'float');
    fwrite(fid,h.meanaccuracy,'float');
    fwrite(fid,h.meanlatency,'float');
    fwrite(fid,h.sortfile,'char');
    fwrite(fid,h.numevents,'int');
    fwrite(fid,h.compoper,'char');
    fwrite(fid,h.avgmode,'char');
    fwrite(fid,h.review,'char');
    fwrite(fid,h.nsweeps,'ushort');
    fwrite(fid,h.compsweeps,'ushort');
    fwrite(fid,h.acceptcnt,'ushort');
    fwrite(fid,h.rejectcnt,'ushort');
    fwrite(fid,h.pnts,'ushort');
    fwrite(fid,h.nchannels,'ushort');
    fwrite(fid,h.avgupdate,'ushort');
    fwrite(fid,h.domain,'char');
    fwrite(fid,h.variance,'char');
    fwrite(fid,h.rate,'ushort');
    fwrite(fid,h.scale,'double');
    fwrite(fid,h.veogcorrect,'char');
    fwrite(fid,h.heogcorrect,'char');
    fwrite(fid,h.aux1correct,'char');
    fwrite(fid,h.aux2correct,'char');
    fwrite(fid,h.veogtrig,'float');
    fwrite(fid,h.heogtrig,'float');
    fwrite(fid,h.aux1trig,'float');
    fwrite(fid,h.aux2trig,'float');
    fwrite(fid,h.heogchnl,'short');
    fwrite(fid,h.veogchnl,'short');
    fwrite(fid,h.aux1chnl,'short');
    fwrite(fid,h.aux2chnl,'short');
    fwrite(fid,h.veogdir,'char');
    fwrite(fid,h.heogdir,'char');
    fwrite(fid,h.aux1dir,'char');
    fwrite(fid,h.aux2dir,'char');
    fwrite(fid,h.veog_n,'short');
    fwrite(fid,h.heog_n,'short');
    fwrite(fid,h.aux1_n,'short');
    fwrite(fid,h.aux2_n,'short');
    fwrite(fid,h.veogmaxcnt,'short');
    fwrite(fid,h.heogmaxcnt,'short');
    fwrite(fid,h.aux1maxcnt,'short');
    fwrite(fid,h.aux2maxcnt,'short');
    fwrite(fid,h.veogmethod,'char');
    fwrite(fid,h.heogmethod,'char');
    fwrite(fid,h.aux1method,'char');
    fwrite(fid,h.aux2method,'char');
    fwrite(fid,h.ampsensitivity,'float');
    fwrite(fid,h.lowpass,'char');
    fwrite(fid,h.highpass,'char');
    fwrite(fid,h.notch,'char');
    fwrite(fid,h.autoclipadd,'char');
    fwrite(fid,h.baseline,'char');
    fwrite(fid,h.offstart,'float');
    fwrite(fid,h.offstop,'float');
    fwrite(fid,h.reject,'char');
    fwrite(fid,h.rejstart,'float');
    fwrite(fid,h.rejstop,'float');
    fwrite(fid,h.rejmin,'float');
    fwrite(fid,h.rejmax,'float');
    fwrite(fid,h.trigtype,'char');
    fwrite(fid,h.trigval,'float');
    fwrite(fid,h.trigchnl,'char');
    fwrite(fid,h.trigmask,'short');
    fwrite(fid,h.trigisi,'float');
    fwrite(fid,h.trigmin,'float');
    fwrite(fid,h.trigmax,'float');
    fwrite(fid,h.trigdir,'char');
    fwrite(fid,h.autoscale,'char');
    fwrite(fid,h.n2,'short');
    fwrite(fid,h.dir,'char');
    fwrite(fid,h.dispmin,'float');
    fwrite(fid,h.dispmax,'float');
    fwrite(fid,h.xmin,'float');
    fwrite(fid,h.xmax,'float');
    fwrite(fid,h.automin,'float');
    fwrite(fid,h.automax,'float');
    fwrite(fid,h.zmin,'float');
    fwrite(fid,h.zmax,'float');
    fwrite(fid,h.lowcut,'float');
    fwrite(fid,h.highcut,'float');
    fwrite(fid,h.common,'char');
    fwrite(fid,h.savemode,'char');
    fwrite(fid,h.manmode,'char');
    fwrite(fid,h.ref,'char');
    fwrite(fid,h.rectify,'char');
    fwrite(fid,h.displayxmin,'float');
    fwrite(fid,h.displayxmax,'float');
    fwrite(fid,h.phase,'char');
    fwrite(fid,h.screen,'char');
    fwrite(fid,h.calmode,'short');
    fwrite(fid,h.calmethod,'short');
    fwrite(fid,h.calupdate,'short');
    fwrite(fid,h.calbaseline,'short');
    fwrite(fid,h.calsweeps,'short');
    fwrite(fid,h.calattenuator,'float');
    fwrite(fid,h.calpulsevolt,'float');
    fwrite(fid,h.calpulsestart,'float');
    fwrite(fid,h.calpulsestop,'float');
    fwrite(fid,h.calfreq,'float');
    fwrite(fid,h.taskfile,'char');
    fwrite(fid,h.seqfile,'char');
    fwrite(fid,h.spectmethod,'char');
    fwrite(fid,h.spectscaling,'char');
    fwrite(fid,h.spectwindow,'char');
    fwrite(fid,h.spectwinlength,'float');
    fwrite(fid,h.spectorder,'char');
    fwrite(fid,h.notchfilter,'char');
    fwrite(fid,h.headgain,'short');
    fwrite(fid,h.additionalfiles,'int');
    fwrite(fid,h.unused,'char');
    fwrite(fid,h.fspstopmethod,'short');
    fwrite(fid,h.fspstopmode,'short');
    fwrite(fid,h.fspfvalue,'float');
    fwrite(fid,h.fsppoint,'short');
    fwrite(fid,h.fspblocksize,'short');
    fwrite(fid,h.fspp1,'ushort');
    fwrite(fid,h.fspp2,'ushort');
    fwrite(fid,h.fspalpha,'float');
    fwrite(fid,h.fspnoise,'float');
    fwrite(fid,h.fspv1,'short');
    fwrite(fid,h.montage,'char');
    fwrite(fid,h.eventfile,'char');
    fwrite(fid,h.fratio,'float');
    fwrite(fid,h.minor_rev,'char');
    fwrite(fid,h.eegupdate,'short');
    fwrite(fid,h.compressed,'char');
    fwrite(fid,h.xscale,'float');
    fwrite(fid,h.yscale,'float');
    fwrite(fid,h.xsize,'float');
    fwrite(fid,h.ysize,'float');
    fwrite(fid,h.acmode,'char');
    fwrite(fid,h.commonchnl,'uchar');
    fwrite(fid,h.xtics,'char');
    fwrite(fid,h.xrange,'char');
    fwrite(fid,h.ytics,'char');
    fwrite(fid,h.yrange,'char');
    fwrite(fid,h.xscalevalue,'float');
    fwrite(fid,h.xscaleinterval,'float');
    fwrite(fid,h.yscalevalue,'float');
    fwrite(fid,h.yscaleinterval,'float');
    fwrite(fid,h.scaletoolx1,'float');
    fwrite(fid,h.scaletooly1,'float');
    fwrite(fid,h.scaletoolx2,'float');
    fwrite(fid,h.scaletooly2,'float');
    fwrite(fid,h.port,'short');
    fwrite(fid,h.numsamples,'ulong');
    fwrite(fid,h.filterflag,'char');
    fwrite(fid,h.lowcutoff,'float');
    fwrite(fid,h.lowpoles,'short');
    fwrite(fid,h.highcutoff,'float');
    fwrite(fid,h.highpoles,'short');
    fwrite(fid,h.filtertype,'char');
    fwrite(fid,h.filterdomain,'char');
    fwrite(fid,h.snrflag,'char');
    fwrite(fid,h.coherenceflag,'char');
    fwrite(fid,h.continuoustype,'char');
    fwrite(fid,h.eventtablepos,'long');
    fwrite(fid,h.continuousseconds,'float');
    fwrite(fid,h.channeloffset,'long');
    fwrite(fid,h.autocorrectflag,'char');
    fwrite(fid,h.dcthreshold,'uchar');  % = 383H
end % if WriteOptions.header

% ELECT.DESCRIPTIONS : 75*n.channels bytes.
% Starts with 384h, finishes with (899+75*nchannels)dec
% 10 channels: 671h   % 12 channels: 707h
% 24 channels: A8Bh   % 32 channels: CE3h
% 64 channels: 1643h  % 128 channels: 2903h

if WriteOptions.electrodes
    for n = 1:h.nchannels
        lablength = fwrite(fid,e(n).lab,'char');
        fwrite(fid,e(n).reference,'char',10-lablength);
        fwrite(fid,e(n).skip,'char');
        fwrite(fid,e(n).reject,'char');
        fwrite(fid,e(n).display,'char');
        fwrite(fid,e(n).bad,'char');
        fwrite(fid,e(n).n,'ushort');
        fwrite(fid,e(n).avg_reference,'char');
        fwrite(fid,e(n).clipadd,'char');
        fwrite(fid,e(n).x_coord,'float');
        fwrite(fid,e(n).y_coord,'float');
        fwrite(fid,e(n).veog_wt,'float');
        fwrite(fid,e(n).veog_std,'float');
        fwrite(fid,e(n).snr,'float');
        fwrite(fid,e(n).heog_wt,'float');
        fwrite(fid,e(n).heog_std,'float');
        fwrite(fid,e(n).baseline,'short');
        fwrite(fid,e(n).filtered,'char');
        fwrite(fid,e(n).fsp,'char');
        fwrite(fid,e(n).aux1_wt,'float');
        fwrite(fid,e(n).aux1_std,'float');
        fwrite(fid,e(n).senstivity,'float');
        fwrite(fid,e(n).gain,'char');
        fwrite(fid,e(n).hipass,'char');
        fwrite(fid,e(n).lopass,'char');
        fwrite(fid,e(n).page,'uchar');
        fwrite(fid,e(n).size,'uchar');
        fwrite(fid,e(n).impedance,'uchar');
        fwrite(fid,e(n).physicalchnl,'uchar');
        fwrite(fid,e(n).rectify,'char');
        fwrite(fid,e(n).calib,'float');
    end % for 
end % if WriteOptions.electrodes

%% SET FILE POINTER TO END OF HEADER INFORMATION
%   Need to advance the pointer in the event that the header and/or
%   electrodes have not been written (or have already been written to
%   file) since the file pointer position is used to determine the
%   beginning of the data. 
%
%   This is a bit clunky, but ensures that 'append' mode works properly
%   when calculated the number of data points to write.
fseek(fid, 900+75*h.nchannels, 'bof'); 
% finding if 32-bits of 16-bits file
% ----------------------------------
begdata = ftell(fid);
enddata = h.eventtablepos;   % after data
if strcmpi(WriteOptions.dataformat, 'int16')
     nums    = (enddata-begdata)/h.nchannels/2;
else nums    = (enddata-begdata)/h.nchannels/4;
end;

% number of sample to write
% -------------------------
if ~isempty(WriteOptions.sample1)
   WriteOptions.t1      = WriteOptions.sample1/h.rate;
else 
   WriteOptions.sample1 = WriteOptions.t1*h.rate;
end;
if strcmpi(WriteOptions.dataformat, 'int16')
     startpos = WriteOptions.t1*h.rate*2*h.nchannels;
else startpos = WriteOptions.t1*h.rate*4*h.nchannels;
end;
if isempty(WriteOptions.ldnsamples)
     if ~isempty(WriteOptions.lddur)
          WriteOptions.ldnsamples = round(WriteOptions.lddur*h.rate); 
     else WriteOptions.ldnsamples = nums;
     end;
end;

% scaling data from microvolts
% ----------------------------
if strcmpi(WriteOptions.scale, 'on')
%     disp('Scaling data .....')
    for i=1:h.nchannels
       bas=e(i).baseline;
       sen=e(i).senstivity;
       cal=e(i).calib;
       mf=sen*(cal/204.8);
       dat(i,:)=(dat(i,:)/mf)+bas;
   end
end


% write data
% ----------

% disp('Writing data .....')
if type == 'cnt' 
    
    if WriteOptions.data
        channel_off = h.channeloffset/2;
  
        fseek(fid, startpos, 0);
        if channel_off <= 1
              for temploop =1:WriteOptions.ldnsamples;
                  fwrite(fid, dat(1:h.nchannels, temploop), WriteOptions.dataformat)';
            end;
        else
              for temploop =1:h.nchannels;
                  fwrite(fid, dat(temploop, 1:channel_off), WriteOptions.dataformat)';
            end;
          
            counter = 1;	
            while counter*channel_off < WriteOptions.ldnsamples

                  for temploop =1:h.nchannels;
                    fwrite(fid, dat(temploop, counter*channel_off+1:counter*channel_off+channel_off), WriteOptions.dataformat)';
                end;
                counter = counter + 1;
            end;
        end;	
    end % if WriteOptions.data
      
    if WriteOptions.eventtable
        % write event table
        % -----------------
        frewind(fid);
        fseek(fid,h.eventtablepos,'bof');

        disp('Writing Event Table...')
        fwrite(fid,eT.teeg,'uchar');
        fwrite(fid,eT.size,'ulong');
        fwrite(fid,eT.offset,'ulong');
    
        if eT.teeg==2
            nevents=eT.size/sizeEvent2;
        elseif eT.teeg==1
            nevents=eT.size/sizeEvent1;
        else
            disp('No !!! teeg <> 2 and teeg <> 1');
            ev2 = [];
        end % if eT.teeg==2
  
        %%%% to change offset in points back to bytes 
        if ~isempty(ev2)
            ev2p=ev2; 
            ioff=900+(h.nchannels*75); %% initial offset : header + electordes desc 
            
            % 140219 CWB: Writing events to file in byte form depends on
            % data precision. Need to select the correct data precision
            % here.
            if strcmpi(WriteOptions.dataformat, 'int16')
                NBYTES=2;
            elseif strcmpi(WriteOptions.dataformat, 'int32')
                NBYTES=4;
            end % if strcmpi
            
            for i=1:nevents 
                ev2p(i).offset=((ev2p(i).offset + WriteOptions.sample1)*NBYTES*h.nchannels) +ioff; %% 2 short int end 
            end     
            ev2 = ev2p;
        end; % if ~isempty(ev2)

        if eT.teeg==2
            nevents=eT.size/sizeEvent2;
            for i=1:nevents
                ev2(i).stimtype      = fwrite(fid,ev2(i).stimtype,'ushort');
                ev2(i).keyboard      = fwrite(fid,ev2(i).keyboard,'char');
                ev2(i).keypad_accept = fwrite(fid,ev2(i).keypad_accept,'char');
                ev2(i).offset        = fwrite(fid,ev2(i).offset,'long');
                ev2(i).type          = fwrite(fid,ev2(i).type,'short'); 
                ev2(i).code          = fwrite(fid,ev2(i).code,'short');
                ev2(i).latency       = fwrite(fid,ev2(i).latency,'float');
                ev2(i).epochevent    = fwrite(fid,ev2(i).epochevent,'char');
                ev2(i).accept        = fwrite(fid,ev2(i).accept,'char');
                ev2(i).accuracy      = fwrite(fid,ev2(i).accuracy,'char');
            end % for i=1:nevents
        elseif eT.teeg==1
            nevents=eT.size/sizeEvent1;
            for i=1:nevents
                ev2(i).stimtype      = fwrite(fid,ev2(i).stimtype,'ushort');
                ev2(i).keyboard      = fwrite(fid,ev2(i).keyboard,'char');
                ev2(i).keypad_accept = fwrite(fid,ev2(i).keypad_accept,'char');
                ev2(i).offset        = fwrite(fid,ev2(i).offset,'long');
            end % for i=1:nevents
        else
            disp('No !!! teeg <> 2 and teeg <> 1');
            ev2 = [];
        end %if/elseif/else
    end % if WriteOptions.eventtable
end % if type=='cnt'

%% WE ARE AT h.nextfile position here
%   So there has to be additional information at the end of the CNT file
%   that is NOT read in by loadcnt.m but we ultimately need in order to
%   rewrite the file. 
%
%   Ah, there IS more information that loadcnt is not actually reading,
%   although it depends on some of it to determine if the data are 16 or 32
%   bit. 
%
%   CWB put this information in a junk field (cntdataset.junk) that can be
%   written if it's present. 
if WriteOptions.endtag
    if size(endtag,1)>1
        fwrite(fid,endtag);
    else
        fwrite(fid,t,'char');
    end
end % if WriteOptions.endtag

fclose(fid);
% disp(['Finished writing file ' filename ' ...'])

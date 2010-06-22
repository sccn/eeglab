% read_erpss() - read an compressed and uncompressed ERPSS file formats 
%                (.RAW or .RDF) 
%                
% Usage: 
%         >> [data,events,header] = read_erpss(filename);
% Inputs:
%   filename - Name of ERPSS data file (including extension) 
%
% Outputs:
%   data     - Data array [nchans samples]. Data is rescaled to microVolt.
%   events   - Event information structure:
%              events.sample_offset[]: Event offsets in samples
%                                      from the first sample (0)
%              events.event_code[]     Event codes (integers: 1-128)
%   header - header information structure:
%            header.nchans   Number of channels
%            header.nframes  Number of data frames (i.e., samples, timepoints)
%
% Notes: ERPSS format was developed by Jonathan Hansen at the Hillyard lab 
%              at UCSD (http://sdepl.ucsd.edu/erpss/).
%
% Authors: Jeng-Ren Duann & Arnaud Delorme, CNL/Salk & INC/UCSD, 2002-12-12
%          with help from Andrey Vankov

% $Log: read_erpss.m,v $
% Revision 1.29  2007/08/07 20:01:08  arno
% unused variables
%
% Revision 1.28  2007/08/07 19:48:28  arno
% fix help (bug 382)
%
% Revision 1.27  2006/12/15 20:09:32  toby
% Reverted to version 1.25. Andre reports 1.26 not working
%
% Revision 1.25  2005/05/20 18:39:53  hilit
% if header is empty not rescaling to uV.
%
% Revision 1.24  2005/04/08 23:23:41  arno
% rescale to microvolt
%
% Revision 1.23  2005/03/11 00:07:03  arno
% message
%
% Revision 1.22  2005/03/10 23:59:37  arno
% message
%
% Revision 1.21  2005/03/10 23:47:15  arno
% edit message
%
% Revision 1.20  2003/12/10 17:00:03  arno
% comment
%
% Revision 1.19  2003/11/19 03:05:22  arno
% debuging calibration
%
% Revision 1.18  2003/06/28 00:28:38  arno
% (implmenting new header info
% and rescaling to mircouV
%
% Revision 1.19  2010/06/22 00:28:38  Christina Karns & Bill Troyer
% Old version assumed little endian
% Current version checks endian first 
% and loads as big endian if needed


function [eeg,ev,header] = read_erpss(filename)
    
    if nargin < 1
      help read_erpss;
      return;
    end;
    
    eeg = [];
    ev = [];
    header = [];

    
    fp = fopen(filename,'rb','ieee-le');
    if fp == -1,
        error('read_erpss(): Cannot open data file...!');
    else
        disp('File opened:');
    end
    
    %Check for endian-ness
    atag = fread(fp,1,'uint32');
    if atag == hex2dec('b0aa55')
        fprintf('Using Little-Endian byte format');
        tmptag = hex2dec('f0aa55');
        endian = 'ieee-le';
    else 
        fprintf('Using Big-Endian byte format');
        tmptag = hex2dec('aa5500f0');
        fclose(fp);
        endian = 'ieee-be';
        fp = fopen(filename,'rb',endian);
    end
    
    fseek(fp,4,-1);
    compressed = fread(fp,1,'uint16');
    if compressed, 
        disp('Reading compressed data: if an error occurs while calling the MEX'); 
        disp('function, try recompiling it for your platform "mex decompresserpss.cc"'); 
    end;
    fseek(fp,6,-1);
    header.nchans = fread(fp,1,'uint16');
    
    cnt = 0;
    ev_cnt = 0;
    firstpass = 1;
    
    % first pass, scan data
    totalsize = 0;
    totblocks = 0;
    disp('finding total number of blocks ...');
    while(~feof(fp)),
        tag = fread(fp,1,'uint32');
        if length(tag) == 0,
            break;
        end
        if tag == tmptag,
            cnt = cnt + 1;

            % Read nchans and block length
            fseek(fp,2,0);
            nchans = fread(fp,1,'uint16');
            fread(fp,1,'uint16');
            block_size = fread(fp,1,'uint16');

            % Read events
            if ~firstpass
                if compressed
                    fseek(fp,10,0);
                    block_size_compress = fread(fp,1,'uint16');
                    fseek(fp,48+110*4+block_size_compress*2,0);
                else
                    fseek(fp,60+110*4+nchans*block_size*2,0);
                end;
            else
                if compressed
                    fseek(fp,10,0);
                    block_size_compress = double(fread(fp,1,'uint16'));
                    header.amplif              = double(fread(fp,1,'uint32'));
                    header.clock_freq          = double(fread(fp,1,'uint32'));
                    header.divider             = double(fread(fp,1,'uint32'));
                    header.ad_range_mv         = double(fread(fp,1,'uint32'));
                    header.ad_bits             = double(fread(fp,1,'uint32'));
                    header.nsteps              = double(fread(fp,1,'uint32'));
                    fseek(fp,24+110*4+block_size_compress*2,0);
                else
                    fseek(fp,12,0);
                    header.amplif              = double(fread(fp,1,'uint32'));
                    header.clock_freq          = double(fread(fp,1,'uint32'));
                    header.divider             = double(fread(fp,1,'uint32'));
                    header.ad_range_mv         = double(fread(fp,1,'uint32'));
                    header.ad_bits             = double(fread(fp,1,'uint32'));
                    header.nsteps              = double(fread(fp,1,'uint32'));
                    fseek(fp,24+110*4+nchans*block_size*2,0);
                end;
                
                if (header.ad_range_mv ~= 0) & (header.amplif ~= 0) | (header.ad_bits ~= 0)
                    header.rescaleuv = header.ad_range_mv*1000 / header.amplif / pow2(header.ad_bits);
                end                
                if header.divider & header.clock_freq
                    header.srate     = header.clock_freq/header.divider/header.nsteps;
                    if round(header.srate)
                        fprintf('Sampling rate is %4.4fHz\n', header.srate);
                    else
                        fprintf('Unknown sampling rate.\n');
                    end;
                end;
                firstpass = 0;
            end; 
            totalsize = totalsize + block_size;
            
        %else
        %    fprintf('.');
        end
        totblocks = totblocks+1;
    end
    eeg = zeros(header.nchans, totalsize);
    
    % second pass, read data
    disp(['Reading blocks (out of ' num2str(cnt) '):']);
    cnt = 0;
    totalsize = 0;
    fclose(fp);
    fp = fopen(filename,'rb',endian);
    fseek(fp,552,-1);
    temp = fread(fp,1,'uint16');
    if ~isfield(header, 'srate')
        header.srate = 1000000.0/temp;
    end;
    fseek(fp,6,-1);
    header.nchans = fread(fp,1,'uint16');

    while(~feof(fp)),
        tag = fread(fp,1,'uint32');
        if length(tag) == 0,
            break;
        end
        if tag == tmptag,
            cnt = cnt + 1;
            if ~mod(cnt,100)
                fprintf('%d ', cnt);
            end;
            if ~mod(cnt,1000)
                fprintf('\n');
            end;
            
            % Read nchans and block length
            fseek(fp,2,0);
            nchans = fread(fp,1,'uint16');
            fread(fp,1,'uint16');
            block_size = fread(fp,1,'uint16');
            header.ndupsamp   = fread(fp,1,'uint16');
            header.nrun       = fread(fp,1,'uint16');
            header.err_detect = fread(fp,1,'uint16');
            header.nlost      = fread(fp,1,'uint16');
            nevents    = fread(fp,1,'uint16');
            block_size_compress = fread(fp,1,'uint16');
            %if cnt == 3, return; end;
            
            % Read events
            fseek(fp,48,0);
            for i=1:nevents,
                samp_off  = fread(fp,1,'uint8');
                cond_code = fread(fp,1,'uint8'); 
                ev_code   = fread(fp,1,'uint16');
                ev_cnt    = ev_cnt + 1;
                ev(ev_cnt).sample_offset = samp_off + (cnt-1)*block_size+1; %+1 for Matlab 
                ev(ev_cnt).event_code    = ev_code;
                ev(ev_cnt).cond_code     = cond_code;
            end
            fseek(fp,4*(110-nevents),0);
            if compressed
                data = uint8(fread(fp,2*block_size_compress,'uchar'));
                try
                    data = decompresserpss(data, nchans, block_size, 0); % 0 is little endian, 1 is big endian
                catch,
                    disp('Error executing the decompresserpss function');
                    disp('Likely this function has not been compiled for your system');
                    disp('Under Matlab, go to eeglabXXXX/functions/sigprocfunc/');
                    disp('and type mex decompresserpss.cc');
                    error('decompresserpss function error (see message above)');
                end;
            else 
                data = fread(fp,nchans*block_size,'int16');
            end;
            try, 
                eeg(:,totalsize+1:totalsize+length(data)/nchans) = reshape(data,nchans,length(data)/nchans); % concatenate data blocks
                totalsize = totalsize + length(data)/nchans;
            catch,
                fprintf('\nWarning: block %d truncated, skipped\n', cnt);
            end;
        end
    end
    fprintf('\n');
    
    % rescale to  uv
    % --------------
    if isfield(header, 'rescaleuv') & header.rescaleuv ~= 0
        disp('Rescaling data to microVolt');
        eeg = eeg*header.rescaleuv;
    end;
    
    fclose(fp);
    header.nframes = size(eeg,2);
        

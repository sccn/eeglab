% read_erpss() - read an compressed and uncompressed ERPSS file formats 
%                (.RAW or .RDF) 
%                
% Usage: 
%         >> [data,events,datasize] = read_erpss(filename);
% Inputs:
%   filename - Name of uncompressed ERPSS data file (including extension) 
%
% Outputs:
%   data     - Data array [nchans samples]
%   events   - Event information structure:
%              events.sample_offset[]: Event offsets in samples
%                                      from the first sample (0)
%              events.event_code[]     Event codes (integers: 1-128)
%   datasize - Data size information structure:
%              datasize.nchans   Number of channels
%              datasize.nframes  Number of data frames (i.e., samples, timepoints)
%
% Notes: ERPSS format was developed by Jonathan Hansen at the Hillyard lab 
%              at UCSD (http://sdepl.ucsd.edu/erpss/).
%
% Authors: Jeng-Ren Duann, CNL/Salk & INC/UCSD, 2002-12-12
%          Arnaud Delorme, debug, spped-up and read compressed data 2003-06
%          with help from Andrey Vankov

function [eeg,ev,header] = read_erpss(filename)
  
    eeg = [];
    ev = [];
    header = [];
    tmptag = hex2dec('f0aa55');
    
    fp = fopen(filename,'rb','ieee-le');
    if fp == -1,
        error('read_erpss(): Cannot open data file...!');
    else
        disp('File opened:');
    end
    
    fseek(fp,4,-1);
    compressed = fread(fp,1,'uint16');
    if compressed, 
        disp('Reading compressed data'); 
        complendian = 0;
        switch computer
         case {'MAC2','SUN4','SOL2','SGI','SGI64'}, complendian = 1;
         case {'PCWIN','LNX86','GLNX86'}, complendian = 0;
         otherwise 
          disp('Disp: endian of computer not known, set to PC (windows/linux) by default');
        end;
    end;
    fseek(fp,6,-1);
    header.nchans = fread(fp,1,'uint16');
    
    cnt = 0;
    ev_cnt = 0;
    
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
            if compressed
                fseek(fp,10,0);
                block_size_compress = fread(fp,1,'uint16');
                fseek(fp,48+110*4+block_size_compress*2,0);
            else
                fseek(fp,60+110*4+nchans*block_size*2,0);
            end;
            totalsize = totalsize + block_size;
        end
        totblocks = totblocks+1;
    end
    eeg = zeros(header.nchans, totalsize);
    
    % second pass, read data
    disp(['Reading blocks (out of ' num2str(cnt) '):']);
    cnt = 0;
    totalsize = 0;
    fclose(fp);
    fp = fopen(filename,'rb','ieee-le');
    fseek(fp,552,-1);
    temp = fread(fp,1,'uint16');
    header.srate = 1000000.0/temp;
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
            ndupsamp = fread(fp,1,'uint16');
            nrun = fread(fp,1,'uint16');
            err_detect = fread(fp,1,'uint16');
            nlost = fread(fp,1,'uint16');
            nevents = fread(fp,1,'uint16');
            block_size_compress = fread(fp,1,'uint16');
	    
            % Read events
            fseek(fp,48,0);
            for i=1:nevents,
                samp_off = fread(fp,1,'uint8');
                cond_code = fread(fp,1,'uint8');
                ev_code = fread(fp,1,'uint16');
                ev_cnt = ev_cnt + 1;
                ev(ev_cnt).sample_offset = samp_off + (cnt-1)*block_size+1; %+1 for Matlab 
                ev(ev_cnt).event_code = ev_code;
            end
            fseek(fp,4*(110-nevents),0);
            if compressed
                data = uint8(fread(fp,2*block_size_compress,'uchar'));
                data = decompresserpss(data, nchans, block_size, complendian);
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
    
    fclose(fp);
    header.nframes = size(eeg,2);
        
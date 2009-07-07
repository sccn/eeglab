function [HDR]=asn1read2(arg1,offset)

% modified by AS
%  cleanup unneed parts: checkForText, checkEncapsulate 
%  add linear and recursive mode
%  return complete HDR-structure

if ischar(arg1)
        HDR.FILE.OPEN = 0; 
        HDR.FileName = arg1;
        HDR.FILE.FID = fopen(HDR.FileName,'r');
        fseek(HDR.FILE.FID,offset,'bof'); 	% skip preamble
else
        HDR = arg1;
end;    

fid=HDR.FILE.FID;
index=0;

level = 0; 
recursive = 1;  % switch between "recursive" and "linear" mode

if recursive,
        HDR=getItem(fid,level,recursive);
else
        K = 1; 
        while ~feof(fid),
                hdr=getItem(fid,level,recursive);
                
                HDR.fpos(K) = hdr.fpos;
                if any(size(hdr.tag)<1), return; end; 
                HDR.Primitive(K) = hdr.Primitive;
                HDR.tag(K) = hdr.tag;
                HDR.len(K) = hdr.len;
                K=K+1;
        end;
end;
if ~fid, 
        fclose(fid);
end; 
return;


%%%% Read tag and length%%%%
function [HDR]=getItem(fid,level,recursive);

HDR.fpos = ftell(fid);
tag = fread(fid,1,'uint8');
tagclass = bitand(tag,192)/64;             % get class 
HDR.Primitive = ~bitand(tag,32);

if bitand(tag,31)==31,
        t  = fread(fid,1,'uint8');
        ac = bitand(t,127);
        while bitand(t,128),
                t  = fread(fid,1,'uint8');
                ac = [ac*128 + bitand(t,127)];
        end;
        HDR.tag = [tag,ac];
else
        HDR.tag = tag;
end;
if (tagclass==2),
        HDR.tag=bitand(HDR.tag,31);
        HDR.class='Context Tag';
elseif(tagclass==0)
        HDR.class='Universal Tag';
elseif(tagclass==1)
        HDR.class='Application Tag';
elseif(tagclass==3)
        HDR.class='Private Tag';
end;


len = fread(fid,1,'uint8');
HDR.len = len;
if bitand(len,128),
        len1 = bitand(len,127);
        if len1==0,
                % unbestimmte form
                HDR.len = NaN; 
        else
                % Langform 
                len1 = fread(fid,len1,'uint8');
                HDR.len = sum(len1.*256.^[length(len1)-1:-1:0]');
        end
        % else % Kurzform
end;           
%fseek(fid,HDR.len,'cof'); 

if isempty(tag),tag=0;end;

if HDR.Primitive,
        %[ret1,HDR.buffer]=checkForText( fid, HDR.len);
        HDR.Value = getdata(fid, HDR.len);
elseif recursive, 
        pos = ftell(fid);
        K=0; val={};
        while (ftell(fid) < (HDR.len+pos))
                K=K+1; 
                val{K} = getItem(fid,level+1,recursive);
        end;
        HDR.Value = val; 
end;	
HDR.level=level; 



function [buffer]=getdata( fid, len)
buffer = [];
if(len <= 2 ),
        fseek(fid,len,'cof');
elseif (len < 4 ),
        [buffer,sampleLength] = fread(fid,[1,len]);
else
        [buffer,sampleLength] = fread(fid,[1,len]);
end;





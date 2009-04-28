function f = eeg_load_scan4_cnt_event(fid)

% eeg_load_scan4_cnt_event - Read Neuroscan 4.x Continuous EEG Events
% 
% eventTable = eeg_load_scan4_cnt_event(fid)
%               
%   fid         -  File identifier
%
%   eventTable  -  structure of event info
%
%   Note: Developed with Scan 4.1 CNT files
%
%   See also: eeg_load_scan4_cnt_data
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no express or implied warranties
% History:  2002, Sean.Fitzgibbon@flinders.edu.au
%           06/2002, Darren.Weber_at_radiology.ucsf.edu
%                    adapted to eeg_toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





fseek(fid,886,'bof');           eventPos = fread(fid,1,'long');


fseek(fid,eventPos,'bof');      f.type   = fread(fid,1,'uchar');

                                f.size   = fread(fid,1,'long');
                                f.offset = fread(fid,1,'long');
                                

fseek(fid,f.offset,'cof');


for i = 1:(f.size/19),
    

    f.event(i).stimType   = fread(fid,1,'short');
    

    f.event(i).keyBoard   = fread(fid,1,'char');
    

    f.event(i).keyPad     = fread(fid,1,'bit4');
    

    f.event(i).Accept     = fread(fid,1,'bit4');
    

    f.event(i).offset     = fread(fid,1,'long');
    

    f.event(i).type       = fread(fid,1,'short');
    

    f.event(i).code       = fread(fid,1,'short');
    

    f.event(i).latency    = fread(fid,1,'float');
    

    f.event(i).epochEvent = fread(fid,1,'char');
    

    f.event(i).accept     = fread(fid,1,'char');
    

    f.event(i).accuracy   = fread(fid,1,'char');
    

end

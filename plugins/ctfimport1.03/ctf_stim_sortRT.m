function ctf = sortrtctf(ctfds)

RESPTHRESH = 10^7;

if nargin < 1
    ctfds = '';
end

ctf = ctf_read(ctfds);

stim_chan = find(strcmp(ctf.sensor.label,'STIM'));

if (length(stim_chan) == 1)
    stim_data = ctf.data(:,stim_chan,:);
else
    error('ERROR: could not properly locate STIM channel');
end

RTs = [];
for i = 1:size(stim_data,3)
    tmp = stim_data(:,1,i);
    RTs(end+1) = min(find(tmp > RESPTHRESH));
end

[RTs_sorted,RTsx] = sort(RTs);

ctf.data = ctf.data(:,:,RTsx);


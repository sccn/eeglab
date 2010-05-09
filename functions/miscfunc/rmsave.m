% rmsave() - return the RMS in each channel, epoch
%
% Usage:
%         >> ave = rmsave(data,frames);

% Scott Makeig, CNL/Salk Institute, La Jolla, 9/98

function ave = rmsave(data,frames)

if nargin<1
  help rmsave
  return
end
if nargin<2
	frames = size(data,2); 
	data = reshape(data, size(data,1), size(data,2)*size(data,3));
end;

chans = size(data,1);
datalength = size(data,2);
if rem(datalength,frames)
   fprintf('frames should divide data length.\n');
   return
end
if frames < 1
   fprintf('frames should be > 1.\n');
   return
end

epochs = datalength/frames;
ave = zeros(chans,epochs);
i=1;
while i<= epochs
  dat = matsel(data,frames,0,0,i);
  dat = dat.*dat;
  ave(:,i) = sqrt(mean(dat'))';
  i = i+1;
end


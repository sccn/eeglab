% readlocsold() - Read electrode locations file in style of topoplot() or headplot().
%              Output channel information is ordered by channel numbers.
%
% Usage: >> [nums labels th r x y] = readlocsold(locfile);% {default, polar 2-D}
%        >> [nums labels th r x y] = readlocsold(locfile,'polar');     % 2-D
%        >> [nums labels x y z]    = readlocsold(locfile,'spherical'); % 3-D
%        >> [nums labels x y z]    = readlocsold(locfile,'cartesian'); % 3-D
% Input: 
%        locfile  = 'filename' of electrode location file.
%        loctype  =  File type:'polar'     2-D {default}. See >> topoplot example
%                              'spherical' 3-D. See >> headplot example
%                              'cartesian' 3-D. See >> headplot cartesian
% Outputs:
%          nums   = ordered channel numbers (from locfile)
%          labels = Matrix of 4-char channel labels (nchannels,4)
% 2-D:     r,th   = polar coordinates of each electrode (th in radians)
%          x,y    = 2-D Cartesian coordinates (from pol2cart())
% 3-D:     x,y,z  = 3-D Cartesian coordinates (normalized, |x,y,z| = 1)

% Scott Makeig, CNL / Salk Institute, La Jolla CA 3/01
%  code from topoplot()

function [channums,labels,o1,o2,o3,o4] = readlocsold(loc_file,loctype)

verbose = 0;

if nargin<2
  loctype = 'polar'; % default file type
end
icadefs
if nargin<1
  loc_file = DEFAULT_ELOC;
end

fid = fopen(loc_file);
if fid<1,
  fprintf('readlocsold(): cannot open electrode location file (%s).\n',loc_file);
  return
end
if strcmp(loctype,'spherical')| strcmp(loctype,'polar')
  A = fscanf(fid,'%d %f %f  %s',[7 Inf]);  
elseif strcmp(loctype,'cartesian')
  A = fscanf(fid,'%d %f %f %f %s',[8 Inf]);  
else
  fprintf('readlocsold(): unknown electrode location file type %s.\n',loctype);
  return
end

fclose(fid);

A = A';
channums = A(:,1);
[channums csi] = sort(channums);

if strcmp(loctype,'cartesian')
  labels = setstr(A(csi,5:8));
else
  labels = setstr(A(csi,4:7));
end
idx = find(labels == '.');                       % some labels have dots
labels(idx) = setstr(abs(' ')*ones(size(idx)));  % replace them with spaces

badchars = find(double(labels)<32|double(labels)>127);
if ~isempty(badchars)
  fprintf(...
  'readlocsold(): Bad label(s) read - Each label must have 4 chars (. => space)\n');
  return
end

for c=1:length(channums)
  if labels(c,3)== ' ' & labels(c,4)== ' '
      labels(c,[2 3]) = labels(c,[1 2]);
      labels(c,1) = ' '; % move 1|2-letter labels to middle of string
  end
end

if strcmp(loctype,'polar')
  th = pi/180*A(csi,2);      % convert degrees to radians
  rad  = A(csi,3);
  o2 = rad; o1 = th;
  [x,y] = pol2cart(th,rad);  % transform from polar to cartesian coordinates
  o3 = x; o4 = y;
elseif strcmp(loctype,'spherical')
  th = pi/180*A(csi,2);
  phi = pi/180*A(csi,3);
  x = sin(th).*cos(phi);
  y = sin(th).*sin(phi);
  z = cos(th);
  o1 = x; o2 =y; o3 = z;
elseif strcmp(loctype,'cartesian')
  x = A(csi,2);
  y = A(csi,3);
  z = A(csi,4);
  dists = sqrt(x.^2+y.^2+z.^2);
  x = y./dists;              % normalize [x y z] vector lengths to 1
  y = y./dists;
  z = z./dists;
  o1 = x; o2 =y; o3 = z;
end

if verbose
  fprintf('Location data for %d electrodes read from file %s.\n',...
                   size(A,1),loc_file);
end

if nargout<1
  fprintf('\n');
  for c=1:length(channums)
     if strcmp(loctype,'polar')
       fprintf('   %d	%s	%4.3f	%4.3f	%4.3f	%4.3f\n',...
             channums(c),labels(c,:),rad(c),th(c),x(c),y(c));
     end
  end
  fprintf('\n');
  o1 = []; o2=[]; o3=[]; o4=[];
  labels = [];
end

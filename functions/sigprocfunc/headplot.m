% headplot() - plot a spherically-splined EEG field map on a semi-realistic 
%              3-D head model. Rotate head using left mouse button.
%
% Example:
%   >> headplot example   - show an example spherical 'eloc_angles' file
%   >> headplot cartesian - show an example cartesian 'eloc_angles' file
%
% Usage (do once):
%   >> headplot('setup',['eloc_angles'],['splinefile'],['comment'],['type'])
%
% Inputs: 
%   'eloc_angles' - file of electrode locations in spherical (or cartesian) coords.
%                    cf. functions: cart2topo() and topo2sph()
%   'splinefile'  - name of spline file to save splining info into
%   'comment'     - optional string vector containing info for spline file
%   'type'        - type of electrode location file ('cartesian' or default
%                     'spherical')
%
% General usage:
%   >> headplot(values,'spline_file','Param','Value',...)
%
% Inputs:
%   values        - vector containing value at each electrode position
%   'spline_file' - spline filename computed and saved by running 'setup'
%
% Optional Parameters:
%   'title'      -  Plot title {default none}
%   'electrodes' - 'on'|'off' -> show electrode positions {default 'on'}
%   'labels'     -  2 -> plot stored electrode labels;
%                   1 -> plot channel numbers; 0 -> no labels {default}
%   'cbar'       -  0 -> Plot colorbar {default: no colorbar}
%                   H -> Colorbar axis handle (e.g., choose location)
%   'view'       - Camera viewpoint in deg. [azimuth elevation]
%                  'back'|'b'=[  0 30]; 'front'|'f'=[180 30] 
%                  'left'|'l'=[-90 30]; 'right'|'r'=[ 90 30];
%                  'frontleft'|'bl','backright'|'br', etc.,
%                  'top'=[0 90]   {default [143 18]}
%   'maplimits'  - 'absmax' -> make limits +/- the absolute-max
%                  'maxmin' -> scale to data range
%                   [min,max] -> user-definined values
%                   {default = 'absmax'; red +, blue -, green 0}
%   'lights'     - (3,N) matrix whose rows give x,y,z pos. of 
%                   N lights {default: 4 lights at corners}
%   'lighting'   - 'off' = show wire frame {default 'on'} 
%   'colormap'   -  3-column colormap matrix {default jet(64)}
%   'verbose'    - 'off' -> no msgs, no rotate3d {default 'on'}
%
% Note: if an error is generated, headplot may close the current figure
%
% Authors: Colin Humphries & Scott Makeig, SCCN/INC/UCSD, La Jolla, 1998 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) Colin Humphries and Scott Makeig, CNL / Salk Institute, Feb. 1998
% Spherical spline method: Perrin et al. (1989) Electroenceph clin Neurophys
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

% $Log: not supported by cvs2svn $
% Revision 1.5  2002/08/27 00:35:39  arno
% closing current figure
%
% Revision 1.4  2002/08/14 16:48:54  arno
% remove ICADIR
%
% Revision 1.3  2002/07/25 18:24:08  arno
% debugging
%
% Revision 1.2  2002/04/17 20:56:14  arno
% changing XYZ coordinate transformation for eloc structure
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 12-12-98 changed electrode label lines to MarkerColor -sm
% 12-12-98 added colorbar option -sm (still graphically marred by tan rect.)
% 12-13-98 implemented colorbar option using enhanced cbar -sm 
% 12-13-98 implemented 'setup' comment option -sm 
% 03-20-00 added cartesian electrode locations option -sm
% 07-14-00 fixed line in calgx() -sm from -ch
% 03-23-01 documented 'cartesian' locfile option -sm
% 01-25-02 reformated help & license, added links -ad 
% 03-21-02 added readlocs and the use of eloc input structure -ad 

function [] = headplot(values,arg1,p1,v1,p2,v2,p3,v3)

if nargin < 1
    help headplot
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Set Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

icadefs   % load definitions
set(gca,'Color',BACKCOLOR);
% mesh_file  = [ICADIR '/newupper.mat']; % whole head model file (183K)
mesh_file  = ['mhead.mat'];      % upper head model file (987K)

Lighting   = 'on';
Maplimits  = 'absmax';
Lights = [-125  125  80; ...
           125  125  80; ...
           125 -125 125; ...
          -125 -125 125];    % default lights at four corners

View       = [143 18];      % default viewpoint
Colorbar   = 0;              % default no colorbar (Note: has tan top - bug)
ColorbarAxes = 0;

HeadCenter = [0 0 30];
Cmap       = jet(64);        % default colormap
FaceColor  = [.8 .55 .35]*1.1; % ~= ruddy Caucasian - pick your complexion!
Electrodes = 'on';           % show electrode positions by default
MAX_ELECTRODES = 1024;
Elecnums   = 0;     % 0 is 'off; 1 is 'on'
Elecnames  = 0;     % 0 is 'off; 1 is 'on'
ElectDFac  = 1.06;  % plot electrode marker dots out from head surface
NamesDFac  = 1.05;  % plot electrode names/numbers out from markers
NamesColor = 'k'; % 'r';
NamesSize  =  10;   % FontSize for electrode names
MarkerColor= 'k';

sqaxis     = 1;     % if non-zero, make head proportions anatomical
title_font = 18;
titl       = [];
verbose    = 'on';
if isstr(values)
  values   = lower(values);
  if strcmp(values,'setup')
%
%%%%%%%%%%%%%%%%%%% Perform splining file setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    if nargin<3|nargin>5
       error(['headplot(): setup requires 3-5 arguments.\n'])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute and save splining matrix for the electrode array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Setting up splining matrix.\n');
    
    eloc_file = arg1;
    if nargin==5
       loctype = p2;
    else
       loctype = 'spherical';
    end
    spline_file = p1
    return
    if nargin >= 4
      comment = v1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open electrode file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isstr(eloc_file)
	    fid = fopen(eloc_file);
	    if fid == -1
	      error(['headplot(): Error opening file: ',eloc_file])
	    end
	    if strcmp(loctype,'spherical')
	      A = fscanf(fid,'%d %f %f  %s',[7 MAX_ELECTRODES]);  
	    elseif strcmp(loctype,'cartesian')
	      A = fscanf(fid,'%d %f %f %f %s',[8 MAX_ELECTRODES]);  
	    else
	      error(['headplot: unknown electrode location file type.\n']);
	    end
	    fclose(fid);
	    fprintf(['Electrode file ',eloc_file,' opened.\n'])
	    A = A';
	    fprintf('Location data for %d electrodes read.\n',size(A,1));
        
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%
	    % Record ElectrodeNames
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%
	    ElectrodeNames = zeros(size(A,1),4);
	    for r = 1:size(ElectrodeNames,1)
	    if strcmp(loctype,'spherical')
	      ElectrodeNames(r,:) = sprintf('%s',A(r,4:7));
	    elseif strcmp(loctype,'cartesian')
	      ElectrodeNames(r,:) = sprintf('%s',A(r,5:8));
	    end
	      for c=1:4
	        if ElectrodeNames(r,c) == '.'
	          ElectrodeNames(r,c) = ' ';
	        end
	      end
	      if ElectrodeNames(r,3)== ' ' & ElectrodeNames(r,4)== ' '
	            ElectrodeNames(r,[2 3]) = ElectrodeNames(r,[1 2]);
	            ElectrodeNames(r,1) = ' '; % move 1|2-letter label to center
	      end
	    end
    
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    % Convert from spherical coords to Cartesian
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    if strcmp(loctype,'spherical')
	      Th  = pi/180*A(:,2);
	      Phi = pi/180*A(:,3);
	      Xe = sin(Th).*cos(Phi);
	      Ye = sin(Th).*sin(Phi);
	      Ze = cos(Th);
	    elseif strcmp(loctype,'cartesian')
	      Xe = A(:,2);
	      Ye = A(:,3);
	      Ze = A(:,4);
	      dists = sqrt(Xe.^2+Ye.^2+Ze.^2);
	      Xe = Xe./dists;
	      Ye = Ye./dists;
	      Ze = Ze./dists;
	    end
    else
	    %%%%%%%%%%%%%%%%%%%%%
	    % Electrode structure
	    %%%%%%%%%%%%%%%%%%%%%
	    if ~isfield(eloc_file, 'X') | isempty(eloc_file(1).X) % no X Y Z coordinates
			[tmp labels Th Rd] = readlocs(eloc_file);
			ElectrodeNames = strvcat(labels); 
			labels = strvcat(labels);
            [Phi Th] = topo2sph( [Th(:) Rd(:)]);
   	        Th  = pi/180*Th;
		    Phi = pi/180*Phi;
            [Xe Ye Ze] = sph2cart( Th, Phi, ones(size(Th)));
		    fprintf('Headplot: generating XYZ coordinates from polar coordinates\n');
	        dists = sqrt(Xe.^2+Ye.^2+Ze.^2);
	        Xe = Xe./dists;
	        Ye = Ye./dists;
	        Ze = Ze./dists;
        else
		    fprintf('Headplot: using existing XYZ coordinates\n');
		    ElectrodeNames = strvcat({ eloc_file.labels });
		    Xe = cell2mat( { eloc_file.X } )';
		    Ye = cell2mat( { eloc_file.Y } )';
		    Ze = cell2mat( { eloc_file.Z } )';
	        dists = sqrt(Xe.^2+Ye.^2+Ze.^2);
	        Xe = Xe./dists;
	        Ye = Ye./dists;
	        Ze = Ze./dists;
        end;
		Xetmp = Xe;
		Xe = -Ye;
		Ye = Xetmp;
    end;
    Xe = Xe(:);
    Ye = Ye(:);
    Ze = Ze(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate g(x) for electrodes 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    enum = length(Xe);
    onemat = ones(enum,1);
    G = zeros(enum,enum);
    for i = 1:enum
      ei = onemat-((Xe(i)*onemat-Xe).^2 + (Ye(i)*onemat-Ye).^2 + ...
        (Ze(i)*onemat-Ze).^2)/2;
      gx = zeros(1,enum);
      for j = 1:enum
        gx(j) = calcgx(ei(j));
      end
      G(i,:) = gx;
    end
    fprintf('Calculating splining matrix...\n')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open mesh file - contains POS and index1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist(mesh_file)
       fprintf('headplot(): mesh_file "%s" not found.\n',mesh_file);
       fprintf('            Change file name in headplot.m.\n');
       return
    end
    eval(['load ',mesh_file,' -mat'])
    fprintf('Loaded mesh file %s\n',mesh_file);
    newPOS = POS(index1,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Project head vertices onto unit sphere
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nPOS = newPOS-ones(size(newPOS,1),1)*HeadCenter;
    spherePOS = sqrt(ones./(sum((nPOS.*nPOS)')))'*ones(1,3).*nPOS;
    x = spherePOS(:,1);
    y = spherePOS(:,2);
    z = spherePOS(:,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate new electrode positions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Computing electrode locations on head...\n');
    for i=1:length(Xe)
      elect = [Xe(i) Ye(i) Ze(i)];
      dists = distance(elect,spherePOS');
      [S,I] = sort(dists);
      npoints = I(1:3);
      diffe = newPOS(npoints,:)-spherePOS(npoints,:);
      newElect(i,:) = elect+mean(diffe)*ElectDFac;
      if Ze(i) < -1 % WAS 0 !!!! % HeadCenter(3) ????
        newElect(i,:) = [0 0 0]; % mark as electrode position not to be plotted
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate g(x) for sphere mesh vertices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Computing %d vertices. Should take a while (see wait bar)\n',...
                      length(x))
    fprintf('            but doesnt have to be done again for this montage...\n');
    hwb = waitbar(0,'Computing spline file (percent done)...');

    hwbend = length(x);
    for j = 1:length(x)
      % fprintf('%d ',j)
      X = x(j);
      Y = y(j);
      Z = z(j);
      ei = onemat-((X*onemat-Xe).^2 + (Y*onemat-Ye).^2 + (Z*onemat-Ze).^2)/2;
      for i = 1:length(ei)
        gx(j,i) = calcgx(ei(i));  
      end
      waitbar(j/hwbend)
    end
    close(hwb)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save spline file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Saving (587k) file ',spline_file])
    if nargin == 3
      eval(['save ',spline_file,' Xe Ye Ze G gx newElect ElectrodeNames'])
    else
      eval(['save ',spline_file,' Xe Ye Ze G gx newElect ElectrodeNames comment'])
    end
    fprintf('\n')
    eval(['! ls -l ',spline_file]);
    fprintf('\n')
    return

  elseif strcmp(values,'example') | strcmp(values,'demo')
%
%%%%%%%%%%%%%%%%%% Show an example electrode angles file  %%%%%%%%%%%%%%%%%%%%%%%%
%
       fprintf(['\nExample of a headplot() electrode angles file (spherical coords.)\n',...
               'Fields:  chan_num cor_deg horiz_deg channel_name\n\n',...
               '           1        -90     -72        Fp1.\n',...
               '           2         90      72        Fp2.\n',...
               '           3        -62     -57        F3..\n',...
               '           4         62      57        F4..\n',...
               '           5        -45       0        C3..\n',...
               '           6         45       0        C4..\n',...
               '           7       -118       2        A1..\n',...
               '           8        118      -2        A2..\n',...
               '           9        -62      57        P3..\n',...
               '           10        62     -57        P4..\n',...
               '           11       -90      72        O1..\n',...
               '           12        90     -72        O2..\n',...
               '           13       -90     -36        F7..\n',...
               '           14        90      36        F8..\n',...
               '           15       -90       0        T3..\n',...
               '           16        90       0        T4..\n',...
               '           17       -90      36        T5..\n',...
               '           18        90     -36        T6..\n',...
               '           19        45      90        Fz..\n',...
               '           20         0       0        Cz..\n',...
               '           21        45     -90        Pz..\n',...
             '\nA 90 deg coronal rotation points to right ear, -90 to left.\n' ,...
               'A positive horizontal rotation is counterclockwise from above.\n',...
               'Use pol2sph() to convert from topoplot() format to spherical.\n',...
               'Channel names should have 4 chars (. = space).\n\n\n']);
       return
  elseif strcmp(values,'cartesian') 
%
%%%%%%%%%%%%%%%%%% Show an example cartesian electrode file  %%%%%%%%%%%%%%%%%%%
%
       fprintf(['\nExample of a headplot() electrode location file (cartesian coords.)\n',...
               'Fields:  chan_num  x        y        z        channel_name\n\n',...
               '           1       0.4528   0.8888  -0.0694          Fp1.\n',...
               'Channel names should have 4 chars (. = space).\n\n\n']);
       return
  else
    fprintf('headplot(): Unknown first argument (%s).\n',values)
    help headplot
  end
else
%
%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
   if nargin < 2
       help headplot
       return
   end
   spline_file = arg1;
   nargs = nargin;
  if nargs > 2
     if ~(floor(nargs/2) == nargs/2)
      error('headplot(): Cannot have an odd number of inputs.')
     end
    for i = 3:2:nargs
      Param = eval(['p',int2str((i-3)/2 +1)]);
      Value = eval(['v',int2str((i-3)/2 +1)]);
      if ~isstr(Param)
         error('headplot(): Parameters must be strings.')
      end
      switch lower(Param)
        case 'cbar'
          if isstr(Value)
            error('headplot(): Colorbar value must be 0 or axis handle.')
          end
          Colorbar = 1;
          if size(Value,1) == 1 & size(Value,2) == 1
              ColorbarAxes = 0;
          else
              ColorbarAxes = Value;
          end
        case 'lighting'
          if ~isstr(Value)
            close; error('headplot(): Lighting value must be on or off.')
          end
          Value = lower(Value);
          if ~strcmp(Value,'on') & ~strcmp(Value,'off')
            close; error('headplot(): Lighting value must be on or off.')
          end
          Lighting = lower(Value);
        case 'maplimits'
          Maplimits = Value;
        case 'title'
          titl = Value;
        case 'lights'
          Lights = Value;
          if size(Lights,2) ~= 3
            close; error('headplot(): Light matrix must be (3,N).')
          end
        case 'view'
          View = Value;
        case 'verbose'
          if ~isstr(Value)
            close; error('headplot(): verbose value must be on or off.')
          end
          Value = lower(Value);
          if ~strcmp(Value,'on') & ~strcmp(Value,'off')
            close; error('headplot(): verbose value must be on or off.')
          end
          verbose = Value;
	    case {'colormap','cmap'}
	      if size(Value,2) ~= 3
	        close; error('Colormap must be an n x 3 matrix.')
	      end
	      Cmap = Value;
        case {'electrodes','elec'}
          if ~isstr(Value)
            close; error('headplot(): electrodes value must be on or off.')
          elseif ~strcmp(Value,'on') & ~strcmp(Value,'off')
            close; error('headplot(): electrodes value must be on or off.')
          end
          Electrodes = Value;
        case 'labels'
          if isstr(Value)
            close; error(['headplot(): labels value must be 0, 1, or 2.'])
          end
          if Value>2 | Value<0 
               close; error(['headplot(): labels value must be 0, 1, or 2.'])
          end
          if Value == 1
              Elecnums = 1;
          elseif Value == 2
              Elecnames = 1;
          end
        otherwise
          fprintf('headplot(): Unknown Parameter %s\n',Param)
          return
      end
    end
  end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Open head mesh and electrode spline files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~exist(spline_file)
       fprintf(...
 'headplot(): spline_file "%s" not found. Run headplot in "setup" mode.\n',...
           spline_file);
       return
  end
  eval(['load ',spline_file, ' -mat'])

  if ~exist(mesh_file)
	  close;
	  error(['headplot(): head mesh file ',meshfile,...
	   ' not found.\n  Change file name in headplot.m.\n'])
  end
  eval(['load ',mesh_file,' -mat'])
  enum = length(values);
  if enum ~= length(Xe)
	  close;
	  error(['headplot(): Number of values in spline file should equal number of electrodes.'])
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Perform interpolation
  %%%%%%%%%%%%%%%%%%%%%%%%%%

  meanval = mean(values);
  values = values - meanval;
  onemat = ones(enum,1);
  lamd = 0.1;
  C = [(G + lamd);ones(1,enum)]\[values(:);0];
  P = zeros(1,size(gx,1));
  for j = 1:size(gx,1)
    P(j) = dot(C,gx(j,:));
  end
  P = P + meanval;

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot surfaces
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  cla % clear axis
  HeadAxes = gca;

  W = zeros(1,size(POS,1));
  m = size(Cmap,1);
  if size(Maplimits) == [1,2]
     amin = Maplimits(1);
     amax = Maplimits(2);
  elseif strcmp(Maplimits,'maxmin') | strcmp(Maplimits,'minmax')
   amin = min(min(abs(P)))*1.02; % 2% shrinkage keeps within color bounds
   amax = max(max(abs(P)))*1.02; 
  elseif strcmp(Maplimits,'absmax')
   amin = -max(max(abs(P)))*1.02; % 2% shrinkage keeps within color bounds
   amax = -amin;
  end
  idx = min(m,round((m-1)*(P-amin)/(amax-amin))+1); % get colormap indices
  W(index1) = idx;
  colormap(Cmap)
  p1 = patch('Vertices',POS,'Faces',TRI1,'FaceVertexCdata',W(:),...
      'FaceColor','interp');    %%%%%%%%% plot scalp map %%%%%%%%%

  FCmap = [Cmap; Cmap(end,:); FaceColor; FaceColor; FaceColor];
  colormap(FCmap)
  W = ones(1,size(POS,1))*(m+2);
  p2 = patch('Vertices',POS,'Faces',TRI2,'FaceColor','interp',...
      'FaceVertexCdata',W(:)); %%%%%%%% plot face and lower head %%%%%%

  axis([-125 125 -125 125 -125 125])
  axis off % hide axis frame

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Draw colorbar - Note: uses enhanced cbar() function by Colin Humphries
  %%%%%%%%%%%%%%%%%%%%%%%%%      
  if Colorbar==1                
    BACKCOLOR = get(gcf,'Color');
    if ColorbarAxes == 0       
	  ColorbarHandle = cbar(0,3,[amin amax]); 
    else
	  ColorbarHandle = cbar(ColorbarAxes,3,[amin amax]); 
    end
    pos = get(ColorbarHandle,'position');  % move left & shrink to match head size
    set(ColorbarHandle,'position',[pos(1)-.05 pos(2)+0.13 pos(3)*0.7 pos(4)-0.26]);
  end
  axes(HeadAxes);

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Turn on lights
  %%%%%%%%%%%%%%%%%%%%%%%%%

  if strcmp(Lighting,'on')
    set([p1 p2],'EdgeColor','none')
    
    for i = 1:size(Lights,1)
      hl(i) = light('Position',Lights(i,:),'Color',[1 1 1],...
      'Style','infinite');
    end
    set(p2,'DiffuseStrength',.6,'SpecularStrength',0,...
    'AmbientStrength',.4,'SpecularExponent',5)
    set(p1,'DiffuseStrength',.6,'SpecularStrength',0,...
    'AmbientStrength',.3,'SpecularExponent',5)
    lighting phong  % all this gives a matte reflectance
  end  

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Set viewpoint
  %%%%%%%%%%%%%%%%%%%%%%%%%

  if isstr(View)
    switch lower(View)
      case {'front','f'}
        view(-180,30)
      case {'back','b'}
        view(0,30)
      case {'left','l'}
        view(-90,30)
      case {'right','r'}
        view(90,30)
      case {'frontright','fr'}
        view(135,30)
      case {'backright','br'}
        view(45,30)
      case {'frontleft','fl'}
        view(-135,30)
      case {'backleft','bl'}
        view(-45,30)
      case 'top'
        view(0,90)
      case 'bottom'    % undocumented option!
        view(0,-90)
        Lights = [-125 125 80;   ...
                   125 125 80;   ...
                   125 -125 125; ...
                  -125 -125 125; ...
                   0 10 -80]; % add light from below!
      otherwise
        close; error(['headplot(): Invalid View value %s',View])
    end
  else
    if ~isstr(View)
      [h,a] = size(View);
      if h~= 1 | a~=2
          close; error('headplot(): View matrix size must be (1,2).')
      end
    end
    view(View)   % set camera viewpoint
  end
  
  if strcmp(Electrodes,'on') % plot the electrode locations
   if exist('newElect')
    newNames = newElect*NamesDFac; % Calculate electrode label positions
    if Elecnames | Elecnums
      LineColor = MarkerColor; % 'y';
    else
      LineColor = MarkerColor;
    end
    for i = 1:size(newElect,1)
      if newElect(i,:) ~= [0 0 0]  % plot radial lines to electrode sites
	    line([newElect(i,1) HeadCenter(1)],[newElect(i,2) HeadCenter(2)],...
	            [newElect(i,3) HeadCenter(3)],'color',LineColor,'linewidth',1);

        if Elecnums        % plot electrode numbers
          t=text(newNames(i,1),newNames(i,2),newNames(i,3),int2str(i)); 
          set(t,'Color',NamesColor,'FontSize',NamesSize,'FontWeight','bold',...
                    'HorizontalAlignment','center');

        elseif Elecnames   % plot electrode names
         if exist('ElectrodeNames')
          name = sprintf('%s',ElectrodeNames(i,:));
          t=text(newNames(i,1),newNames(i,2),newNames(i,3),name);
          set(t,'Color',NamesColor,'FontSize',NamesSize,'FontWeight','bold',...
                    'HorizontalAlignment','center'); 
         else
           fprintf('Variable ElectrodeNames not read from spline file.\n');
         end
        else               % plot electrode markers
          line(newElect(:,1),newElect(:,2),newElect(:,3),'marker',...
                    '.','markersize',20,'color',MarkerColor,'linestyle','none');
        end
      end
    end
   else
    fprintf('Variable newElect not read from spline file.\n');
   end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Turn on rotate3d, allowing rotation of the plot using the mouse
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if strcmp(verbose,'on')
    rotate3d on;   % Allow 3-D rotation of the plot by dragging the
  else             % left mouse button while cursor is on the plot
    rotate3d off
  end              
  % Make axis square
  if sqaxis
    axis image    % keep the head proportions human and as large as possible
  end
  % Add a plot title
  if ~isempty(titl);
    % title(['\n' titl],'fontsize',title_font);
    title([titl],'fontsize',title_font); % Note: \n not interpreted by matlab-5.2
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calcgx() - function used in 'setup'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = calcgx(in)

out = 0;
m = 4;       % 4th degree Legendre polynomial
for n = 1:7  % compute 7 terms
  L = legendre(n,in);
  % out = out + ((2*n+1)/n^m*(n+1)^m)*L(1) ;
  out = out + ((2*n+1)/(n^m*(n+1)^m))*L(1) ; % bug fix by Colin 7/00
end
out = out/(4*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  distance() - function used in 'setup'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = distance(w,p)
% w is a matrix of row vectors
% p is a matrix of column vectors

l1 = size(w,1);
l2 = size(p,2);
out = zeros(l1,l2);

for i = 1:l1
  x = w(i,:)'*ones(1,l2);
  out(i,:) = sum((x-p).^2).^.5;
end

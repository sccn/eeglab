% headplot() - plot a spherically-splined EEG field map on a semi-realistic 
%              3-D head model. Can 3-D rotate the head image using the left 
%              mouse button.
% Example:
%   >> headplot example   % show an example spherical 'eloc_angles' file
%   >> headplot cartesian % show an example cartesian 'eloc_angles' file
%
% Setup usage (do only once for each scalp montage):
%
%   >> headplot('setup', elocs, splinefile, 'Param','Value',...);
%      %
%      % NOTE: previous call format below is still supported
%      %       >> headplot('setup', elocs, splinefile, comment, type);
%
% Required Setup-mode Inputs: 
%
%   elocs         - file of electrode locations (compatible with readlocs()),
%                   or EEG.chanlocs channel location structure. If the channel 
%                   file extension is not standard, use readlocs() to load the 
%                   data file, e.g.
%                      >> headplot('setup', ...
%                            readlocs('myfile.xxx', 'filetype', 'besa'),...
%                               'splinefile');
%   splinefile    - name of spline file to save spline info into. It is saved as a
%                   *.mat file and should be given the extension .spl .
%
% Optional Setup-mode Inputs:
%
%   'meshfile'    - ['string'] Matlab files containing at least two variables:
%                      POS    - 3-D positions of vertices: 
%                               x=left-right; y=back-front; z=up-down
%                      TRI1   - faces on which the scalp map should be computed
%                     plus possible optional variables:
%                     center (optional) - 3-D center of head mesh
%                     TRI2   (optional) - faces in skin color
%                     NORM   (optional) - normal for each vertex (better shading)
%   'orilocs'     - ['off'|'on'] use original electrode locations on the head
%                   {default: 'off'} (extrapolated to spherical). Note that these
%                   electrode locations must be coregisted with the head mesh.
%   'transform'   - [real array] Talairach-model transformation matrix to co-register 
%                   the electrode locations with the head mesh:
%                    [shiftX shiftY shiftZ pitch roll yaw scaleX scaleY scaleZ]
%                   The transform is applied in the order shift(rotate(scale(elocs)))
%                   by the dipfit2.* plugin function traditionaldipfit.m
%                   This array is returned by coregister().
%  'plotmeshonly' - [string] plot only mesh and electrode positions. Options are
%                   'head' to plot the standard head mesh; 'sphere' to plot the
%                   texture of the head on a sphere; 'off' not to plot anything.
%                   {default: 'off'}
%   'comment'     - ['string'] optional string containing comments for spline file
%                   {default: []}
%
% Standard-mode Usage thereafter:
%
%       >> headplot(values,'spline_file','Param','Value',...)
%
% Required Standard-mode Inputs:
%
%   values        - vector containing a data value at each electrode position
%   'spline_file' - spline filename, computed and saved in 'setup' mode (above)
%
% Optional Standard-mode Inputs:
%
%   'meshfile'   - [string] mesh file name. See file content in the setup-mode
%                  description above. {default: the EEGLAB head template file}.
%   'electrodes' - ['on'|'off'] -> show electrode positions {default 'on'}
%   'title'      -  Plot title {default: none}
%   'labels'     -  2 -> plot stored electrode labels;
%                   1 -> plot channel numbers; 0 -> no labels {default 0}
%   'cbar'       -  0 -> Plot colorbar {default: no colorbar} 
%                        Note: standard jet colormap) red = +;blue = -;green=0
%                   h -> Colorbar axis handle (to specify headplot location)
%   'view'       - Camera viewpoint in deg. [azimuth elevation]
%                  'back'|'b'=[  0 30]; 'front'|'f'=[180 30] 
%                  'left'|'l'=[-90 30]; 'right'|'r'=[ 90 30];
%                  'frontleft'|'bl','backright'|'br', etc.,
%                  'top'=[0 90],  Can rotate with mouse {default [143 18]}
%   'maplimits'  - 'absmax' -> make limits +/- the absolute-max
%                  'maxmin' -> scale to data range
%                   [min,max] -> user-definined values
%                      {default = 'absmax'}
%   'lights'     - (3,N) matrix whose rows give [x y z] pos. of each of
%                   N lights {default: four lights at corners}
%   'lighting'   - 'off' = show wire frame head {default 'on'} 
%   'colormap'   -  3-column colormap matrix {default: jet(64)}
%   'verbose'    - 'off' -> no msgs, no rotate3d {default: 'on'}
%   'orilocs'    - [channel structure or channel file name] Use original 
%                  channel locations instead of the one extrapolated from 
%                  spherical locations. Note that if you use 'orilocs'
%                  during setup, this is not necessary here since the 
%                  original channel location have already been saved.
%                  This option might be useful to show more channels than
%                  the ones actually used for interpolating (e.g., fiducials).
%   'transform'  - [real array] homogeneous transformation matrix to apply
%                  to the original locations ('orilocs') before plotting them.
%
% Note: if an error is generated, headplot() may close the current figure
%
% Authors: Arnaud Delorme, Colin Humphries, Scott Makeig, SCCN/INC/UCSD, 
%          La Jolla, 1998-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) Arnaud Delorme, Colin Humphries and Scott Makeig, 
%               CNL / Salk Institute, Feb. 1998
%
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

% 12-12-98 changed electrode label lines to MarkerColor -sm
% 12-12-98 added colorbar option -sm (still graphically marred by tan rect.)
% 12-13-98 implemented colorbar option using enhanced cbar -sm 
% 12-13-98 implemented 'setup' comment option -sm 
% 03-20-00 added cartesian electrode locations option -sm
% 07-14-00 fixed line in calgx() -sm from -ch
% 03-23-01 documented 'cartesian' locfile option -sm
% 01-25-02 reformated help & license, added links -ad 
% 03-21-02 added readlocs and the use of eloc input structure -ad 

function [HeadAxes, ColorbarHandle] = headplot(values, arg1, varargin)

if nargin < 1
    help headplot
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Set Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

icadefs   % load definitions
set(gca,'Color',BACKCOLOR);
DEFAULT_MESH      = ['mheadnew.mat'];      % upper head model file (987K)
DEFAULT_TRANSFORM = [0 -5 0 -0.1000 0 -1.5700 1040 800 950]; % stretching in different dimensions
DEFAULT_TRANSFORM = [0 -10 0 -0.1000 0 -1.600 1100 1100 1100]; % keep spherical shape.
%DEFAULT_MESH  = '/home/arno/matlab/juliehiresmesh.mat';
%DEFAULT_MESH  = ['/home/scott/matlab/old' '/newupper.mat']; % whole head model file (183K)

DEFAULT_LIGHTS = [-125  125  80; ...
                  125  125  80; ...
                  125 -125 125; ...
                  -125 -125 125];    % default lights at four corners

HeadCenter = [0 0 30];
FaceColor  = [.8 .55 .35]*1.1; % ~= ruddy Caucasian - pick your complexion!
MAX_ELECTRODES = 1024;
ElectDFac  = 1.06;  % plot electrode marker dots out from head surface
plotelecopt.NamesDFac  = 1.05;  % plot electrode names/numbers out from markers
plotelecopt.NamesColor = 'k'; % 'r';
plotelecopt.NamesSize  =  10;   % FontSize for electrode names
plotelecopt.MarkerColor= 'k';

sqaxis     = 1;     % if non-zero, make head proportions anatomical
title_font = 18;
if isstr(values)
    values   = lower(values);
    if strcmp(values,'setup')
        
%
%%%%%%%%%%%%%%%%%%% Perform splining file setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    if nargin < 3
        help headplot;
        return;
    end;
    eloc_file = arg1;
    spline_file = varargin{1};
        
    g = finputcheck(varargin(2:end), { 'orilocs'      'string'  { 'on','off' }             'off';
                                       'plotmeshonly' 'string'  { 'head','off','sphere' }  'off';
                                       'meshfile'     'string'  []                         DEFAULT_MESH;
                                       'chaninfo'     'struct'  []                         struct([]);
                                       'plotchans'    'integer' []                         [];
                                       'ica'          'string'  { 'on','off' }             'off';
                                       'transform'    'real'    []                         DEFAULT_TRANSFORM;
                                       'comment'      'string'  []                         '' });
    if isstr(g), 
        fprintf(g);
        clear g; 
        g.comment   = varargin{2}; 
        g.orilocs   = 'off';
        g.meshfile  = DEFAULT_MESH;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open electrode file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [eloc_file labels Th Rd indices] = readlocs(eloc_file);
    indices = find(~cellfun('isempty', { eloc_file.X }));
    
    % channels to plot
    % ----------------
    if isempty(g.plotchans), g.plotchans = [1:length(eloc_file)]; end;
    if ~isfield(g.chaninfo, 'nosedir'),     g.chaninfo(1).nosedir     = '+x'; end;
    indices = intersect(g.plotchans, indices);
    
    % if ICA select subset of channels if necessary
    % ---------------------------------------------
    if ~isfield(g.chaninfo, 'icachansind'), g.chaninfo(1).icachansind = 1:length(eloc_file); end;
    if strcmpi(g.ica, 'on'), 
        rmchans2 = setdiff( g.chaninfo.icachansind, indices ); % channels to remove (non-plotted) 
        newinds = 1:length(g.chaninfo.icachansind);     
        allrm = [];       
        % remove non-plotted channels from indices
        for index = 1:length(rmchans2)
            chanind = find(g.chaninfo.icachansind == rmchans2(index));
            allrm   = [ allrm chanind ];
        end;

        newinds(allrm) = [];        
        indices   = newinds;
        eloc_file = eloc_file(g.chaninfo.icachansind); 
    end;
    
    fprintf('Headplot: using existing XYZ coordinates\n');
    ElectrodeNames = strvcat({ eloc_file.labels });
    ElectrodeNames = ElectrodeNames(indices,:);
        
    Xeori = [ eloc_file(indices).X ]';
    Yeori = [ eloc_file(indices).Y ]';
    Zeori = [ eloc_file(indices).Z ]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rotate channel coordinates if necessary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(lower(g.chaninfo.nosedir), '+x')
        rotate = 0;
    else
        if strcmpi(lower(g.chaninfo.nosedir), '+y')
            rotate = 3*pi/2;
        elseif strcmpi(lower(g.chaninfo.nosedir), '-x')
            rotate = pi;
        else rotate = pi/2;
        end;
        allcoords = (Yeori + Xeori*sqrt(-1))*exp(sqrt(-1)*rotate);
        Xeori     = imag(allcoords);
        Yeori     = real(allcoords);
    end;
    newcoords = [ Xeori Yeori Zeori ];
    
    %newcoords = transformcoords( [ Xe Ye Ze ], [0 -pi/16 -1.57], 100, -[6 0 46]);
    %newcoords = transformcoords( [ Xeori Yeori Zeori ], g.transform(4:6), g.transform(7:9), g.transform(1:3));
    % same performed below with homogenous transformation matrix
    
    transmat  = traditionaldipfit( g.transform ); % arno
    newcoords = transmat*[ newcoords ones(size(newcoords,1),1)]';
    newcoords = newcoords(1:3,:)';
    
    % original center was [6 0 16] but the center of the sphere is [0 0 30] 
    % which compensate (see variable Headcenter)
    %newcoords = transformcoords( [ Xe Ye Ze ], -[0 0 -pi/6]);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % normalize with respect to head center
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newcoordsnorm      = newcoords - ones(size(newcoords,1),1)*HeadCenter;
    tmpnorm            = sqrt(sum(newcoordsnorm.^2,2));
    Xe = newcoordsnorm(:,1)./tmpnorm;
    Ye = newcoordsnorm(:,2)./tmpnorm;
    Ze = newcoordsnorm(:,3)./tmpnorm;
    %plotchans3d([ Xe Ye Ze], cellstr(ElectrodeNames)); return;
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate g(x) for electrodes 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(g.plotmeshonly, 'off')
        fprintf('Setting up splining matrix.\n');    
        enum = length(Xe);
        onemat = ones(enum,1);
        G = zeros(enum,enum);
        for i = 1:enum
            ei = onemat-sqrt((Xe(i)*onemat-Xe).^2 + (Ye(i)*onemat-Ye).^2 + ...
                             (Ze(i)*onemat-Ze).^2); % default was /2 and no sqrt
            gx = zeros(1,enum);
            for j = 1:enum
                gx(j) = calcgx(ei(j));
            end
            G(i,:) = gx;
        end
    end;
    fprintf('Calculating splining matrix...\n')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open mesh file - contains POS and index1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist(g.meshfile)
       error(sprintf('headplot(): mesh file "%s" not found\n',g.meshfile));
    end
    try, load(g.meshfile,'-mat');
    catch,
        POS  = load('mheadnewpos.txt', '-ascii');
        TRI1 = load('mheadnewtri1.txt', '-ascii'); % upper head
        %try, TRI2 = load('mheadnewtri2.txt', '-ascii'); catch, end; % lower head
        %index1 = load('mheadnewindex1.txt', '-ascii');
        center = load('mheadnewcenter.txt', '-ascii');
    end;
    if exist('index1') ~= 1, index1 = sort(unique(TRI1(:))); end;
    if exist('TRI2')   ~= 1, TRI2 = []; end;
    if exist('NORM')   ~= 1, NORM = []; end;
    if exist('TRI1')   ~= 1, error('Variable ''TRI1'' not defined in mesh file'); end;
    if exist('POS')    ~= 1, error('Variable ''POS'' not defined in mesh file'); end;
    if exist('center') ~= 1, center = [0 0 0]; disp('Using [0 0 0] for center of head mesh'); end;
    
    fprintf('Loaded mesh file %s\n',g.meshfile);
    newPOS = POS(index1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Project head vertices onto unit sphere
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spherePOS      = newPOS-ones(size(newPOS,1),1)*HeadCenter; % recenter
    nPOSnorm       = sqrt(sum(spherePOS.^2,2));
    spherePOS(:,1) = spherePOS(:,1)./nPOSnorm;
    spherePOS(:,2) = spherePOS(:,2)./nPOSnorm;
    spherePOS(:,3) = spherePOS(:,3)./nPOSnorm;
    x = spherePOS(:,1);
    y = spherePOS(:,2);
    z = spherePOS(:,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate new electrode positions on head
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(g.orilocs, 'off')
        fprintf('Computing electrode locations on head...\n');
        for i=1:length(Xe)
            elect = [Xe(i) Ye(i) Ze(i)];
            dists = distance(elect,spherePOS');
            [S,I] = sort(dists);
            npoints = I(1:3); % closest 3 points
            diffe = newPOS(npoints,:)-spherePOS(npoints,:);
            newElect(i,:) = elect+mean(diffe)*ElectDFac;
            %if Ze(i) < 0               % Plot superior electrodes only.
            %    newElect(i,:) = [0 0 0]; % Mark lower electrodes  as having
            %end                        % an electrode position not to be plotted
        end
    else 
        fprintf('Using original electrode locations on head...\n');
        newElect = newcoords;
    end;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot mesh and electrodes only
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~strcmpi(g.plotmeshonly, 'off')
        if strcmpi(g.plotmeshonly, 'sphere')
            newElect(:,1) = Xe;
            newElect(:,2) = Ye;
            newElect(:,3) = Ze;        
            POS(index1,:) = spherePOS; HeadCenter = [ 0 0 0 ];
        end;
        plotmesh(TRI1, POS, NORM);
        plotelecopt.labelflag = 0;
        plotelec(newElect, ElectrodeNames, HeadCenter, plotelecopt);
        rotate3d;
        return;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate g(x) for sphere mesh vertices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Computing %d vertices. Should take a while (see wait bar)\n',...
                      length(x))
    fprintf('            but doesnt have to be done again for this montage...\n');
    icadefs;

    gx = fastcalcgx(x,y,z,Xe,Ye,Ze);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save spline file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    comment          = g.comment;
    headplot_version = 2;
    transform        = g.transform;
    try, save(spline_file, '-V6', '-mat', 'Xe', 'Ye', 'Ze', 'G', 'gx', 'newElect', ...
              'ElectrodeNames', 'indices', 'comment', 'headplot_version', 'transform');   
    catch,
        try,  save(spline_file, '-mat', 'Xe', 'Ye', 'Ze', 'G', 'gx', 'newElect', ...
                   'ElectrodeNames', 'indices', 'comment', 'headplot_version', 'transform');
        catch, error('headplot: save spline file error, out of space or file permission problem');
        end;
    end;
    tmpinfo = dir(spline_file);
    fprintf('Saving (%dk) file %s\n',round(tmpinfo.bytes/1000), spline_file);
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
               'Channel names should have 4 chars (. = space).\n',...
               'See also >> headplot cartesian\n\n\n']);
       return
  elseif strcmp(values,'cartesian') 
%
%%%%%%%%%%%%%%%%%% Show an example cartesian electrode file  %%%%%%%%%%%%%%%%%%%
%
       fprintf(['\nExample of a headplot() electrode location file (cartesian coords.)\n',...
               'Fields:  chan_num  x        y        z      channel_name\n\n',...
               '           1       0.4528   0.8888  -0.0694       Fp1.\n',...
               'Channel names should have 4 chars (. = space).\n',...
               'See also >> headplot example\n\n\n']);
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
   
   g = finputcheck( varargin, { ...
       'cbar'       'real'   [0 Inf]         []; % Colorbar value must be 0 or axis handle.'
       'lighting'   'string' { 'on','off' }  'on';
       'verbose'    'string' { 'on','off' }  'on';
       'maplimits'  { 'string','real' }  []  'absmax'; 
       'title'      'string' []              '';
       'lights'     'real'   []              DEFAULT_LIGHTS;
       'view'       'real'   []              [143 18];
       'colormap'   'real'   []              jet(64);
       'transform'  'real'   []              [];
       'meshfile'   'string' []              DEFAULT_MESH;
       'electrodes' 'string' { 'on','off' }  'on';            
       'orilocs'    { 'string','struct' } [] '';            
       'labels'     'integer' [0 1 2]        0 }, 'headplot');
   if isstr(g) error(g); end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Open head mesh and electrode spline files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~exist(spline_file)
       error(sprintf('headplot(): spline_file "%s" not found. Run headplot in "setup" mode\n',...
           spline_file));
  end
  load(spline_file, '-mat');
  if exist('indices'), 
      try,
          values = values(indices);
      catch, error('problem of index or electrode number with splinefile'); end;
  end;
  enum = length(values);
  if enum ~= length(Xe)
	  close;
	  error('headplot(): Number of values in spline file should equal number of electrodes')
  end
  
  % change electrode if necessary
  % -----------------------------
  if ~isempty(g.orilocs)
      eloc_file = readlocs( g.orilocs );
      fprintf('Using original electrode locations on head...\n');
      indices = find(~cellfun('isempty', { eloc_file.X } ));
      newElect(:,1) = [ eloc_file(indices).X ]'; % attention inversion before
      newElect(:,2) = [ eloc_file(indices).Y ]';
      newElect(:,3) = [ eloc_file(indices).Z ]';        
      
      % optional transformation
      % -----------------------
      if ~isempty(g.transform)
          transmat  = traditionaldipfit( g.transform ); % arno
          newElect  = transmat*[ newElect ones(size(newElect,1),1)]';
          newElect  = newElect(1:3,:)';
      end;
  end;
  
  % load mesh file
  % --------------
  if ~exist(g.meshfile)
      error(sprintf('headplot(): mesh file "%s" not found\n',g.meshfile));
  end
  try, load(g.meshfile,'-mat');
  catch,
      POS  = load('mheadnewpos.txt', '-ascii');
      TRI1 = load('mheadnewtri1.txt', '-ascii'); % upper head
      % TRI2 = load('mheadnewtri2.txt', '-ascii'); % lower head
      % index1 = load('mheadnewindex1.txt', '-ascii');
      center = load('mheadnewcenter.txt', '-ascii');
  end;
  if exist('index1') ~= 1, index1 = sort(unique(TRI1(:))); end;
  if exist('TRI2')   ~= 1, TRI2 = []; end;
  if exist('NORM')   ~= 1, NORM = []; end;
  if exist('TRI1')   ~= 1, error('Variable ''TRI1'' not defined in mesh file'); end;
  if exist('POS')    ~= 1, error('Variable ''POS'' not defined in mesh file'); end;
  if exist('center') ~= 1, center = [0 0 0]; disp('Using [0 0 0] for center of head mesh'); end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Perform interpolation
  %%%%%%%%%%%%%%%%%%%%%%%%%%

  meanval = mean(values); values = values - meanval; % make mean zero
  onemat = ones(enum,1);
  lamd = 0.1;
  C = pinv([(G + lamd);ones(1,enum)]) * [values(:);0]; % fixing division error
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
  m = size(g.colormap,1);
  if size(g.maplimits) == [1,2]
      amin = g.maplimits(1);
      amax = g.maplimits(2);
  elseif strcmp(g.maplimits,'maxmin') | strcmp(g.maplimits,'minmax')
      amin = min(min(abs(P)))*1.02; % 2% shrinkage keeps within color bounds
      amax = max(max(abs(P)))*1.02; 
  elseif strcmp(g.maplimits,'absmax')
      amin = min(min(abs(P)))*1.02; % 2% shrinkage keeps within color bounds
      amax = max(max(abs(P)))*1.02; 
      amax = max(-amin, amax);
      amin = -amax;
      %amin = -max(max(abs(P)))*1.02; % 2% shrinkage keeps within color bounds
      %amax = -amin;    
  end 
  
  idx = min(m,round((m-1)*(P-amin)/(amax-amin))+1); % get colormap indices
  %subplot(1,2,1); hist(P(:));
  %idx = round((m-1)*P/(amax-amin))+m/2;
  %idx = max(1,min(m,idx)); % get colormap indices
  %subplot(1,2,2); hist(idx(:)); 
  %return;

  W(index1) = idx;
  colormap(g.colormap)
  p1 = patch('Vertices',POS,'Faces',TRI1,'FaceVertexCdata',W(:),...
      'FaceColor','interp', 'cdatamapping', 'direct', 'tag', 'mesh');    %%%%%%%%% Plot scalp map %%%%%%%%%
  if exist('NORM') == 1 & ~isempty(NORM)
      set(p1, 'vertexnormals', NORM);
  end;
  
  if ~isempty(TRI2)
      FCmap = [g.colormap; g.colormap(end,:); FaceColor; FaceColor; FaceColor];
      colormap(FCmap)
      W = ones(1,size(POS,1))*(m+2);
      p2 = patch('Vertices',POS,'Faces',TRI2,'FaceColor','interp',...
                 'FaceVertexCdata',W(:)); %%%%%%%% Plot face and lower head %%%%%%
  else 
      p2 = [];
  end;

  axis([-125 125 -125 125 -125 125])
  axis off % hide axis frame
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Draw colorbar - Note: uses enhanced cbar() function by Colin Humphries
  %%%%%%%%%%%%%%%%%%%%%%%%%      
  if ~isempty(g.cbar)
      BACKCOLOR = get(gcf,'Color');
      if g.cbar == 0       
          ColorbarHandle = cbar(0,3,[amin amax]); 
          pos = get(ColorbarHandle,'position');  % move left & shrink to match head size
          set(ColorbarHandle,'position',[pos(1)-.05 pos(2)+0.13 pos(3)*0.7 pos(4)-0.26]);
      else
          ColorbarHandle = cbar(g.cbar,3,[amin amax]); 
      end
  end
  axes(HeadAxes);

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Turn on lights
  %%%%%%%%%%%%%%%%%%%%%%%%%

  if strcmp(g.lighting,'on')
    set([p1 p2],'EdgeColor','none')
    
    for i = 1:size(g.lights,1)
      hl(i) = light('Position',g.lights(i,:),'Color',[1 1 1],...
      'Style','infinite');
    end
    if ~isempty(p2)
        set(p2,'DiffuseStrength',.6,'SpecularStrength',0,...
               'AmbientStrength',.4,'SpecularExponent',5)
    end;
    set(p1,'DiffuseStrength',.6,'SpecularStrength',0,...
    'AmbientStrength',.3,'SpecularExponent',5)
    lighting phong  % all this gives a matte reflectance
  end  

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Set viewpoint
  %%%%%%%%%%%%%%%%%%%%%%%%%

  if isstr(g.view)
    switch lower(g.view)
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
        close; error(['headplot(): Invalid View value %s',g.view])
    end
  else
      if ~isstr(g.view)
          [h,a] = size(g.view);
          if h~= 1 | a~=2
              close; error('headplot(): View matrix size must be (1,2).')
          end
      end
      view(g.view)   % set camera viewpoint
  end
  
  if strcmp(g.electrodes,'on') % plot the electrode locations
      if exist('newElect')
          plotelecopt.labelflag = g.labels;
          plotelec(newElect, ElectrodeNames, HeadCenter, plotelecopt);
      else
          fprintf('Variable newElect not read from spline file.\n');
      end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Turn on rotate3d, allowing rotation of the plot using the mouse
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if strcmp(g.verbose,'on')
    rotate3d on;   % Allow 3-D rotation of the plot by dragging the
  else             % left mouse button while cursor is on the plot
    rotate3d off
  end              
  % Make axis square
  if sqaxis
    axis image    % keep the head proportions human and as large as possible
  end
  % Add a plot title
  if ~isempty(g.title);
    % title(['\n' g.title],'fontsize',title_font);
    title([g.title],'fontsize',title_font); % Note: \n not interpreted by matlab-5.2
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
    out = out + ((2*n+1)/(n^m*(n+1)^m))*L(1);
end
out = out/(4*pi);

    
%%%%%%%%%%%%%%%%%%%
function gx = fastcalcgx(x,y,z,Xe,Ye,Ze)

onemat = ones(length(x),length(Xe));
EI = onemat - sqrt((repmat(x,1,length(Xe)) - repmat(Xe',length(x),1)).^2 +... 
                    (repmat(y,1,length(Xe)) - repmat(Ye',length(x),1)).^2 +...
                    (repmat(z,1,length(Xe)) - repmat(Ze',length(x),1)).^2);
%
gx = zeros(length(x),length(Xe));
m = 4;
icadefs;
hwb = waitbar(0,'Computing spline file (only done once)...', 'color', BACKEEGLABCOLOR);
hwbend = 7;
for n = 1:7
    L = legendre(n,EI);
    gx = gx + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
    waitbar(n/hwbend,hwb);
end
gx = gx/(4*pi);    
close(hwb);

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

% %%%%%%%%%%%%%%%
% plot electrodes
% %%%%%%%%%%%%%%%
function plotelec(newElect, ElectrodeNames, HeadCenter, opt);
    
    newNames = newElect*opt.NamesDFac; % Calculate electrode label positions
    for i = 1:size(newElect,1)
        if newElect(i,:) ~= [0 0 0]  % plot radial lines to electrode sites
            line([newElect(i,1) HeadCenter(1)],[newElect(i,2) HeadCenter(2)],...
                 [newElect(i,3) HeadCenter(3)],'color',opt.MarkerColor,'linewidth',1);
            
            if opt.labelflag == 1        % plot electrode numbers
                t=text(newNames(i,1),newNames(i,2),newNames(i,3),int2str(i)); 
                set(t,'Color',opt.NamesColor,'FontSize',opt.NamesSize,'FontWeight','bold',...
                      'HorizontalAlignment','center');
                
            elseif opt.labelflag == 2   % plot electrode names
                if exist('ElectrodeNames')
                    name = sprintf('%s',ElectrodeNames(i,:));
                    t=text(newNames(i,1),newNames(i,2),newNames(i,3),name);
                    set(t,'Color',opt.NamesColor,'FontSize',opt.NamesSize,'FontWeight','bold',...
                          'HorizontalAlignment','center'); 
                else
                    fprintf('Variable ElectrodeNames not read from spline file.\n');
                end
            else               % plot electrode markers
                line(newElect(:,1),newElect(:,2),newElect(:,3),'marker',...
                     '.','markersize',20,'color',opt.MarkerColor,'linestyle','none');
            end
        end
    end;


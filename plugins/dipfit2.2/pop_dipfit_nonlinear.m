% pop_dipfit_nonlinear() - interactively do dipole fit of selected ICA components
%
% Usage: 
%  >> EEGOUT = pop_dipfit_nonlinear( EEGIN )
%
% Inputs:
%   EEGIN       input dataset
%
% Outputs:
%   EEGOUT      output dataset
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003
%         Arnaud Delorme, SCCN, La Jolla 2003
%         Thanks to Nicolas Robitaille for his help on the CTF MEG
%         implementation

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl/

% Copyright (C) 2003 Robert Oostenveld, SMI/FCDC roberto@smi.auc.dk
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EEGOUT, com] = pop_dipfit_nonlinear( EEG, subfunction, parent, dipnum )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the code for this interactive dialog has 4 major parts
% - draw the graphical user interface
% - synchronize the gui with the data
% - synchronize the data with the gui
% - execute the actual dipole analysis

% the subfunctions that perform handling of the gui are
% - dialog_selectcomponent
% - dialog_checkinput
% - dialog_setvalue
% - dialog_getvalue
% - dialog_plotmap
% - dialog_plotcomponent
% - dialog_flip
% the subfunctions that perform the fitting are
% - dipfit_position
% - dipfit_moment

if ~plugin_askinstall('Fieldtrip-lite', 'ft_sourceanalysis'), return; end;

if nargin<1
  help pop_dipfit_nonlinear;
  return
elseif nargin==1
  
  EEGOUT = EEG;
  com = '';

  if ~isfield(EEG, 'chanlocs')
    error('No electrodes present');
  end
  
  if ~isfield(EEG, 'icawinv')
    error('No ICA components to fit');
  end
  
  if ~isfield(EEG, 'dipfit')
    error('General dipolefit settings not specified');
  end
  
  if ~isfield(EEG.dipfit, 'vol') & ~isfield(EEG.dipfit, 'hdmfile')
    error('Dipolefit volume conductor model not specified');
  end
  
  % select all ICA components as 'fitable'
  select = 1:size(EEG.icawinv,2);
  if ~isfield(EEG.dipfit, 'current')
    % select the first component as the current component
    EEG.dipfit.current = 1;
  end
  
  % verify the presence of a dipole model
  if ~isfield(EEG.dipfit, 'model')
    % create empty dipole model for each component
    for i=select
      EEG.dipfit.model(i).posxyz = zeros(2,3);
      EEG.dipfit.model(i).momxyz = zeros(2,3);
      EEG.dipfit.model(i).rv     = 1;
      EEG.dipfit.model(i).select = [1];
    end
  end
  
  % verify the size of each dipole model
  for i=select
    if ~isfield(EEG.dipfit.model, 'posxyz') | length(EEG.dipfit.model) < i | isempty(EEG.dipfit.model(i).posxyz)
      % replace all empty dipole models with a two dipole model, of which one is active
      EEG.dipfit.model(i).select = [1];
      EEG.dipfit.model(i).rv = 1;
      EEG.dipfit.model(i).posxyz = zeros(2,3);
      EEG.dipfit.model(i).momxyz = zeros(2,3);
    elseif size(EEG.dipfit.model(i).posxyz,1)==1
      % replace all one dipole models with a two dipole model
      EEG.dipfit.model(i).select = [1];
      EEG.dipfit.model(i).posxyz = [EEG.dipfit.model(i).posxyz; [0 0 0]];
      EEG.dipfit.model(i).momxyz = [EEG.dipfit.model(i).momxyz; [0 0 0]];
    elseif size(EEG.dipfit.model(i).posxyz,1)>2
      % replace all more-than-two dipole models with a two dipole model
      warning('pruning dipole model to two dipoles');
      EEG.dipfit.model(i).select = [1];
      EEG.dipfit.model(i).posxyz = EEG.dipfit.model(i).posxyz(1:2,:);
      EEG.dipfit.model(i).momxyz = EEG.dipfit.model(i).momxyz(1:2,:);
    end
  end
  
  % default is not to use symmetry constraint
  constr = [];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % construct the graphical user interface
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % define the callback functions for the interface elements
  cb_plotmap         = 'pop_dipfit_nonlinear(EEG, ''dialog_plotmap'', gcbf);';
  cb_selectcomponent = 'pop_dipfit_nonlinear(EEG, ''dialog_selectcomponent'', gcbf);';
  cb_checkinput      = 'pop_dipfit_nonlinear(EEG, ''dialog_checkinput'', gcbf);';
  cb_fitposition     = 'pop_dipfit_nonlinear(EEG, ''dialog_getvalue'', gcbf); pop_dipfit_nonlinear(EEG, ''dipfit_position'', gcbf); pop_dipfit_nonlinear(EEG, ''dialog_setvalue'', gcbf);';
  cb_fitmoment       = 'pop_dipfit_nonlinear(EEG, ''dialog_getvalue'', gcbf); pop_dipfit_nonlinear(EEG, ''dipfit_moment''  , gcbf); pop_dipfit_nonlinear(EEG, ''dialog_setvalue'', gcbf);';
  cb_close           = 'close(gcbf)';
  cb_help            = 'pophelp(''pop_dipfit_nonlinear'');';
  cb_ok              = 'uiresume(gcbf);'; 
  cb_plotdip         = 'pop_dipfit_nonlinear(EEG, ''dialog_plotcomponent'', gcbf);';
  cb_flip1           = 'pop_dipfit_nonlinear(EEG, ''dialog_flip'', gcbf, 1);';
  cb_flip2           = 'pop_dipfit_nonlinear(EEG, ''dialog_flip'', gcbf, 2);';
  cb_sym             = [ 'set(findobj(gcbf, ''tag'', ''dip2sel''), ''value'', 1);' cb_checkinput ];
  
  % vertical layout for each line 
  geomvert  =   [1 1 1 1 1 1 1 1 1]; 
  
  % horizontal layout for each line 
  geomhoriz = {
    [0.8 0.5 0.8 1 1]
    [1]
    [0.7  0.7     2 2 1]
    [0.7 0.5 0.2 2 2 1]
    [0.7 0.5 0.2 2 2 1]
    [1] 
    [1 1 1]
    [1]
    [1 1 1]
  };
  
  % define each individual graphical user element
  elements  = { ...
      { 'style' 'text'       'string' 'Component to fit'    } ...
      { 'style' 'edit'       'string' 'dummy'    'tag' 'component' 'callback' cb_selectcomponent } ...
      { 'style' 'pushbutton' 'string' 'Plot map'                   'callback' cb_plotmap } ...
      { 'style' 'text'       'string' 'Residual variance = '                                             } ...
      { 'style' 'text'       'string' 'dummy'         'tag' 'relvar'                                  } ...
      { } ...
      { 'style' 'text'    'string' 'dipole'     } ...
      { 'style' 'text'    'string' 'fit'        } ...
      { 'style' 'text'    'string' 'position'   } ...
      { 'style' 'text'    'string' 'moment'     } ...
      { } ...
      ...
      { 'style' 'text'        'string' '#1' 'tag' 'dip1'                                 } ...
      { 'style' 'checkbox'    'string' ''   'tag' 'dip1sel'    'callback' cb_checkinput  } { } ...
      { 'style' 'edit'        'string' ''   'tag' 'dip1pos'    'callback' cb_checkinput  } ...
      { 'style' 'edit'        'string' ''   'tag' 'dip1mom'    'callback' cb_checkinput  } ...
      { 'style' 'pushbutton'  'string' 'Flip (in|out)'         'callback' cb_flip1       } ...   
      ...
      { 'style' 'text'        'string' '#2' 'tag' 'dip2'                                 } ...
      { 'style' 'checkbox'    'string' ''   'tag' 'dip2sel'    'callback' cb_checkinput  } { } ...
      { 'style' 'edit'        'string' ''   'tag' 'dip2pos'    'callback' cb_checkinput  } ...
      { 'style' 'edit'        'string' ''   'tag' 'dip2mom'    'callback' cb_checkinput  } ...
      { 'style' 'pushbutton'  'string' 'Flip (in|out)'         'callback' cb_flip2       } ...   
      ...
      { } { 'style' 'checkbox' 'string' 'Symmetry constrain for dipole #2' 'tag' 'dip2sym' 'callback' cb_sym  'value' 1 } ...
      { } { } { } ...
      { 'style' 'pushbutton'  'string' 'Fit dipole(s)'' position & moment' 'callback' cb_fitposition } ...   
      { 'style' 'pushbutton'  'string' 'OR fit only dipole(s)'' moment'   'callback' cb_fitmoment   } ...   
      { 'style' 'pushbutton'  'string' 'Plot dipole(s)'                  'callback' cb_plotdip     } ...   
    };
  
  % add the cancel, help and ok buttons at the bottom
  
  geomvert  = [geomvert 1 1];
  
  geomhoriz = {geomhoriz{:} [1] [1 1 1]};
  
  elements  = { elements{:} ...
      { } ...
      { 'Style', 'pushbutton', 'string', 'Cancel', 'callback', cb_close } ...
      { 'Style', 'pushbutton', 'string', 'Help',   'callback', cb_help  } ...
      { 'Style', 'pushbutton', 'string', 'OK',     'callback', cb_ok    } ...
    };
  
  % activate the graphical interface
  supergui(0, geomhoriz, geomvert, elements{:});
  dlg = gcf;
  set(gcf, 'name', 'Manual dipole fit -- pop_dipfit_nonlinear()');
  set(gcf, 'userdata', EEG);
  pop_dipfit_nonlinear(EEG, 'dialog_setvalue', dlg);
  uiwait(dlg);
  if ishandle(dlg)
    pop_dipfit_nonlinear(EEG, 'dialog_getvalue', dlg);
    % FIXME, rv is undefined since the user may have changed dipole parameters
    % FIXME, see also dialog_getvalue subfucntion
    EEGOUT = get(dlg, 'userdata');
    close(dlg);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % implement all subfunctions through a switch-yard
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif nargin>=3
  
  %disp(subfunction)
  EEG = get(parent, 'userdata');
  
  switch subfunction
    
    case 'dialog_selectcomponent'
      current = get(findobj(parent, 'tag', 'component'), 'string');
      current = str2num(current);
      current = current(1);
      current = min(current, size(EEG.icaweights,1));
      current = max(current, 1);
      set(findobj(parent, 'tag', 'component'), 'string', int2str(current));
      EEG.dipfit.current = current;
      % reassign the global EEG object back to the dialogs userdata
      set(parent, 'userdata', EEG);
      % redraw the dialog with the current model
      pop_dipfit_nonlinear(EEG, 'dialog_setvalue', parent);
      
    case 'dialog_plotmap'
      current = str2num(get(findobj(parent, 'tag', 'component'), 'string'));
      figure; pop_topoplot(EEG, 0, current, [ 'IC ' num2str(current) ], [1 1], 1); 
      title([ 'IC ' int2str(current) ]);
      
    case 'dialog_plotcomponent'
      current = get(findobj(parent, 'tag', 'component'), 'string');
      EEG.dipfit.current = str2num(current);
      if ~isempty( EEG.dipfit.current )
        pop_dipplot(EEG, 'DIPFIT',  EEG.dipfit.current, 'normlen', 'on', 'projlines', 'on', 'mri', EEG.dipfit.mrifile);
      end;
      
    case 'dialog_checkinput'
      if get(findobj(parent, 'tag', 'dip1sel'), 'value') & ~get(findobj(parent, 'tag', 'dip1act'), 'value')
        set(findobj(parent, 'tag', 'dip1act'), 'value', 1);
      end
      if get(findobj(parent, 'tag', 'dip2sel'), 'value') & ~get(findobj(parent, 'tag', 'dip2act'), 'value')
        set(findobj(parent, 'tag', 'dip2act'), 'value', 1);
      end
      if ~all(size(str2num(get(findobj(parent, 'tag', 'dip1pos'), 'string')))==[1 3])
        set(findobj(parent, 'tag', 'dip1pos'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).posxyz(1,:)));
      else
        EEG.dipfit.model(EEG.dipfit.current).posxyz(1,:) = str2num(get(findobj(parent, 'tag', 'dip1pos'), 'string'));
      end
      if ~all(size(str2num(get(findobj(parent, 'tag', 'dip2pos'), 'string')))==[1 3])
        set(findobj(parent, 'tag', 'dip2pos'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).posxyz(2,:)));
      else
        EEG.dipfit.model(EEG.dipfit.current).posxyz(2,:) = str2num(get(findobj(parent, 'tag', 'dip2pos'), 'string'));
      end
      if ~all(size(str2num(get(findobj(parent, 'tag', 'dip1mom'), 'string')))==[1 3])
        set(findobj(parent, 'tag', 'dip1mom'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(1,:)));
      else
        EEG.dipfit.model(EEG.dipfit.current).momxyz(1,:) = str2num(get(findobj(parent, 'tag', 'dip1mom'), 'string'));
      end
      if ~all(size(str2num(get(findobj(parent, 'tag', 'dip2mom'), 'string')))==[1 3])
        set(findobj(parent, 'tag', 'dip2mom'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(2,:)));
      else
        EEG.dipfit.model(EEG.dipfit.current).momxyz(2,:) = str2num(get(findobj(parent, 'tag', 'dip2mom'), 'string'));
      end
      if get(findobj(parent, 'tag', 'dip2sel'), 'value') & get(findobj(parent, 'tag', 'dip2sym'), 'value') & ~get(findobj(parent, 'tag', 'dip1sel'), 'value')
        set(findobj(parent, 'tag', 'dip2sel'), 'value', 0);
      end
      set(parent, 'userdata', EEG);
      
    case 'dialog_setvalue'
      % synchronize the gui with the data
      set(findobj(parent, 'tag', 'component'), 'string', EEG.dipfit.current);
      set(findobj(parent, 'tag', 'relvar' ), 'string', sprintf('%0.2f%%', EEG.dipfit.model(EEG.dipfit.current).rv * 100));
      set(findobj(parent, 'tag', 'dip1sel'), 'value', ismember(1, EEG.dipfit.model(EEG.dipfit.current).select));
      set(findobj(parent, 'tag', 'dip2sel'), 'value', ismember(2, EEG.dipfit.model(EEG.dipfit.current).select));
      set(findobj(parent, 'tag', 'dip1pos'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).posxyz(1,:)));
      if strcmpi(EEG.dipfit.coordformat, 'CTF')
           set(findobj(parent, 'tag', 'dip1mom'), 'string', sprintf('%f %f %f', EEG.dipfit.model(EEG.dipfit.current).momxyz(1,:)));
      else set(findobj(parent, 'tag', 'dip1mom'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(1,:)));
      end;
      Ndipoles = size(EEG.dipfit.model(EEG.dipfit.current).posxyz, 1);
      if Ndipoles>=2
          set(findobj(parent, 'tag', 'dip2pos'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).posxyz(2,:)));
          if strcmpi(EEG.dipfit.coordformat, 'CTF')
               set(findobj(parent, 'tag', 'dip2mom'), 'string', sprintf('%f %f %f', EEG.dipfit.model(EEG.dipfit.current).momxyz(2,:)));
          else set(findobj(parent, 'tag', 'dip2mom'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(2,:)));
          end;
      end
      
    case 'dialog_getvalue'
      % synchronize the data with the gui
      if get(findobj(parent, 'tag', 'dip1sel'), 'value'); select = [1]; else select = []; end; 
      if get(findobj(parent, 'tag', 'dip2sel'), 'value'); select = [select 2]; end; 
      posxyz(1,:) = str2num(get(findobj(parent, 'tag', 'dip1pos'), 'string')); 
      posxyz(2,:) = str2num(get(findobj(parent, 'tag', 'dip2pos'), 'string')); 
      momxyz(1,:) = str2num(get(findobj(parent, 'tag', 'dip1mom'), 'string')); 
      momxyz(2,:) = str2num(get(findobj(parent, 'tag', 'dip2mom'), 'string')); 
      % assign the local values to the global EEG object
      EEG.dipfit.model(EEG.dipfit.current).posxyz = posxyz;
      EEG.dipfit.model(EEG.dipfit.current).momxyz = momxyz;
      EEG.dipfit.model(EEG.dipfit.current).select = select;
      % FIXME, rv is undefined after a manual change of parameters
      % FIXME, this should either be undated continuously or upon OK buttonpress
      % EEG.dipfit.model(EEG.dipfit.current).rv = nan;
      
      % reassign the global EEG object back to the dialogs userdata
      set(parent, 'userdata', EEG);
      
    case 'dialog_flip'
      % flip the orientation of the dipole
      current = EEG.dipfit.current;
      moment  = EEG.dipfit.model(current).momxyz;
      EEG.dipfit.model(current).momxyz(dipnum,:) = [ -moment(dipnum,1) -moment(dipnum,2) -moment(dipnum,3)];
      set(findobj(parent, 'tag', ['dip' int2str(dipnum) 'mom']), 'string', ...
        sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(current).momxyz(dipnum,:)));
      set(parent, 'userdata', EEG);
      
    case {'dipfit_moment', 'dipfit_position'}
      % determine the selected dipoles and components
      current  = EEG.dipfit.current;
      select   = find([get(findobj(parent, 'tag', 'dip1sel'), 'value') get(findobj(parent, 'tag', 'dip2sel'), 'value')]);
      if isempty(select)
        warning('no dipoles selected for fitting');
        return
      end
      % remove the dipoles from the model that are not selected, but keep
      % the original dipole model (to keep the GUI consistent)
      model_before_fitting = EEG.dipfit.model(current);
      EEG.dipfit.model(current).posxyz = EEG.dipfit.model(current).posxyz(select,:);
      EEG.dipfit.model(current).momxyz = EEG.dipfit.model(current).momxyz(select,:);
      if strcmp(subfunction, 'dipfit_moment')
        % the default is 'yes' which should only be overruled for fitting dipole moment
        cfg.nonlinear  = 'no';
      end
      dipfitdefs;
      if get(findobj(parent, 'tag', 'dip2sym'), 'value') & get(findobj(parent, 'tag', 'dip2sel'), 'value')
          if strcmpi(EEG.dipfit.coordformat,'MNI')
              cfg.symmetry = 'x';
          else
              cfg.symmetry = 'y';
          end;
      else
          cfg.symmetry = [];
      end
      
      cfg.component  = current;
      % convert structure into list of input arguments
      arg = [fieldnames(cfg)' ; struct2cell(cfg)']; 
      arg = arg(:)';
      
      % make a dialog to interrupt the fitting procedure
      fig = figure('visible', 'off');
      supergui( fig, {1 1}, [], ...
        {'style' 'text' 'string' 'Press button below to stop fitting' }, ...
        {'style' 'pushbutton' 'string' 'Interupt' 'callback' 'figure(gcbf); set(gcbf, ''tag'', ''stop'');' } );
      drawnow;
      % start the dipole fitting
      try
          warning backtrace off;
          EEG = dipfit_nonlinear(EEG, arg{:});
          warning backtrace on;
      catch,
          disp('Dipole localization failed');
      end;
      
      % should the following string be put into com? ->NOT SUPPORTED
      % --------------------------------------------------------
      com = sprintf('%s = dipfit_nonlinear(%s,%s)\n', inputname(1), inputname(1), vararg2str(arg));
      
      % this GUI always requires two sources in the dipole model
      % first put the original model back in and then replace the dipole parameters that have been fitted
      model_after_fitting = EEG.dipfit.model(current);
      newfields = fieldnames( EEG.dipfit.model );
      for index = 1:length(newfields)
      	  eval( ['EEG.dipfit.model(' int2str(current) ').' newfields{index} ' = model_after_fitting.' newfields{index} ';' ]);
      end;
      EEG.dipfit.model(current).posxyz(select,:) = model_after_fitting.posxyz;
      EEG.dipfit.model(current).momxyz(select,:) = model_after_fitting.momxyz;
      EEG.dipfit.model(current).rv               = model_after_fitting.rv;
      %EEG.dipfit.model(current).diffmap          = model_after_fitting.diffmap;

      % reassign the global EEG object back to the dialogs userdata
      set(parent, 'userdata', EEG);
      % close the interrupt dialog
      if ishandle(fig)
        close(fig);
      end
            
    otherwise
      error('unknown subfunction for pop_dipfit_nonlinear');
    end		% switch subfunction
    
  end		% if nargin
  
  
